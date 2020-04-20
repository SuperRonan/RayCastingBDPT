#pragma once

#include <mutex>
#include <vector>
#include <Math/Vectorf.h>
#include <Geometry/BoundingBox.h>
#include <functional>
#include <limits>

namespace Integrator
{
	// Guideline for the photons
	template <class Float>
	class PhotonBase
	{
	public:
		Math::Vector<Float, 3> point()const;
	};
	struct PhotonId
	{
		unsigned int cell = 0;
		unsigned int id = 0;

		static constexpr PhotonId invalid()
		{
			PhotonId res;
			res.cell = res.id = std::numeric_limits<unsigned int>::max();
			return res;
		}

		bool isInvalid()const
		{
			return cell == std::numeric_limits<unsigned int>::max() && id == std::numeric_limits<unsigned int>::max();
		}

		bool isValid()const
		{
			return cell != std::numeric_limits<unsigned int>::max() || id != std::numeric_limits<unsigned int>::max();
		}
	};
	template <class PhotonType, class Float=double>
	class PhotonMap
	{
	public:
		
		template <class T>
		using Collection = std::vector<T>;

		std::mutex m_mutex;
		
		Math::Vector<Float, 3> m_cell_size;
		Math::Vector<int, 3> m_size;
		std::vector<Collection<PhotonType>> m_map;
		Geometry::BoundingBox m_bb;

		bool m_built = false;

	public:

		PhotonMap()
		{}

		void setParams(Geometry::BoundingBox const& bb, Math::Vector<int, 3> const& resolution)
		{
			m_bb = bb;
			m_size = resolution;

			Math::Vector<Float, 3> diag = bb.diag();
			m_cell_size = diag.simdDiv(Math::Vector<Float, 3>(resolution));
		}

		void init()
		{
			int N = m_size.prod();
			m_map.clear();
			m_map.resize(N);
			m_map.shrink_to_fit();
			m_built = false;
		}

		void init(Geometry::BoundingBox const& bb, Math::Vector<int, 3> const& resolution)
		{
			setParams(bb, resolution);
			init();
		}

		void freeMemory()
		{
			clear();
			m_map.shrink_to_fit();
		}

		void clear()
		{
			m_map.clear();
			m_built = false;
		}

		int index_int(Math::Vector<int, 3> const& ijk)const
		{
			return ijk[2] + m_size[2] * (ijk[1] + ijk[0] * m_size[1]);
		}

		Math::Vector<Float, 3> cell_index(Math::Vector<Float, 3> const& xyz)const
		{
			return m_size.simdMin(Math::Vector<int, 3>(0, 0, 0).simdMax((Math::Vector<int, 3>((xyz - m_bb[0]).simdDiv(m_cell_size)))));
		}

		int index(Math::Vector<Float, 3> const& xyz)const
		{
			return index_int(cell_index(xyz));
		}

		bool insideMap(Math::Vector<int, 3> const& cell_index)const
		{
			return cell_index[0] >= 0 && cell_index[0] < m_size[0] &&
				cell_index[1] >= 0 && cell_index[1] < m_size[1] &&
				cell_index[2] >= 0 && cell_index[2] < m_size[2];
		}

		PhotonId addPhoton(PhotonType const& photon)
		{
			PhotonId res;
			res.cell = index(photon.point());
			Collection<PhotonType>& list = m_map[res.cell];
			m_mutex.lock();
			list.push_back(std::move(photon));
			res.id = list.size() - 1;
			m_mutex.unlock();
			return res;
		}

		PhotonType const& operator[](PhotonId const& id)const
		{
			return m_map[id.cell][id.id];
		}

		PhotonType & operator[](PhotonId const& id)
		{
			return m_map[id.cell][id.id];
		}

		void buildDone()
		{
			m_built = true;
		}

		bool built()const
		{
			return m_built;
		}

		size_t size()const
		{
			return m_map.size();
		}

		template <class Function>
		__forceinline void loopThroughPhotons(const Function& function, Math::Vector<Float, 3> const& point)const
		{
			const Math::Vector<int, 3> ijk = cell_index(point);
			for (int i = -1; i <= 1; ++i)
			{
				for (int j = -1; j <= 1; ++j)
				{
					for (int k = -1; k <= 1; ++k)
					{
						const Math::Vector<int, 3> cell = ijk + Math::Vector<int, 3>(i, j, k);
						if (insideMap(cell))
						{
							const Collection<PhotonType>& photons = m_map[index_int(cell)];
							for (PhotonType const& photon : photons)
							{
								function(photon);
							}
						}
					}
				}
			}
		}

		template <class Function>
		__forceinline void loopThroughPhotons(const Function& function, Math::Vector<Float, 3> const& point)
		{
			const Math::Vector<int, 3> ijk = cell_index(point);
			for (int i = -1; i <= 1; ++i)
			{
				for (int j = -1; j <= 1; ++j)
				{
					for (int k = -1; k <= 1; ++k)
					{
						const Math::Vector<int, 3> cell = ijk + Math::Vector<int, 3>(i, j, k);
						if (insideMap(cell))
						{
							Collection<PhotonType>& photons = m_map[index_int(cell)];
							for (PhotonType & photon : photons)
							{
								function(photon);
							}
						}
					}
				}
			}
		}

		template <class Function>
		__forceinline void loopThroughPhotons(const Function& function, Math::Vector<Float, 3> const& point, Float custom_radius)const
		{
			Math::Vector<int, 3> range = Math::Vector<Float, 3>(custom_radius).simdDiv(m_cell_size).ceil();
			const Math::Vector<int, 3> ijk = cell_index(point);
			for (int i = -range[0]; i <= range[0]; ++i)
			{
				for (int j = -range[1]; j <= range[1]; ++j)
				{
					for (int k = -range[2]; k <= range[2]; ++k)
					{
						const Math::Vector<int, 3> cell = ijk + Math::Vector<int, 3>(i, j, k);
						if (insideMap(cell))
						{
							const Collection<PhotonType>& photons = m_map[index_int(cell)];
							for (PhotonType const& photon : photons)
							{
								function(photon);
							}
						}
					}
				}
			}
		}

		template <class Function>
		__forceinline void loopThroughPhotons(const Function& function, Math::Vector<Float, 3> const& point, Float custom_radius)
		{
			Math::Vector<int, 3> range = Math::Vector<Float, 3>(custom_radius).simdDiv(m_cell_size).ceil();
			const Math::Vector<int, 3> ijk = cell_index(point);
			for (int i = -range[0]; i <= range[0]; ++i)
			{
				for (int j = -range[1]; j <= range[1]; ++j)
				{
					for (int k = -range[2]; k <= range[2]; ++k)
					{
						const Math::Vector<int, 3> cell = ijk + Math::Vector<int, 3>(i, j, k);
						if (insideMap(cell))
						{
							const Collection<PhotonType>& photons = m_map[index_int(cell)];
							for (PhotonType & photon : photons)
							{
								function(photon);
							}
						}
					}
				}
			}
		}

	};
}