#pragma once

#include <mutex>
#include <vector>
#include <Math/Vectorf.h>
#include <Geometry/BoundingBox.h>
#include <functional>

namespace Integrator
{
	// Guideline for the photons
	template <class Float>
	class PhotonBase
	{
	public:
		Math::Vector<Float, 3> point()const;
	};

	template <class PhotonType, class Float=double>
	class PhotonMap
	{
	protected:

		template <class T>
		using Collection = std::vector<T>;

		std::mutex m_mutex;
		
		Math::Vector<Float, 3> m_cell_size;
		Math::Vector<int, 3> m_size;
		std::vector<Collection<PhotonType>> m_map;
		Geometry::BoundingBox m_bb;

	public:

		PhotonMap()
		{}

		void init(Geometry::BoundingBox const& bb, Math::Vector<int, 3> const& resolution)
		{
			m_bb = bb;
			m_size = resolution;

			Math::Vector<Float, 3> diag = bb.diag();
			m_cell_size = diag.simdDiv(Math::Vector<Float, 3>(resolution));

			int N = m_size.prod();
			m_map.resize(N);
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

		void addPhoton(PhotonType const& photon)
		{
			int id = index(photon.point());
			Collection<PhotonType>& list = m_map[id];
			m_mutex.lock();
			list.push_back(photon);
			m_mutex.unlock();
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

	};
}