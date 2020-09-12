#pragma once

#include <vector>
#include <cassert>
#include <Math/Vectorf.h>
#include <Geometry/BoundingBox.h>

namespace Integrator
{
	template <class Float, class Int=int64_t>
	class VisibilityCache
	{
	protected:

		struct CacheInfo
		{
			Int sucess, failed;
		};

		using Vector3f = Math::Vector<Float, 3>;
		using Vector3i = Math::Vector<Int, 3>;

		// Cells are doubled
		std::vector<CacheInfo> m_data;
		Vector3i m_size;
		Geometry::BoundingBox m_bb;
		Vector3f m_cell_size;

		Int m_3D_size;

	public:

		VisibilityCache()
		{}

		void init(Geometry::BoundingBox const& bb, Vector3i const& resolution)
		{
			m_bb = bb;
			m_size = resolution;
			m_cell_size = bb.diag().simdDiv(Vector3f(resolution));

			m_3D_size = m_size.prod();
			m_data.resize(m_3D_size * m_3D_size);
			m_data.shrink_to_fit();
			reset();
		}

		bool isInitialized()const
		{
			return m_data.size();
		}

		void reset()
		{
			std::fill(m_data.begin(), m_data.end(), CacheInfo{ 1, 0 });
		}

		void free()
		{
			m_size = 0;
			m_cell_size = 0;
			m_data.resize(0);
			m_data.shrink_to_fit();
		}

		Int index_int(Vector3i const& ijk)const
		{
			return ijk[2] + m_size[2] * (ijk[1] + ijk[0] * m_size[1]);
		}

		Vector3i cell_index(Vector3f const& xyz)const
		{
			assert(m_bb.isInside(xyz));
			return m_size.simdMin(
				Vector3i(0, 0, 0).simdMax(
					Vector3i(
						(xyz - m_bb[0]).simdDiv(m_cell_size)
					)
				)
			);
		}

		Int index(Vector3f const& xyz)const
		{
			return index_int(cell_index(xyz));
		}

		bool inside(Vector3i const& ijk)const
		{
			return ijk[0] >= 0 && ijk[0] < m_size[0] &&
				ijk[1] >= 0 && ijk[1] < m_size[1] &&
				ijk[2] >= 0 && ijk[2] < m_size[2];
		}

		Float estimateSucess(Vector3f const& a, Vector3f const& b)const
		{
			Int a_index = index(a), b_index = index(b);
			assert(a_index >= 0 && a_index < m_3D_size);
			assert(b_index >= 0 && b_index < m_3D_size);
			Int sucess = 0, failed = 0;
			const CacheInfo& cell1 = m_data[a_index * m_3D_size + b_index];
			sucess += cell1.sucess;
			failed += cell1.failed;
			const CacheInfo& cell2 = m_data[b_index * m_3D_size + a_index];
			sucess += cell2.sucess;
			failed += cell2.failed;
			return Float(sucess) / Float(sucess + failed);
		}

		void report(Vector3f const& a, Vector3f const& b, Int sucess, Int failed)
		{
			Int a_index = index(a), b_index = index(b);
			assert(a_index >= 0 && a_index < m_3D_size);
			assert(b_index >= 0 && b_index < m_3D_size);
			// Only update one of the two cells
			Int cell_index = a_index * m_3D_size + b_index;
			assert(cell_index < m_3D_size* m_3D_size);
			CacheInfo& cell1 = m_data[cell_index];
			cell1.sucess += sucess;
			cell1.failed += failed;
		}

		void report(Vector3f const& a, Vector3f const& b, bool sucess)
		{
			report(a, b, sucess ? 1 : 0, sucess ? 0 : 1);
		}

		void reportSucess(Vector3f const& a, Vector3f const& b)
		{
			report(a, b, 1, 0);
		}

		void reportFailure(Vector3f const& a, Vector3f const& b)
		{
			report(a, b, 0, 1);
		}

	};
}