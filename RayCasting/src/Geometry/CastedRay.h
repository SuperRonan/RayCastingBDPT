#ifndef _Geometry_CastedRay
#define _Geometry_CastedRay

#pragma warning (push)
#pragma warning (disable: 4355)

#include <assert.h>
#include <Geometry/Ray.h>
#include <Geometry/Intersection.h>

namespace Geometry
{
	template <class Primitive>
	class CastedRay : public Ray
	{
	protected:
		Intersection<Primitive> m_intersection;

	public:


		CastedRay(Math::Vector3f const & source, Math::Vector3f const & direction)
			: Ray(source, direction)
		{}

		CastedRay(Ray const & ray)
			: Ray(ray)
		{}


		void intersect(const Primitive & p)
		{
			m_intersection.update(*this, p);
		}

		Intersection<Primitive> const & intersectionFound() const
		{ 
			return m_intersection ;
		}


		bool validIntersectionFound() const
		{ 
			return m_intersection.valid() ;
		}
	};
}

#pragma warning (pop)

#endif
