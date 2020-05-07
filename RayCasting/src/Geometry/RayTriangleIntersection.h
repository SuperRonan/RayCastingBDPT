#ifndef _Geometry_RayTriangleIntersection_H
#define _Geometry_RayTriangleIntersection_H

#include <Geometry/Ray.h>
#include "Triangle.h"
#include <Spy/Spy.h>
#include <assert.h>
#include "Intersection.h"
#include <Geometry/Hit.h>

namespace Geometry
{

	class RayTriangleIntersection: public Intersection<double>
	{
	protected:

		double m_u ;

		double m_v ;

		const Triangle * m_triangle;

	public:

		RayTriangleIntersection(const Triangle * triangle, const Ray & ray)
			: m_triangle(triangle)
		{
			m_valid=triangle->intersection(ray, m_t, m_u, m_v) ;
		}


		RayTriangleIntersection():
			Intersection(),
			m_triangle(nullptr)
		{

		}




		double uTriangleValue() const
		{ return m_u ; }

		double vTriangleValue() const
		{ return m_v ; }


		const Triangle * triangle() const
		{ return m_triangle ; }


		Math::Vector3f intersection() const
		{
			return m_triangle->samplePoint(m_u, m_v);
		}





		virtual void fill_hit(Hit& hit, Ray const& ray)const
		{
			assert(m_valid);
			hit.z = m_t;
			hit.geometry = m_triangle->geometry();
			hit.primitive = (const void*)m_triangle;
			hit.point = ray.sample_point(m_t);
			hit.to_view = -ray.direction();
			hit.facing = m_triangle->facing(hit.to_view);
			hit.normal = m_triangle->sampleNormal(m_u, m_v, hit.facing);
			hit.primitive_normal = m_triangle->normal(hit.facing);
			hit.primitive_uv = { m_u, m_v };
			hit.tex_uv = m_triangle->interpolateTextureCoordinate(m_u, m_v);
			hit.reflected = hit.normal * (2 * (hit.normal * hit.to_view)) - hit.to_view;
		}
	} ;
}

#endif
