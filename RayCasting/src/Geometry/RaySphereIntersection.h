#pragma once

#include "Ray.h"
#include "Sphere.h"
#include <Math/Vector.h>
#include <limits>
#include "Intersection.h"
namespace Geometry
{
	class RaySphereIntersection : public Intersection
	{
	protected:

		const Sphere * m_sphere;

	public:

		RaySphereIntersection():
			Intersection(false),
			m_sphere(nullptr)
		{

		}

		RaySphereIntersection(RaySphereIntersection const& other):
			Intersection(other.valid(), other.t()),
			m_sphere(other.m_sphere)
		{}

		RaySphereIntersection(Ray const& ray, Sphere const& sphere):
			m_sphere(&sphere)
		{
			Math::Vector3f oc = ray.source() - sphere.center();
			const double a = 1.0; 
			const double b = (oc * ray.direction())*2;
			const double c = oc * oc - sphere.radius2();
			const double delta = b * b - 4 * a * c; 
			if (std::abs(delta) < std::numeric_limits<double>::epsilon())//delta ~= 0
			{
				
				m_t = (-b) / (2*a);
				m_valid = m_t > 0.00000001;
				//m_valid = false;
			}
			else if(delta > 0)
			{
				const double left = -b / (2*a);
				const double right = std::sqrt(delta) / (2*a);
				m_t = left - right;
				if (m_t > 0.00000001)
				{
					m_valid = true;
				}
				else
				{
					m_t = left + right;
					m_valid = m_t > 0.00000001;
				}
			}
			else
			{
				m_valid = false;
			}
		}


		const Sphere * sphere()const noexcept
		{
			return m_sphere;
		}


		RaySphereIntersection & operator=(RaySphereIntersection const& other) = default;


		virtual void fill_hit(Hit& hit, Ray const& ray)const
		{
			assert(m_valid);
			hit.z = m_t;
			hit.geometry = m_sphere;
			hit.primitive = (const void*)m_sphere;
			hit.point = ray.sample_point(m_t);
			hit.to_view = -ray.direction();
			hit.facing = m_sphere->facing(hit.point - m_sphere->center(), hit.to_view);
			hit.normal = m_sphere->normal(hit.point, hit.facing);
			hit.primitive_normal =hit.normal;
			hit.primitive_uv = m_sphere->uv(hit.normal, hit.facing);
			hit.tex_uv = m_sphere->texture_uv(hit.primitive_uv);
			hit.reflected = hit.normal * (2 * (hit.normal * hit.to_view)) - hit.to_view;
		}

		

	};
}