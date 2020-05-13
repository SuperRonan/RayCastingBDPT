#pragma once

#include <Geometry/Ray.h>
#include <Geometry/Hit.h>
#include <cassert>
#include <Geometry/Shapes/Sphere.h>
#include <Geometry/Triangle.h>
#include <Geometry/Shapes/Disk.h>

namespace Geometry
{
	template <class Primitive>
	class Intersection
	{
	public:
		double m_t = std::numeric_limits<double>::max();
		bool m_valid;

		const Primitive* m_primitive=nullptr;

		Math::Vector2f m_uv = { 0, 0 };

	public:

		Intersection(bool v=false, double t=std::numeric_limits<double>::max()):
			m_t(t),
			m_valid(v)
		{}

		template <class Other>
		bool operator<(Intersection<Other> const& i)const noexcept
		{
			return (m_valid & i.m_valid & (m_t < i.m_t)) | (!i.m_valid);
		}

		double t()const noexcept
		{
			return m_t;
		}

		bool valid()const noexcept
		{
			return m_valid;
		}

		__forceinline void update(Ray const& ray, Primitive const& p)
		{}

		__forceinline void fill(Hit & hit, Ray const& ray)const
		{}
	};

	///////////////////////////////////////////////////////////////////////
	// Sphere specialization
	///////////////////////////////////////////////////////////////////////

	template <>
	__forceinline void Intersection<Sphere>::update(Ray const& ray, Sphere const& sphere)
	{
		const Math::Vector3f oc = ray.source() - sphere.center();
		const double a = 1.0;
		const double b = (oc * ray.direction()) * 2;
		const double c = oc * oc - sphere.radius2();
		const double delta = b * b - 4 * a * c;
		if (std::abs(delta) < std::numeric_limits<double>::epsilon())
		{
			const double t = -b / (2 * a);
			if (t > 0.000000001 && t < m_t)
			{
				m_valid = true;
				m_t = t;
				m_primitive = &sphere;
			}
		}
		else if(delta > 0)
		{
			const double left = -b / (2 * a);
			const double right = std::sqrt(delta) / (2 * a);
			double t = left - right;
			if (t > 0.000000001 && t < m_t)
			{
				m_valid = true;
				m_t = t;
				m_primitive = &sphere;
				return;
			}
			t = left + right;
			if (t > 0.0000000001 && t < m_t)
			{
				m_valid = true;
				m_t = t;
				m_primitive = &sphere;
			}
		}
		return;
	}

	template <>
	__forceinline void Intersection<Sphere>::fill(Hit& hit, Ray const& ray)const
	{
		assert(m_valid);
		hit.z = m_t;
		hit.geometry = m_primitive;
		hit.primitive = m_primitive;
		hit.point = ray.sample_point(m_t);
		hit.to_view = -ray.direction();
		hit.facing = m_primitive->facing(hit.point - m_primitive->center(), hit.to_view);
		hit.normal = m_primitive->normal(hit.point);
		hit.primitive_normal = hit.normal;
		hit.primitive_uv = m_primitive->uv(hit.normal, hit.facing);
		hit.tex_uv = m_primitive->texture_uv(hit.primitive_uv);
		//hit.reflected = hit.normal * (2 * (hit.normal * hit.to_view)) - hit.to_view;
	}







	/////////////////////////////////////////////////////////////////
	// Triangle Specialization
	/////////////////////////////////////////////////////////////////

	template<>
	__forceinline void Intersection<Triangle>::update(Ray const& ray, Triangle const& tri)
	{
		const Math::Vector3f pvec = ray.direction() ^ tri.vAxis();
		const double det = tri.uAxis() * pvec;

		if (std::abs(det) < std::numeric_limits<double>::epsilon() * 1000)
		{
			return;
		}
		const double inv_det = 1.0 / det;

		const Math::Vector3f tvec = ray.source() - tri.vertex0();

		const double u = (tvec * pvec) * inv_det;
		if (std::abs(u - 0.5) > 0.5)
		{
			return;
		}

		Math::Vector3f qvec = tvec ^ tri.uAxis();

		const double v = (ray.direction() * qvec) * inv_det;
		if (v <= 0.0 || u + v > 1.0)
		{
			return;
		}
		const double t = (tri.vAxis() * qvec) * inv_det;
		if (t < 0.000000001 || t > m_t)
		{
			return;
		}
		m_valid = true;
		m_t = t;
		m_uv = { u, v };
		m_primitive = &tri;
		return;
	}

	template <>
	__forceinline void Intersection<Triangle>::fill(Hit& hit, Ray const& ray)const
	{
		assert(m_valid);
		hit.z = m_t;
		hit.geometry = m_primitive->geometry();
		hit.primitive = m_primitive;
		hit.point = ray.sample_point(m_t);
		hit.to_view = -ray.direction();
		hit.facing = m_primitive->facing(hit.to_view);
		hit.normal = m_primitive->sampleNormal(m_uv[0], m_uv[1]);
		hit.primitive_normal = m_primitive->normal();
		hit.primitive_uv = m_uv;
		hit.tex_uv = m_primitive->interpolateTextureCoordinate(m_uv[0], m_uv[1]);
		//hit.reflected = hit.normal * (2 * (hit.normal * hit.to_view)) - hit.to_view;
	}



	//////////////////////////////////////////////////////////////////////////
	// Disk specialization
	//////////////////////////////////////////////////////////////////////////
	template<>
	__forceinline void Intersection<Disk>::update(Ray const& ray, Disk const& disk)
	{
		double denum = disk.normal() * ray.direction();
		double num = disk.normal() * (disk.center() - ray.source());
		double t = num / denum;
		if (t < 0.0000000001 || t > m_t)
		{
			return;
		}
		Math::Vector3f inter_point = ray.sample_point(t);
		Math::Vector3f point_to_center = inter_point - disk.center();
		double u = point_to_center * disk.tg();
		double v = point_to_center * disk.ctg();
		if (u * u + v * v > disk.radius2())
		{
			return;
		}

		m_uv = { u, v };
		m_t = t;
		m_valid = true;
		m_primitive = &disk;
	}

	template <>
	__forceinline void Intersection<Disk>::fill(Hit& hit, Ray const& ray)const
	{
		assert(m_valid);
		hit.z = m_t;
		hit.geometry = m_primitive;
		hit.primitive = m_primitive;
		hit.point = ray.sample_point(m_t);
		hit.to_view = -ray.direction();
		hit.facing = m_primitive->facing(hit.to_view);
		hit.normal = m_primitive->normal();
		hit.primitive_normal = m_primitive->normal(hit.facing);
		hit.primitive_uv = m_uv;
		hit.tex_uv = m_uv;
		//hit.reflected = hit.normal * (2 * (hit.normal * hit.to_view)) - hit.to_view;
	}
}