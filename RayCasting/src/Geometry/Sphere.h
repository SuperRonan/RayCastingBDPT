#pragma once

#include <Math/Vector.h>
#include <Geometry/GeometryBase.h>
#include <utils.h>
#include <Math/SolidAngleSampler.h>
#include <Geometry/Lambert.h>
#include <Geometry/Primitive.h>

namespace Geometry
{
	class Sphere : public GeometryBase, public Primitive
	{
	protected:
		Math::Vector3f m_center;
		double m_radius, m_radius_2;

		unsigned int m_inclination_div = 1;
		unsigned int m_azimuth_div = 1;

	public:

		Sphere(Math::Vector3f pos=0.0, double radius=0, Material * mat= new Lambertian(0.5)) :
			GeometryBase(mat),
			m_center(pos),
			m_radius(radius),
			m_radius_2(radius*radius)
		{
			m_surface = 4 * Math::pi * m_radius_2;
		}

		inline Math::Vector3f center()const final override
		{
			return m_center;
		}

		void setCenter(Math::Vector3f const& vec)noexcept
		{
			m_center = vec;
		}

		BoundingBox box()const override final
		{
			Math::Vector3f unit = 1;
			return BoundingBox(center() - unit * radius(), center() + unit * radius());
		}

		void setRadius(double rad)noexcept
		{
			m_radius = rad;
			m_radius_2 = rad * rad;
			m_surface = 4 * Math::pi * m_radius_2;
		}

		inline double radius()const noexcept
		{
			return m_radius;
		}

		inline double radius2()const noexcept
		{
			return m_radius_2;
		}

		inline bool facing(Math::Vector3f const& normal, Math::Vector3f const& to_view)const noexcept
		{
			return normal * to_view >= 0;
		}

		inline Math::Vector3f normal(Math::Vector3f const& point, bool facing=true)const
		{
			return (point - m_center) / (m_radius * (facing ? 1 : -1));
		}

		//assuming the normal is indeed normalized
		inline Math::Vector2f uv(Math::Vector3f normal, bool facing)const
		{
			normal *= (facing ? 1 : -1);
			return Math::Vector2f(acos(normal[2]) / Math::pi , (atan2(normal[1], normal[0]) + Math::pi) / Math::twoPi);
		}

		inline Math::Vector2f texture_uv(Math::Vector2f const& uv)const noexcept
		{
			//TODO
			return uv;
		}

		virtual bool build_lights()override
		{
			m_surface = 4 * Math::pi * m_radius_2;
			return m_material->is_emissive();
		}

		virtual void divide(unsigned int div)override
		{
			m_inclination_div = std::max(1u, unsigned int(std::sqrt(div)));
			m_azimuth_div = std::max(1u, unsigned int(div / m_inclination_div));
			m_divisions = m_inclination_div * m_azimuth_div;

			m_surface = 4 * Math::pi * m_radius_2;
			check_capacity(div);
		}

		virtual void sampleLights(LightSampleStack& res, Math::Sampler & sampler, unsigned int n=1)const final override
		{
			for (unsigned int i = 0; i < n; ++i)
			{
				Sphere::sampleLight(*(res.end()+i), sampler, i+n);
			}
			res.grow(n);
		}

		virtual void sampleLight(SurfaceLightSample& res, Math::Sampler & sampler, unsigned int i = 0)const final override
		{
			unsigned int offset = (m_offset + i) % m_divisions;
			//offset /= 2;
			unsigned int m_sub_inclination = offset / m_azimuth_div;
			unsigned int m_sub_azimuth = offset % m_azimuth_div;
			double xi1 = (sampler.generateContinuous<double>(m_sub_inclination, m_sub_inclination + 1)) / double(m_inclination_div);
			double inclination = std::acos(1 - 2 * xi1);
			double xi2 = (sampler.generateContinuous<double>(m_sub_azimuth, m_sub_azimuth + 1)) / double(m_azimuth_div);
			double azimuth = xi2 * Math::twoPi;
			Math::Vector3f normal = Math::Vector3f::make_sphere(inclination, azimuth);
			Math::Vector3f point = m_center + (normal * m_radius);
			double pdf = 1.0 / surface();
			res = { pdf, this, uv(normal, true), normal, point };
		}

	protected:

		//sample a point on the sphere that is visible from the point of he hit, according to the solid angle sampler
		void _sampleLight(SurfaceLightSample& res, Math::SolidAngleSampler const& point_sampler, Math::Sampler & sampler, double _pdf, unsigned int i=0)const
		{
			unsigned int offset = (m_offset + i) % m_divisions;
			unsigned int m_sub_inclination = offset % m_azimuth_div;
			unsigned int m_sub_azimuth = offset / m_azimuth_div;
			double xi1 = (sampler.generateContinuous<double>(m_sub_inclination, m_sub_inclination + 1)) / double(m_inclination_div);
			double xi2 = (sampler.generateContinuous<double>(m_sub_azimuth, m_sub_azimuth + 1)) / double(m_azimuth_div);
			Math::Vector3f normal = point_sampler.generate(xi1, xi2);
			Math::Vector3f point = m_center + (normal * m_radius);
			res = { _pdf, this, uv(normal, true), normal, point };
		}

	public: 


		virtual void sampleLights(LightSampleStack& res, Hit const& hit, Math::Sampler& sampler, unsigned int n)const final override
		{
			const Math::Vector3f d = (hit.point - m_center);
			const double dist = d.norm();
			double cost = m_radius / dist;
			const double theta = acos(cost);
			Math::SolidAngleSampler sasampler(d / dist, theta);
			double pdf = 1.0 / (Math::twoPi * (1 - cost) * m_radius_2);
			for (unsigned int i = 0; i < n; ++i)
			{
				SurfaceLightSample& sls = *(res.end() + i);
				Sphere::_sampleLight(sls, sasampler, sampler, pdf, i + n);
				assert(sls.normal * d >= 0);
			}
			res.grow(n);
		}

		//the hit must not be in the sphere
		virtual void sampleLight(SurfaceLightSample& res, Hit const& hit, Math::Sampler& sampler, unsigned int i = 0)const final override
		{
			const Math::Vector3f d = (hit.point - m_center);
			const double dist = d.norm();
			double cost = m_radius / dist;
			const double theta = acos(cost);
			double pdf = 1.0 / (Math::twoPi * (1 - cost) * m_radius_2);
			Math::SolidAngleSampler sasampler(d / dist, theta);

			Sphere::_sampleLight(res, sasampler, sampler, pdf, i);
		}

		virtual double pdfSamplingPoint(Hit const& hit, Math::Vector3f const& point)const override
		{
			double d = (hit.point - point).norm();
			double cost = m_radius / d;
			return 1.0 / (Math::twoPi * (1 - cost) * m_radius_2);
		}
		
	};
}