#pragma once

#include <Math/Vector.h>
#include <Geometry/GeometryBase.h>
#include <utils.h>
#include <Math/SolidAngleSampler.h>
#include <Geometry/Materials/Lambert.h>
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

		Sphere(Math::Vector3f pos=0.0, double radius=0, Material * mat= new Lambertian<Geometry::REFLECT>(0.5)) :
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

		const GeometryBase* geometry()const final override
		{
			return this;
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


		Math::Vector2f uv(Math::Vector3f const& point)const final override
		{
			Math::Vector3f normal = (point - m_center).normalized();
			return uv(normal, true);
		}

		Math::Vector2f tuv(Math::Vector2f const& uv)const final override
		{
			return uv;
		}

		Math::Vector3f point(Math::Vector3f const& uv)const final override
		{
			return Math::Vector3f::make_sphere(uv[0] * Math::pi, uv[1] * Math::twoPi);
		}

		Math::Vector3f normal(Math::Vector3f const& point, Math::Vector2f const& uv)const final override
		{
			return (point - m_center).normalized();
		}

		Math::Vector3f shading_normal(Math::Vector3f const& point, Math::Vector2f const& uv)const final override
		{
			return normal(point, uv);
		}



		virtual bool build_lights()override
		{
			m_surface = 4 * Math::pi * m_radius_2;
			return m_material->is_emissive();
		}

		virtual void divide(unsigned int div)override
		{
			unsigned int sqrt_div = std::max(1.0, std::sqrt(div));
			if (sqrt_div * sqrt_div == div)
			{
				m_inclination_div = sqrt_div;
				m_azimuth_div = sqrt_div;
			}
			else
			{
				m_inclination_div = div;
				m_azimuth_div = 1;
			}
			
			m_divisions = m_inclination_div * m_azimuth_div;
			m_surface = 4 * Math::pi * m_radius_2;
			check_capacity(div);
		}

		virtual void sampleLight(SurfaceSample& res, Math::Sampler & sampler, unsigned int i = 0)const final override
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
			res = { pdf, this, this, uv(normal, true), normal, point };
		}

		virtual void sampleLight(SurfaceSample& res, double xi1, double xi2)const final override
		{
			double inclination = std::acos(1 - 2 * xi1);
			double azimuth = xi2 * Math::twoPi;
			Math::Vector3f normal = Math::Vector3f::make_sphere(inclination, azimuth);
			Math::Vector3f point = m_center + (normal * m_radius);
			double pdf = 1.0 / surface();
			res = { pdf, this, this, uv(normal, true), normal, point };
		}

		
		virtual void sampleLight(SurfaceSample& res, Hit const& hit, Math::Sampler& sampler, unsigned int i = 0)const final override
		{
			const Math::Vector3f cp = (m_center - hit.point);
			const double dist2 = cp.norm2();
			const double dist = std::sqrt(dist2);
			if (dist2 <= m_radius_2)
			{
				return Sphere::sampleLight(res, sampler, i);
			}
			const double sint = m_radius / dist;
			const double cost = std::sqrt(1 - sint * sint);
			const double theta = acos(cost);

			const double pdf_solid_angle = 1.0 / (Math::twoPi * (1 - cost));
			const Math::SolidAngleSampler sasampler(cp / dist, theta, cost);

			const unsigned int offset = (m_offset + i) % m_divisions;
			const unsigned int m_sub_inclination = offset / m_azimuth_div;
			const unsigned int m_sub_azimuth = offset % m_azimuth_div;
			const double xi1 = (sampler.generateContinuous<double>(m_sub_inclination, m_sub_inclination + 1)) / double(m_inclination_div);
			const double xi2 = (sampler.generateContinuous<double>(m_sub_azimuth, m_sub_azimuth + 1)) / double(m_azimuth_div);
			const Math::Vector3f sampled_dir = sasampler.generate(xi1, xi2);
			const double sampled_cos_theta = sampled_dir * cp.normalized();
			const double sampled_theta = std::acos(sampled_cos_theta);
			assert(sampled_theta <= theta);
			//Now get the point on the sphere
			double t;
			{
				const double a = 1.0;
				const double b = -2.0 * (cp * sampled_dir);
				const double c = (dist2 - m_radius_2);
				const double delta = b * b - 4.0 * a * c;
				assert(delta >= std::numeric_limits<double>::epsilon());
				if (std::abs(delta) < std::numeric_limits<double>::epsilon())
				{
					t = -b / (2.0 * a);
				}
				else
				{
					const double left = -b / (2.0 * a);
					const double right = std::sqrt(delta) / (2.0 * a);
					t = left - right;
				}
				assert(t >= 0);
			}
			const Math::Vector3f sampled_point = hit.point + sampled_dir * t;
			const Math::Vector3f sampled_normal = (sampled_point - m_center) / m_radius;
			const double area_pdf = pdf_solid_angle * std::abs(sampled_dir * sampled_normal) / (t*t);

			res = { area_pdf, this, this, uv(sampled_normal, true), sampled_normal, sampled_point};
		}

		virtual void sampleLight(SurfaceSample& res, Hit const& hit, double xi1, double xi2)const final override
		{
			const Math::Vector3f cp = (m_center - hit.point);
			const double dist2 = cp.norm2();
			const double dist = std::sqrt(dist2);
			if (dist2 <= m_radius_2)
			{
				return Sphere::sampleLight(res, xi1, xi2);
			}
			const double sint = m_radius / dist;
			const double cost = std::sqrt(1 - sint * sint);
			const double theta = acos(cost);

			const double pdf_solid_angle = 1.0 / (Math::twoPi * (1 - cost));
			const Math::SolidAngleSampler sasampler(cp / dist, theta, cost);

			const Math::Vector3f sampled_dir = sasampler.generate(xi1, xi2);
			const double sampled_cos_theta = sampled_dir * cp.normalized();
			const double sampled_theta = std::acos(sampled_cos_theta);
			assert(sampled_theta <= theta);
			//Now get the point on the sphere
			double t;
			{
				const double a = 1.0;
				const double b = -2.0 * (cp * sampled_dir);
				const double c = (dist2 - m_radius_2);
				const double delta = b * b - 4.0 * a * c;
				assert(delta >= std::numeric_limits<double>::epsilon());
				if (std::abs(delta) < std::numeric_limits<double>::epsilon())
				{
					t = -b / (2.0 * a);
				}
				else
				{
					const double left = -b / (2.0 * a);
					const double right = std::sqrt(delta) / (2.0 * a);
					t = left - right;
				}
				assert(t >= 0);
			}
			const Math::Vector3f sampled_point = hit.point + sampled_dir * t;
			const Math::Vector3f sampled_normal = (sampled_point - m_center) / m_radius;
			const double area_pdf = pdf_solid_angle * std::abs(sampled_dir * sampled_normal) / (t * t);

			res = { area_pdf, this, this, uv(sampled_normal, true), sampled_normal, sampled_point };
		}

		virtual double pdfSamplingPoint(Hit const& hit, Math::Vector3f const& point)const override
		{
			const Math::Vector3f cp = (m_center - hit.point);
			const double dist_to_center2 = cp.norm2();
			const double dist_to_center = std::sqrt(dist_to_center2);
			if (dist_to_center2 <= m_radius_2)
			{
				return 1.0 / surface();
			}
			const Math::Vector3f normal = (point - m_center) / m_radius;
			Math::Vector3f to_shaded_point = hit.point - point;
			const double dist2 = to_shaded_point.norm2();
			const double dist = std::sqrt(dist2);
			to_shaded_point /= dist;
			
			const double cos_theta_on_light = normal * to_shaded_point;
			if (cos_theta_on_light < 0) // Not visible
			{
				return 0;
			}

			// Compute the solid angle
			const double sint = m_radius / dist_to_center;
			const double cost = std::sqrt(1 - sint * sint);
			const double theta = acos(cost);
			const double phi = Math::piDiv2 - theta;
			const double cosphi = std::cos(phi);
			const double pdf_solid_angle = 1.0 / (Math::twoPi * (1 - cost));

			return pdf_solid_angle * cos_theta_on_light / dist2;
		}
		
	};
}