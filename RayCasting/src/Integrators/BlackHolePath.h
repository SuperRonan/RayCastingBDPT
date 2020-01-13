#pragma once

#include <Integrators/RayTracingBaseIntegrator.h>

namespace Integrator
{
	class BlackHolePath final : public RayTracingBaseIntegrator
	{

	public:

		Math::Vector3f normal = Math::Vector3f(0, 0, 1).normalized();

		double delta = 0.1;

		double radius = 1;
		double intersecting_radius = 1.5 * radius;
		double min_ad_radius = 1.6 * radius;
		double max_ad_radius = 2.2 * radius;

		double mass = 1e10;

		double scene_rad = 10;

		int u_div = 10;
		int v_div = 5;
	
	protected:

		Math::Vector3f tg, ctg;

		bool intersectDisk(Ray const& ray, Math::Vector3f & point, Math::Vector2f & uv, double& t)const
		{
			double denum = normal * ray.direction();
			double num = normal * ( - ray.source());
			t = num / denum;
			if (t < 0.0000000001)
			{
				return false;
			}
			point = ray.sample_point(t);
			Math::Vector3f point_to_center = point;
			double u = point_to_center * tg;
			double v = point_to_center * ctg;
			double uv2 = u * u + v * v;
			if (uv2 > (max_ad_radius * max_ad_radius) || uv2 < (min_ad_radius * min_ad_radius))
			{
				return false;
			}
			uv[1] = (std::sqrt(uv2) - min_ad_radius) / (max_ad_radius - min_ad_radius);
			uv[0] = my_atan(v, u) / Math::twoPi;
			return true;
		}

		RGBColor diskColor(Math::Vector2f& uv)const
		{
			//return { uv[0], uv[1], 0 };
			int i = uv[0] * u_div;
			int j = uv[1] * v_div;
			if ((i % 2) ^ (j % 2))
				return RGBColor(0, 0, 1);
			else
				return 1;
		}

	public:

		BlackHolePath(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{
			tg = Math::Vector3f(1, 0, 0);
			tg = tg - normal * (normal * tg);
			std::cout << tg << std::endl;
			if (tg.norm() < std::numeric_limits<double>::epsilon() * 10)
			{
				tg = { 0, 1, 0 };
				tg = tg - normal * (normal * tg);
				std::cout << tg << std::endl;
				if (tg.norm() < std::numeric_limits<double>::epsilon() * 10)
				{
					tg = { 0, 0, 1 };
					tg = tg - normal * (normal * tg);
					std::cout << tg << std::endl;
					assert(tg.norm() > std::numeric_limits<double>::epsilon() * 10);
				}
			}
			tg = tg.normalized();
			ctg = normal ^ tg;
			ctg = ctg.normalized();
			tg = ctg ^ normal;
			tg = tg.normalized();
			std::cout <<"normal: "<< normal << std::endl;
			std::cout << tg << std::endl;
			std::cout << ctg << std::endl;
		}

		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			Ray ray = pray;
			{
				Math::Vector3f point;
				Math::Vector2f uv;
				double t;
				if (intersectDisk(ray, point, uv, t))
				{
					return diskColor(uv);
				}
			}
			return 0;
		}
	};
}