#pragma once

#include <Integrators/RayTracingBaseIntegrator.h>


namespace Integrator
{
	
	class RayTracingIntegrator final : public RayTracingBaseIntegrator
	{

	public:

		RayTracingIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height):
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{}

		
		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final
		{
			Ray ray = pray;
			unsigned int len = 1;
			RGBColor res=0;
			RGBColor T=1;
			while (len <= m_max_len)
			{
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					++len;
					Material const& material = *hit.geometry->getMaterial();
					
					res += T * material.Le(hit.primitive_normal, hit.tex_uv, hit.to_view);

					if (!material.spicky())
					{
						res += T * RayTracingBaseIntegrator::addOneDirectIllumination(scene, hit, sampler);
						break;
					}
					else
					{
						DirectionSample next_dir;
						material.sampleBSDF(hit, next_dir, sampler);

						ray = Ray(hit.point, next_dir.direction);
						T *= next_dir.bsdf * std::abs(next_dir.direction * hit.primitive_normal) / next_dir.pdf;

						if (T.isBlack())
						{
							break;
						}
					}
				}
				else
				{
					res = T * scene.getBackgroundColor(ray.direction());
					break;
				}
			}
			return res;
		}
		
	};
	
}