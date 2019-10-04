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
			unsigned int depth = 0;
			RGBColor res=0;
			RGBColor T=1;
			while (depth <= m_max_depth)
			{
				++depth;
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					Material const& material = *hit.geometry->getMaterial();
					
					res += T * material.Le(hit.facing, hit.tex_uv);

					if (material.use_direct())
					{
						res += T * RayTracingBaseIntegrator::addDirectIllumination(scene, hit, sampler);
						break;
					}
					else
					{
						DirectionSample next_dir;
						material.sampleBSDF(hit, m_diffuse_samples, m_specular_samples, next_dir, sampler);

						ray = Ray(hit.point, next_dir.direction);
						T *= next_dir.bsdf * (next_dir.direction * hit.primitive_normal) / next_dir.pdf;

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