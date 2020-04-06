#pragma once

#include <Integrators/RayTracingBaseIntegrator.h>


namespace Integrator
{
	class IterativePathTracingIntegrator final: public RayTracingBaseIntegrator
	{

	public:

		IterativePathTracingIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			Ray ray = pray;
			bool use_emissive = true;
			double cost = ray.direction() * scene.m_camera.m_front;
			RGBColor prod_color = scene.m_camera.We<true>(ray.direction());
			double prod_pdf = scene.m_camera.pdfWeSolidAngle<true>(pray.direction());

			RGBColor res = 0;
			for (int len = 2; len <= m_max_len; ++len)
			{
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					const Material& material = *hit.geometry->getMaterial();
					if (use_emissive && material.is_emissive())
					{
						res += prod_color * material.Le(hit.facing, hit.tex_uv) / prod_pdf;
					}

					use_emissive = hit.geometry->getMaterial()->spicky();
					if (!use_emissive && len < m_max_len)
					{
						res += prod_color * addOneDirectIllumination(scene, hit, sampler) / prod_pdf;
					}


#ifdef LATE_RUSSIAN
					double alpha = len < 4 ? 1 : m_alpha;
#else
					double alpha = m_alpha;
#endif
					double xi = sampler.generateContinuous<double>();
					if (xi < alpha)
					{
						DirectionSample next_dir;
						material.sampleBSDF(hit, m_diffuse_samples, m_specular_samples, next_dir, sampler);
						prod_color *= next_dir.bsdf * std::abs(next_dir.direction * hit.primitive_normal);
						prod_pdf *= next_dir.pdf * alpha;
						ray = Ray(hit.point, next_dir.direction);
					}
					else
					{
						break;
					}

					if (prod_color.isBlack() || prod_pdf == 0)
					{
						break;
					}
				}
				else
				{
					if(prod_pdf != 0)
						res += prod_color * scene.getBackgroundColor(ray.direction()) / prod_pdf;
					break;
				}
			}
			return res;
		}
	};



	class NaivePathTracingIntegrator final : public RayTracingBaseIntegrator
	{

	public:

		NaivePathTracingIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			Ray ray = pray;
			RGBColor prod_color = scene.m_camera.We<true>(ray.direction());
			double prod_pdf = scene.m_camera.pdfWeSolidAngle<true>(ray.direction());
			RGBColor res = 0;
			for (int len = 2; len <= m_max_len; ++len)
			{
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					const Material& material = *hit.geometry->getMaterial();
					
					res += prod_color * material.Le(hit.facing, hit.tex_uv) / prod_pdf;
					

#ifdef LATE_RUSSIAN
					double alpha = len < 4  ? 1 : m_alpha;
#else
					double alpha = m_alpha;
#endif
					double xi = sampler.generateContinuous<double>();
					if (xi < alpha)
					{
						DirectionSample next_dir;
						material.sampleBSDF(hit, m_diffuse_samples, m_specular_samples, next_dir, sampler);
						prod_color *= next_dir.bsdf * std::abs(next_dir.direction * hit.primitive_normal);
						prod_pdf *= next_dir.pdf * alpha;
						ray = Ray(hit.point, next_dir.direction);
					}
					else
					{
						break;
					}

					if (prod_color.isBlack() || prod_pdf == 0)
					{
						break;
					}
				}
				else
				{
					res += prod_color * scene.getBackgroundColor(ray.direction()) / prod_pdf;
					break;
				}
			}
			return res;
		}
	};


}