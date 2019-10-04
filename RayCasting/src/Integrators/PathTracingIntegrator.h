#pragma once

#include <Integrators/RayTracingBaseIntegrator.h>


namespace Integrator
{
	//just less efficient than the iterative path integrator
	//but you can more than one samle (never happens)
	// I think I don't use it anymore
	class RecursivePathTracingIntegrator final: public RayTracingBaseIntegrator
	{
	protected:


		RGBColor addIndirectIllumination(Scene const& scene, Hit const& hit, bool next_em, int depth, Math::Sampler& sampler)const
		{
			RGBColor res = 0;
			DirectionStack ds;
			hit.geometry->getMaterial()->sampleBSDF(hit, m_diffuse_samples, m_specular_samples, ds, sampler);

			for (DirectionSample const& dir : ds)
			{
				const RGBColor& bsdf = dir.bsdf;
				const double& pdf = dir.pdf;
				const Math::Vector3f& direction = dir.direction;
				//const double& weight = dir.weight;
				assert(pdf > 0);

				RGBColor prod = bsdf * ((abs(direction * hit.primitive_normal)) / pdf);
				if (!prod.isBlack())
				{
					//TODO check for the separation of lights

					RGBColor Li = tracePath(scene, Ray(hit.point, direction), sampler, next_em, hit.primitve, depth + 1);

					res += prod * Li;

				}
			}
			return res;
		}



		RGBColor tracePath(Scene const& scene, Ray const& ray, Math::Sampler& sampler, bool use_emissive = true, const void* current = nullptr, unsigned int depth = 0)const 
		{
			if (depth > m_max_depth)
			{
				return 0;
			}
			Hit hit;
			if (scene.full_intersection(ray, hit))
			{

				RGBColor res = use_emissive ? hit.geometry->getMaterial()->Le(hit.facing, hit.tex_uv) : 0;


#ifdef SAMPLE_DIRECT
				const bool use_direct = hit.geometry->getMaterial()->use_direct();
#else
				const bool use_direct = false;
#endif
				if (use_direct)
				{
					RGBColor direct = addDirectIllumination(scene, hit, sampler);
					res += direct;
				}


#ifdef SHORT_RUSSIAN
				double alpha = depth < 2 ? 1 : m_alpha;
#else
				double alpha = m_alpha;
#endif
				double xi = sampler.generateContinuous<double>();
				if (xi < alpha)
				{
					const bool use_next_em = !use_direct;



					RGBColor indirect = addIndirectIllumination(scene, hit, use_next_em, depth, sampler) / alpha;

					res += indirect;

				}

				return res;
			}
			else
			{
				return scene.getBackgroundColor(ray.direction());
			}
		}

	public: 

		RecursivePathTracingIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& ray, Math::Sampler& sampler)const final override
		{
			return tracePath(scene, ray, sampler);
		}
	};



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
			RGBColor prod_color = scene.m_camera.We(ray.direction()) * (1 / (cost * cost * cost));
			double prod_pdf = scene.m_camera.pdfWeSolidAngle(pray.direction());
			RGBColor res = 0;
			for (int depth = 0; depth < m_max_depth + 1; ++depth)
			{
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					const Material& material = *hit.geometry->getMaterial();
					if (use_emissive && material.is_emissive())
					{
						res += prod_color * material.Le(hit.facing, hit.tex_uv) / prod_pdf;
					}

					use_emissive = !hit.geometry->getMaterial()->use_direct();
					if (!use_emissive)
					{
						res += prod_color * addDirectIllumination(scene, hit, sampler) / prod_pdf;
					}


#ifdef SHORT_RUSSIAN
					double alpha = depth < 3 ? 1 : m_alpha;
#else
					double alpha = m_alpha;
#endif
					double xi = sampler.generateContinuous<double>();
					if (xi < alpha)
					{
						DirectionSample next_dir;
						material.sampleBSDF(hit, m_diffuse_samples, m_specular_samples, next_dir, sampler);
						prod_color *= next_dir.bsdf * (next_dir.direction * hit.primitive_normal);
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



	class NaivePathTracingIntegrator final : public RayTracingBaseIntegrator
	{

	public:

		NaivePathTracingIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			Ray ray = pray;
			RGBColor prod_color = 1;
			double prod_pdf = 1;
			RGBColor res = 0;
			for (int depth = 0; depth <= m_max_depth + 1; ++depth)
			{
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					const Material& material = *hit.geometry->getMaterial();
					
					res += prod_color * material.Le(hit.facing, hit.tex_uv) / prod_pdf;
					

#ifdef SHORT_RUSSIAN
					double alpha = depth < 3  ? 1 : m_alpha;
#else
					double alpha = m_alpha;
#endif
					double xi = sampler.generateContinuous<double>();
					if (xi < alpha)
					{
						DirectionSample next_dir;
						material.sampleBSDF(hit, m_diffuse_samples, m_specular_samples, next_dir, sampler);
						prod_color *= next_dir.bsdf * (next_dir.direction * hit.primitive_normal);
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