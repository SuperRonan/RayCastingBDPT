#pragma once

#include <Integrators/RayTracingBaseIntegrator.h>
#include <cassert>


namespace Integrator
{
	class MISPathTracingIntegrator : public RayTracingBaseIntegrator
	{
	protected:

	public:

		MISPathTracingIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{}


		static __forceinline void balance(__out double& w0, __out double& w1, const double p00, const double p01, const double p10, const double p11)
		{
			w0 = p00 / (p00 + p10);
			w1 = p11 / (p01 + p11);
		}

		static __forceinline void power(__out double& w0, __out double& w1, const double p00, const double p01, const double p10, const double p11, double beta = 2)
		{
			w0 = pow(p00, beta) / (pow(p00, beta) + pow(p10, beta));
			w1 = pow(p11, beta) / (pow(p01, beta) + pow(p11, beta));
		}


		RGBColor MISAddDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler)const
		{
			RGBColor res, light_contribution, surface_contribution;
			double light_pdf = 1, surface_pdf, total_pdf;
			bool multi_light = false; // To sample all the lights at a time
			//switching from point to direction from the hit

			// Light sampling ///////////////////////////////////////////
			StackN<double> distStack;
			SurfaceLightSample sls;

			LightSampleStack ps;
			if (multi_light)
			{
				for (const GeometryBase* const sl : scene.m_surface_lights)
				{
					if (sl != hit.geometry)
					{
						sl->sampleLights(ps, sampler, m_direct_samples); //old
						//sl->sampleLights(ps, hit, sampler, m_direct_samples); //new clever but biased (need to correct it)
					}
				}
			}
			else
			{
				sampleOneLight(scene, hit, sampler, sls);
				ps.push(sls);
			}

			for (SurfaceLightSample& sls : ps)
			{
				Math::Vector3f dir = (sls.vector - hit.point);
				double dist = dir.norm();
				distStack.push(dist);
				sls.vector = dir / dist;
			}

			ColorStack bsdfs;
			//bsdfs.push(material.BSDF(hit, sls.vector));
			hit.geometry->getMaterial()->BSDF(hit, bsdfs, ps);
			for (size_t i = 0; i < bsdfs.size(); ++i)
			{
				RGBColor const& bsdf = bsdfs[i];
				const double cosi = std::abs(hit.primitive_normal * ps[i].vector);
				RGBColor prod = bsdf * cosi;
				if (!prod.isBlack()) //or almost black ???
				{
					Hit light_hit;
					Ray ray(hit.point, ps[i].vector);
					if (scene.full_intersection(ray, light_hit) && samePoint(light_hit, distStack[i], ps[i].geo))
					{
						RGBColor color = light_hit.geometry->getMaterial()->Le(light_hit.facing, light_hit.tex_uv);
						double dist2 = light_hit.z * light_hit.z;
						double cosl = light_hit.primitive_normal * (-ray.direction().normalized());
						light_pdf *= ps[i].pdf;
						light_contribution += prod * cosl * color / (ps[i].pdf * m_direct_samples * dist2);
					}
				}
			}

			// Surface sampling ///////////////////////////////////////////
			//for (size_t i = 0; i < bsdfs.size(); ++i)
			{
				//RGBColor const& bsdf = bsdfs[i];
				// sample direction
				Geometry::DirectionSample dir;
				hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, dir, sampler);
				RGBColor const& bsdf = hit.geometry->getMaterial()->BSDF(hit, dir.direction.normalized());
				const double cosi = std::abs(hit.primitive_normal * dir.direction.normalized());
				RGBColor prod = bsdf * cosi;

				Math::Vector3f dir2 = (dir.direction - hit.point);
				double dist = dir2.norm();

				if (!prod.isBlack()) //or almost black ???
				{
					Hit surface_hit;
					Ray ray(hit.point, dir.direction.normalized());
					if (scene.full_intersection(ray, surface_hit) && dir.direction * ray.direction() > 0 /*&& samePoint(surface_hit, dist, hit.geometry)*/)
					{
						double dist2 = surface_hit.z * surface_hit.z;
						double cosl = surface_hit.primitive_normal * (-ray.direction());
						surface_pdf = dir.pdf;
						RGBColor color;
						if (surface_hit.geometry->getMaterial()->is_emissive())
						{
							color = surface_hit.geometry->getMaterial()->Le(surface_hit.facing, surface_hit.tex_uv);
						}
						surface_contribution = prod * cosl * color / (dir.pdf * m_direct_samples * dist2) * dist2 / cosl;
					}
				}
			}
			total_pdf = light_pdf + surface_pdf;
			//light_pdf = surface_pdf = 0.5; total_pdf = 1;
			//light_pdf = surface_pdf = 1; total_pdf = 2;
			res = light_contribution * light_pdf / total_pdf + surface_contribution * surface_pdf / total_pdf;
			if (res.anythingWrong())
			{
				return 0;
			}
			return res;
		}



		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			Ray ray = pray;
			//bool use_emissive = true;
			double cost = ray.direction() * scene.m_camera.m_front;
			RGBColor prod_color = scene.m_camera.We(ray.direction()) * (1 / (cost * cost * cost));
			double prod_pdf = scene.m_camera.pdfWeSolidAngle(pray.direction());
			RGBColor res = 0;
			bool first_bounce = true;
			for (int depth = 0; depth < m_max_depth + 1; ++depth)
			{
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					const Material& material = *hit.geometry->getMaterial();
					RGBColor light_contribution;
					double pdf_light = 0;
					if (first_bounce && material.is_emissive())
					{
						//pdf_light = material.pdfLight(hit, ray.direction().normalized());
						light_contribution = prod_color * material.Le(hit.facing, hit.tex_uv) / prod_pdf; //emissif
					}
					first_bounce = false;

					//double pdf_surface = material.pdf(hit, ray.direction().normalized());
					RGBColor direct_contribution = prod_color * MISAddDirectIllumination(scene, hit, sampler) / prod_pdf;

					res += light_contribution + direct_contribution;

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
}