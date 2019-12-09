#pragma once

#include <Integrators/LightIntegratorBase.h>
#include <cassert>

namespace Integrator
{
	class LightIntegrator final : public LightIntegratorBase
	{
	protected:

	public:

		LightIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			LightIntegratorBase(sample_per_pixel, width, height)
		{

		}

		template <bool LIGHT_VERTEX>
		inline bool connectVertexToCamera(Scene const& scene, RGBColor const& beta, Hit const& hit, LightVertexStack & res)const
		{
			Math::Vector3f vertex_to_camera = scene.m_camera.getPosition() - hit.point;
			const double dist2 = vertex_to_camera.norm2();
			const double dist = std::sqrt(dist2);
			vertex_to_camera /= dist;

			LightVertex vertex;
			vertex.uv = scene.m_camera.raster(-vertex_to_camera);

			double We = scene.m_camera.We(-vertex_to_camera);
			if (We == 0)
			{
				return false;
			}

			double G = std::abs(vertex_to_camera * hit.primitive_normal) / dist2;

			if constexpr (LIGHT_VERTEX)
			{
				vertex.light = beta * hit.geometry->getMaterial()->Le(vertex_to_camera * hit.primitive_normal > 0, hit.tex_uv);
			}
			else
			{
				vertex.light = beta * hit.geometry->getMaterial()->BSDF(hit, vertex_to_camera);
			}

			vertex.light = vertex.light * G * We;

			if (vertex.light.isBlack() || !cameraVisibility(scene, hit.point))
			{
				return false;
			}

			res.push(vertex);
			return true;
		}

		void traceLight(Scene const& scene, LightVertexStack& lvs, Math::Sampler& sampler)const final override
		{
			//first, sample a point a light
			SurfaceLightSample light_point;
			sampleLight(scene, light_point, sampler);


			RGBColor beta = light_point.geo->getMaterial()->Le(true, light_point.uv) / light_point.pdf;

			Hit hit;
			hit.point = light_point.vector;
			hit.normal = hit.primitive_normal = light_point.normal;
			hit.tex_uv = light_point.uv;
			hit.geometry = light_point.geo;

			connectVertexToCamera<true>(scene, 1.0 / light_point.pdf, hit, lvs);

			Math::RandomDirection Le_sampler(&sampler, light_point.normal, 1);
			Math::Vector3f next_dir = Le_sampler.generate();
			double cosl = light_point.normal * next_dir;
			if (cosl < 0)
			{
				cosl = -cosl;
				next_dir = -next_dir;
			}
			double next_dir_pdf = cosl / Math::pi;


			beta = beta * (cosl / next_dir_pdf);

			unsigned int depth = 0;
			while(!beta.isBlack())
			{
				++depth;
				Ray ray(hit.point, next_dir);
				if (scene.full_intersection(ray, hit))
				{
					connectVertexToCamera<false>(scene, beta, hit, lvs);
					
					//sample next dir

					if (depth <= m_max_depth)
					{
						DirectionSample next_dir_sample;
						hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, next_dir_sample, sampler, true);

						beta = beta * next_dir_sample.bsdf * std::abs(next_dir_sample.direction * hit.primitive_normal) / next_dir_sample.pdf;
						next_dir = next_dir_sample.direction;
					}
					else
						break;
				}
				else
					break;
			}
			
		}

	};
}