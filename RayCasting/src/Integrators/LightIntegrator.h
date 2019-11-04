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

		inline bool connectVertexToCamera(Scene const& scene, RGBColor const& beta, Hit const& hit, LightVertexStack& res, bool isLight = false)const
		{
			Math::Vector3f vertex_to_camera = scene.m_camera.getPosition() - hit.point;
			const double dist2 = vertex_to_camera.norm2();
			const double dist = std::sqrt(dist2);
			//vertex_to_camera /= dist; //normalize
			vertex_to_camera = vertex_to_camera.normalized();

			LightVertex vertex;
			vertex.uv = scene.m_camera.raster(-vertex_to_camera);

			double We = scene.m_camera.We(-vertex_to_camera);
			if (We == 0)
			{
				return false;
			}

			double G = std::abs(vertex_to_camera * hit.primitive_normal) / dist2; // Geometric term

			if (isLight)
			{
				vertex.light = beta;
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
			SurfaceLightSample sls;
			sampleLight(scene, sls, sampler);

			RGBColor beta = sls.geo->getMaterial()->Le(true, sls.uv) / sls.pdf;

			Hit hit;
			hit.point = sls.vector;
			hit.normal = hit.primitive_normal = sls.normal;
			hit.tex_uv = sls.uv;
			hit.geometry = sls.geo;

			connectVertexToCamera(scene, beta, hit, lvs, true);

			DirectionSample dir_from_light = sls.geo->getMaterial()->sampleLightDirection(sls, sampler); // Uses dedicated function

			double cosl = std::abs(sls.normal * dir_from_light.direction);

			beta = beta * cosl / dir_from_light.pdf;

			unsigned int depth = 0;
			Math::Vector3f next_dir(dir_from_light.direction);

			while (!beta.isBlack())
			{
				++depth;
				Ray ray(hit.point, next_dir);
				if (scene.full_intersection(ray, hit))
				{
					connectVertexToCamera(scene, beta, hit, lvs);

					//sample next dir

					if (depth <= m_max_depth)
					{
						DirectionSample next_dir_sample;
						hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, next_dir_sample, sampler);

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