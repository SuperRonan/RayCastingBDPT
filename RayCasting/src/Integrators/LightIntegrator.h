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
				vertex.light = beta * hit.geometry->getMaterial()->Le(vertex_to_camera, hit.tex_uv, hit.primitive_normal);
			}
			else
			{
				vertex.light = beta * hit.geometry->getMaterial()->BSDF(hit, vertex_to_camera, true);
			}

			vertex.light = vertex.light * G * We;

			if (vertex.light.isBlack() || !scene.cameraVisibility(hit.point))
			{
				return false;
			}

			res.push(vertex);
			return true;
		}

		void traceLight(Scene const& scene, LightVertexStack& lvs, Math::Sampler& sampler)const final override
		{
			// The fist sample can be drawn with the information of the camera
			{
				SurfaceSample light_point;
				Hit camera_hit;
				camera_hit.point = scene.m_camera.m_position;
				scene.sampleLi(sampler, light_point, camera_hit);
				//sampleOneLight(scene, sampler, light_point);
				Hit hit;
				hit.point = light_point.vector;
				hit.normal = hit.primitive_normal = light_point.normal;
				hit.tex_uv = light_point.uv;
				hit.geometry = light_point.geo;

				connectVertexToCamera<true>(scene, 1.0 / light_point.pdf, hit, lvs);
			}

			if (m_max_len == 2)
			{
				return;
			}

			//first, sample a point a light
			SurfaceSample light_point;
			scene.sampleLe(sampler, light_point);


			RGBColor beta = 1.0 / light_point.pdf;

			Hit hit;
			hit.point = light_point.vector;
			hit.normal = hit.primitive_normal = light_point.normal;
			hit.tex_uv = light_point.uv;
			hit.geometry = light_point.geo;

			//connectVertexToCamera<true>(scene, 1.0 / light_point.pdf, hit, lvs);

			DirectionSample ds = hit.geometry->getMaterial()->sampleLightDirection(light_point, sampler);
						
			Math::Vector3f next_dir = ds.direction;
			beta *= ds.bsdf / ds.pdf * std::abs(hit.primitive_normal * ds.direction);
			

			unsigned int len = 2;
			while(!beta.isBlack())
			{
				Ray ray(hit.point, next_dir);
				if (scene.full_intersection(ray, hit))
				{
					++len;
					connectVertexToCamera<false>(scene, beta, hit, lvs);
					
					//sample next dir

					if (len < m_max_len)
					{
						DirectionSample next_dir_sample;
						hit.geometry->getMaterial()->sampleBSDF(hit, next_dir_sample, sampler, true);

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