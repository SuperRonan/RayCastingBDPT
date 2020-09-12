#pragma once

#include <vector>
#include <Integrators/Integrator.h>

namespace Integrator
{
	class StratifiedPT: public Integrator
	{
	protected:

		std::vector<unsigned int> m_strata;

		unsigned int m_u_strata, m_v_strata;

		Image::Image<Image::MultiSample<RGBColor>> m_frame_buffer;
		
		void resizeFrameBuffer(size_t w, size_t h)
		{
			m_frame_buffer.resize(w, h);
		}


		bool computeStrata(Scene const& scene, unsigned int width, unsigned int height)
		{
			unsigned int N = width * height;
			if (N != m_strata.size())
			{
				m_strata.resize(N);
				std::iota(m_strata.begin(), m_strata.end(), 0);
				// use the sampler
				std::shuffle(m_strata.begin(), m_strata.end(), std::mt19937_64());

				m_u_strata = width;
				m_v_strata = height;
				return true;
			}
			return false;
		}

	public:


		StratifiedPT(unsigned int spp, unsigned int width, unsigned int height) :
			m_frame_buffer(width, height)
		{
			m_sample_per_pixel = (spp);
		}

		RGBColor directIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler, unsigned int stratum_index)const
		{
			Geometry::SurfaceSample sample;
			unsigned int u_stratum = stratum_index % m_u_strata;
			unsigned int v_stratum = stratum_index / m_u_strata;
			double x1 = sampler.generateStratified<double>(u_stratum, m_u_strata);
			double x2 = sampler.generateStratified<double>(v_stratum, m_v_strata);
			{// SampleLi
				
				x2 *= scene.m_surface_lights.size();
				int light_index = x2;
				sample.geo = scene.m_surface_lights[light_index];
				x2 = (x2 - light_index) / scene.m_surface_lights.size();
				sample.geo->sampleLight(sample, hit, x1, x2);
				sample.pdf /= scene.m_surface_lights.size();
			}
			//scene.sampleLe(sampler, sample);
			Math::Vector3f to_light = sample.vector - hit.point;
			const double dist2 = to_light.norm2();
			const double dist = std::sqrt(dist2);
			to_light /= dist;
			const RGBColor bsdf = hit.geometry->getMaterial()->BSDF(hit, to_light, hit.to_view);
			const double cos_on_light = -(sample.normal * to_light);
			const double cos_theta = std::abs(to_light * hit.normal);
			const RGBColor Le = sample.geo->getMaterial()->Le(sample.normal, sample.uv, -to_light);
			const RGBColor contrib = bsdf * Le * std::abs(cos_on_light) * cos_theta / dist2;
			if (contrib.isBlack() || contrib.anythingWrong())
				return 0;
			Hit light_hit;
			Ray ray(hit.point, to_light);
			if (scene.full_intersection(ray, light_hit) && samePoint(light_hit, dist, sample.geo))
			{
				return contrib / sample.pdf;
			}
			return 0;
		}





		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler, unsigned int stratum_index)const
		{
			Ray ray = pray;
			unsigned int len = 1;
			RGBColor res = 0;
			RGBColor T = 1;
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
						res += T * directIllumination(scene, hit, sampler, stratum_index);
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






		virtual void render(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			resizeFrameBuffer(visu.width(), visu.height());
			m_frame_buffer.fill(Image::MultiSample<RGBColor>());


			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;

			ProgressReporter reporter;
			reporter.start(m_sample_per_pixel);
			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				OMP_PARALLEL_FOR
					for (long y = 0; y < m_frame_buffer.height(); y++)
					{
						int tid = omp_get_thread_num();

						for (size_t x = 0; x < visu.width(); x++)
						{

							size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), pass);
							Math::Sampler sampler(seed);

							double xp = sampler.generateContinuous<double>();
							double yp = sampler.generateContinuous<double>();

							double v = ((double)y + yp) / m_frame_buffer.height();
							double u = ((double)x + xp) / visu.width();

							Ray ray = scene.m_camera.getRay(u, v);


							RGBColor result = sendRay(scene, ray, sampler, 0);



							Image::MultiSample<RGBColor>& pixel = m_frame_buffer(x, y);
							pixel.add(result);

							visu.plot(x, y, pixel.mean());
						}//pixel x
					}//pixel y
					//the pass has been computed
				reporter.report(pass + 1, -1);

				kbr = visu.update();

				if (kbr == Visualizer::Visualizer::KeyboardRequest::done)
				{
					goto __render__end__loop__;
				}
				else if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
				{
					Image::ImWrite::write(m_frame_buffer);
				}

			}//pass per pixel
		__render__end__loop__:

			reporter.finish();
			scene.reset_surface_lights();

			while (kbr != Visualizer::Visualizer::KeyboardRequest::done)
			{
				kbr = visu.update();

				if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
				{
					Image::ImWrite::write(m_frame_buffer);
				}
			}
			m_frame_buffer.resize(1, 1);
		}


		virtual void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult& res)final override
		{
			resizeFrameBuffer(width, height);
			m_frame_buffer.fill(Image::MultiSample<RGBColor>());
			ProgressReporter reporter;
			reporter.start(m_sample_per_pixel);
			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				OMP_PARALLEL_FOR
					for (long y = 0; y < m_frame_buffer.height(); y++)
					{
						int tid = omp_get_thread_num();

						for (size_t x = 0; x < width; x++)
						{
							size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), pass);
							Math::Sampler sampler(seed);

							double xp = sampler.generateContinuous<double>();
							double yp = sampler.generateContinuous<double>();

							double u = ((double)x + xp) / width;
							double v = ((double)y + yp) / height;

							Ray ray = scene.m_camera.getRay(u, v);


							RGBColor result = sendRay(scene, ray, sampler, 0);

							Image::MultiSample<RGBColor>& pixel = m_frame_buffer(x, y);
							pixel.add(result);

						}//pixel x
					}//pixel y
					//the pass has been computed

				reporter.report(pass + 1, -1);
			}//pass per pixel
			reporter.finish();
			scene.reset_surface_lights();

			//fill the result
			{
				res = Auto::RenderResult();
				res.time = reporter.time();
				res.image.resize(width, height);
				OMP_PARALLEL_FOR
					for (long i = 0; i < m_frame_buffer.size(); ++i)
					{
						res.image.m_data[i] = m_frame_buffer.m_data[i].mean();
					}
			}

		}



		void fastRender(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
#ifdef TIME_SEED
			size_t time_seed = nano();
#endif

			computeStrata(scene, visu.width(), visu.height());

			OMP_PARALLEL_FOR
				for (long y = 0; y < visu.height(); ++y)
				{
					for (size_t x = 0; x < visu.width(); ++x)
					{
						size_t seed = pixelSeed(x, y, visu.width(), visu.height(), 0);
#ifdef TIME_SEED
						seed += time_seed;
#endif
						Math::Sampler sampler(seed);
						Ray ray(scene.m_camera.getRay(((double)x + 0.5) / visu.width(), ((double)y + 0.5) / visu.height()));
						RGBColor result = sendRay(scene, ray, sampler, m_strata[x * visu.width() + y]);

						visu.plot(x, y, result);
					}
				}
			visu.update();
		}


		virtual void debug(Scene const&, Visualizer::Visualizer&) {}

	};
}