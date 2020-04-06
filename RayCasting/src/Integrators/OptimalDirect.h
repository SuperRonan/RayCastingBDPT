#pragma once

#include "Integrator.h"
#include <omp.h>
#include <Geometry/RGBColor.h>
#include <Image/MultiSample.h>
#include <Image/ImWrite.h>
#include <armadillo>


namespace Integrator
{

	class OptimalDirect : public Integrator
	{
	protected:

		Image::Image<Image::MultiSample<RGBColor>> m_frame_buffer;

		unsigned int m_direct_samples = 1;

		void resizeFrameBuffer(size_t w, size_t h)
		{
			m_frame_buffer.resize(w, h);
		}

	public:

		OptimalDirect(unsigned int sample_per_pixel, unsigned int width, unsigned int height):
			m_frame_buffer(width, height),
			m_direct_samples(sample_per_pixel)
		{
			m_sample_per_pixel = 1;
		}


		RGBColor L(Scene const& scene, Ray const& ray, Math::Sampler& sampler, Hit & hit)const
		{
			if (scene.full_intersection(ray, hit))
			{
				if (hit.geometry->getMaterial()->is_emissive())
				{
					return hit.geometry->getMaterial()->Le(hit.facing, hit.tex_uv);
				}

				
			}
			else
			{
				return scene.getBackgroundColor(ray.direction());
			}
		}

		RGBColor directLighting(Scene const& scene, Hit const& hit, Math::Sampler& sampler)const
		{
			double mat[3] = { 0, 0, 0 };
			RGBColor vec[2] = { 0, 0 };

			const int bsdf_samples = std::max(1u, m_direct_samples / 2);
			const int surface_samples = std::max(1u, m_direct_samples / 2);

			const auto updateSystem = [&](double p0, double p1, RGBColor const& f)
			{
				double q0 = bsdf_samples * p0;
				double q1 = surface_samples * p1;
				double sum_pdf = p0 + p1;
				double sum_q = q0 + q1;
				double S = 1.0 / (sum_q);
				double w0 = p0 * S;
				double w1 = p1 * S;

				//RGBColor balance_estimate = estimated_f * (tech == 0 ? w0 : w1);

				mat[0] += w0 * w0;
				mat[1] += w0 * w1;
				mat[2] += w1 * w1;

				for (int k = 0; k < 3; ++k)
				{
					vec[0][k] += f[k] * w0 * S;
					vec[1][k] += f[k] * w1 * S;
				}
			};

			Material const& material = *hit.geometry->getMaterial();

			RGBColor res;

			

			//tech 0: BSDF
			for (int j = 0; j < bsdf_samples; ++j)
			{
				Hit light_hit;
				DirectionSample dir;
				material.sampleBSDF(hit, 1, 1, dir, sampler);

				RGBColor f = L(scene, Ray(hit.point, dir.direction), sampler, light_hit) * dir.bsdf * std::abs(dir.direction * hit.primitive_normal);

				double dist2 = light_hit.z * light_hit.z;
				double cos_light = light_hit.valid() ? std::abs(dir.direction * light_hit.primitive_normal) : 1;

				double convert = (cos_light / dist2);
				double bsdf_pdf = convert * dir.pdf;

				f *= convert;

				double surface_pdf = 0;
				if (light_hit.valid() && light_hit.geometry->getMaterial()->is_emissive())
				{
					surface_pdf = scene.pdfSamplingLight(light_hit.geometry);
				}

				res += f / bsdf_pdf * 0.5 / bsdf_samples;

				updateSystem(bsdf_pdf, surface_pdf, f);
			}

			//tech 1: Sampling the light
			for (int j = 0; j < surface_samples; ++j)
			{
				SurfaceSample sls;
				sampleOneLight(scene, sampler, sls, j);
				Math::Vector3f dir = sls.vector - hit.point;
				const double dist2 = dir.norm2();
				const double dist = std::sqrt(dist2);
				dir /= dist;
				double cos_light = std::abs(sls.normal * dir);
				double convert = dist2 == 0 ? 1 : (cos_light / dist2);

				RGBColor f = sls.geo->getMaterial()->Le(sls.normal * dir < 0, sls.uv) * hit.geometry->getMaterial()->BSDF(hit, dir) * std::abs(dir * hit.primitive_normal) * convert;
				Hit light_hit;
				double bsdf_pdf = hit.geometry->getMaterial()->pdf(hit, dir) * convert;
				if (!(scene.full_intersection(Ray(hit.point, dir), light_hit) && samePoint(light_hit, dist)))
				{
					f = 0;
					bsdf_pdf = 0;
				}

				res += f / sls.pdf * 0.5 / surface_samples;

				updateSystem(bsdf_pdf, sls.pdf, f);
			}
			return res;

			//solve the system
			arma::mat22 armat;
			armat(0, 0) = mat[0];
			armat(0, 1) = mat[1];
			armat(1, 0) = mat[1];
			armat(1, 1) = mat[2];
			for (int i = 0; i < 3; ++i)
				if (std::isnan(mat[i]) || std::isinf(mat[i]))
					__debugbreak();
			const auto mat_inv = arma::pinv(armat);
			
			
			for (int k = 0; k < 3; ++k)
			{
				arma::vec2 contrib_vec;
				contrib_vec[0] = vec[0][k];
				contrib_vec[1] = vec[1][k];
				
				contrib_vec = mat_inv * contrib_vec;

				res[k] = arma::sum(contrib_vec);
			}
			
			return res;
		}


		RGBColor sendRay(Scene const& scene, Ray const& ray, Math::Sampler& sampler)const 
		{
			RGBColor beta = scene.m_camera.We(ray.direction()) / scene.m_camera.pdfWeSolidAngle(ray.direction());
			RGBColor res = 0;

			Hit hit;
			if (scene.full_intersection(ray, hit))
			{
				const Material& material = *hit.geometry->getMaterial();
				res += beta * material.Le(hit.facing, hit.tex_uv);

				res += beta * directLighting(scene, hit, sampler);
			}
			else
			{
				res += beta * scene.getBackgroundColor(ray.direction());
			}

			return res;
		}

		

		void render(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			resizeFrameBuffer(visu.width(), visu.height());
			m_frame_buffer.fill(Image::MultiSample<RGBColor>());

			// 1 - Rendering time
			LARGE_INTEGER frequency;        // ticks per second
			LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			QueryPerformanceFrequency(&frequency);
			// start timer
			QueryPerformanceCounter(&t1);

			double pixelArea = scene.m_camera.m_down.norm() * scene.m_camera.m_right.norm() / (m_frame_buffer.size());

			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;
			const size_t number_of_pixels = m_frame_buffer.size();
			const size_t sample_pass = number_of_pixels;
			size_t pass = 0;
			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_sample_per_pixel; ++passPerPixelCounter)
			{

				::std::cout << "Pass: " << pass << "/" << Integrator::m_sample_per_pixel << ::std::endl;


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


							RGBColor result = sendRay(scene, ray, sampler);


							Image::MultiSample<RGBColor>& pixel = m_frame_buffer(x, y);
							pixel.add(result);

							visu.plot(x, y, pixel.mean());
						}//pixel x
					}//pixel y
					//the pass has been computed
				++pass;

				scene.update_lights_offset(m_direct_samples);
				kbr = visu.update();
				QueryPerformanceCounter(&t2);
				elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
				double remainingTime = (elapsedTime / pass) * (Integrator::m_sample_per_pixel - pass);
				::std::cout << "time: " << elapsedTime << "s. " << ", remaining time: " << remainingTime << "s. " << ", total time: " << elapsedTime + remainingTime << ::std::endl;

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
			// stop timer
			QueryPerformanceCounter(&t2);
			elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
			::std::cout << "time: " << elapsedTime << "s. " << ::std::endl;

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

			// 1 - Rendering time
			LARGE_INTEGER frequency;        // ticks per second
			LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			QueryPerformanceFrequency(&frequency);
			// start timer
			QueryPerformanceCounter(&t1);

			double pixelArea = scene.m_camera.m_down.norm() * scene.m_camera.m_right.norm() / (m_frame_buffer.size());



			const size_t number_of_pixels = m_frame_buffer.size();
			const size_t sample_pass = number_of_pixels;
			size_t pass = 0;
			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_sample_per_pixel; ++passPerPixelCounter)
			{

				std::cout << '\r' + progession_bar(pass, m_sample_per_pixel, 100) << std::flush;


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


							RGBColor result = sendRay(scene, ray, sampler);

							Image::MultiSample<RGBColor>& pixel = m_frame_buffer(x, y);
							pixel.add(result);

						}//pixel x
					}//pixel y
					//the pass has been computed
				++pass;

				scene.update_lights_offset(m_direct_samples);

			}//pass per pixel
			std::cout << '\r' + progession_bar(100, 100, 100) << std::endl;
			// stop timer
			QueryPerformanceCounter(&t2);
			elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;


			scene.reset_surface_lights();

			//fill the result
			{
				res = Auto::RenderResult();
				res.time = elapsedTime;
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
			int tmp = m_direct_samples;
			m_direct_samples = 1;
#ifdef TIME_SEED
			size_t time_seed = nano();
#endif
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

						RGBColor result = sendRay(scene, ray, sampler);

						visu.plot(x, y, result);
					}
				}
			visu.update();
			m_direct_samples = tmp;
		}




		void debug(Scene const& scene, Visualizer::Visualizer& visu)final override
		{
			size_t x, y;

			bool print_zero = false;

			while (true)
			{
				std::cout << "x= ";
				std::cin >> x;
				std::cout << "y= ";
				std::cin >> y;

				std::vector<RGBColor> samples(m_sample_per_pixel);
				size_t width = visu.width(), height = visu.height();

				///OMP_PARALLEL_FOR
				for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_sample_per_pixel; ++passPerPixelCounter)
				{
					size_t pass = passPerPixelCounter;
					size_t seed = pixelSeed(x, y, width, height, pass);
					Math::Sampler sampler(seed);

					double xp = sampler.generateContinuous<double>();
					double yp = sampler.generateContinuous<double>();

					double u = ((double)x + xp) / width;
					double v = ((double)y + yp) / height;

					Ray ray = scene.m_camera.getRay(u, v);

					RGBColor sample = sendRay(scene, ray, sampler);

					samples[pass] = sample;
				}

				RGBColor sum = 0;
				for (size_t i = 0; i < samples.size(); ++i)
				{
					RGBColor const& sample = samples[i];
					if (!print_zero)
					{
						if (sample.isBlack())
						{
							continue;
						}
					}
					std::cout << "sample " << i << ": " << sample << std::endl;
					sum += sample;
				}
				RGBColor mean = sum / double(samples.size());
				std::cout << "mean: " << mean << std::endl;

				visu.plot(x, y, mean);
				visu.update();

			}
		}
	};

}