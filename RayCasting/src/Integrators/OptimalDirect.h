#pragma once

#include "Integrator.h"
#include <omp.h>
#include <Geometry/RGBColor.h>
#include <Image/ImWrite.h>
#include <MIS/Estimators.h>

namespace Integrator
{

	class OptimalDirect : public Integrator
	{
	protected:

		using Estimator = MIS::DirectEstimator<RGBColor>;
		//using Estimator = MIS::BalanceEstimator<RGBColor>;

		Image::Image<RGBColor> m_frame_buffer;

		void resizeFrameBuffer(size_t w, size_t h)
		{
			m_frame_buffer.resize(w, h);
		}

	public:

		OptimalDirect(unsigned int sample_per_pixel, unsigned int width, unsigned int height):
			m_frame_buffer(width, height)
		{
			m_sample_per_pixel = sample_per_pixel;

		}


		// returns the balance estimate
		RGBColor MISAddDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler, double * weights)const
		{
			RGBColor res = 0;
			weights[0] = 0;// by default
			weights[1] = 1;
			//sample the surface
			SurfaceSample light_sample;
			//scene.sampleLe(sampler, light_sample);
			scene.sampleLi(sampler, light_sample, hit);
			Math::Vector3f to_light = (light_sample.vector - hit.point);
			const double dist2 = to_light.norm2();
			const double dist = sqrt(dist2);
			to_light /= dist;
			RGBColor bsdf = hit.geometry->getMaterial()->BSDF(hit, to_light);
			
			if (!bsdf.isBlack())
			{
				Ray ray(hit.point, to_light);
				Hit light_hit;
				if (scene.full_intersection(ray, light_hit) && samePoint(light_hit, dist))
				{
					RGBColor contribution = bsdf * std::abs(to_light * hit.primitive_normal) * light_hit.geometry->getMaterial()->Le(light_hit.primitive_normal, light_hit.tex_uv, light_hit.to_view);
					const double surface_pdf = light_sample.pdf * dist2 / std::abs(light_hit.primitive_normal * to_light);
					const double bsdf_pdf = hit.geometry->getMaterial()->pdf(hit, to_light);
					const double sum = bsdf_pdf + surface_pdf;
					weights[0] = bsdf_pdf / sum;
					weights[1] = surface_pdf / sum;
					res = contribution / (sum);
				}
			}
			return res;
		}

		RGBColor MISAddRISDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler)const
		{
			RGBColor res = 0;

			//sample the surface
			SurfaceSample light_sample;
			RGBColor contribution;
			scene.sampleLiRIS(sampler, light_sample, hit, &contribution);
			if (contribution.isBlack())
				return contribution;
			Math::Vector3f to_light = (light_sample.vector - hit.point);
			const double dist2 = to_light.norm2();
			const double dist = sqrt(dist2);
			to_light /= dist;
			Ray ray(hit.point, to_light);
			Hit light_hit;
			if (scene.full_intersection(ray, light_hit) && samePoint(light_hit, dist))
			{
				const double surface_pdf = light_sample.pdf * dist2 / std::abs(light_hit.primitive_normal * to_light);
				const double bsdf_pdf = hit.geometry->getMaterial()->pdf(hit, to_light);
				double weight = surface_pdf / (surface_pdf + bsdf_pdf);
				res = contribution / light_sample.pdf * weight;
			}
			return res;
		}


		// Returns the emissive lighting (paths with 2 vertices)
		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler, std::vector<Estimator> & estimators)const 
		{
			Ray ray(pray);
			RGBColor beta = scene.m_camera.We(ray.direction()) / scene.m_camera.pdfWeSolidAngle(ray.direction());
			RGBColor emissive = 0;
			int len = 1;
			Hit hit, prev_hit;
			bool prev_delta = true;
			double dir_pdf;
			while (len < m_max_len)
			{
				prev_hit = hit;
				if (scene.full_intersection(ray, hit))
				{
					++len;
					if (len == 2)
					{
						emissive = beta * hit.geometry->getMaterial()->Le(hit.primitive_normal, hit.tex_uv, hit.to_view);
					}
					else if (prev_delta)
					{
						RGBColor contribution = hit.geometry->getMaterial()->Le(hit.primitive_normal, hit.tex_uv, hit.to_view);
					}
					else
					{
						
						RGBColor contribution = beta * hit.geometry->getMaterial()->Le(hit.primitive_normal, hit.tex_uv, hit.to_view);
						if (!contribution.isBlack())
						{
							Estimator& estimator = estimators[len - 2 - 1];
							double surface_pdf_area;
							/*if constexpr (USE_RIS)
								surface_pdf_area = scene.pdfRISEstimate(prev_hit, hit, sampler, contribution);
							else*/
							surface_pdf_area = scene.pdfSampleLi(hit.geometry, prev_hit, hit.point);
							const double conversion = hit.z * hit.z / (std::abs(hit.primitive_normal * ray.direction()));
							const double surface_pdf = surface_pdf_area * conversion;
							const double sum = dir_pdf + surface_pdf;
							const double weights[2] = { dir_pdf / sum, surface_pdf / sum };
							RGBColor balance_estimate = contribution * weights[0];
							estimator.addEstimate(balance_estimate, weights, 0);
						}
					}

					prev_delta = hit.geometry->getMaterial()->delta();

					if (!prev_delta && len < m_max_len)
					{
						/*if constexpr (USE_RIS)
							res += beta * MISAddRISDirectIllumination(scene, hit, sampler);
						else*/
						double weights[2];
						RGBColor balance_estimate = beta * MISAddDirectIllumination(scene, hit, sampler, weights);
						Estimator& estimator = estimators[len - 2];
						estimator.addEstimate(balance_estimate, weights, 1);
					}

					double xi = sampler.generateContinuous<double>();

					DirectionSample next_dir;
					hit.geometry->getMaterial()->sampleBSDF(hit, next_dir, sampler);
					beta *= next_dir.bsdf * std::abs(hit.primitive_normal * next_dir.direction) / next_dir.pdf;
					ray = Ray(hit.point, next_dir.direction);
					dir_pdf = next_dir.pdf;
					
					if (beta.isBlack())
						break;
				}
				else
				{
					emissive = beta * scene.getBackgroundColor(ray.direction());
					break;
				}
			}
			return emissive;
		}

		

		void render(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			std::vector<std::vector<Estimator>> estimators_buffer;
			// There are two techniques, 0 -> bsdf sampling, 1 -> light sampling
			{
				std::vector<Estimator> estimators(m_max_len - 2, Estimator(2));
				estimators_buffer = Parallel::preAllocate(estimators);
			}
			resizeFrameBuffer(visu.width(), visu.height());
			m_frame_buffer.fill(Image::MultiSample<RGBColor>());
			ProgressReporter reporter;
			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;
			reporter.start(visu.height());
			OMP_PARALLEL_FOR
			for (long y = 0; y < m_frame_buffer.height(); y++)
			{
				int tid = omp_get_thread_num();
				std::vector<Estimator>& estimators = estimators_buffer[tid];
				for (size_t x = 0; x < visu.width(); x++)
				{
					// Start the pixel
					for (Estimator& estimator : estimators)		estimator.reset();
					size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), 0);
					Math::Sampler sampler(seed);
					RGBColor emissive = 0;
					for (int n = 0; n < m_sample_per_pixel; ++n)
					{
						double xp = sampler.generateContinuous<double>();
						double yp = sampler.generateContinuous<double>();

						double v = ((double)y + yp) / m_frame_buffer.height();
						double u = ((double)x + xp) / visu.width();

						Ray ray = scene.m_camera.getRay(u, v);

						RGBColor result = sendRay(scene, ray, sampler, estimators);
						emissive += result;
					}

					RGBColor estimate = 0;
					for (Estimator estimator : estimators)
					{
						estimate += estimator.solve(m_sample_per_pixel);
					}
					RGBColor result = emissive / m_sample_per_pixel + estimate;
					visu.plot(x, y, result);
					m_frame_buffer(x, y) = result;
				}//pixel x
				if(tid == 0)
					reporter.report(y + 1, -1);
			}//pixel y
			reporter.finish();
			kbr = visu.update();
			visu.show();
			if (kbr == Visualizer::Visualizer::KeyboardRequest::done)
			{
				goto __render__end__loop__;
			}
			else if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
			{
				Image::ImWrite::write(m_frame_buffer);
			}
		__render__end__loop__:
			
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


							RGBColor result = 0;


						}//pixel x
					}//pixel y
					//the pass has been computed
				++pass;


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
						res.image.m_data[i] = m_frame_buffer.m_data[i];
					}
			}

		}



		void fastRender(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
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

						RGBColor result = 0;// sendRay(scene, ray, sampler);

						visu.plot(x, y, result);
					}
				}
			visu.update();
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

					RGBColor sample = 0;// sendRay(scene, ray, sampler);

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