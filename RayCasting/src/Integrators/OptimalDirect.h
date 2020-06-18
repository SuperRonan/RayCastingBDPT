#pragma once

#include "Integrator.h"
#include <omp.h>
#include <Geometry/RGBColor.h>
#include <Image/ImWrite.h>
#include <MIS/Estimators.h>

namespace Integrator
{
	template <int MODE>
	class OptimalDirect : public Integrator
	{
	protected:

		// MODE: 
		// 0 -> bsdf + Li
		// 1 -> bsdf + RIS
		// 2 -> bsdf + Li + RIS

		static_assert(MODE >= 0 && MODE <= 2);

#define USE_RIS (MODE >= 1)
#define USE_Li (MODE == 0 || MODE == 2)

#define numTechs ((MODE == 2) ? 3 : 2)

#define Li_ID 1
#define RIS_ID (MODE == 1 ? 1 : 2)
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
			weights[0] = 0;// By default
			weights[Li_ID] = 1;
			if (USE_RIS)
				weights[RIS_ID] = 0;
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
					const double conversion = dist2 / std::abs(light_hit.primitive_normal * to_light);
					const double cos_theta = std::abs(to_light * hit.primitive_normal);
					const RGBColor Le = light_hit.geometry->getMaterial()->Le(light_hit.primitive_normal, light_hit.tex_uv, light_hit.to_view);
					// in Solid angle
					const RGBColor contribution = bsdf * cos_theta * Le;
					
					const double surface_pdf = light_sample.pdf * conversion;
					const double bsdf_pdf = hit.geometry->getMaterial()->pdf(hit, to_light);
					double sum = bsdf_pdf + surface_pdf;
					if constexpr(USE_RIS)
					{
						const double ris_pdf = scene.pdfRISEstimate(hit, light_hit, sampler, contribution / conversion) * conversion;
						sum += ris_pdf;
						weights[RIS_ID] = ris_pdf / sum;
					}
					weights[0] = bsdf_pdf / sum;
					weights[Li_ID] = surface_pdf / sum;
					res = contribution / surface_pdf * weights[Li_ID];
				}
			}
			return res;
		}

		RGBColor MISAddRISDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler, double * weights)const
		{
			RGBColor res = 0;
			weights[0];
			weights[RIS_ID] = 1;
			if constexpr (USE_Li)
				weights[Li_ID] = 0; 
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
				const double conversion = dist2 / std::abs(light_hit.primitive_normal * to_light);
				const double ris_pdf = light_sample.pdf * conversion;
				const double bsdf_pdf = hit.geometry->getMaterial()->pdf(hit, to_light);
				double sum = bsdf_pdf + ris_pdf;
				if constexpr (USE_Li)
				{
					const double Li_pdf = scene.pdfSampleLi(light_hit.geometry, hit, light_hit.point) * conversion;
					sum += Li_pdf;
					weights[Li_ID] = Li_pdf / sum;
				}
				weights[0] = bsdf_pdf / sum;
				weights[RIS_ID] = ris_pdf / sum;
				res = contribution / light_sample.pdf * weights[RIS_ID];
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
			RGBColor prev_bsdf;
			double prev_cos_theta;
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
						RGBColor Le = hit.geometry->getMaterial()->Le(hit.primitive_normal, hit.tex_uv, hit.to_view);
						RGBColor estimate = beta * Le;
						Estimator& estimator = estimators[len - 2 - 1];
						estimator.addOneTechniqueEstimate(estimate, 0);
					}
					else
					{
						RGBColor Le = hit.geometry->getMaterial()->Le(hit.primitive_normal, hit.tex_uv, hit.to_view);
						if (!Le.isBlack())
						{
							Estimator& estimator = estimators[len - 2 - 1];
							const double conversion = hit.z * hit.z / (std::abs(hit.primitive_normal * ray.direction()));
							double sum = dir_pdf;
							double weights[numTechs];
							double* const& pdf = weights;
							pdf[0] = dir_pdf;
							if constexpr (USE_RIS)
							{
								double G = prev_cos_theta / conversion;
								RGBColor contribution = prev_bsdf * Le * G;
								const double surface_pdf_ris_area = scene.pdfRISEstimate(prev_hit, hit, sampler, contribution);
								const double pdf_ris = surface_pdf_ris_area * conversion;
								sum += pdf_ris;
								pdf[RIS_ID] = pdf_ris;
							}
							if constexpr (USE_Li)
							{
								const double surface_pdf_area = scene.pdfSampleLi(hit.geometry, prev_hit, hit.point);
								const double surface_pdf = surface_pdf_area * conversion;
								pdf[Li_ID] = surface_pdf;
								sum += surface_pdf;
							}
							for (int i = 0; i < numTechs; ++i)	weights[i] = pdf[i] / sum;
							
							RGBColor balance_estimate = beta * Le * weights[0];
							estimator.addEstimate(balance_estimate, weights, 0);
						}
					}

					prev_delta = hit.geometry->getMaterial()->delta();

					if (!prev_delta && len < m_max_len)
					{
						double weights[numTechs];
						RGBColor balance_estimate;
						Estimator& estimator = estimators[len - 2];
						
						if constexpr (USE_Li)
						{
							balance_estimate = MISAddDirectIllumination(scene, hit, sampler, weights);
							estimator.addEstimate(beta * balance_estimate, weights, Li_ID);
						}
						if constexpr(USE_RIS)
						{
							balance_estimate = MISAddRISDirectIllumination(scene, hit, sampler, weights);
							estimator.addEstimate(beta * balance_estimate, weights, RIS_ID);
						}
					}

					DirectionSample next_dir;
					hit.geometry->getMaterial()->sampleBSDF(hit, next_dir, sampler);
					prev_cos_theta = std::abs(hit.primitive_normal * next_dir.direction);
					beta *= next_dir.bsdf * prev_cos_theta / next_dir.pdf;
					ray = Ray(hit.point, next_dir.direction);
					dir_pdf = next_dir.pdf;
					prev_bsdf = next_dir.bsdf;
					if (beta.isBlack())
						break;
				}
				else
				{
					break;
				}
			}
			return emissive;
		}

		

		void render(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			std::vector<std::vector<Estimator>> estimators_buffer;
			// There are two / three techniques, 0 -> bsdf sampling, 1 -> light sampling, 2 -> RISLi
			{
				std::vector<Estimator> estimators(m_max_len - 2, Estimator(numTechs));
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
#undef numTechs
#undef USE_RIS
#undef USE_Li
#undef Li_ID
#undef RIS_ID
	};
}