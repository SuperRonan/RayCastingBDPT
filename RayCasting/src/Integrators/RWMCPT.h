#pragma once

#include <Integrators/Integrator.h>
#include <Image/ImWrite.h>

namespace Integrator
{
	// Reweighted MC PT
	class RWMCPT : public Integrator
	{
	protected:

		Image::Image<RGBColor> m_frame_buffer;

		struct Sample
		{
			RGBColor f;
			double pdf;

			Sample():
				f(0),
				pdf(1)
			{}

			Sample(RGBColor const& f, double pdf):
				f(f),
				pdf(pdf)
			{}

			double w()const
			{
				return pdf;
			}

			RGBColor estimate()const
			{
				return f / pdf;
			}

			bool isZero()const
			{
				return f.isBlack();
			}
		};

		RGBColor addOneDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler, double & pdf)const
		{
			Geometry::SurfaceSample sample;
			scene.sampleLi(sampler, sample, hit);
			pdf = sample.pdf;
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
				return contrib;
			}
			return 0;
		}

		void sendRay(Scene const& scene, Ray ray, Math::Sampler& sampler, Sample * samples, double * sums)const
		{
			bool use_emissive = true;
			double cost = ray.direction() * scene.m_camera.m_front;
			RGBColor f = scene.m_camera.We<true>(ray.direction());
			double pdf = scene.m_camera.pdfWeSolidAngle<true>(ray.direction());
			RGBColor res = 0;
			for (int len = 2; len <= m_max_len; ++len)
			{
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					
					const Material& material = *hit.geometry->getMaterial();
					if (use_emissive)
					{
						const double conversion = len == 2 ? 1.0 : (std::abs(ray.direction() * hit.primitive_normal) / (hit.z * hit.z));
						const int iid = len - 2;
						Sample& sample = samples[iid];
						RGBColor c = f * material.Le(hit.primitive_normal, hit.tex_uv, hit.to_view);
						sample = Sample(c, pdf * conversion);
						sums[iid] += sample.w();
					}

					//use_emissive = hit.geometry->getMaterial()->spicky();
					use_emissive = hit.geometry->getMaterial()->delta();
					if (!use_emissive && len < m_max_len)
					{
						const int iid = len - 1;
						Sample& sample = samples[iid];
						double pdf_l;
						RGBColor c = f * addOneDirectIllumination(scene, hit, sampler, pdf_l);
						sample = Sample(c, pdf * pdf_l);
						sums[iid] += sample.w();
					}

					DirectionSample next_dir;
					material.sampleBSDF(hit, next_dir, sampler);
					f *= next_dir.bsdf * std::abs(next_dir.direction * hit.primitive_normal);
					pdf *= next_dir.pdf;
					ray = Ray(hit.point, next_dir.direction);
					
					if (f.isBlack() || pdf == 0)
					{
						break;
					}
				}
				else
				{
					const int iid = len - 2;
					//RGBColor estimate = beta * scene.getBackgroundColor(ray.direction());
					Sample& sample = samples[iid];
					sample = Sample(0, pdf);
					sums[iid] += sample.w();
					break;
				}
			}
		}

	public:

		RWMCPT(unsigned int spp, unsigned int width, unsigned int height) :
			m_frame_buffer(width, height)
		{
			m_sample_per_pixel = spp;
		}

		void render(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			const int Nintegrals = m_max_len - 1;
			std::vector<std::vector<Sample>> samples_buffers = Parallel::preAllocate(std::vector<Sample>(m_sample_per_pixel * Nintegrals));
			std::vector<std::vector<double>> sums_buffers = Parallel::preAllocate(std::vector<double>(Nintegrals));
			m_frame_buffer.resize(visu.width(), visu.height());
			m_frame_buffer.fill(Image::MultiSample<RGBColor>());
			ProgressReporter reporter;
			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;
			reporter.start(visu.height());
			OMP_PARALLEL_FOR
				for (long y = 0; y < m_frame_buffer.height(); y++)
				{
					int tid = omp_get_thread_num();
					std::vector<Sample>& samples_buffer = samples_buffers[tid];
					std::vector<double>& sums_buffer = sums_buffers[tid];
					for (size_t x = 0; x < visu.width(); x++)
					{
						// Start the pixel
						size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), 0);
						Math::Sampler sampler(seed);
						// Reset the buffers
						std::fill(samples_buffer.begin(), samples_buffer.end(), Sample());
						std::fill(sums_buffer.begin(), sums_buffer.end(), 0.0);
						for (int n = 0; n < m_sample_per_pixel; ++n)
						{
							double xp = sampler.generateContinuous<double>();
							double yp = sampler.generateContinuous<double>();

							double v = ((double)y + yp) / m_frame_buffer.height();
							double u = ((double)x + xp) / visu.width();

							Ray ray = scene.m_camera.getRay(u, v);

							sendRay(scene, ray, sampler, samples_buffer.data() + n * Nintegrals, sums_buffer.data());
						}
						RGBColor result = 0;
						for (int n = 0; n < m_sample_per_pixel; ++n)
						{
							for (int iid = 0; iid < Nintegrals; ++iid)
							{
								if (sums_buffer[iid] > 0)
								{
									const Sample& sample = samples_buffer[n * Nintegrals + iid];
									//double weight = 1.0 / (double)m_sample_per_pixel;
									double weight = sample.w() / sums_buffer[iid];
									result += sample.estimate() * weight;
								}
							}
						}

						visu.plot(x, y, result);
						m_frame_buffer(x, y) = result;
					}//pixel x
					if (tid == 0)
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
			m_frame_buffer.resize(width, height);
			m_frame_buffer.fill(0);

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

		virtual void debug(Scene const&, Visualizer::Visualizer&) override
		{

		}

	};
}