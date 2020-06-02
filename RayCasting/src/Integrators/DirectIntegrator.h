#pragma once

#include <Integrators/Integrator.h>
#include <Image/MultiSample.h>
#include <settings.h>
#include <utils.h>
#include <Image/ImWrite.h>
#include <System/ProgressReporter.h>

namespace Integrator
{
	class DirectIntegrator: public Integrator
	{
	protected:

		Image::Image<Image::MultiSample<RGBColor>> m_frame_buffer;

		unsigned int m_direct_samples = 1;

	public:

		DirectIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			m_frame_buffer(width, height)
		{
			m_sample_per_pixel = sample_per_pixel;
		}

		void resizeFrameBuffer(size_t w, size_t h)
		{
			m_frame_buffer.resize(w, h);
		}

		virtual RGBColor sendRay(Scene const& scene, Ray const& ray, Math::Sampler& sampler)const = 0;

		

		void render(Scene const& scene, Visualizer::Visualizer& visu) final override
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
						

						RGBColor result = sendRay(scene, ray, sampler);



						Image::MultiSample<RGBColor>& pixel = m_frame_buffer(x, y);
						pixel.add(result);

						visu.plot(x, y, pixel.mean());
					}//pixel x
				}//pixel y
				//the pass has been computed
				reporter.report(pass+1, -1);

				scene.update_lights_offset(m_direct_samples);
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


		virtual void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult & res)final override
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
							

							RGBColor result = sendRay(scene, ray, sampler);

							Image::MultiSample<RGBColor>& pixel = m_frame_buffer(x, y);
							pixel.add(result);

						}//pixel x
					}//pixel y
					//the pass has been computed

				scene.update_lights_offset(m_direct_samples);
				reporter.report(pass+1, -1);
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
					Ray ray(scene.m_camera.getRay(((double)x +0.5) / visu.width(), ((double)y +0.5) / visu.height()));
					RGBColor result = sendRay(scene, ray, sampler);

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