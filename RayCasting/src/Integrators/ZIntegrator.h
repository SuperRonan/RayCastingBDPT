#pragma once

#include <Integrators/Integrator.h>
#include <limits>
#include <Image/ImWrite.h>
namespace Integrator
{
	class ZIntegrator : public Integrator
	{
	protected:

		Image::Image<Image::MultiSample<RGBColor>> m_frame_buffer;

		double m_z_near = std::numeric_limits<double>::max();
		double m_z_far = std::numeric_limits<double>::epsilon();

		static constexpr double Z_BACKGROUND = std::numeric_limits<double>::max();

	public:

		ZIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			m_frame_buffer(width, height)
		{
			m_sample_per_pixel = sample_per_pixel;
		}

		void resizeFrameBuffer(size_t w, size_t h)
		{
			m_frame_buffer.resize(w, h);
		}

		double depth(Scene const& scene, Ray const& ray)const
		{
			Hit hit;

			if (scene.full_intersection(ray, hit))
			{
				return hit.z;
			}
			return Z_BACKGROUND;
		}

		void render(Scene const& scene, Visualizer::Visualizer& visu)
		{
			resizeFrameBuffer(visu.width(), visu.height());
			m_frame_buffer.fill(RGBColor());

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

				for(unsigned int sample = 0; sample < m_sample_per_pixel; ++sample)
				{
					::std::cout << "Pass: " << pass << "/" << Integrator::m_sample_per_pixel << ::std::endl;
					++pass;

					OMP_PARALLEL_FOR
						for (long y = 0; y < m_frame_buffer.height(); y++)
						{
							int tid = omp_get_thread_num();
							
							for (size_t x = 0; x < visu.width(); x++)
							{
								size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), sample);
								Math::Sampler sampler(seed);

								double xp = sampler.generateContinuous<double>();
								double yp = sampler.generateContinuous<double>();

								double u = ((double)x + xp) / m_frame_buffer.width();
								double v = ((double)y + yp) / m_frame_buffer.height();

								Ray ray = scene.m_camera.getRay(u, v);

								double result = depth(scene, ray);

								Image::MultiSample<RGBColor>& pixel = m_frame_buffer(x, y);
								pixel.add(result);

								visu.plot(x, y, pixel.mean());
							}//pixel x
						}//pixel y
						//the pass has been computed
						
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
				}// sample
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
		}

		virtual void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult& res)final override
		{
			//AQUA
		}



		void fastRender(Scene const& scene, Visualizer::Visualizer& visu)
		{
			double next_near = std::numeric_limits<double>::max();
			double next_far = std::numeric_limits<double>::min();

			OMP_PARALLEL_FOR
			for (long y = 0; y < visu.height(); ++y)
			{
				double min_near = std::numeric_limits<double>::max();
				double max_far = std::numeric_limits<double>::min();

				for (size_t x = 0; x < visu.width(); ++x)
				{
					Ray ray(scene.m_camera.getRay(((double)x) / visu.width(), ((double)y) / visu.height()));
					double z = depth(scene, ray);

					if (z != Z_BACKGROUND)
					{
						if (z < min_near)
						{
							min_near = z;
						}
						if (z > max_far)
						{
							max_far = z;
						}
					}
					//z buffer normalisation
					z = (std::max(z - m_z_near, 0.0)) / (m_z_far - m_z_near);

					visu.plot(x, y, z);
				}
				double tmp = std::min(min_near, next_near);

				next_near = tmp;

				tmp = std::max(max_far, next_far);

				next_far = tmp;
			}
			m_z_near = next_near;
			m_z_far = next_far;
			visu.update();
		}

		void debug(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			//AQUA
		}
	};
}