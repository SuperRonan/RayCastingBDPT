#pragma once

#include <Integrators/BidirectionalBase.h>

namespace Integrator
{
	class LightIntegratorBase : public BidirectionalBase
	{
	public:

		LightIntegratorBase(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			BidirectionalBase(sample_per_pixel, width, height)
		{
			
		}



		virtual void traceLight(Scene const& scene, LightVertexStack& lvs, Math::Sampler& sampler)const = 0;


		static __forceinline size_t seed(size_t sample, size_t pass, size_t sample_per_pass)
		{
#ifdef SAMPLER_BIAS
			return pass * sample_per_pass;
#else
			return pass * sample_per_pass + sample;
#endif
		}

		void render(Scene const& scene, Visualizer::Visualizer& visu)
		{
			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;

			LARGE_INTEGER frequency;        // ticks per second
			LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			QueryPerformanceFrequency(&frequency);
			// start timer
			QueryPerformanceCounter(&t1);

			m_frame_buffer.resize(visu.width(), visu.height());
			m_frame_buffer.fill();
			visu.clean();
			const double pixel_area = scene.m_camera.m_down.norm() * scene.m_camera.m_right.norm() / (m_frame_buffer.size());
			const size_t npixels = m_frame_buffer.size();
			const size_t sample_pass = npixels;
			size_t total = 0;

			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				::std::cout << "Pass: " << pass << "/" << m_sample_per_pixel << ::std::endl;

				OMP_PARALLEL_FOR
					for (long sample = 0; sample < sample_pass; ++sample)
					{
						LightVertexStack lvs;
						const size_t seed = LightIntegratorBase::seed(sample, pass, sample_pass);
						Math::Sampler sampler(seed);

						traceLight(scene, lvs, sampler);

						for (LightVertex const& lv : lvs)
						{
							unsigned int u = (lv.uv[0]) * m_frame_buffer.width(), v = (lv.uv[1]) * m_frame_buffer.height();
							RGBColor& pixel = m_frame_buffer(u, v);
							RGBColor Ct = (lv.light);
							//TODO manage omp atomic
							pixel += (Ct);
						}
					}
				total += sample_pass;

				showFrame(visu, total);

				scene.update_lights_offset(1);

				kbr = visu.update();

				QueryPerformanceCounter(&t2);
				elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
				double remainingTime = (elapsedTime / pass) * (m_sample_per_pixel - pass);
				::std::cout << "time: " << elapsedTime << "s. " << ", remaining time: " << remainingTime << "s. " << ", total time: " << elapsedTime + remainingTime << ::std::endl;

				if (kbr == Visualizer::Visualizer::KeyboardRequest::done)
				{
					goto __render__end__loop__;
				}

				if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
				{
					Image::ImWrite::write(m_frame_buffer, 1.0 / (double)total);
				}
			}

		__render__end__loop__:
			// stop timer
			QueryPerformanceCounter(&t2);
			elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
			::std::cout<<"time: "<<elapsedTime<<"s. "<<::std::endl ;

			scene.reset_surface_lights();

			while (kbr != Visualizer::Visualizer::KeyboardRequest::done)
			{
				kbr = visu.update();

				if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
				{
					Image::ImWrite::write(m_frame_buffer, 1.0 / (double)total);
				}
			}
		}







		virtual void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult& res)final override
		{

			LARGE_INTEGER frequency;        // ticks per second
			LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			QueryPerformanceFrequency(&frequency);
			// start timer
			QueryPerformanceCounter(&t1);

			m_frame_buffer.resize(width, height);
			m_frame_buffer.fill();

			const double pixel_area = scene.m_camera.m_down.norm() * scene.m_camera.m_right.norm() / (m_frame_buffer.size());
			const size_t npixels = m_frame_buffer.size();
			const size_t sample_pass = npixels;
			size_t total = 0;

			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				std::cout << '\r' + progession_bar(pass, m_sample_per_pixel, 100) << std::flush;
				OMP_PARALLEL_FOR
					for (long sample = 0; sample < sample_pass; ++sample)
					{
						LightVertexStack lvs;
						const size_t seed = LightIntegratorBase::seed(sample, pass, sample_pass);
						Math::Sampler sampler(seed);

						traceLight(scene, lvs, sampler);

						for (LightVertex const& lv : lvs)
						{
							unsigned int u = (lv.uv[0]) * m_frame_buffer.width(), v = (lv.uv[1]) * m_frame_buffer.height();
							RGBColor& pixel = m_frame_buffer(u, v);
							RGBColor Ct = (lv.light);
							//TODO manage omp atomic
							pixel += Ct;
						}
					}
				total += sample_pass;

				scene.update_lights_offset(1);

			}
			std::cout << '\r' + progession_bar(100, 100, 100) << std::endl;
			// stop timer
			QueryPerformanceCounter(&t2);
			elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;

			scene.reset_surface_lights();

			//fill the result
			{
				res.time = elapsedTime;
				res.image.resize(width, height);
				OMP_PARALLEL_FOR
					for (long i = 0; i < m_frame_buffer.size(); ++i)
					{
						res.image.m_data[i] = m_frame_buffer.m_data[i] / total;
					}
			}

		}








		void fastRender(Scene const& scene, Visualizer::Visualizer& visu)
		{
			m_frame_buffer.resize(visu.width(), visu.height());
			m_frame_buffer.fill();
			visu.clean();
			const double pixel_area = scene.m_camera.m_down.norm() * scene.m_camera.m_right.norm() / (m_frame_buffer.size());
			const size_t npixels = m_frame_buffer.size();
			const size_t sample_pass = npixels / 2;

			OMP_PARALLEL_FOR
				for (long sample = 0; sample < sample_pass; ++sample)
				{
					LightVertexStack lvs;
					size_t seed = LightIntegratorBase::seed(sample, 0, sample_pass);
#ifdef TIME_SEED
					seed += nano();
#endif
					Math::Sampler sampler(seed);

					traceLight(scene, lvs, sampler);

					for (LightVertex const& lv : lvs)
					{
						unsigned int u = (lv.uv[0]) * m_frame_buffer.width(), v = (lv.uv[1]) * m_frame_buffer.height();
						RGBColor& pixel = m_frame_buffer(u, v);
						RGBColor Ct = (lv.light);
						
						//TODO manage omp atomic
						pixel += (Ct);
					}
				}

			showFrame(visu, sample_pass);

			visu.update();
		}

		void debug(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			//TODO
		}

	};
}