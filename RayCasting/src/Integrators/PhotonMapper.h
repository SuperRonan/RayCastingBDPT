#pragma once

#include <settings.h>
#include <Integrators/Integrator.h>
#include <Math/Vectorf.h>
#include <Integrators/PhotonMap.h>
#include <Image/Image.h>
#include <Image/MultiSample.h>

namespace Integrator
{
	class PhotonMapper: public Integrator
	{
	protected:

		Image::Image<Image::MultiSample<RGBColor>> m_frame;

		void resizeFrameBuffer(size_t w, size_t h)
		{
			m_frame.resize(w, h);
		}

		using Vector3i = Math::Vector<int, 3>;
		using Vector3f = Math::Vector3f;
		using Vector2f = Math::Vector2f;

		using PhotonFloat = float;

		template <class Float>
		class Photon
		{
		public:
			const Primitive* m_primitive;
			Math::Vector<Float, 3> m_point;
			uint8_t m_depth;
			Math::Vector<Float, 3> m_dir;
			Float m_beta[3];

			Photon(Hit const& hit, uint8_t depth, RGBColor const& beta):
				m_primitive(hit.primitve),
				m_point(hit.point),
				m_depth(depth),
				m_dir(hit.to_view)
			{
				m_beta[0] = beta[0];
				m_beta[1] = beta[1];
				m_beta[2] = beta[2];
			}

			void fillHit(Hit& res)const
			{
				res.primitve = m_primitive;
				res.geometry = m_primitive->geometry();
				res.point = m_point;
				res.to_view = m_dir;

				res.primitive_uv = m_primitive->uv(m_point);
				res.tex_uv = m_primitive->tuv(res.primitive_uv);

				res.primitive_normal = m_primitive->normal(m_point, res.primitive_uv);
				res.normal = m_primitive->shading_normal(m_point, res.primitive_uv);

				res.z = -1;
				res.facing = res.normal * res.to_view > 0;
			}

			RGBColor beta()const
			{
				return RGBColor(m_beta[0], m_beta[1], m_beta[2]);
			}

			Math::Vector<Float, 3> point()const
			{
				return m_point;
			}

			Vector3f pnormal()const
			{
				return m_primitive->normal(m_point, m_primitive->tuv(m_primitive->uv(m_point)));
			}
		};

		using Photonf = Photon<PhotonFloat>;

	public:

		static int photon_size()
		{
			return sizeof(Photonf);
		}

	protected:

		PhotonMap<Photonf> m_map;

		double m_relative_radius;
		double m_radius, m_radius2;

		unsigned int m_number_of_photons;

		unsigned int m_spicky_samples;

		double m_alpha, m_beta;

		double kernel(double dist2)const
		{
			return 1.0 / (Math::pi * m_radius2);
			double dist = std::sqrt(dist2);
			return m_alpha - dist * m_beta;
		}

	public:

		PhotonMapper(unsigned int sample_per_pixel, unsigned int width, unsigned int height):
			m_frame(width, height)
		{
			m_sample_per_pixel = sample_per_pixel;
		}

		void setParams(Scene const& scene, double relative_radius, int pcount)
		{
			m_relative_radius = relative_radius;
			BoundingBox m_bb = scene.m_sceneBoundingBox;
			m_bb[0] -= Vector3f(0.001, 0.001, 0.001);
			m_bb[1] += Vector3f(0.001, 0.001, 0.001);
			Vector3f dim = m_bb.diag();
			double max_dir = dim.simdAbs().max();

			m_radius = max_dir * m_relative_radius;
			m_radius2 = m_radius * m_radius;
			Math::Vector3f m_pixel_size = m_radius;
			Vector3f sizef = dim.simdDiv(m_pixel_size);
			Vector3i m_size = sizef.ceil();

			m_map.init(m_bb, m_size);

			m_alpha = 3.0 / (Math::pi * m_radius2);
			m_beta = 3.0 / (Math::pi * m_radius2 * m_radius);

			m_number_of_photons = pcount;

			m_spicky_samples = 1;
		}

		void buildMap(Scene const& scene, int offset = 0)
		{
			//tic();
			OMP_PARALLEL_FOR
				for (int sample = 0; sample < m_number_of_photons; ++sample)
				{
					Math::Sampler sampler(sample+offset);
					SurfaceSample sls;
					sampleOneLight(scene, sampler, sls);
					
					DirectionSample dirSample = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);
					Hit hit;
					RGBColor beta = dirSample.bsdf / (sls.pdf * dirSample.pdf) * std::abs(dirSample.direction * sls.normal);
					Ray ray(sls.vector, dirSample.direction);
					for (int len = 2; len <= m_max_len-1; ++len)
					{
						if (scene.full_intersection(ray, hit))
						{
							if (!hit.geometry->getMaterial()->spicky())
							{
								Photonf photon = { hit, (uint8_t)(len - 2), beta / ((double)m_number_of_photons) };
								m_map.addPhoton(photon);
							}

							hit.geometry->getMaterial()->sampleBSDF(hit, dirSample, sampler, true);
							beta *= dirSample.bsdf / dirSample.pdf * std::abs(dirSample.direction * hit.primitive_normal);
							ray = { hit.point, dirSample.direction };
						}
						else
						{
							break;
						}
					}
				}
			//toc();
			m_map.buildDone();
		}

		RGBColor addPhotons(Hit const& hit, int len)const
		{
			RGBColor photons_contrib;
			m_map.loopThroughPhotons([&](Photonf const& photon) {
				int path_len = photon.m_depth + 2 + len - 1;
				if (path_len <= m_max_len)
				{
					const Vector3f d = hit.point - photon.m_point;
					const double dist2 = d.norm2();
					if (dist2 < m_radius2 && std::abs(hit.primitive_normal * photon.pnormal()) > 0.8)
					{
						const double k = kernel(dist2);
						photons_contrib += photon.beta() * hit.geometry->getMaterial()->BSDF(hit, photon.m_dir, hit.to_view) * k;
					}
				}
				}, hit.point);
			return photons_contrib;
		}

		template <bool FAST>
		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const 
		{
			Ray ray = pray;
			RGBColor beta = scene.m_camera.We<true>(ray.direction()) / scene.m_camera.pdfWeSolidAngle<true>(ray.direction());
			RGBColor res = 0;
			Hit hit;
			if (scene.full_intersection(ray, hit))
			{
				int len = 2;
				const Material& material = *hit.geometry->getMaterial();
				res += beta * material.Le(hit.facing, hit.tex_uv);
				if (material.spicky())
				{
					const int N = FAST ? 1 : m_spicky_samples;
					RGBColor subres = 0;
					for (int subsample = 0; subsample < N; ++subsample)
					{
						Hit subhit = hit;
						RGBColor subbeta = beta;
						DirectionSample next_dir;
						material.sampleBSDF(subhit, next_dir, sampler, false);
						ray = { subhit.point, next_dir.direction };
						subbeta *= next_dir.bsdf / next_dir.pdf * std::abs(next_dir.direction * subhit.primitive_normal);
						for (int sublen = 3; sublen <= m_max_len; ++sublen)
						{
							if (scene.full_intersection(ray, subhit))
							{ 
								const Material& submaterial = *subhit.geometry->getMaterial();
								subres += subbeta * submaterial.Le(subhit.facing, subhit.tex_uv);
								if (submaterial.spicky())
								{
									submaterial.sampleBSDF(subhit, next_dir, sampler, false);
									ray = { subhit.point, next_dir.direction };
									subbeta *= next_dir.bsdf / next_dir.pdf * std::abs(next_dir.direction * subhit.primitive_normal);
								}
								else
								{
									subres += subbeta * addPhotons(subhit, sublen);
									break;
								}
							}
							else
							{
								subres += subbeta * scene.getBackgroundColor(ray.direction());
								break;
							}
						}
					}

					res += subres / N;
				}
				else
				{
					res += beta * addPhotons(hit, len);
				}
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
			m_frame.fill(Image::MultiSample<RGBColor>());

			// 1 - Rendering time
			LARGE_INTEGER frequency;        // ticks per second
			LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			QueryPerformanceFrequency(&frequency);
			// start timer
			QueryPerformanceCounter(&t1);

			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;
			const size_t number_of_pixels = m_frame.size();
			const size_t sample_pass = number_of_pixels;
			size_t pass = 0;
			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_sample_per_pixel; ++passPerPixelCounter)
			{

				::std::cout << "Pass: " << pass << "/" << Integrator::m_sample_per_pixel << ::std::endl;

				m_map.init();
				buildMap(scene, m_frame.size() * passPerPixelCounter);

				OMP_PARALLEL_FOR
					for (long y = 0; y < m_frame.height(); y++)
					{
						int tid = omp_get_thread_num();

						for (size_t x = 0; x < visu.width(); x++)
						{

							size_t seed = pixelSeed(x, y, m_frame.width(), m_frame.height(), pass);
							Math::Sampler sampler(seed);

							double xp = sampler.generateContinuous<double>();
							double yp = sampler.generateContinuous<double>();

							double v = ((double)y + yp) / m_frame.height();
							double u = ((double)x + xp) / visu.width();

							Ray ray = scene.m_camera.getRay(u, v);


							RGBColor result = sendRay<false>(scene, ray, sampler);



							Image::MultiSample<RGBColor>& pixel = m_frame(x, y);
							pixel.add(result);

							visu.plot(x, y, pixel.mean());
						}//pixel x
					}//pixel y
					//the pass has been computed
				++pass;

				scene.update_lights_offset(1);
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
					Image::ImWrite::write(m_frame);
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
					Image::ImWrite::write(m_frame);
				}
			}
			m_frame.resize(1, 1);
		}


		virtual void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult& res)final override
		{
			if (!m_map.built())
			{
				buildMap(scene);
			}
			resizeFrameBuffer(width, height);
			m_frame.fill(Image::MultiSample<RGBColor>());

			// 1 - Rendering time
			LARGE_INTEGER frequency;        // ticks per second
			LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			QueryPerformanceFrequency(&frequency);
			// start timer
			QueryPerformanceCounter(&t1);

			const size_t number_of_pixels = m_frame.size();
			const size_t sample_pass = number_of_pixels;
			size_t pass = 0;
			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_sample_per_pixel; ++passPerPixelCounter)
			{

				std::cout << '\r' + progession_bar(pass, m_sample_per_pixel, 100) << std::flush;


				OMP_PARALLEL_FOR
					for (long y = 0; y < m_frame.height(); y++)
					{
						int tid = omp_get_thread_num();

						for (size_t x = 0; x < width; x++)
						{
							size_t seed = pixelSeed(x, y, m_frame.width(), m_frame.height(), pass);
							Math::Sampler sampler(seed);

							double xp = sampler.generateContinuous<double>();
							double yp = sampler.generateContinuous<double>();

							double u = ((double)x + xp) / width;
							double v = ((double)y + yp) / height;

							Ray ray = scene.m_camera.getRay(u, v);


							RGBColor result = sendRay<false>(scene, ray, sampler);

							Image::MultiSample<RGBColor>& pixel = m_frame(x, y);
							pixel.add(result);

						}//pixel x
					}//pixel y
					//the pass has been computed
				++pass;

				scene.update_lights_offset(1);

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
					for (long i = 0; i < m_frame.size(); ++i)
					{
						res.image.m_data[i] = m_frame.m_data[i].mean();
					}
			}

		}



		void fastRender(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			if (!m_map.built())
			{
				buildMap(scene);
			}
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
						RGBColor result = sendRay<true>(scene, ray, sampler);

						visu.plot(x, y, result);
					}
				}
			visu.update();
		}




		void debug(Scene const& scene, Visualizer::Visualizer& visu)final override
		{
			if (!m_map.built())
			{
				buildMap(scene);
			}
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

					RGBColor sample = sendRay<false>(scene, ray, sampler);

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