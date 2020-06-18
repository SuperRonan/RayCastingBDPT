#pragma once

#include <Integrators/Integrator.h>
#include <Math/Vectorf.h>
#include <Integrators/PhotonMap.h>
#include <Image/Image.h>
#include <Image/MultiSample.h>
#include <settings.h>
#include <omp.h>
#include <Image/ImWrite.h>

namespace Integrator
{
	class ProgressivePhotonMapper : public Integrator
	{
	protected:

		Image::Image<RGBColor> m_frame;

		Image::Image<RGBColor> m_emissive;

		Image::Image<PhotonId> m_importons_index;



		void showFrame(Visualizer::Visualizer& visu, size_t total)const
		{
			if (visu.visible())
			{
				OMP_PARALLEL_FOR
					for (long x = 0; x < m_frame.width(); ++x)
					{
						for (size_t y = 0; y < m_frame.height(); ++y)
						{
							visu.plot(x, y, m_frame(x, y) / double(total));
						}
					}
			}
		}

		void resizeFrameBuffer(size_t w, size_t h)
		{
			m_frame.resize(w, h);
		}

		using Vector3i = Math::Vector<int, 3>;
		using Vector3f = Math::Vector3f;
		using Vector2f = Math::Vector2f;

		using PhotonFloat = float;

		template <class Float>
		class Importon
		{
		public:
			const Primitive* m_primitive;
			Math::Vector<Float, 3> m_point;
			uint8_t m_depth;
			Math::Vector<Float, 3> m_dir;
			Float m_beta[3];
			
			int m_pixel;

			Float m_radius2;
			int NdivbyAlpha;
			int M = 0;
			Float n_tau[3];
			Float m_tau[3];

			Importon(Hit const& hit, uint8_t depth, RGBColor const& beta, int pixel) :
				m_primitive(hit.primitive),
				m_point(hit.point),
				m_depth(depth),
				m_dir(hit.to_view),
				m_pixel(pixel)
			{
				m_beta[0] = beta[0];
				m_beta[1] = beta[1];
				m_beta[2] = beta[2];
			}

			void fillHit(Hit& res)const
			{
				res.primitive = m_primitive;
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

			RGBColor Ntau()const
			{
				return RGBColor(n_tau[0], n_tau[1], n_tau[2]);
			}

			void setNTau(RGBColor const& t)
			{
				for (int i = 0; i < 3; ++i)
					n_tau[i] = t[i];
			}

			RGBColor Mtau()const
			{
				return RGBColor(m_tau[0], m_tau[1], m_tau[2]);
			}

			void setMTau(RGBColor const& t)
			{
				for (int i = 0; i < 3; ++i)
					m_tau[i] = t[i];
			}
		};

		using Importonf = Importon<PhotonFloat>;

	protected:

		PhotonMap<Importonf> m_map;

		double m_relative_radius;
		double m_radius, m_radius2;

		unsigned int m_number_of_photons;

		unsigned int m_spicky_samples;


		double alpha = 0.9999;

	public:

		ProgressivePhotonMapper(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
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

			m_number_of_photons = pcount;

			m_spicky_samples = 1;
		}

		void buildMap(Scene const& scene, int offset = 0)
		{
			//tic();
			OMP_PARALLEL_FOR
				for (int x = 0; x < m_frame.width(); ++x)
				{
					for (int y = 0; y < m_frame.height(); ++y)
					{
						int index = m_frame.index(x, y);
						Math::Sampler sampler(index + offset);

						double xp = sampler.generateContinuous<double>()*0+0.5;
						double yp = sampler.generateContinuous<double>()*0+0.5;

						double u = ((double)x + xp) / m_frame.width();
						double v = ((double)y + yp) / m_frame.height();

						Ray ray = scene.m_camera.getRay(u, v);

						Hit hit;
						RGBColor beta = scene.m_camera.We(ray.direction()) / scene.m_camera.pdfWeSolidAngle(ray.direction());
						bool broke = false;
						for (int len = 2; len <= m_max_len - 1; ++len)
						{
							if (scene.full_intersection(ray, hit))
							{
								if (hit.geometry->getMaterial()->is_emissive())
								{
									m_emissive[index] += beta * hit.geometry->getMaterial()->Le(hit.primitive_normal, hit.tex_uv, hit.to_view);
								}
								if (!hit.geometry->getMaterial()->delta())
								{
									Importonf imp = Importonf(hit, len - 2, beta, index);
									imp.m_radius2 = m_radius2;
									imp.setNTau(0);
									imp.NdivbyAlpha = 0;
									imp.setMTau(0);
									imp.M = 0;
									
									PhotonId map_index = m_map.addPhoton(imp);
									m_importons_index[index] = map_index;
									broke = true;
									break;
								}
								DirectionSample dirSample;
								hit.geometry->getMaterial()->sampleBSDF(hit, dirSample, sampler, true);
								beta *= dirSample.bsdf / dirSample.pdf * std::abs(dirSample.direction * hit.primitive_normal);
								ray = { hit.point, dirSample.direction };
							}
						}
						if (!broke)
						{
							m_importons_index[index] = PhotonId::invalid();
						}
					}
				}
					
				
			//toc();
			m_map.buildDone();
		}

		__forceinline bool accept(Hit const& hit, Importonf const& imp)const
		{
			Vector3f d = hit.point - imp.m_point;
			const double dist2 = d.norm2();
			d /= std::sqrt(dist2);
			return dist2 < imp.m_radius2 && 
				std::abs(hit.primitive_normal * imp.pnormal() > 0.8) && 
				std::abs(hit.primitive_normal * d) < 0.3;
		}

		void updateImportons(RGBColor const& beta, Hit const& hit, int len)
		{
			m_map.loopThroughPhotons([&](Importonf & imp) {
				int path_len = imp.m_depth + 2 + len - 1;
				if (path_len <= m_max_len)
				{
					double radius2 = imp.m_radius2;
					
					if (accept(hit, imp))
					{
						
						//m_mutex[imp.m_pixel].lock();
						RGBColor tau = beta * hit.geometry->getMaterial()->BSDF(hit, imp.m_dir, hit.to_view);

						double g = (imp.NdivbyAlpha*m_alpha+m_alpha) / (imp.NdivbyAlpha*m_alpha+1.0);
						RGBColor new_tau = (imp.Ntau() + tau) * g;
						double new_radius2 = radius2 * g;

						imp.NdivbyAlpha++;
						imp.setNTau(new_tau);
						imp.m_radius2 = new_radius2;

						//imp.M++;
						//imp.setMTau(imp.Mtau() + tau);
							
						//m_mutex[imp.m_pixel].unlock();
					}
				}
				}, hit.point);
		}


		void sendPhoton(Scene const& scene, Math::Sampler& sampler, int offset)
		{
			SurfaceSample sls;
			scene.sampleLe(sampler, sls, offset);

			
			Hit hit;
			DirectionSample ds = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);
			RGBColor beta = sls.geo->getMaterial()->Le(sls.normal, sls.uv, ds.direction) / (sls.pdf * ds.pdf) * std::abs(ds.direction * sls.normal);
			Ray ray(sls.vector, ds.direction);

			for (int len = 2; len < m_max_len; ++len)
			{
				if (scene.full_intersection(ray, hit))
				{
					if (!hit.geometry->getMaterial()->spicky())
					{
						updateImportons(beta, hit, len);
					}

					hit.geometry->getMaterial()->sampleBSDF(hit, ds, sampler, true);
					beta *= ds.bsdf / ds.pdf * std::abs(ds.direction * hit.primitive_normal);
					ray = { hit.point, ds.direction };
				}
			}
		}

		void render(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			resizeFrameBuffer(visu.width(), visu.height());
			m_emissive.resize(visu.width(), visu.height());
			m_importons_index.resize(visu.width(), visu.height());
			m_frame.fill(0);
			m_emissive.fill(0);


			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;

			m_map.init();
			buildMap(scene, m_frame.size() * 0);
			ProgressReporter reporter;
			reporter.start(m_sample_per_pixel);
			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				OMP_PARALLEL_FOR
					for (int p = 0; p < m_number_of_photons; ++p)
					{
						Math::Sampler sampler(p + m_number_of_photons * pass);
						sendPhoton(scene, sampler, p + pass * m_number_of_photons);
					}

				OMP_PARALLEL_FOR
				for (int i = 0; i < m_importons_index.size(); ++i)
				{
					if (!m_importons_index[i].isInvalid())
					{
						Importonf& imp = m_map[m_importons_index[i]];
						//double g = (imp.NdivbyAlpha * alpha + imp.M * alpha) / std::max(1.0, (imp.NdivbyAlpha * alpha + imp.M));
						//if (g == 0)	g = 1;
						//imp.NdivbyAlpha += imp.M;
						//imp.setNTau(imp.Ntau() + imp.Mtau() * g);
						//imp.m_radius2 *= g;
						//m_frame[i] = imp.beta() * imp.Ntau() / (Math::pi * imp.m_radius2 * m_number_of_photons);
						//imp.M = 0;
						//imp.setMTau(0);

						m_frame[i] = imp.beta() * imp.Ntau() / (Math::pi * imp.m_radius2 * m_number_of_photons) + m_emissive[i] * (pass + 1);
					}
				}
				
				showFrame(visu, pass+1);
				
				reporter.report(pass + 1, -1);
				scene.update_lights_offset(1);
				kbr = visu.update();

				if (kbr == Visualizer::Visualizer::KeyboardRequest::done)
				{
					goto __render__end__loop__;
				}
				else if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
				{
					Image::ImWrite::write(m_frame, 1.0 / (double)pass);
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
					Image::ImWrite::write(m_frame, 1.0 / (double)m_sample_per_pixel);
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
						res.image.m_data[i] = m_frame.m_data[i];
					}
			}

		}



		void debug(Scene const& scene, Visualizer::Visualizer& visu)final override
		{
			
		}





	};
}