#pragma once

#include <Integrators/BidirectionalBase.h>
#include <Integrators/PhotonMap.h>
#include <System/ScopedAssignment.h>
#include <omp.h>
#include <Image/ImWrite.h>

namespace Integrator
{
	class SimpleVCM : public BidirectionalBase
	{
	protected:

		using Vector3i = Math::Vector<int, 3>;
		using Vector3f = Math::Vector3f;
		using Vector2f = Math::Vector2f;

		using PhotonFloat = float;

		template <class Float>
		class Photon
		{
		public:

			const Primitive* m_primitive;
			const GeometryBase* m_geometry;
			Math::Vector<Float, 3> m_point;
			uint8_t m_len;
			Math::Vector<Float, 3> m_dir;
			Float m_beta[3];
			
			Float pdfProd = 0.0;
			Float pt_pdf, npt_pdf;
			PhotonId prev;

			Photon(Hit const& hit, uint8_t len, RGBColor const& beta) :
				m_geometry(hit.geometry),
				m_primitive(hit.primitve),
				m_point(hit.point),
				m_len(len),
				m_dir(hit.to_view)
			{
				m_beta[0] = beta[0];
				m_beta[1] = beta[1];
				m_beta[2] = beta[2];
			}

			Photon(const GeometryBase* geo, const Primitive* primitive, Vector3f const& point, Vector3f const& dir, uint8_t len, RGBColor const& beta) :
				m_geometry(geo),
				m_primitive(primitive),
				m_point(point),
				m_len(len),
				m_dir(dir)
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

			Vector3f pnormal()const
			{
				return m_primitive->normal(m_point, m_primitive->tuv(m_primitive->uv(m_point)));
			}

			RGBColor beta()const
			{
				return RGBColor(m_beta[0], m_beta[1], m_beta[2]);
			}

			Math::Vector<Float, 3> point()const
			{
				return m_point;
			}

		};

		using Photonf = Photon<PhotonFloat>;

		PhotonMap<Photonf> m_map;

		double m_relative_radius;
		double m_radius, m_radius2;

		unsigned int m_number_of_photons;

		double kernel(double dist2)const
		{
			return dist2 <= m_radius2 ? 1.0 / (Math::pi * m_radius2) : 0.0;
		}

	public:

		SimpleVCM(unsigned int sample_per_pixel, unsigned int width, unsigned int height):
			BidirectionalBase(sample_per_pixel, width, height)
		{}

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
			Math::Vector3f m_pixel_size = m_radius * 2.0;
			Vector3f sizef = dim.simdDiv(m_pixel_size);
			Vector3i m_size = sizef.ceil();

			m_map.init(m_bb, m_size);

			m_number_of_photons = pcount;
		}


		void buildMap(Scene const& scene, int offset = 0)
		{
			//tic();
			OMP_PARALLEL_FOR
				for (int sample = 0; sample < m_number_of_photons; ++sample)
				{
					Math::Sampler sampler(sample + offset);
					SurfaceSample sls;
					sampleOneLight(scene, sampler, sls);
					double prod_pdf = sls.pdf;
					
					Photonf light_photon = Photonf(sls.geo, sls.primitive, sls.vector, 0.0, 1, 1);
					light_photon.pdfProd = prod_pdf;
					
					
					DirectionSample dirSample = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);
					double pdf_solid_angle = sls.pdf;
					double pt_pdf = 1, npt_pdf = 1;
					RGBColor beta = dirSample.bsdf / (sls.pdf * dirSample.pdf) * std::abs(dirSample.direction * sls.normal);
					Hit hit;
					Ray ray(sls.vector, dirSample.direction);

					PhotonId prev_photon;
					double prev_pdf_solid_angle;
					for (int len = 2; len <= m_max_len - 1; ++len)
					{
						Hit prev_hit = hit;
						if (scene.full_intersection(ray, hit))
						{
							if (len == 2)
							{
								prev_photon = m_map.addPhoton(light_photon);
								pt_pdf = scene.pdfSamplingLight(hit.geometry, prev_hit, hit.point);
							}
							else
							{
								
							}
							double pdf_area = pdf_solid_angle * std::abs(hit.primitive_normal * ray.direction()) / (hit.z * hit.z);
							prod_pdf *= pdf_area;

							hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, dirSample, sampler, true);
							{
								Photonf photon = { hit, (uint8_t)(len), beta / ((double)m_number_of_photons) };
								photon.pdfProd = prod_pdf;
								photon.prev = prev_photon;
								prev_photon = m_map.addPhoton(std::move(photon));
							}

							

							double pdf_rev_sa = hit.geometry->getMaterial()->delta() ? dirSample.pdf : hit.geometry->getMaterial()->pdf(hit, -ray.direction(), dirSample.direction);

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






		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const
		{
			Ray ray = pray;
			RGBColor beta = scene.m_camera.We(ray.direction()) / scene.m_camera.pdfWeSolidAngle<true>(ray.direction());
			RGBColor res = 0;
			Hit hit;
			int len;
			for (len = 2; len <= m_max_len; ++len)
			{
				if (scene.full_intersection(ray, hit))
				{
					res += hit.geometry->getMaterial()->Le(hit.facing, hit.tex_uv);
					if (hit.geometry->getMaterial()->spicky())
					{
						DirectionSample ds;
						hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, ds, sampler);
						ray = { hit.point, ds.direction };
						beta *= ds.bsdf / ds.pdf * std::abs(ds.direction * hit.primitive_normal);
					}
					else
						return res + beta * evaluateRadiance(scene, hit, len, sampler);
				}
				else
				{
					return res + beta * scene.getBackgroundColor(ray.direction());
				}
			}
			return res;
		}


		RGBColor connectToLight(Scene const& scene, Hit const& hit, Math::Sampler& sampler, double path_pdf)const
		{
			if (hit.geometry->getMaterial()->delta())
				return 0;
			SurfaceSample sls;
			sampleOneLight(scene, sampler, sls);
			Vector3f to_light = (sls.vector - hit.point);
			const double dist2 = to_light.norm2();
			const double dist = std::sqrt(dist2);
			to_light /= dist;
			RGBColor bsdf = hit.geometry->getMaterial()->BSDF(hit, to_light);
			if (!bsdf.isBlack())
			{
				Ray ray(hit.point, to_light);
				Hit light_hit;
				if (scene.full_intersection(ray, light_hit) && samePoint(light_hit, dist))
				{
					RGBColor contribution = bsdf * light_hit.geometry->getMaterial()->Le(light_hit.facing, light_hit.tex_uv) * 
						std::abs(to_light * hit.primitive_normal) * std::abs(to_light * light_hit.primitive_normal) / dist2;
					double bsdf_pdf = hit.geometry->getMaterial()->pdf(hit, to_light) * std::abs(to_light * light_hit.primitive_normal) / dist2;
					double weight = sls.pdf / (sls.pdf + bsdf_pdf) * 0.5;
					return contribution / (sls.pdf) * weight;
				}
			}
			return 0;
		}
		
		

		RGBColor evaluateRadiance(Scene const& scene, Hit const& base, int base_len, Math::Sampler & sampler)const
		{
			Vector3f const& wo = base.to_view;
			RGBColor pt_res = 0, npt_res = 0, pm_res=0;
			// Path tracings
			{
				double path_pdf = 1;
				Hit hit = base;
				RGBColor beta = 1;
				Ray ray;
				for (int len = base_len; len < m_max_len; ++len)
				{
					pt_res += beta * connectToLight(scene, hit, sampler, path_pdf);
					
					DirectionSample ds;
					hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, ds, sampler);
					beta *= ds.bsdf / ds.pdf * std::abs(hit.normal * ds.direction);
					ray = { hit.point, ds.direction };
					double pdf_solid_angle = ds.pdf;
					bool prev_delta = hit.geometry->getMaterial()->delta();
					Hit prev_hit = hit;
					if (scene.full_intersection(ray, hit))
					{
						RGBColor contrib = beta * hit.geometry->getMaterial()->Le(hit.facing, hit.tex_uv);
						
						double light_pdf = prev_delta ? 0.0 : (path_pdf * scene.pdfSamplingLight(hit.geometry, prev_hit, hit.point));
						path_pdf *= pdf_solid_angle * std::abs(ray.direction() * hit.primitive_normal) / (hit.z * hit.z);
						double weight = path_pdf / (path_pdf + light_pdf) * 0.5;
						npt_res += contrib * weight;
					}
					else
						break;
				}
			}


			// VM
			m_map.loopThroughPhotons([&](Photonf const& photon){
				const int path_len = photon.m_len + base_len - 1;
				if (path_len <= m_max_len && photon.m_len != 1)
				{
					const Vector3f d = base.point - photon.point();
					const double dist2 = d.norm2();
					const double k = kernel(dist2);
					const Material& mat = *photon.m_primitive->geometry()->getMaterial();
					if (kernel(dist2) > 0 && base.geometry == photon.m_primitive->geometry())
					{
						RGBColor contrib = photon.beta() * mat.BSDF(base, photon.m_dir, wo) * k;
						if (!contrib.isBlack())
						{
							double npt_pdf, pt_pdf;
							double path_pdf = 1;
							const Photonf* current_photon = &photon;
							Hit current_hit;
							current_photon->fillHit(current_hit);
							Vector3f prev_wi = wo;
							for (int len = photon.m_len; len > 2; --len)
							{
								const Photonf* next_photon = &m_map[current_photon->prev];
								const Material& photon_mat = *current_hit.geometry->getMaterial();
								Math::Vector3f dir = next_photon->point() - current_photon->point();
								const double dist2 = dir.norm2();
								dir /= std::sqrt(dist2);
								double pdf_solid_angle = photon_mat.pdf(current_hit, dir, prev_wi);
								next_photon->fillHit(current_hit);
								double pdf_area = pdf_solid_angle * std::abs(current_hit.primitive_normal * dir) / dist2;
								path_pdf *= pdf_area;
								prev_wi = -dir;
								current_photon = next_photon;
							}
							assert(current_photon->m_len == 2);
							const Photonf* light_photon = &m_map[current_photon->prev];
							Vector3f dir = light_photon->point() - current_hit.point;
							const double dist2 = dir.norm2();
							dir /= std::sqrt(dist2);
							const Vector3f light_normal = light_photon->pnormal();
							const double cos_light = std::abs(dir * light_normal);

							const double pdf_dir_sa = current_hit.geometry->getMaterial()->pdf(current_hit, dir, prev_wi);
							const double pdf_dir = pdf_dir_sa * cos_light / dist2;
							const double pdf_light = scene.pdfSamplingLight(light_photon->m_geometry, current_hit, light_photon->point());

							const double epdf_pt = path_pdf * pdf_light;
							const double epdf_npt = path_pdf * pdf_dir;
							const double epdf_pm = photon.pdfProd * Math::pi * m_radius2 * m_number_of_photons;


							double weight = epdf_pm / (epdf_npt + epdf_pm + epdf_pt);
							//double weight = 0.5;
							pm_res += contrib * weight;
						}
					}
				}
			}, base.point);

			return (pt_res + npt_res + pm_res);
		}










		void render(Scene const& scene, Visualizer::Visualizer& visu) final override
		{

			resizeFrameBuffer(visu.width(), visu.height());
			m_frame_buffer.fill(0);

			// 1 - Rendering time
			LARGE_INTEGER frequency;        // ticks per second
			LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			QueryPerformanceFrequency(&frequency);
			// start timer
			QueryPerformanceCounter(&t1);

			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;
			const size_t number_of_pixels = m_frame_buffer.size();
			const size_t sample_pass = number_of_pixels;
			size_t pass = 0;
			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_sample_per_pixel; ++passPerPixelCounter)
			{

				::std::cout << "Pass: " << pass << "/" << Integrator::m_sample_per_pixel << ::std::endl;

				m_map.init();
				buildMap(scene, m_frame_buffer.size() * passPerPixelCounter);

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


							

							m_frame_buffer(x, y) += sendRay(scene, ray, sampler);

							visu.plot(x, y, m_frame_buffer(x, y) / (passPerPixelCounter+1));
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
					Image::ImWrite::write(m_frame_buffer, 1.0 / (double)(passPerPixelCounter+1));
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
					Image::ImWrite::write(m_frame_buffer, 1.0 / (double)m_sample_per_pixel);
				}
			}
			m_frame_buffer.resize(1, 1);
		}


		virtual void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult& res)final override
		{
			if (!m_map.built())
			{
				buildMap(scene);
			}
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



							RGBColor& pixel = m_frame_buffer(x, y);
							pixel += sendRay(scene, ray, sampler);

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
					for (long i = 0; i < m_frame_buffer.size(); ++i)
					{
						res.image.m_data[i] = m_frame_buffer.m_data[i];
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
						

						visu.plot(x, y, sendRay(scene, ray, sampler));
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

					//RGBColor sample = sendRay<false>(scene, ray, sampler);

					//samples[pass] = sample;
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