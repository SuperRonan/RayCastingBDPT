#pragma once

#include <Integrators/BidirectionalBase.h>
#include <Integrators/PhotonMap.h>
#include <omp.h>
#include <Image/ImWrite.h>
#include <utils.h>
#include <System/ScopedAssignment.h>

namespace Integrator
{
	class VCM : public BidirectionalBase
	{
	protected:

		struct Photon
		{
			Math::Vector3f m_point;

			Math::Vector3f point()const
			{
				return m_point;
			}

			unsigned int pixel;
			unsigned char len;
		};

	public:


		

		mutable PhotonMap<Photon> m_map;

		double m_relative_radius;
		double m_radius, m_radius2;
		double m_photon_emitted;

		mutable Image::Image<Path> m_light_paths;
		bool m_light_paths_built = false;



		VCM(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			BidirectionalBase(sample_per_pixel, width, height)
		{

		}


		virtual void setLen(unsigned int len)override
		{
			Integrator::setLen(len);
		}

		// Relative readius is relative to the scene radius
		void setParams(Scene const& scene, double relative_radius, double photon_relative_radius = 1)
		{
			m_relative_radius = relative_radius;
			BoundingBox m_bb = scene.m_sceneBoundingBox;
			m_bb[0] -= Math::Vector3f(0.001, 0.001, 0.001);
			m_bb[1] += Math::Vector3f(0.001, 0.001, 0.001);
			Math::Vector3f dim = m_bb.diag();
			double max_dir = dim.simdAbs().max();

			m_radius = max_dir * m_relative_radius;
			
			Math::Vector3f m_pixel_size = m_radius;
			Math::Vector3f sizef = dim.simdDiv(m_pixel_size);
			Math::Vector<int, 3> m_size = sizef.ceil();
			
			m_radius *= photon_relative_radius;
			m_radius2 = m_radius * m_radius;

			m_map.init(m_bb, m_size);
		}


		void buildLightPaths(Scene const& scene, size_t width, size_t height, size_t pass)
		{
			m_light_paths.resize(width, height);
			m_map.dumpPhotons();
			m_light_paths_built = false;
			OMP_PARALLEL_FOR
			for (int x = 0; x < width; ++x)
			{
				for (size_t y = 0; y < height; ++y)
				{
					Path& path = m_light_paths(x, y);
					path.reset();
					Math::Sampler sampler(pixelSeed(x, y, width, height, pass) * 2);
					traceLightSubPath(scene, sampler, path);
					
					// Fill the photon map
					for (int i = 1; i < path.size(); ++i)
					{
						if (!path[i].hit.geometry->getMaterial()->spicky())
						{
							Photon photon;
							photon.m_point = path[i].hit.point;
							photon.pixel = m_light_paths.index(x, y);
							photon.len = i + 1;
							m_map.addPhoton(photon);
						}
					}
				}
			}
			m_photon_emitted = width * height;
			m_light_paths_built = true;
			m_map.buildDone();
		}


		
		void VertexConnection(Scene const& scene, Path& cameraSubPath, Path& lightSubPath, Math::Sampler& sampler, RGBColor& pixel_res, LightVertexStack& lt_vertices)const
		{
			int first_t_not_spicky=-1;
			double s1_pdf;
			for (int t = 1; t <= cameraSubPath.size(); ++t)
			{
				Vertex& camera_top = cameraSubPath[t - 1];
				if (first_t_not_spicky == -1 && t >= 2 && !camera_top.hit.geometry->getMaterial()->spicky())
					first_t_not_spicky = t;
				int last_s_not_spicky = -1;
				for (int s = 0; s <= lightSubPath.size(); ++s)
				{
					if (s + t > m_max_len)
					{
						break;
					}
					RGBColor L = 0;

					// special cases of connections strategies
					if (s + t == 1)
						continue;
					if (s == 0)
					{
						// naive path tracing
						if (camera_top.hit.geometry->getMaterial()->is_emissive())
						{
							L = camera_top.beta * camera_top.hit.geometry->getMaterial()->Le(camera_top.pNormal(), camera_top.hit.tex_uv, camera_top.omega_o());
							double pdf = scene.pdfSampleLe(camera_top.hit.geometry);
							s1_pdf = scene.pdfSampleLi(camera_top.hit.geometry, cameraSubPath[t - 2].hit, camera_top.hit.point);
							double weight = VCWeight(cameraSubPath, lightSubPath, s, t, first_t_not_spicky, last_s_not_spicky, s1_pdf, pdf);
							pixel_res += L * weight;
						}
					}
					else
					{
						ScopedAssignment<Vertex> resampled_vertex_sa;
						if (s == 1) 
						{
							SurfaceSample sls;
							scene.sampleLi(sampler, sls, camera_top.hit);
							s1_pdf = sls.pdf;
							Vertex light_resampled;
							light_resampled.delta = false;
							light_resampled.type = Vertex::Type::Light;
							light_resampled.hit.geometry = sls.geo;
							light_resampled.hit.normal = light_resampled.hit.primitive_normal = sls.normal;
							light_resampled.hit.tex_uv = sls.uv;
							light_resampled.hit.point = sls.vector;
							double Le_pdf = scene.pdfSampleLe(sls.geo);
							light_resampled.fwd_pdf = Le_pdf;
							light_resampled.beta = 1.0 / sls.pdf;
							resampled_vertex_sa = { lightSubPath.begin(), light_resampled };
						}
						else if (s == 2) // BTW this assumes that the connecting loops are in the order t then s
						{
							s1_pdf = scene.pdfSampleLi(lightSubPath[0].hit.geometry, lightSubPath[1].hit, lightSubPath[0].hit.point);
						}

						Vertex& light_top = lightSubPath[s - 1];
						if (camera_top.delta || light_top.delta)
							continue;
						
						if (s >= 2 && !light_top.hit.geometry->getMaterial()->spicky())
							last_s_not_spicky = s;
						
						Math::Vector3f dir = camera_top.hit.point - light_top.hit.point;
						const double dist2 = dir.norm2();
						const double dist = sqrt(dist2);
						dir /= dist;

						double G = std::abs(dir * light_top.hit.primitive_normal) / dist2;

						RGBColor camera_connection, light_connection;

						if (t == 1)
						{
							//light tracing
							camera_connection = scene.m_camera.We(-dir);
						}
						else
						{
							//general case
							camera_connection = camera_top.hit.geometry->getMaterial()->BSDF(camera_top.hit, -dir, false);
							G *= std::abs(camera_top.hit.primitive_normal * dir);
						}

						if (s == 1)
						{
							//direct lighting estimation
							light_connection = light_top.hit.geometry->getMaterial()->Le(light_top.pNormal(), light_top.hit.tex_uv, dir);
						}
						else
						{
							//general case
							light_connection = light_top.hit.geometry->getMaterial()->BSDF(light_top.hit, dir, true);
						}

						L = camera_top.beta * camera_connection * G * light_connection * light_top.beta;
						

						if (!L.isBlack() && scene.visibility(light_top.hit.point, camera_top.hit.point))
						{
							L *= VCWeight(cameraSubPath, lightSubPath, s, t, first_t_not_spicky, last_s_not_spicky, s1_pdf);
							if (t == 1)
							{
								LightVertex lv;
								lv.light = L / scene.m_camera.resolution;
								lv.uv = scene.m_camera.raster(-dir);
								lt_vertices.push(lv);
							}
							else
							{
								pixel_res += L;
							}
						}

					}
				}
			}
		}

		__forceinline bool accept(Hit const& hit, Vertex const& photon)const
		{
			Math::Vector3f d = hit.point - photon.hit.point;
			const double dist2 = d.norm2();
			d /= std::sqrt(dist2);
			return dist2 < m_radius2 && 
				std::abs(hit.primitive_normal * photon.hit.primitive_normal > 0.8) &&
				std::abs(hit.primitive_normal * d) < 0.3;
		}
		
		RGBColor VertexMerging(Scene const& scene, Path& cameraSubPath, Math::Sampler& sampler)const
		{
			RGBColor res = 0;
			for (int t = 2; t <= cameraSubPath.size(); ++t)
			{
				const Vertex& pt = cameraSubPath[t-1];
				if (!pt.material()->spicky())
				{
					m_map.loopThroughPhotons([&](Photon const& photon)
						{
							const int len = photon.len + t - 1;
							if (len <= m_max_len)
							{
								Path& lightSubPath = m_light_paths[photon.pixel];
								const int s = photon.len - 1;
								assert(s >= 1);
								const Vertex& qs_plus = lightSubPath[photon.len - 1];
								const Vertex& qs = lightSubPath[s - 1];
								if (accept(pt.hit, qs_plus))
								{
									const double k = 1.0 / (Math::pi * m_radius2);
									RGBColor L = qs_plus.beta * pt.hit.geometry->getMaterial()->BSDF(pt.hit, qs_plus.hit.to_view, pt.hit.to_view) * pt.beta;
									
									const double s1_pdf = scene.pdfSampleLi(lightSubPath[0].hit.geometry, s == 1 ? pt.hit : lightSubPath[1].hit, lightSubPath[0].hit.point);
									const double w = VMWeight(cameraSubPath, lightSubPath, s, t, s1_pdf);

									res += L * k / m_photon_emitted * w;
								}
							}
						}, pt.hit.point);
					break;
				}
			}
			return res;
		}

		void computeSample(Scene const& scene, int pix_index, double u, double v, __in Math::Sampler& sampler, __out RGBColor& pixel_res, LightVertexStack& lt_vertices)const
		{
			pixel_res = 0;
			Ray ray = scene.m_camera.getRay(u, v);
			Path cameraSubPath;
			traceCameraSubPath(scene, sampler, cameraSubPath, ray);
			Path& lightSubPath = m_light_paths[pix_index];


			VertexConnection(scene, cameraSubPath, lightSubPath, sampler, pixel_res, lt_vertices);

			pixel_res += VertexMerging(scene, cameraSubPath, sampler);
		}


		////////////////////////////////////////////////////////////////
		//Computes te MIS weights for Vertex Connection
		// - the last parameter if the probability of sampling the last point on the camera sub path if s == 0 (pure path tracing), else it is not necessary 
		////////////////////////////////////////////////////////////////
		double VCWeight(
			Path& cameras, Path& lights,
			const int main_s, const int main_t,
			const int first_t_not_spicky, const int last_s_not_spicky,
			double s1_pdf, double pdf_sampling_point = -1)const
		{
			Vertex* xt = cameras.begin() + (main_t - 1);
			Vertex* ys = main_s == 0 ? nullptr : lights.begin() + (main_s - 1);
			Vertex* xtm = main_t == 1 ? nullptr : cameras.begin() + (main_t - 2);
			Vertex* ysm = main_s < 2 ? nullptr : lights.begin() + (main_s - 2);

			ScopedAssignment<double> xt_pdf_rev_sa;
			ScopedAssignment<double> ys_pdf_rev_sa;
			ScopedAssignment<double> xtm_pdf_rev_sa;
			ScopedAssignment<double> ysm_pdf_rev_sa;

			xt_pdf_rev_sa = { &xt->rev_pdf,[&]() {
				if (main_t == 1)
					return 0.0;
				else if (ys)
					return ys->pdf<TransportMode::Radiance, true>(*xt, ys->hit.to_view);
				else
					return pdf_sampling_point;
			}() };
			assert(xt->rev_pdf >= 0);

			if (ys)
			{
				ys_pdf_rev_sa = { &ys->rev_pdf, xt->pdf<TransportMode::Importance, true>(*ys, xt->hit.to_view) };
			}

			if (xtm)
			{
				xtm_pdf_rev_sa = { &xtm->rev_pdf, xt->pdf<TransportMode::Importance, true>(*xtm, xt->dir_to_vertex(ys)) };
			}

			if (ysm)
			{
				ysm_pdf_rev_sa = { &ysm->rev_pdf, ys->pdf<TransportMode::Radiance, true>(*ysm, ys->dir_to_vertex(xt)) };
			}

			const double actual_ni = main_t == 1 ? cameras[0].hit.camera->resolution : 1;
			const double actual_main_ri = (main_s == 1 ? s1_pdf / lights[0].fwd_pdf : 1);

			
			double sum = actual_main_ri * actual_ni;

			sum += sumRatioVC(cameras, lights, main_s, main_t, s1_pdf);

			// Find the VM pdf
			if (main_s + main_t > 2)
			{
				double vm_ri = 1.0;
				if (first_t_not_spicky != -1 && !(first_t_not_spicky == main_t && main_s == 0)) // A merge is possible on the camera subpath
				{
					const Vertex* camera_end = xt;
					const Vertex* light_end = ys; // the photon is merged with camera_end
					for (int t = main_t; t > first_t_not_spicky; --t)
					{
						light_end = camera_end;
						camera_end = &cameras[t - 2];

						vm_ri *= light_end->rev_pdf / light_end->fwd_pdf;
					}
					double pacc = Math::pi * m_radius2 * camera_end->rev_pdf;
					vm_ri *= pacc;
					assert(vm_ri >= 0);
					sum += vm_ri * m_photon_emitted;
				}
				else if (first_t_not_spicky == -1 && last_s_not_spicky != -1) // The merge is on the light subpath
				{
					const Vertex* camera_end = nullptr;
					const Vertex* light_end;
					assert(main_s >= last_s_not_spicky);
					for (int s = main_s; s >= last_s_not_spicky; --s)
					{
						camera_end = &lights[s - 1];
						light_end = &lights[s - 2];

						vm_ri *= camera_end->rev_pdf / camera_end->fwd_pdf;
					}
					double pacc = Math::pi * m_radius2 * camera_end->fwd_pdf;
					vm_ri *= pacc;
					sum += vm_ri * m_photon_emitted;
				}
			}
			

			double weight = (actual_main_ri*actual_ni) / sum;
			return weight;
		}


		////////////////////////////////////////////////////////////////
		//Computes te MIS weights for the Vertex Merging
		////////////////////////////////////////////////////////////////
		double VMWeight(
			Path& cameras, Path& lights,
			const int main_s, const int main_t,
			double s1_pdf)const
		{
			assert(lights.size() >= main_s + 1);
			assert(main_s > 0);
			assert(main_t >= 2);
		
			Vertex* xt = &cameras[main_t - 1];
			Vertex* ys = &lights[main_s - 1];
			Vertex* xtm = main_t == 1 ? nullptr : cameras.begin() + (main_t - 2);
			Vertex* ysm = main_s < 2 ? nullptr : lights.begin() + (main_s - 2);
			const Vertex* photon = &lights[main_s];

			ScopedAssignment<double> xt_pdf_rev_sa;
			ScopedAssignment<double> ys_pdf_rev_sa;
			ScopedAssignment<double> xtm_pdf_rev_sa;
			ScopedAssignment<double> ysm_pdf_rev_sa;

			//xt_pdf_rev_sa = { &xt->pdf_rev, ys->pdf<TransportMode::Radiance, true>(*xt, ys->hit.to_view) };
			xt_pdf_rev_sa = { &xt->rev_pdf, photon->fwd_pdf }; // Merging

			ys_pdf_rev_sa = { &ys->rev_pdf, xt->pdf<TransportMode::Importance, true>(*ys, xt->hit.to_view) };

			if (xtm)
			{
				xtm_pdf_rev_sa = { &xtm->rev_pdf, xt->pdf<TransportMode::Importance, true>(*xtm, xt->dir_to_vertex(ys)) };
			}

			if (ysm)
			{
				// Since it is already computed
				//ysm_pdf_rev_sa = { &ysm->pdf_rev, ys->pdf<TransportMode::Radiance, true>(*ysm, ys->dir_to_vertex(xt)) };
			}

			const double actual_ni = m_photon_emitted;
			const double actual_main_ri = photon->fwd_pdf * Math::pi * m_radius2;


			double sum = actual_main_ri * actual_ni;

			sum += sumRatioVC(cameras, lights, main_s, main_t, s1_pdf);
			sum += (main_s == 1) ? (s1_pdf / ys->fwd_pdf) : 1.0;

			double weight = (actual_main_ri * actual_ni) / sum;
			return weight;
		}










		static __forceinline size_t pixelSeed(size_t x, size_t y, size_t width, size_t height, size_t pass)
		{
#ifdef SAMPLER_BIAS
			return (pass) * (width * height);
#else
			return (y * width + x) + (pass) * (width * height);
#endif
		}


		void render(Scene const& scene, Visualizer::Visualizer& visu)final override
		{
			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;

			m_frame_buffer.resize(visu.width(), visu.height());
			m_frame_buffer.fill();
			visu.clean();
			ProgressReporter reporter;
			reporter.start(m_sample_per_pixel);
			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				buildLightPaths(scene, visu.width(), visu.height(), pass);
				const size_t sample_pass = m_frame_buffer.size();
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

							double v = ((double)y + yp) / visu.height();
							double u = ((double)x + xp) / visu.width();

							RGBColor pixel = 0;
							LightVertexStack lvs;

							computeSample(scene, m_light_paths.index(x, y), u, v, sampler, pixel, lvs);

							m_frame_buffer(x, y) += pixel;
							for (LightVertex const& lv : lvs)
							{
								int lx = lv.uv[0] * visu.width();
								int ly = lv.uv[1] * visu.height();
								m_frame_buffer(lx, ly) += lv.light;
							}
						}//pixel x
					}//pixel y
					//the pass has been computed
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
					Image::ImWrite::write(m_frame_buffer, 1.0 / double(pass+1));
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
					Image::ImWrite::write(m_frame_buffer, 1.0 / (double)m_sample_per_pixel);
				}
			}
			m_light_paths_built = false;
		}










		void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult& res)final override
		{
			m_frame_buffer.resize(width, height);
			m_frame_buffer.fill();
			ProgressReporter reporter;
			reporter.start(m_sample_per_pixel);
			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				const size_t sample_pass = m_frame_buffer.size();
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

							double v = ((double)y + yp) / m_frame_buffer.height();
							double u = ((double)x + xp) / m_frame_buffer.width();

							RGBColor pixel = 0;
							LightVertexStack lvs;

							computeSample(scene, m_light_paths.index(x, y), u, v, sampler, pixel, lvs);

							m_frame_buffer(x, y) += pixel;
							for (LightVertex const& lv : lvs)
							{
								int lx = lv.uv[0] * m_frame_buffer.width();
								int ly = lv.uv[1] * m_frame_buffer.height();
								m_frame_buffer(lx, ly) += lv.light;
							}

							//visu.plot(x, y, pixel);
						}//pixel x
					}//pixel y
					//the pass has been computed
				scene.update_lights_offset(1);
				reporter.report(pass + 1, -1);
			}//pass per pixel

			reporter.finish();
			scene.reset_surface_lights();

			//fill the result
			{
				res.time = reporter.time();
				res.image.resize(width, height);
				OMP_PARALLEL_FOR
					for (long i = 0; i < m_frame_buffer.size(); ++i)
					{
						res.image.m_data[i] = m_frame_buffer.m_data[i] / m_sample_per_pixel;
					}
			}
		}




		void fastRender(Scene const& scene, Visualizer::Visualizer& visu)final override
		{
			m_frame_buffer.resize(visu.width(), visu.height());
			m_frame_buffer.fill();
			visu.clean();
			if (!m_light_paths_built)
				buildLightPaths(scene, visu.width(), visu.height(), 0);
			OMP_PARALLEL_FOR
				for (long y = 0; y < m_frame_buffer.height(); y++)
				{
					int tid = omp_get_thread_num();
					double v = ((double)y + 0.5) / m_frame_buffer.height();
					for (size_t x = 0; x < visu.width(); x++)
					{
						double u = ((double)x + 0.5) / visu.width();

						Ray ray = scene.m_camera.getRay(u, v);
						size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), 0);
						Math::Sampler sampler(seed);

						RGBColor pixel = 0;
						LightVertexStack lvs;

						computeSample(scene, m_light_paths.index(x, y), u, v, sampler, pixel, lvs);

						m_frame_buffer(x, y) += pixel;
						for (LightVertex const& lv : lvs)
						{
							int lx = lv.uv[0] * visu.width();
							int ly = lv.uv[1] * visu.height();
							m_frame_buffer(lx, ly) += lv.light;
						}

					}//pixel x
				}//pixel y
				//the pass has been computed
			
			showFrame(visu, 1);

			visu.update();
		}

		void debug(Scene const& scene, Visualizer::Visualizer& visu) final override
		{

		}




	};
}