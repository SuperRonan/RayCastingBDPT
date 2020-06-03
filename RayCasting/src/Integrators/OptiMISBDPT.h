#pragma once

#include <Integrators/BidirectionalBase.h>
#include <System/ScopedAssignment.h>
#include <atomic>
#include <MIS/ImageEstimators.h>

namespace Integrator
{
	class OptiMISBDPT : public BidirectionalBase
	{
	protected:

		void computeSample(Scene const& scene, double u, double v, __in Math::Sampler& sampler, double * weights)const
		{
			Path cameraSubPath, LightSubPath;

			Ray ray = scene.m_camera.getRay(u, v);

			traceCameraSubPath(scene, sampler, cameraSubPath, ray);
			traceLightSubPath(scene, sampler, LightSubPath);
			double s1_pdf;
			for (int t = 1; t <= cameraSubPath.size(); ++t)
			{
				Vertex& camera_top = cameraSubPath[t - 1];

				for (int s = 0; s <= LightSubPath.size(); ++s)
				{
					if (s + t > m_max_len)
					{
						break;
					}
					// The balance estimate
					RGBColor L = 0;
					Math::Vector2f p = { u, v };
					// special cases of connections strategies
					if (s + t == 1)
						continue;
					auto& solver = solvers[s + t - 2];
					if (s == 0)
					{
						// naive path tracing
						if (camera_top.hit.geometry->getMaterial()->is_emissive())
						{
							L = camera_top.beta * camera_top.hit.geometry->getMaterial()->Le(camera_top.pNormal(), camera_top.hit.tex_uv, camera_top.omega_o());
							double pdf_light = scene.pdfSampleLe(camera_top.hit.geometry);
							s1_pdf = scene.pdfSampleLi(camera_top.hit.geometry, cameraSubPath[t - 2].hit, camera_top.hit.point);
							double weight = computeWeights(weights, cameraSubPath, LightSubPath, s, t, s1_pdf, pdf_light);
							if (weight == -1)
							{
								solver.addOneTechniqueEstimate(L, s, p);
							}
							else
							{
								L *= weight;
								solver.addEstimate(L, weights, s, p);
							}
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
							resampled_vertex_sa = { LightSubPath.begin(), light_resampled };
						}
						else if (s == 2)
						{
							// BTW this assumes that the connecting loops are in the order t then s
							s1_pdf = scene.pdfSampleLi(LightSubPath[0].hit.geometry, LightSubPath[1].hit, LightSubPath[0].hit.point);
						}

						Vertex& light_top = LightSubPath[s - 1];

						if (camera_top.delta || light_top.delta)
							continue;

						Math::Vector3f dir = camera_top.hit.point - light_top.hit.point;
						const double dist2 = dir.norm2();
						const double dist = sqrt(dist2);
						dir /= dist;

						double G = std::abs(dir * light_top.hit.primitive_normal) / dist2;

						RGBColor camera_connection, light_connection;
						double ni = 1;
						if (t == 1)
						{
							//light tracing
							camera_connection = scene.m_camera.We(-dir);
							p = scene.m_camera.raster(-dir);
							ni = scene.m_camera.resolution;
							if (!scene.m_camera.validRaster(p))
							{
								continue;
							}
						}
						else
						{
							//general case
							camera_connection = camera_top.hit.geometry->getMaterial()->BSDF(camera_top.hit, -dir, false);
							G *= std::abs(camera_top.hit.primitive_normal * dir);
						}

						if (G == 0)
							continue;

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

						L = camera_top.beta * camera_connection * G * light_connection * light_top.beta / ni;
						

						if (scene.visibility(light_top.hit.point, camera_top.hit.point))
						{
							double weight = computeWeights(weights, cameraSubPath, LightSubPath, s, t, s1_pdf);
							if (weight == -1)
							{
								solver.addOneTechniqueEstimate(L, s, p);
							}
							else
							{
								solver.addEstimate(L* weight, weights, s, p);
							}
						}
					}
				}
			}
		}


		////////////////////////////////////////////////////////////////
		//Computes pdf_light of all technique of the bdpt
		// - the last parameter if the probability of sampling the last point on the camera sub path if s == 0 (pure path tracing), else it is not necessary 
		// returns the sum of all the qi (= ni * pi)
		// returns -1 if only the actual tech could produce the sample
		////////////////////////////////////////////////////////////////
		double computeWeights(
			double * weights,
			Path& cameras, Path& lights,
			const int main_s, const int main_t,
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
					return ys->pdf<TransportMode::Radiance, true>(*xt, ys->omega_o());
				else
					return pdf_sampling_point;
			}() };
			assert(xt->rev_pdf >= 0);

			if (ys)
			{
				ys_pdf_rev_sa = { &ys->rev_pdf, xt->pdf<TransportMode::Importance, true>(*ys, xt->omega_o()) };
			}

			if (xt->rev_pdf == 0 && (!ys || ys->rev_pdf == 0))
				return -1;

			if (xtm)
			{
				xtm_pdf_rev_sa = { &xtm->rev_pdf, xt->pdf<TransportMode::Importance, true>(*xtm, xt->dir_to_vertex(ys)) };
			}

			if (ysm)
			{
				ysm_pdf_rev_sa = { &ysm->rev_pdf, ys->pdf<TransportMode::Radiance, true>(*ysm, ys->dir_to_vertex(xt)) };
			}

			double*& ratios = weights;
			
			const double actual_ni = main_t == 1 ? cameras[0].camera().resolution : 1;
			const double actual_main_ri = (main_s == 1 ? s1_pdf / lights[0].fwd_pdf : 1);
			ratios[main_s] = actual_main_ri * actual_ni;

			double sum = actual_main_ri * actual_ni;

			sum += ComputeSumRatioVC(cameras, lights, main_s, main_t, s1_pdf, ratios);
			
			for (int i = 0; i < main_s + main_t; ++i)
			{
				weights[i] = ratios[i] / sum;
			}

			return weights[main_s];
		}

	public:

		OptiMISBDPT(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			BidirectionalBase(sample_per_pixel, width, height)
		{

		}


		virtual void setLen(unsigned int d)override
		{
			Integrator::setLen(d);
		}


		static __forceinline size_t pixelSeed(size_t x, size_t y, size_t width, size_t height, size_t pass)
		{
#ifdef SAMPLER_BIAS
			return (pass) * (width * height);
#else
			return (y * width + x) + (pass) * (width * height);
#endif
		}

		mutable std::vector<MIS::DirectEstimatorImage<Image::IMAGE_ROW_MAJOR>> solvers;
		//mutable std::vector<MIS::BalanceEstimatorImage<Image::IMAGE_ROW_MAJOR>> solvers;

		const bool DEBUG = false;

		void render(Scene const& scene, Visualizer::Visualizer& visu)final override
		{
			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;
			m_frame_buffer.resize(visu.width(), visu.height());
			m_frame_buffer.fill();
			visu.clean();
			ProgressReporter reporter;
			std::cout << "allocating..." << std::endl;
			reporter.start(m_max_len-1);
			solvers.reserve(m_max_len);
			for (int len = 2; len <= m_max_len; ++len)
			{
				int num_tech = len;
				solvers.emplace_back(num_tech, visu.width(), visu.height());
				solvers.back().setSampleForTechnique(len - 1, m_frame_buffer.size());
				reporter.report(len - 1);
			}
			std::vector<std::vector<double>> weights_buffer = Parallel::preAllocate(std::vector<double>(m_max_len));
			reporter.finish();
			std::cout << "Sampling..." << std::endl;
			
			reporter.start(m_sample_per_pixel);
			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				const size_t sample_pass = m_frame_buffer.size();
				OMP_PARALLEL_FOR
					for (long y = 0; y < m_frame_buffer.height(); y++)
					{
						int tid = Parallel::tid();
						double* buffer = weights_buffer[tid].data();

						for (size_t x = 0; x < visu.width(); x++)
						{
							size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), pass);
							Math::Sampler sampler(seed);

							double xp = sampler.generateContinuous<double>();
							double yp = sampler.generateContinuous<double>();

							double v = ((double)y + yp) / visu.height();
							double u = ((double)x + xp) / visu.width();

							computeSample(scene, u, v, sampler, buffer);
						}//pixel x
					}//pixel y
					//the pass has been computed
				reporter.report(pass + 1, -1);
				scene.update_lights_offset(1);
				kbr = visu.update();
				
				if (kbr == Visualizer::Visualizer::KeyboardRequest::done)
				{
					goto __render__end__loop__;
				}
				else if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
				{
					reporter.start(m_max_len - 1);
					std::cout << "Solving..." << std::endl;
					for (int len = 2; len <= m_max_len; ++len)
					{
						int d = len - 2;
						solvers[d].solve(m_frame_buffer, pass+1);
						reporter.report(len - 1);
					}
					showFrame(visu, 1);
					reporter.finish();
					visu.show();
					Image::ImWrite::write(m_frame_buffer);
				}
			}//pass per pixel
			reporter.finish();
			std::cout << "Solving..." << std::endl;
			reporter.start(m_max_len - 1);
			for (int len = 2; len <= m_max_len; ++len)
			{
				int d = len - 2;
				if (DEBUG)
					solvers[d].debug(m_sample_per_pixel, 1, 1, 1, 1);
				solvers[d].solve(m_frame_buffer, m_sample_per_pixel);
				reporter.report(len - 1);
			}
			reporter.finish();
			showFrame(visu, 1);

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

			solvers.clear();
			
		}










		void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult& res)final override
		{
			ProgressReporter reporter;
			m_frame_buffer.resize(width, height);
			m_frame_buffer.fill();
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

							double v = ((double)y + yp) / m_frame_buffer.height();
							double u = ((double)x + xp) / m_frame_buffer.width();

							RGBColor pixel = 0;
							LightVertexStack lvs;

							//computeSample(scene, u, v, sampler, pixel, lvs);

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
				reporter.report(pass + 1, -1);
				scene.update_lights_offset(1);

			}//pass per pixel
			reporter.finish();
			scene.reset_surface_lights();

			//fill the result
			{
				res.time = reporter.time();;
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
			//m_frame_buffer.resize(visu.width(), visu.height());
			//m_frame_buffer.fill();
			visu.clean();
			OMP_PARALLEL_FOR
				for (long y = 0; y < visu.height(); y++)
				{
					int tid = omp_get_thread_num();
					double v = ((double)y + 0.5) / visu.height();
					for (size_t x = 0; x < visu.width(); x++)
					{
						double u = ((double)x + 0.5) / visu.width();

						Ray ray = scene.m_camera.getRay(u, v);
						size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), 0);
						Math::Sampler sampler(seed);

						RGBColor pixel = 0;

						Hit hit;
						if (scene.full_intersection(ray, hit))
						{
							pixel = hit.geometry->getMaterial()->ID_COLOR() * std::abs(ray.direction() * hit.primitive_normal);
						}
						visu.plot(x, y, pixel);
					}//pixel x
				}//pixel y
				//the pass has been computed
			visu.update();
		}

		void debug(Scene const& scene, Visualizer::Visualizer& visu) final override
		{

		}


	};
}