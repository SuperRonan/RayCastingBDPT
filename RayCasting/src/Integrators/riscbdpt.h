#pragma once

#include <Integrators/BidirectionalBase.h>
#include <Image/Image.h>
#include <Image/ImWrite.h>

namespace Integrator
{
	class RISCBDPT: public BidirectionalBase
	{
	protected:

		mutable Image::Image<Path> m_light_paths;
		bool m_light_paths_built = false;

		int m_num_candidates = 4;

		mutable std::vector<std::vector<const Path*>> m_candidates_buffer;

	public:


		void buildLightPaths(Scene const& scene, size_t width, size_t height, size_t pass)
		{
			m_light_paths_built = false;
			m_light_paths.resize(width, height);
			OMP_PARALLEL_FOR
				for (int x = 0; x < width; ++x)
				{
					for (size_t y = 0; y < height; ++y)
					{
						Path& path = m_light_paths(x, y);
						path.reset();
						Math::Sampler sampler(pixelSeed(x, y, width, height, pass));
						traceLightSubPath(scene, sampler, path);
					}
				}
			m_light_paths_built = true;
		}

		double computeWeight(Path& csp, Path& lsp, int S, int T, double s1_pdf)const
		{
			return 0;
		}


		void computeSample(Scene const& scene, int light_index, double u, double v, __in Math::Sampler& sampler, __out RGBColor& pixel_res, LightVertexStack& lt_vertices)const
		{
			Path cameraSubPath;
			const Path& pixel_LightSubPath = m_light_paths[light_index];
			Path LightSubPath;
			std::memcpy(&LightSubPath, &pixel_LightSubPath, sizeof(Path));
			Ray ray = scene.m_camera.getRay(u, v);
			traceCameraSubPath(scene, sampler, cameraSubPath, ray);
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
					RGBColor L = 0;

					if (s + t == 1)
						continue;

					if (s == 0)
					{
						// naive path tracing
						// The light subpath is the one associated with the pixel
						if (camera_top.hit.geometry->getMaterial()->is_emissive())
						{
							const RGBColor Le = camera_top.hit.geometry->getMaterial()->Le(camera_top.pNormal(), camera_top.hit.tex_uv, camera_top.omega_o());
							L = camera_top.beta * Le;
							double pdf = scene.pdfSampleLe(camera_top.hit.geometry);
							
							{
								const RGBColor bsdf = cameraSubPath[t - 2].type == Vertex::Type::Camera ?
									cameraSubPath[t - 2].hit.camera->We(-camera_top.omega_o()) :
									cameraSubPath[t - 2].material()->BSDF(cameraSubPath[t - 2].hit, -camera_top.omega_o(), cameraSubPath[t - 2].omega_o(), false);
								const RGBColor contrib = bsdf * Le * Vertex::G(camera_top, cameraSubPath[t - 2]);
								s1_pdf = scene.pdfRISEstimate(cameraSubPath[t - 2].hit, camera_top.hit, sampler, contrib);
							}

							double weight = 1.0;// VCbalanceWeight(cameraSubPath, LightSubPath, s, t, s1_pdf, pdf);
							pixel_res += L * weight;
						}
					}
					else
					{
						if (t > 1 && s > 1)
						{
							// Use RIS to select the light sub path

						}
						else
						{
							// The light subpath is the one associated with the pixel

							Vertex& light_top = LightSubPath[s - 1];
							if (camera_top.delta || light_top.delta)
								continue;

							ScopedAssignment<Vertex> resampled_vertex_sa;
							if (s == 1)
							{
								// simply resample the vertex on the light source
								SurfaceSample sls;
								// I could get the contribution from the call to this function
								// I could also give the beta of the camera_top as common, but it makes the other call to the pdf a bit complex for now
								scene.sampleLiRIS(sampler, sls, camera_top.hit);
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
							
							Math::Vector3f y_to_x = camera_top.hit.point - light_top.hit.point;
							const double dist2 = y_to_x.norm2();
							const double dist = sqrt(dist2);
							y_to_x /= dist;

							double G = std::abs(y_to_x * light_top.hit.primitive_normal) / dist2;
						}
					}
				}
			}
		}








		RISCBDPT(int sample_per_pixel, int width, int height) :
			BidirectionalBase(sample_per_pixel, width, height),
			m_light_paths(width, height)
		{}

		void render(Scene const& scene, Visualizer::Visualizer& visu)final override
		{
			Visualizer::Visualizer::KeyboardRequest kbr = Visualizer::Visualizer::KeyboardRequest::none;

			m_frame_buffer.resize(visu.width(), visu.height());
			m_frame_buffer.fill();
			ProgressReporter reporter;
			reporter.start(m_sample_per_pixel);
			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				buildLightPaths(scene, visu.width(), visu.height(), pass);
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
				showFrame(visu, pass + 1);
				reporter.report(pass + 1, -1);

				scene.update_lights_offset(1);
				kbr = visu.update();

				if (kbr == Visualizer::Visualizer::KeyboardRequest::done)
				{
					goto __render__end__loop__;
				}
				else if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
				{
					Image::ImWrite::write(m_frame_buffer, 1.0 / double(pass + 1));
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
						}//pixel x
					}//pixel y
					//the pass has been computed
				scene.update_lights_offset(1);
				reporter.report(pass + 1, pass - 1);
			}//pass per pixel

			scene.reset_surface_lights();
			reporter.finish();
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
			{
				buildLightPaths(scene, visu.width(), visu.height(), 0);
			}
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