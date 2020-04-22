#pragma once

#include <Integrators/BidirectionalBase.h>
#include <Integrators/PhotonMap.h>
#include <omp.h>
#include <Image/ImWrite.h>

namespace Integrator
{
	class VCM : public BidirectionalBase
	{
	protected:

		enum TransportMode { Importance, Radiance };
		struct Vertex
		{
			enum Type { Camera, Light, Surface } type;

			RGBColor beta;

			//maybe use something lighter than a full hit?
			Hit hit;

			//wether the bsdf of this has a delta distribution, which would make it unconnectable
			bool delta;

			double pdf_fwd, pdf_rev;


			Vertex(Type t, RGBColor b, Hit const& h, bool d) :
				type(t),
				beta(b),
				hit(h),
				delta(d)
			{}

			Vertex() {}

			//////////////////////////////////
			// returns the probability of sampling the direction from this to next, knowing this has been sampled from prev
			// the probability returned is in area density
			// the function should handle most cases (all for now)
			// delta_works: if true, the pdf returned by the delta pdf will be assumed to be 1, 
			// delta works should be true when the connection has beed sampled by the bsdf (like during the random walk), for deterministic connection, it should be false
			//////////////////////////////////
			template <TransportMode MODE, bool DENSITY_AREA = true>
			double pdf(Vertex const& next, const Vertex* prev)const
			{
				double pdf_solid_angle;
				Math::Vector3f to_vertex = next.hit.point - hit.point;
				const double dist2 = to_vertex.norm2();
				const double dist = std::sqrt(dist2);
				to_vertex /= dist;
				if (type == Type::Camera)
				{
					pdf_solid_angle = hit.camera->pdfWeSolidAngle(to_vertex);
				}
				else if (type == Type::Light || prev == nullptr)
				{
					pdf_solid_angle = hit.geometry->getMaterial()->pdfLight(hit, to_vertex);
				}
				else
				{
					if (hit.geometry->getMaterial()->delta())
						return 0;
					else
						pdf_solid_angle = hit.geometry->getMaterial()->pdf(hit, to_vertex, (prev->hit.point - hit.point).normalized(), MODE == TransportMode::Radiance);
				}
				if constexpr (DENSITY_AREA)
				{
					double res = pdf_solid_angle / dist2;

					if (next.type != Type::Camera)
					{
						res *= std::abs(next.hit.primitive_normal * to_vertex);
					}
					return res;
				}
				else
				{
					return pdf_solid_angle;
				}
			}
		};

		template <class T>
		using Stack = StackN<T>;
		using Path = Stack<Vertex>;

		///////////////////////////////////////////////
		//Takes a random walk through the scene (draws a path and record it in res)
		// -Starts at ray
		// -beta is the throughput of the path
		// -pdf is the solid angle probability of sampling ray.dir
		// -type: true >> importance transport aka camera subpath / false >> luminance transport aka light subpath
		///////////////////////////////////////////////
		template <TransportMode MODE>
		__forceinline unsigned int randomWalk(Scene const& scene, Math::Sampler& sampler, Path& res, Ray ray, RGBColor beta, const double pdf, const unsigned int max_len)const
		{
			Hit hit;
			double pdf_solid_angle = pdf;
			Vertex* prev = &res.top();
			double cos_prev = std::abs(prev->hit.primitive_normal * ray.direction());
			int nv = 0;
			for (nv = 0; nv < max_len; ++nv)
			{
				if (scene.full_intersection(ray, hit))
				{
					const double dist = hit.z;
					const double dist2 = dist * dist;

					res.grow();
					Vertex& vertex = res.top();
					vertex = Vertex(Vertex::Type::Surface, beta, hit, hit.geometry->getMaterial()->delta());
					const double cos_vertex = std::abs(hit.primitive_normal * ray.direction());
					vertex.pdf_fwd = pdf_solid_angle * cos_vertex / dist2;

					//sample next direction
					DirectionSample next_dir;
					double pdf_rev;
					vertex.hit.geometry->getMaterial()->sampleBSDF(vertex.hit, next_dir, sampler, MODE == TransportMode::Radiance, &pdf_rev);

					//prev->setPdfReverse<MODE>(vertex.hit.geometry->getMaterial()->pdf(vertex.hit, -ray.direction(), next_dir.direction) * cos_prev / dist2);
					prev->pdf_rev = pdf_rev * cos_prev / dist2;

					//update info for the next loop
					ray = Ray(hit.point, next_dir.direction);
					prev = &vertex;

					cos_prev = std::abs(next_dir.direction * hit.primitive_normal);
					beta = beta * next_dir.bsdf * cos_prev;
					if (beta.isBlack())
						break;
					beta = beta / next_dir.pdf;
					pdf_solid_angle = next_dir.pdf;
				}
				else
					break;
			}
			return nv;
		}

		unsigned int traceCameraSubPath(Scene const& scene, Math::Sampler& sampler, Path& res, Ray const& ray)const
		{
			res.grow();
			Vertex& camera_vertex = res.top();
			camera_vertex.beta = 1;
			camera_vertex.delta = false;
			camera_vertex.type = Vertex::Type::Camera;
			camera_vertex.pdf_fwd = 1;
			camera_vertex.hit.normal = camera_vertex.hit.primitive_normal = ray.direction();
			camera_vertex.hit.point = scene.m_camera.getPosition();
			camera_vertex.hit.camera = &scene.m_camera;
			double pdf_sa = scene.m_camera.pdfWeSolidAngle<true>(ray.direction());
			RGBColor beta = scene.m_camera.We<true>(ray.direction()) / pdf_sa;
			int nv = randomWalk<TransportMode::Importance>(scene, sampler, res, ray, beta, pdf_sa, m_max_len-1);

			camera_vertex.pdf_rev = 0;
			return nv + 1;
		}

		unsigned int traceLightSubPath(Scene const& scene, Math::Sampler& sampler, Path& res)const
		{
			res.grow();
			Vertex& light_vertex = res.top();
			light_vertex.delta = false;
			light_vertex.type = Vertex::Type::Light;

			SurfaceSample sls;
			sampleOneLight(scene, sampler, sls);

			light_vertex.hit.geometry = sls.geo;
			light_vertex.hit.normal = light_vertex.hit.primitive_normal = sls.normal;
			light_vertex.hit.tex_uv = sls.uv;
			light_vertex.hit.point = sls.vector;

			light_vertex.pdf_fwd = sls.pdf;
			light_vertex.beta = 1.0 / sls.pdf;

			//generate a direction
			Math::RandomDirection Le_sampler(&sampler, sls.normal, 1);
			DirectionSample next_dir = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);

			Ray ray(sls.vector, next_dir.direction);

			RGBColor beta = next_dir.bsdf * std::abs(next_dir.direction * sls.normal) / (sls.pdf * next_dir.pdf);

			int nv = randomWalk<TransportMode::Radiance>(scene, sampler, res, ray, beta, (next_dir.pdf), m_max_len-2);

			return 1 + nv;
		}
		

	public:


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



		void buildLightPaths(Scene const& scene, size_t width, size_t height, size_t pass)
		{
			m_light_paths.resize(width, height);
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
				}
			}
			m_light_paths_built = true;
		}









		void computeSample(Scene const& scene, int pix_index, double u, double v, __in Math::Sampler& sampler, __out RGBColor& pixel_res, LightVertexStack& lt_vertices)const
		{
			Path cameraSubPath;

			Ray ray = scene.m_camera.getRay(u, v);

			traceCameraSubPath(scene, sampler, cameraSubPath, ray);
			Path& lightSubPath = m_light_paths[pix_index];

			double s1_pdf;
			for (int t = 1; t <= cameraSubPath.size(); ++t)
			{
				Vertex& camera_top = cameraSubPath[t - 1];

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
							L = camera_top.beta * camera_top.hit.geometry->getMaterial()->Le(camera_top.hit.facing, camera_top.hit.tex_uv);
							double pdf = scene.pdfSamplingLight(camera_top.hit.geometry);
							s1_pdf = scene.pdfSamplingLight(camera_top.hit.geometry, cameraSubPath[t - 2].hit, camera_top.hit.point);
							double weight = MISWeight(cameraSubPath, lightSubPath, scene.m_camera, s, t, s1_pdf, pdf);
							pixel_res += L * weight;
						}
					}
					else
					{
						ScopedAssignment<Vertex> resampled_vertex_sa;
						if (s == 1) {
							double prev_pdf = lightSubPath[0].pdf_fwd;
							SurfaceSample sls;
							sampleOneLight(scene, camera_top.hit, sampler, sls);
							s1_pdf = sls.pdf;
							Vertex light_resampled;
							light_resampled.delta = false;
							light_resampled.type = Vertex::Type::Light;
							light_resampled.hit.geometry = sls.geo;
							light_resampled.hit.normal = light_resampled.hit.primitive_normal = sls.normal;
							light_resampled.hit.tex_uv = sls.uv;
							light_resampled.hit.point = sls.vector;

							light_resampled.pdf_fwd = prev_pdf;
							light_resampled.beta = 1.0 / sls.pdf;
							resampled_vertex_sa = { lightSubPath.begin(), light_resampled };
						}
						else if (s == 2) // BTW this assumes that the connecting loops are in the order t then s
						{
							s1_pdf = scene.pdfSamplingLight(lightSubPath[0].hit.geometry, lightSubPath[1].hit, lightSubPath[0].hit.point);
						}

						Vertex& light_top = lightSubPath[s - 1];
						if (camera_top.delta || light_top.delta)
							continue;

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
							light_connection = light_top.hit.geometry->getMaterial()->Le((light_top.hit.primitive_normal * dir) > 0, light_top.hit.tex_uv);
						}
						else
						{
							//general case
							light_connection = light_top.hit.geometry->getMaterial()->BSDF(light_top.hit, dir, true);
						}

						L = camera_top.beta * camera_connection * G * light_connection * light_top.beta;
						L *= MISWeight(cameraSubPath, lightSubPath, scene.m_camera, s, t, s1_pdf);

						if (!L.isBlack() && visibility(scene, light_top.hit.point, camera_top.hit.point))
						{
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


		////////////////////////////////////////////////////////////////
		//Computes te MIS weights for the bidirectional path tracer
		// - the last parameter if the probability of sampling the last point on the camera sub path if s == 0 (pure path tracing), else it is not necessary 
		////////////////////////////////////////////////////////////////
		double MISWeight(
			Path& cameras, Path& lights,
			const Camera& camera,
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

			xt_pdf_rev_sa = { &xt->pdf_rev,[&]() {
				if (main_t == 1)
					return 0.0;
				else if (ys)
					return ys->pdf<TransportMode::Radiance, true>(*xt, ysm);
				else
					return pdf_sampling_point;
			}() };

			if (ys)
			{
				ys_pdf_rev_sa = { &ys->pdf_rev, xt->pdf<TransportMode::Importance, true>(*ys, xtm) };
			}

			if (xtm)
			{
				xtm_pdf_rev_sa = { &xtm->pdf_rev, xt->pdf<TransportMode::Importance, true>(*xtm, ys) };
			}

			if (ysm)
			{
				ysm_pdf_rev_sa = { &ysm->pdf_rev, ys->pdf<TransportMode::Radiance, true>(*ysm, xt) };
			}

			const double actual_ni = main_t == 1 ? camera.resolution : 1;
			const double actual_main_ri = (main_s == 1 ? s1_pdf / lights[0].pdf_fwd : 1);

			
			double sum = actual_main_ri * actual_ni;


			//expand the camera sub path
			{
				double ri = 1.0;
				for (int s = main_s; s >= 1; --s)
				{
					const Vertex& camera_end = lights[s - 1];
					const Vertex* light_end = s == 1 ? nullptr : &lights[s - 2];
					ri *= camera_end.pdf_rev / (camera_end.pdf_fwd);
					
					const double actual_ri = s != 2 ? ri :
						ri * s1_pdf / light_end->pdf_fwd;

					if (!(camera_end.delta || (light_end && light_end->delta)))
					{
						sum += actual_ri;
					}
				}
			}

			//expand the light subpath
			{
				double ri = 1.0;
				for (int t = main_t; t >= 2; --t)
				{
					const Vertex& light_end = cameras[t - 1];
					const Vertex& camera_end = cameras[t - 2];
					ri *= light_end.pdf_rev / light_end.pdf_fwd;

					double actual_ri = ri;
					if (main_s == 0 && t == main_t)
						actual_ri = ri * s1_pdf / light_end.pdf_rev;

					if (!(light_end.delta || camera_end.delta))
					{
						double ni = (t == 2 ? camera.resolution : 1); // account for the extra samples of the light tracer
						sum += (actual_ri * ni);
					}
				}
			}

			double weight = (actual_main_ri*actual_ni) / sum;
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
			size_t pass = 0;
			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_sample_per_pixel; ++passPerPixelCounter)
			{
				::std::cout << "Pass: " << pass << "/" << Integrator::m_sample_per_pixel << ::std::endl;
				buildLightPaths(scene, visu.width(), visu.height(), passPerPixelCounter);
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
				total += sample_pass;
				++pass;

				showFrame(visu, pass);

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
					Image::ImWrite::write(m_frame_buffer, 1.0 / double(m_sample_per_pixel));
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
			m_light_paths_built = false;
		}










		void render(Scene const& scene, size_t width, size_t height, Auto::RenderResult& res)final override
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
			size_t pass = 0;
			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_sample_per_pixel; ++passPerPixelCounter)
			{

				std::cout << '\r' + progession_bar(pass, m_sample_per_pixel, 100) << std::flush;

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
				total += sample_pass;
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
				res.time = elapsedTime;
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
			const double pixel_area = scene.m_camera.m_down.norm() * scene.m_camera.m_right.norm() / (m_frame_buffer.size());
			const size_t npixels = m_frame_buffer.size();
			const size_t sample_pass = npixels;
			size_t total = 0;
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
			total = npixels;
			showFrame(visu, 1);

			visu.update();
		}

		void debug(Scene const& scene, Visualizer::Visualizer& visu) final override
		{

		}




	};
}