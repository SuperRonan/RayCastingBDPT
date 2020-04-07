#pragma once

#include <Integrators/BidirectionalBase.h>
#include <System/ScopedAssignment.h>
#include <armadillo>
#include <atomic>
#include <Integrators/optimalmissolver.h>

namespace Integrator
{
	class OptiMISBDPT : public BidirectionalBase
	{
	protected:

		using CameraType = Camera;

		unsigned int max_camera_len, max_light_len;

		enum TransportMode { Importance, Radiance };

		struct Vertex
		{
			enum Type { Camera, Light, Surface } type;

			// in this integrator, beta is the raw throughput (not the estimated)
			RGBColor beta;

			//maybe use something lighter than a full hit?
			Hit hit;

			//wether the bsdf of this has a delta distribution, which would make it unconnectable
			bool delta;

		protected:
			double pdf_importance, pdf_radiance;


		public:

			Vertex(Type t, RGBColor b, Hit const& h, bool d) :
				type(t),
				beta(b),
				hit(h),
				delta(d)
			{

			}


			Vertex()
			{}

			template <TransportMode MODE>
			double pdfForward()const
			{
				if constexpr (MODE == TransportMode::Importance)
					return pdf_importance;
				else
					return pdf_radiance;
			}

			template <TransportMode MODE>
			double pdfReverse()const
			{
				if constexpr (MODE == TransportMode::Importance)
					return pdf_radiance;
				else
					return pdf_importance;
			}

			template <TransportMode MODE>
			double& pdfForward()
			{
				if constexpr (MODE == TransportMode::Importance)
					return pdf_importance;
				else
					return pdf_radiance;
			}

			template <TransportMode MODE>
			double& pdfReverse()
			{
				if constexpr (MODE == TransportMode::Importance)
					return pdf_radiance;
				else
					return pdf_importance;
			}

			template <TransportMode MODE>
			void setPdfForward(double pdf)
			{
				if constexpr (MODE == TransportMode::Importance)
					pdf_importance = pdf;
				else
					pdf_radiance = pdf;
			}

			template <TransportMode MODE>
			void setPdfReverse(double pdf)
			{
				if constexpr (MODE == TransportMode::Importance)
					pdf_radiance = pdf;
				else
					pdf_importance = pdf;
			}

			//////////////////////////////////
			// returns the probability of sampling the direction from this to next, knowing this has been sampled from prev
			// the probability returned is in area density
			// the function should handle most cases (all for now)
			// delta_works: if true, the pdf returned by the delta pdf will be assumed to be 1, 
			// delta works should be true when the connection has beed sampled by the bsdf (like during the random walk), for deterministic connection, it should be false
			//////////////////////////////////
			template <bool DENSITY_AREA = true>
			double pdf(Geometry::Camera const& camera, Vertex const& next, const Vertex* prev, bool delta_works)const
			{
				double pdf_solid_angle;
				Math::Vector3f to_vertex = next.hit.point - hit.point;
				const double dist2 = to_vertex.norm2();
				const double dist = std::sqrt(dist2);
				to_vertex /= dist;
				if (type == Type::Camera)
				{
					pdf_solid_angle = camera.pdfWeSolidAngle(to_vertex);
				}
				else if (type == Type::Light)
				{
					pdf_solid_angle = (hit.primitive_normal * to_vertex) / Math::pi;
					if (pdf_solid_angle < 0)
						pdf_solid_angle = 0;
				}
				else if (prev == nullptr)
				{
					Math::Vector3f normal = hit.primitive_normal * (hit.facing ? 1 : -1);
					pdf_solid_angle = (normal * to_vertex) / Math::pi;
					if (pdf_solid_angle < 0)
						pdf_solid_angle = 0;
				}
				else
				{
					if (hit.geometry->getMaterial()->delta() && delta_works)
						pdf_solid_angle = 1;
					else
						pdf_solid_angle = hit.geometry->getMaterial()->pdf(hit, to_vertex, (prev->hit.point - hit.point).normalized());
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

		using VertexStack = StackN<Vertex>;



		///////////////////////////////////////////////
		//Takes a random walk through the scene (draws a path and record it in res)
		// -Starts at ray
		// -beta is the throughput of the path
		// -pdf is the solid angle probability of sampling ray.dir
		// -type: true >> importance transport aka camera subpath / false >> luminance transport aka light subpath
		///////////////////////////////////////////////
		template <TransportMode MODE>
		__forceinline unsigned int randomWalk(Scene const& scene, Math::Sampler& sampler, VertexStack& res, Ray ray, RGBColor beta, const double pdf, const unsigned int max_len)const
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
					const double cos_vertex = std::abs(hit.primitive_normal * ray.direction());
					vertex = Vertex(Vertex::Type::Surface, beta * cos_vertex / dist2, hit, hit.geometry->getMaterial()->delta());
					vertex.setPdfForward<MODE>(pdf_solid_angle * cos_vertex / dist2);

					//sample next direction

					DirectionSample next_dir;
					vertex.hit.geometry->getMaterial()->sampleBSDF(vertex.hit, next_dir, sampler);

					//prev->setPdfReverse<MODE>(vertex.hit.geometry->getMaterial()->pdf(vertex.hit, -ray.direction(), next_dir.direction) * cos_prev / dist2);
					prev->pdfReverse<MODE>() = next_dir.pdf * cos_prev / dist2;

					//update info for the next loop
					ray = Ray(hit.point, next_dir.direction);
					prev = &vertex;

					cos_prev = std::abs(next_dir.direction * hit.primitive_normal);
					beta = prev->beta * next_dir.bsdf * cos_prev;
					if (beta.isBlack())
						break;
					//beta = beta / next_dir.pdf;
					pdf_solid_angle = next_dir.pdf;
				}
				else
					break;
			}
			return nv;
		}

		unsigned int traceCameraSubPath(Scene const& scene, Math::Sampler& sampler, VertexStack& res, Ray const& ray)const
		{
			res.grow();
			Vertex& camera_vertex = res.top();
			camera_vertex.beta = 1;
			camera_vertex.delta = false;
			camera_vertex.type = Vertex::Type::Camera;
			camera_vertex.setPdfForward<TransportMode::Importance>(1);
			camera_vertex.hit.normal = camera_vertex.hit.primitive_normal = ray.direction();
			camera_vertex.hit.point = scene.m_camera.getPosition();

			int nv = randomWalk<TransportMode::Importance>(scene, sampler, res, ray, scene.m_camera.We<true>(ray.direction()), scene.m_camera.pdfWeSolidAngle<true>(ray.direction()), max_camera_len - 1);

			camera_vertex.setPdfReverse<TransportMode::Importance>(0);
			return nv + 1;
		}

		unsigned int traceLightSubPath(Scene const& scene, Math::Sampler& sampler, VertexStack& res)const
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

			light_vertex.setPdfForward<TransportMode::Radiance>(sls.pdf);
			light_vertex.beta = 1.0;// / sls.pdf;

			//generate a direction
			Math::RandomDirection Le_sampler(&sampler, sls.normal, 1);
			DirectionSample next_dir = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);

			Ray ray(sls.vector, next_dir.direction);

			RGBColor beta = next_dir.bsdf * std::abs(next_dir.direction * sls.normal);// / (sls.pdf * next_dir.pdf);

			int nv = randomWalk<TransportMode::Radiance>(scene, sampler, res, ray, beta, (next_dir.pdf), max_light_len-1);

			return 1 + nv;
		}


		template <class Solver>
		void computeSample(Scene const& scene, double u, double v, __in Math::Sampler& sampler, std::vector<Solver> & solvers, double * pdfs)const
		{
			VertexStack cameraSubPath, LightSubPath;

			Ray ray = scene.m_camera.getRay(u, v);

			traceCameraSubPath(scene, sampler, cameraSubPath, ray);
			traceLightSubPath(scene, sampler, LightSubPath);
			double Pt = useLT ? 1 : cameraSubPath[0].pdfForward<TransportMode::Importance>();
			int t_min = (useLT ? 1 : 2);
			double s1_pdf;
			for (int t = t_min; t <= max_camera_len; ++t)
			{
				Vertex* camera_top;
				
				if (t <= cameraSubPath.size())
				{
					camera_top = cameraSubPath.begin() + (t - 1);
					Pt *= camera_top->pdfForward<TransportMode::Importance>();
				}
				double Ps = 1;
				for (int s = 0; s <= max_light_len; ++s)
				{
					if (s + t > m_max_len)
					{
						break;
					}
					Math::Vector2f p = { u, v };
					const int depth = s + t - 2;
					bool zero = false;
					double sumqi = 0;
					//the contribution (not estimated)
					RGBColor L = 0;
					//handle unsampled samples
					if (s > LightSubPath.size() || t > cameraSubPath.size())
					{
						zero = true;
						if (t == 1)
						{
							continue;
						}
					}
					if(!zero)
					{

						if (s + t > m_max_len)
						{
							break;
						}

						// special cases of connections strategies
						if (s + t == 1)
							continue;
						if (s == 0)
						{
							// naive path tracing
							if (camera_top->hit.geometry->getMaterial()->is_emissive())
							{
								L = camera_top->beta * camera_top->hit.geometry->getMaterial()->Le(camera_top->hit.facing, camera_top->hit.tex_uv);
								double pdf = scene.pdfSamplingLight(camera_top->hit.geometry);
								s1_pdf = scene.pdfSamplingLight(camera_top->hit.geometry, cameraSubPath[t - 2].hit, camera_top->hit.point);
								sumqi = computepdfs(pdfs, cameraSubPath, LightSubPath, scene.m_camera, s, t, Pt, Ps, s1_pdf, pdf);
							}
							else
							{
								zero = true;
							}
						}
						else
						{
							ScopedAssignment<Vertex> resampled_vertex_sa;
							if (s == 1) {
								double prev_pdf = LightSubPath[0].pdfForward<TransportMode::Radiance>();
								SurfaceSample sls;
								sampleOneLight(scene, camera_top->hit, sampler, sls);
								Vertex light_resampled;
								light_resampled.delta = false;
								light_resampled.type = Vertex::Type::Light;
								light_resampled.hit.geometry = sls.geo;
								light_resampled.hit.normal = light_resampled.hit.primitive_normal = sls.normal;
								light_resampled.hit.tex_uv = sls.uv;
								light_resampled.hit.point = sls.vector;

								s1_pdf = sls.pdf;
								light_resampled.setPdfForward<TransportMode::Radiance>(prev_pdf);
								light_resampled.beta = 1.0;
								resampled_vertex_sa = { LightSubPath.begin(), light_resampled };
							}
							else if(s == 2)
							{
								// BTW this assumes that the connecting loops are in the order t then s
								s1_pdf = scene.pdfSamplingLight(LightSubPath[0].hit.geometry, LightSubPath[1].hit, LightSubPath[0].hit.point); 
							}

							Vertex * light_top = LightSubPath.begin() + s - 1;
							Ps *= light_top->pdfForward<Radiance>();
							if (camera_top->delta || light_top->delta)
							{
								zero = true;
							}

							Math::Vector3f dir = camera_top->hit.point - light_top->hit.point;
							const double dist2 = dir.norm2();
							const double dist = sqrt(dist2);
							dir /= dist;

							double G = std::abs(dir * light_top->hit.primitive_normal) / dist2;

							RGBColor camera_connection, light_connection;

							if (t == 1)
							{
								//light tracing
								camera_connection = scene.m_camera.We(-dir);
								p = scene.m_camera.raster(-dir);
								if (!scene.m_camera.validRaster(p))
								{
									zero = true;
									p = { -1, -1 };
									continue;
								}
							}
							else
							{
								//general case
								camera_connection = camera_top->hit.geometry->getMaterial()->BSDF(camera_top->hit, -dir, false);
								G *= std::abs(camera_top->hit.primitive_normal * dir);
							}

							if (s == 1)
							{
								//direct lighting estimation
								light_connection = light_top->hit.geometry->getMaterial()->Le((light_top->hit.primitive_normal * dir) > 0, light_top->hit.tex_uv);
							}
							else
							{
								//general case
								light_connection = light_top->hit.geometry->getMaterial()->BSDF(light_top->hit, dir, true);
							}

							if (G == 0)
							{
								zero = true;
							}
							else
							{
								L = camera_top->beta * camera_connection * G * light_connection * light_top->beta;
								sumqi = computepdfs(pdfs, cameraSubPath, LightSubPath, scene.m_camera, s, t, Pt, Ps, s1_pdf);
							}


							if (!zero)
							{
								//test the visibility
								if (!visibility(scene, light_top->hit.point, camera_top->hit.point))
								{
									zero = true;
								}
							}
						}
					}

					Solver& solver = solvers[depth];
					if (zero)
					{
						solver.AddZeroEstimate(p, s);
					}
					else
					{
						solver.AddEstimate(p, L, pdfs, sumqi, s);
					}
					int _ = 0;
				} // for s
			} // for t
		}


		////////////////////////////////////////////////////////////////
		//Computes pdf of all technique of the bdpt
		// - the last parameter if the probability of sampling the last point on the camera sub path if s == 0 (pure path tracing), else it is not necessary 
		// returns the sum of all the qi (= ni * pi)
		////////////////////////////////////////////////////////////////
		double computepdfs(
			double * pdfs,
			VertexStack& cameras, VertexStack& lights,
			const Camera& camera,
			const int main_s, const int main_t,
			const double Pt, const double Ps,
			double s1_pdf, double pdf_sampling_point = -1)const
		{
			//return 1.0 / double(main_t + main_s);
			Vertex* xt = cameras.begin() + (main_t - 1);
			Vertex* ys = main_s == 0 ? nullptr : lights.begin() + (main_s - 1);
			Vertex* xtm = main_t == 1 ? nullptr : cameras.begin() + (main_t - 2);
			Vertex* ysm = main_s < 2 ? nullptr : lights.begin() + (main_s - 2);

			ScopedAssignment<double> xt_pdf_rev_sa;
			ScopedAssignment<double> ys_pdf_rev_sa;
			ScopedAssignment<double> xtm_pdf_rev_sa;
			ScopedAssignment<double> ysm_pdf_rev_sa;

			xt_pdf_rev_sa = { &xt->pdfReverse<TransportMode::Importance>(),[&]() {
				if (main_t == 1)
					return 0.0;
				if (ys)
					return ys->pdf<true>(camera, *xt, ysm, false);
				else
					return pdf_sampling_point;
			}() };

			if (ys)
			{
				ys_pdf_rev_sa = { &ys->pdfReverse<TransportMode::Radiance>(), xt->pdf<true>(camera, *ys, xtm, false) };
			}

			if (xtm)
			{
				xtm_pdf_rev_sa = { &xtm->pdfReverse<TransportMode::Importance>(), xt->pdf<true>(camera, *xtm, ys, false) };
			}

			if (ysm)
			{
				ysm_pdf_rev_sa = { &ysm->pdfReverse<TransportMode::Radiance>(), ys->pdf<true>(camera, *ysm, xt, false) };
			}

			const double main_p = Ps * Pt;
			const double actual_main_p = (main_s == 1 ? Pt * s1_pdf : main_p);
			pdfs[main_s] = actual_main_p;

			double sum_qi = actual_main_p * (main_t == 1 ? camera.resolution : 1);


			//expand the camera sub path
			{
				double Ph = main_p;
				for (int s = main_s; s >= 1; --s)
				{
					const Vertex& camera_end = lights[s - 1];
					const Vertex* light_end = s == 1 ? nullptr : &lights[s - 2];
					Ph *= camera_end.pdfReverse<TransportMode::Radiance>() / camera_end.pdfForward<TransportMode::Radiance>();

					const double actual_Ph = s == 2 ? 
						Ph * s1_pdf / light_end->pdfForward<TransportMode::Radiance>() : Ph;

					if (!(camera_end.delta || (light_end && light_end->delta)))
					{
						sum_qi += actual_Ph;
						pdfs[s - 1] = actual_Ph;
					}
					else
					{
						pdfs[s - 1] = 0;
					}
				}
			}

			//expand the light subpath
			{
				double Ph = main_p;
				int s = main_s+1;
				int t_min = (useLT ? 2 : 3);
				for (int t = main_t; t >= t_min; --t)
				{
					const Vertex& light_end = cameras[t - 1];
					const Vertex& camera_end = cameras[t - 2];
					Ph *= light_end.pdfReverse<TransportMode::Importance>() / light_end.pdfForward<TransportMode::Importance>();

					double actual_Ph = Ph;
					if (main_s == 0 && t == main_t)
						actual_Ph = Pt * s1_pdf / light_end.pdfForward<TransportMode::Importance>();

					if (!(light_end.delta || camera_end.delta))
					{
						sum_qi += (actual_Ph * (t == 2 ? camera.resolution : 1));
						pdfs[s] = actual_Ph;
					}
					else
					{
						pdfs[s] = 0;
					}
					++s;
				}
			}

			return sum_qi;
		}

	public:

		const bool useLT = true;

		OptiMISBDPT(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			BidirectionalBase(sample_per_pixel, width, height),
			max_camera_len(1),
			max_light_len(1)
		{

		}


		virtual void setLen(unsigned int d)override
		{
			Integrator::setLen(d);
			max_camera_len = d;
			max_light_len = d-1;
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


			std::vector<OptimalSolverImage> solvers;
			//std::vector<OptimalSolverImagePBRT> solvers;
			//std::vector<BalanceSolverImage> solvers;
			
			solvers.reserve(m_max_len-1);
			for (int len = 2; len <= m_max_len; ++len)
			{
				int num_tech = useLT ? len : len - 1;
				solvers.emplace_back(num_tech, visu.width(), visu.height(), m_sample_per_pixel, useLT);
			}
			std::cout << omp_get_num_threads() << " threads" << std::endl;
			std::vector<std::vector<double>> pdfs_buffers(omp_get_num_threads()+16, std::vector<double>(m_max_len, 0.0));


			visu.clean();
			const double pixel_area = scene.m_camera.m_down.norm() * scene.m_camera.m_right.norm() / (m_frame_buffer.size());
			const size_t npixels = m_frame_buffer.size();
			const size_t sample_pass = npixels;
			size_t total = 0;
			size_t pass = 0;
			for (size_t pass = 0; pass < m_sample_per_pixel; ++pass)
			{
				::std::cout << "Pass: " << pass << "/" << Integrator::m_sample_per_pixel << ::std::endl;

				const size_t sample_pass = m_frame_buffer.size();
				OMP_PARALLEL_FOR
					for (long y = 0; y < m_frame_buffer.height(); y++)
					{
						int tid = omp_get_thread_num();
						double* pdfs_buffer = pdfs_buffers[tid].data();

						for (size_t x = 0; x < visu.width(); x++)
						{
							size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), pass);
							Math::Sampler sampler(seed);

							double xp = sampler.generateContinuous<double>();
							double yp = sampler.generateContinuous<double>();

							double v = ((double)y + yp) / visu.height();
							double u = ((double)x + xp) / visu.width();

							computeSample(scene, u, v, sampler, solvers, pdfs_buffer);

							
						}//pixel x
					}//pixel y

					//the pass has been computed
				total += sample_pass;

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
					std::cout << "Solving..." << std::endl;
					for (int len = 2; len <= m_max_len; ++len)
					{
						int d = len - 2;
						std::cout << d << " / " << m_max_len - 2 << std::endl;
						solvers[d].DevelopFilm(&m_frame_buffer, m_sample_per_pixel);
					}
					showFrame(visu, 1);
					std::cout << "Solved!" << std::endl;
					visu.show();
					Image::ImWrite::write(m_frame_buffer);
				}
			}//pass per pixel

			std::cout << "Solving..." << std::endl;
			for (int len = 2; len <= m_max_len; ++len)
			{
				int d = len - 2;
				std::cout << d << " / " << m_max_len-2 << std::endl;
				solvers[d].DevelopFilm(&m_frame_buffer, m_sample_per_pixel);
			}
			std::cout << "Solved!" << std::endl;
			showFrame(visu, 1);

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
			//m_frame_buffer.resize(visu.width(), visu.height());
			//m_frame_buffer.fill();
			visu.clean();
			const size_t npixels = m_frame_buffer.size();
			size_t total = 0;
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