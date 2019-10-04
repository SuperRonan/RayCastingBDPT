#pragma once


/*
#include <Integrators/BidirectionalBase.h>
#include <omp.h>
#include <Image/Image.h>
#include <Image/ImWrite.h>
#include <System/ScopedAssignment.h>
#include <Math/OptiMISSolver.h>
#include <Geometry/BoundedStack.h>


namespace Integrator
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Bidirectional path tracer using the optimal weights
	// for now, the light tracer contribution is not considered
	///////////////////////////////////////////////////////////////////////////////////////////////
	class OptiMISBDPT : public BidirectionalBase
	{
	protected:

		unsigned int max_camera_depth, max_light_depth;

		enum TransportMode { Importance, Radiance };

		struct Vertex
		{
			enum Type { Camera, Light, Surface } type;

			RGBColor beta;

			//maybe use something lighter than a full hit?
			Hit hit;

			//wether the bsdf of this has a delta distribution, which would make it unconnectable
			bool delta;

			//Maybe the same concept of pdf fwd and rev as PBRT is actually better for the engine
			double pdf_inportance, pdf_radiance;

			Vertex(Type t, RGBColor b, Hit const& h, bool d) :
				type(t),
				beta(b),
				hit(h),
				delta(d)
			{

			}

			Vertex()
			{}

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
			// returns the probability of sampling the direction from this to next, knowing this has been sampled from prev													   \\
			// the probability returned is in area density																													   //
			// the function should handle most cases (all for now)																											  //
			// delta_works: if true, the pdf returned by the delta pdf will be assumed to be 1,																				 //
			// delta works should be true when the connection has beed sampled by the bsdf (like during the random walk), for deterministic connection, it should be false  //
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			double pdf(Vertex const& next, const Vertex* prev, bool delta_works)const
			{
				double pdf_solid_angle;
				Math::Vector3f wi = next.hit.point - hit.point;
				const double dist2 = wi.norm2();
				assert(dist2 != 0);
				const double dist = sqrt(dist2);
				wi /= dist;
				if (type == Camera)
				{
					assert(prev == nullptr);
					pdf_solid_angle = 1.0;
				}
				else if (type == Light)
				{
					assert(prev == nullptr);
					pdf_solid_angle = pdfLi(next);
				}
				else
				{
					assert(prev != nullptr);
					const Math::Vector3f wo = (hit.point - prev->hit.point).normalized();

					double pdf_dir;
					if (delta && delta_works)
					{
						pdf_dir = 1;
					}
					else
					{
						pdf_dir = hit.geometry->getMaterial()->pdf(hit, wi, wo);
					}
					pdf_solid_angle = pdf_dir;
				}

				//conversion to area measure
				double pdf_area = pdf_solid_angle / dist2;
				{
					pdf_area *= std::abs(next.hit.primitive_normal * wi);
				}
				return pdf_area;
			}

			double pdfLi(Vertex const& next)const
			{
				Math::Vector3f dir = next.hit.point - hit.point;
				return hit.geometry->getMaterial()->pdfLight(hit, dir.normalized());
			}
		};

		using VertexStack = StackN<Vertex>;

		using OptiMISSolver = Math::OptiMISSolver;
		
		//returns the number of techniques for a given length (the number of vertices in the path)
		//TODO add light tracing option
		__forceinline unsigned int numTech(unsigned int len)const
		{
			if (len == 2)	return 1; //direct visibility of the lights -> only path tracing
			return len - 1;
		}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
		//Takes a random walk through the scene (draws a path and record it in res)									  \\
		// -Starts at ray																							  //
		// -beta is the throughput of the path																		 //
		// -pdf is the solid angle probability of sampling ray.dir													//
		// -type: true >> importance transport aka camera subpath / false >> luminance transport aka light subpath //
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline unsigned int randomWalk(Scene const& scene, Math::Sampler& sampler, VertexStack& res, Ray ray, RGBColor beta, const double pdf, const unsigned int max_depth, TransportMode mode)const
		{
			double fwd = pdf, bck = 0;
			Hit hit;
			for (unsigned int depth = 0; depth < max_depth; ++depth)
			{
				if (scene.full_intersection(ray, hit))
				{
					//add the vertex
					double cos_hit = (hit.primitive_normal * hit.to_view);
					double dist2 = (hit.z * hit.z);
					double area_pdf = fwd;
					area_pdf *= cos_hit / dist2;

					beta *= cos_hit / dist2;

					Vertex vertex(Vertex::Type::Surface, beta, hit, hit.geometry->getMaterial()->delta());

					if (mode == TransportMode::Importance)
					{
						vertex.pdf_inportance = area_pdf;
					}
					else // radiance
					{
						vertex.pdf_radiance = area_pdf;
					}

					res.push(vertex);


					//sample next direction
					DirectionSample next_dir;
					hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, next_dir, sampler);

					beta *= next_dir.bsdf * (next_dir.direction * hit.primitive_normal);
					fwd = next_dir.pdf;

					ray = Ray(hit.point, next_dir.direction);

					if (beta.isBlack() || fwd <= 0)
					{
						break;
					}
				}
				else if (mode == TransportMode::Importance)
				{
					//if this is a camera subpath, we keep the ray (for the skybox I guess)
				}
			}
			return max_depth;
		}

		__forceinline void computeReverseProbabilities(VertexStack& sub_path, TransportMode const& mode)const
		{
			//the pdf rev of the two vertices at the end of the sub path will not be computed, because it is impossible
			//i represents the number of vertices in the subpath
			for (int i = sub_path.size() - 2; i > 0; --i)
			{
				Vertex& vertex = sub_path[i - 1];
				Vertex const& sampler = sub_path[i];
				Vertex const& prev = sub_path[i + 1];

				if (mode == TransportMode::Importance && i == 1)
				{
					//probability for real particle tracing, 0 here since we discard this technique.
					vertex.pdf_radiance = 0;
					continue;
				}

				double pdf = sampler.pdf(vertex, &prev, true);

				if (mode == TransportMode::Importance)
				{
					vertex.pdf_radiance = pdf;
				}
				else
				{
					vertex.pdf_inportance = pdf;
				}
			}
		}

		unsigned int traceCameraSubPath(Scene const& scene, Math::Sampler& sampler, VertexStack& res, Ray const& ray, double dir_pdf = 1)const
		{
			//But what if there is no camera!
			Hit camera_stub_hit;
			camera_stub_hit.point = scene.m_camera.m_position;
			camera_stub_hit.primitive_normal = camera_stub_hit.normal = scene.m_camera.m_front;
			Vertex camera_vertex = Vertex(Vertex::Type::Camera, 1, camera_stub_hit, false);
			camera_vertex.pdf_inportance = 1;
			res.push(camera_vertex);
			unsigned int len = 1 + randomWalk(scene, sampler, res, ray, 1, dir_pdf, max_camera_depth, TransportMode::Importance);
			//compute the reverse: luminace pdf
			computeReverseProbabilities(res, TransportMode::Importance);
			return len;
		}

		unsigned int traceLightSubPath(Scene const& scene, Math::Sampler& sampler, VertexStack& res)const
		{
			//I dont care about the skybox for now
			SurfaceLightSample light;
			if (sampleLight(scene, light, sampler))
			{
				//make the fist vertex on the light
				RGBColor Le = light.geo->getMaterial()->Le(true, light.uv);
				Hit light_stub_hit;
				light_stub_hit.point = light.vector;
				light_stub_hit.normal = light_stub_hit.primitive_normal = light.normal;
				light_stub_hit.geometry = light.geo;
				light_stub_hit.tex_uv = light.uv;
				light_stub_hit.facing = true;
				double pdf = light.pdf;

				Vertex vertex(Vertex::Type::Light, Le, light_stub_hit, false);
				vertex.pdf_radiance = pdf;
				res.push(vertex);

				//choose a direction to start the path
				DirectionSample dir = light.geo->getMaterial()->sampleLightDirection(light, sampler);
				Ray ray(light.vector, dir.direction);
				RGBColor beta = Le * (ray.direction() * light_stub_hit.primitive_normal);
				unsigned int len = 1 + randomWalk(scene, sampler, res, ray, beta, dir.pdf, max_light_depth + 1, TransportMode::Radiance);

				computeReverseProbabilities(res, TransportMode::Radiance);

				return len;
			}
			else
			{
				return 0;
			}
		}



		void computeSample(Scene const& scene, double u, double v, __in Math::Sampler& sampler, std::vector<Math::ImageSolver<OptiMISSolver>>& solvers)const
		{
			Ray ray = scene.m_camera.getRay(u, v);

			VertexStack camera_subpath, light_subpath;

			traceCameraSubPath(scene, sampler, camera_subpath, ray, 1);
			traceLightSubPath(scene, sampler, light_subpath);

			//compute the reverse probabilities of the subpaths


#ifdef MIS_SHOW_WEIGHTS
			const RGBColor __light_color = { 0, 1, 0 };
			const RGBColor __camera_color = { 1, 0, 0 };
			auto weightColor = [&__light_color, &__camera_color](double w, unsigned int s, unsigned int t)
			{
				double total = s + t - 1;
				RGBColor res = __light_color * (double(s)) + __camera_color * double(t);
				return res / total * w;
			};
#endif

			//make the connections between the sub paths
			double Ps = 1;
			for (int s = 0; s <= light_subpath.size(); ++s)//s represents the length of the sub path of light
			{
				const Vertex* vl = s == 0 ? nullptr : &light_subpath[s - 1];
				if (s != 0)
				{
					Ps *= vl->pdf_radiance;
				}
				double Pt = 1;
				for (int t = 1; t <= camera_subpath.size(); ++t)//t represents the length of the camera sub path
				{

					//return argument of the inlined function in the loop
					Math::Vector2f pRaster = { u, v };
					RGBColor L = 0;
					double pdf_sampling_point = -1;
					double path_pdf=0;

					const Vertex* vc = &camera_subpath[t - 1];
					Pt *= vc->pdf_inportance;

					path_pdf = Pt * Ps;

#ifdef BDPT_SINGLE_TECHNIQUE
					const std::vector<int> desired_s = { 1 };
					const std::vector<int> desired_t = { 1, 2, 3, 4, 5, 6 };

					if (std::find(desired_s.cbegin(), desired_s.cend(), s) == (desired_s.cend()) ||
						std::find(desired_t.cbegin(), desired_t.cend(), t) == (desired_t.cend()))
					{
						continue;
					}
#endif

					if (t == 1) //light tracing
					{
						if (s <= 1)
						{
							//either s == 0, then only on point on the path, so this is not really a path
							//either s == 1, so direct visibility between the light and the camera. Since this is BDPT, 
							// we also consider pure path tracing, and those would give a noise free contribution, so discard these light tracing short paths
							continue;
						}

						//for now, light tracing is discarded, 
						continue;

						
						
						pRaster = scene.m_camera.screen_position(vl->hit.point);
						if (valid_uv(pRaster) && !vl->delta)
						{
							Math::Vector3f to_camera = scene.m_camera.m_position - vl->hit.point;
							double dist2 = to_camera.norm2();
							double dist = sqrt(dist2);
							to_camera /= dist;
							double cos_camera = std::abs(scene.m_camera.m_front * -to_camera);
							double cos_light = std::abs(to_camera * vl->hit.primitive_normal);
							double G = cos_camera * cos_light / dist2;
							L *= vl->beta * G;

							if (s != 1)
							{
								L *= vl->hit.geometry->getMaterial()->BSDF(vl->hit, to_camera);
							}
							//connect to the camera
							if (!L.isBlack() && (cameraVisibility(scene, vl->hit.point)))
							{
								
							}
							else
							{
								L = 0;
							}
						}
					}
					else if (s != 0) //most general case
					{
						//try to make the connection between both end of the sub paths
						if (!vl->delta && !vc->delta)//if any is a delta, then the bsdf of the connection will be 0 (unless a miracle)
						{
							Math::Vector3f camera_to_light = vl->hit.point - vc->hit.point;
							double dist2 = camera_to_light.norm2();
							double dist = sqrt(dist2);
							camera_to_light /= dist;
							//We need BSDF from both sides, and the geometric term of the connection
							RGBColor bsdf_camera = vc->hit.geometry->getMaterial()->BSDF(vc->hit, camera_to_light);
							RGBColor bsdf_light;
							if (s == 1)//the light vertex is on the light, so the bsdf is 1, only the Geometric term is considered
							{
								bsdf_light = 1;
							}
							else
							{
								bsdf_light = vl->hit.geometry->getMaterial()->BSDF(vl->hit, -camera_to_light);
							}
							const double cos_camera = std::abs(vc->hit.primitive_normal * camera_to_light);
							const double cos_light = std::abs(vl->hit.primitive_normal * camera_to_light);
							const double G = cos_camera * cos_light / dist2;
							L = vl->beta * bsdf_light * G * bsdf_camera * vc->beta;
							if (!L.isBlack() && (visibility(scene, vc->hit.point, vl->hit.point)))
							{

							}
							else
							{
								L = 0;
							}
						}
					}
					else //pure path tracing
					{
						if (vc->hit.geometry->getMaterial()->is_emissive())
						{
							//since this is pure path tracing, we need to compute the probability of sampling the end point of the path
							pdf_sampling_point = scene.pdfSamplingLight(vc->hit.geometry);
							
							RGBColor Le = vc->hit.geometry->getMaterial()->Le(vc->hit.facing, vc->hit.tex_uv);
							 L = vc->beta * Le;
						}
					}
					
					// Now L and pRaster should be computed 

					// Compute the balance weights
					BoundedStack<20, double> weights;
					const unsigned int len = s + t;
					weights.grow(len);
					double sum_pdf = 0;

					MISWeight(camera_subpath, light_subpath, s, t, Ps, Pt, pdf_sampling_point, weights.begin(), sum_pdf);

					solvers[len - 2].addEstimate(pRaster, L, weights.cbegin(), sum_pdf);
					
				}
			}
		}

		////////////////////////////////////////////////////////////|\\\
		// Computes the weight of each technique					   \\
		//															   //
		////////////////////////////////////////////////////////////////
		void MISWeight(VertexStack& cameras, VertexStack& lights, const int this_s, const int this_t, const double Ps, const double Pt, double pdf_sampling_point,
							__out double * weights, __out double & sum_pdf)const
		{
			assert(Ps * Pt > 0);
			if (this_s + this_t == 2)
			{
				assert(numTech(this_s + this_t) == 1);
				if (this_s == 0)
				{
					weights[0] = 1.0;
				}
				else
				{
					assert(false);
					weights[0] = 0;
				}
			}
			assert(this_t > 0);
			Vertex* pt = &cameras[this_t - 1];
			Vertex* qs = this_s == 0 ? nullptr : &lights[this_s - 1];
			Vertex* pt_minus = this_t <= 1 ? nullptr : &cameras[this_t - 2];
			Vertex* qs_minus = this_s <= 1 ? nullptr : &lights[this_s - 2];

			//temporarly update the pdf of the connecting vertices
			ScopedAssignment<double> a1;
			{
				if (qs)
				{
					a1 = ScopedAssignment<double>(&pt->pdf_radiance, qs->pdf(*pt, qs_minus, false));
				}
				else
				{
					assert(pdf_sampling_point > 0);
					a1 = ScopedAssignment<double>(&pt->pdf_radiance, pdf_sampling_point);
				}
			}

			ScopedAssignment<double> a2;
			if (qs)
			{
				a2 = { &qs->pdf_inportance, pt->pdf(*qs, pt_minus, false) };
			}

			ScopedAssignment<double> a3;
			if (pt_minus)
			{
				if (qs)
				{
					a3 = { &pt_minus->pdf_radiance, pt->pdf(*pt_minus, qs, false) };
				}
				else
				{
					a3 = { &pt_minus->pdf_radiance, pt->pdfLi(*pt_minus) };
				}

			}

			ScopedAssignment<double> a4;
			if (qs_minus)
			{
				assert(qs != nullptr);
				a4 = { &qs_minus->pdf_inportance, qs->pdf(*qs_minus, pt, false) };
			}

			double * const pdfs = weights;

			//pdf of the path by the "main" technique
			double this_pdf = Ps * Pt;
			pdfs[this_s] = this_pdf;


			//expand the camera subpath
			{
				const Vertex* camera_end = &cameras[this_t - 1];
				double _ps = Ps, _pt = Pt;
				for (int i = this_s - 1; i >= 0; --i)
				{
					const unsigned int tech_index = i;
					const Vertex& camera_end = lights[i];
					const Vertex* light_end = i > 0 ? &lights[i - 1] : nullptr;
					bool delta_connection = camera_end.delta || (light_end && light_end->delta);

					_pt *= camera_end.pdf_inportance;
					_ps /= camera_end.pdf_radiance;
					double tech_pdf = 0;
					if (!delta_connection)
					{
						double tech_pdf = _pt * _ps;
					}
					sum_pdf += tech_pdf;
					pdfs[tech_index] = tech_pdf;
					
				}
				//here _ps should be close to 1
			}

			//expand the light subpath
			{
				double _ps = Ps, _pt = Pt;
				unsigned int min_t = 1;//0 if use light tracing
				unsigned int tech_index = this_s;
				for (int i = this_t - 1; i > min_t; --i)
				{
					++tech_index;
					const Vertex& light_end = cameras[i];
					const Vertex& camera_end = cameras[i - 1];
					bool delta_connection = camera_end.delta || (light_end.delta);

					_ps *= light_end.pdf_radiance;
					_pt /= light_end.pdf_inportance;

					double tech_pdf = 0;
					if (!delta_connection)
					{
						double tech_pdf = _ps * _pt;
					}
					pdfs[tech_index] = tech_pdf;
					sum_pdf += tech_pdf;
				}
				//here _pt should be arround 1
			}

			// Now all the pdfs are conputed
			assert(sum_pdf > 0);
			for (int i = 0; i < numTech(this_s + this_t); ++i)
			{
				weights[i] = pdfs[i] / sum_pdf;
			}
		}

	public:

		OptiMISBDPT(unsigned int sub_pixel_samples, unsigned int pass_per_pixel, unsigned int width, unsigned int height) :
			BidirectionalBase(sub_pixel_samples, pass_per_pixel, width, height),
			max_camera_depth(1),
			max_light_depth(1)
		{

		}


		virtual void setDepth(unsigned int d)override
		{
			Integrator::setDepth(d);
			max_camera_depth = d / 2;
			max_light_depth = d - max_camera_depth;
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

			std::vector<Math::ImageSolver<OptiMISSolver>> solvers;
			solvers.reserve(m_max_depth);
			for (int i = 0; i < m_max_depth; ++i)
			{
				int lenght = i + 2; 
				solvers.emplace_back(&m_frame_buffer, numTech(lenght));
			}
			

			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_pass_per_pixel; ++passPerPixelCounter)
			{
				for (size_t sub_x = 0; sub_x < m_sub_pixel_samples; ++sub_x)
				{
					double xp = subPixel(sub_x);
					for (size_t sub_y = 0; sub_y < m_sub_pixel_samples; ++sub_y)
					{
						double yp = subPixel(sub_y);
						::std::cout << "Pass: " << pass << "/" << Integrator::m_maximum_pass << ::std::endl;

						const size_t sample_pass = m_frame_buffer.size();
						OMP_PARALLEL_FOR
							for (size_t y = 0; y < m_frame_buffer.height(); y++)
							{
								int tid = omp_get_thread_num();
								double v = ((double)y + yp) / m_frame_buffer.height();
								for (size_t x = 0; x < visu.width(); x++)
								{
									double u = ((double)x + xp) / visu.width();

									Ray ray = scene.m_camera.getRay(u, v);
									size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), pass);
									Math::Sampler sampler(seed);

									computeSample(scene, u, v, sampler, solvers);

								}//pixel x
							}//pixel y
							//the pass has been computed
						total += sample_pass;
						++pass;
						m_frame_buffer.fill();
						for (Math::ImageSolver<OptiMISSolver>& solver : solvers)
						{
							solver.loop();
							solver.solveAll();
						}

						showFrame(visu, total);

						scene.update_lights_offset(1);
						kbr = visu.update();
						QueryPerformanceCounter(&t2);
						elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
						double remainingTime = (elapsedTime / pass) * (Integrator::m_maximum_pass - pass);
						::std::cout << "time: " << elapsedTime << "s. " << ", remaining time: " << remainingTime << "s. " << ", total time: " << elapsedTime + remainingTime << ::std::endl;

						if (kbr == Visualizer::Visualizer::KeyboardRequest::done)
						{
							goto __render__end__loop__;
						}
						else if (kbr == Visualizer::Visualizer::KeyboardRequest::save)
						{
							Image::ImWrite::write(m_frame_buffer);
						}
					}//sub y
				}//sub x
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
					Image::ImWrite::write(m_frame_buffer, 1.0 / (double)total);
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
			for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_pass_per_pixel; ++passPerPixelCounter)
			{
				for (size_t sub_x = 0; sub_x < m_sub_pixel_samples; ++sub_x)
				{

					double xp = subPixel(sub_x);
					for (size_t sub_y = 0; sub_y < m_sub_pixel_samples; ++sub_y)
					{
						double yp = subPixel(sub_y);
						std::cout << '\r' + progession_bar(pass, m_maximum_pass, 100) << std::flush;

						const size_t sample_pass = m_frame_buffer.size();
						OMP_PARALLEL_FOR
							for (size_t y = 0; y < m_frame_buffer.height(); y++)
							{
								int tid = omp_get_thread_num();
								double v = ((double)y + yp) / m_frame_buffer.height();
								for (size_t x = 0; x < width; x++)
								{
									double u = ((double)x + xp) / width;

									Ray ray = scene.m_camera.getRay(u, v);
									size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), pass);
									Math::Sampler sampler(seed);

									RGBColor result=0;
									
									
								}//pixel x
							}//pixel y
							//the pass has been computed
						total += sample_pass;
						++pass;
						scene.update_lights_offset(1);

					}//sub y
				}//sub x
			}//pass per pixel

			std::cout << '\r' + progession_bar(m_pass_per_pixel, m_pass_per_pixel, 100) << std::endl;
			// stop timer
			QueryPerformanceCounter(&t2);
			elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;

			scene.reset_surface_lights();

			//fill the result
			{
				res.time = elapsedTime;
				res.image.resize(width, height);
				OMP_PARALLEL_FOR
					for (size_t i = 0; i < m_frame_buffer.size(); ++i)
					{
						res.image.m_data[i] = m_frame_buffer.m_data[i] / total;
					}
			}
		}




		void fastRender(Scene const& scene, Visualizer::Visualizer& visu)final override
		{
			//for now it is a stub showing the normals
			m_frame_buffer.resize(visu.width(), visu.height());
			m_frame_buffer.fill();
			visu.clean();
			const double pixel_area = scene.m_camera.m_down.norm() * scene.m_camera.m_right.norm() / (m_frame_buffer.size());
			const size_t npixels = m_frame_buffer.size();
			const size_t sample_pass = npixels;
			size_t total = 0;
			OMP_PARALLEL_FOR
				for (size_t y = 0; y < m_frame_buffer.height(); y++)
				{
					int tid = omp_get_thread_num();
					double v = ((double)y + 0.5) / m_frame_buffer.height();
					for (size_t x = 0; x < visu.width(); x++)
					{
						double u = ((double)x + 0.5) / visu.width();

						Ray ray = scene.m_camera.getRay(u, v);
						size_t seed = pixelSeed(x, y, m_frame_buffer.width(), m_frame_buffer.height(), 0);
						Math::Sampler sampler(seed);

						RGBColor result=0;

						Hit hit;
						if (scene.full_intersection(ray, hit))
						{
							for (int s = 0; s < 3; ++s)
							{
								result[s] = (hit.normal[s] * 0.5 + 0.5) * (hit.primitive_uv[s%2]) * (hit.to_view * hit.normal);
							}
						}
						visu.plot(x, y, result);
					}//pixel x
				}//pixel y
				//the pass has been computed
			total = npixels;
			//showFrame(visu, total);

			visu.update();
		}

		void debug(Scene const& scene, Visualizer::Visualizer& visu) final override
		{
			/*
			size_t x, y;

			bool print_zero = false;

			while (true)
			{
				std::cout << "x= ";
				std::cin >> x;
				std::cout << "y= ";
				std::cin >> y;

				std::vector<RGBColor> direct_samples(m_maximum_pass);
				std::vector<LightVertexStack> lt_samples(m_maximum_pass);
				size_t width = visu.width(), height = visu.height();

				///OMP_PARALLEL_FOR
				for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_pass_per_pixel; ++passPerPixelCounter)
				{
					for (size_t sub_x = 0; sub_x < m_sub_pixel_samples; ++sub_x)
					{
						double xp = subPixel(sub_x);
						double u = (x + xp) / double(width);
						for (size_t sub_y = 0; sub_y < m_sub_pixel_samples; ++sub_y)
						{
							double yp = subPixel(sub_y);
							double v = (y + yp) / double(height);
							size_t pass = passPerPixelCounter * m_sub_pixel_samples * m_sub_pixel_samples + sub_x * m_sub_pixel_samples + sub_y;
							size_t seed = pixelSeed(x, y, width, height, pass);
							Math::Sampler sampler(seed);

							RGBColor& direct_sample = direct_samples[pass];
							LightVertexStack& lt_sample = lt_samples[pass];
							computeSample(scene, u, v, sampler, direct_sample, lt_sample);
							int chunchunmaru = 0;

						}
					}
				}

				RGBColor sum = 0;
				std::cout << "--------------------------------------------------------------------------" << std::endl;
				for (size_t i = 0; i < direct_samples.size(); ++i)
				{
					RGBColor const& direct_sample = direct_samples[i];
					LightVertexStack const& lt_sample = lt_samples[i];
					bool printed = false;
					if (!(print_zero && direct_sample.isBlack()))
					{
						std::cout << "direct sample " << i << ":\t" << direct_sample << std::endl;
						printed = true;
					}

					int j = 0;
					for (LightVertex const& lv : lt_sample)
					{
						if (!(print_zero && lv.light.isBlack()))
						{
							std::cout << "light  sample " << i << ", " << j << " @ " << lv.uv << ":\t" << lv.light << std::endl;
							printed = true;
						}
						++j;
					}

					if (printed)
					{
						std::cout << "--------------------------------------------------------------------------" << std::endl;
					}
				}

				//visu.plot(x, y, mean);
				visu.update();

			}
			
		}


	};
}

*/