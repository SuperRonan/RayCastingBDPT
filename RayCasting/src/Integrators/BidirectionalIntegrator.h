#pragma once

#include <Integrators/BidirectionalBase.h>
#include <System/ScopedAssignment.h>

namespace Integrator
{
	class BidirectionalIntegrator: public BidirectionalBase
	{
	protected:

		using CameraType = Camera;

		unsigned int max_camera_depth, max_light_depth;

		enum TransportMode { Importance, Radiance };

		struct Vertex
		{
			enum Type {Camera, Light, Surface} type;

			RGBColor beta;
			
			//maybe use something lighter than a full hit?
			Hit hit;

			//wether the bsdf of this has a delta distribution, which would make it unconnectable
			bool delta;

			//Maybe the same concept of pdf fwd and rev as PBRT is actually better for the engine
			double pdf_importance, pdf_radiance; // proba d'avoir été échant par la cam ou la lumière

			//I am too lazy to do a union of hit for the different types of vertices
			const CameraType* camera=nullptr;

			Vertex(Type t, RGBColor b, Hit const& h, bool d):
				type(t),
				beta(b),
				hit(h),
				delta(d)
			{
				
			}


			Vertex()
			{}

			//////////////////////////////////
			// returns the probability of sampling the direction from this to next, knowing this has been sampled from prev
			// the probability returned is in area density
			// the function should handle most cases (all for now)
			// delta_works: if true, the pdf returned by the delta pdf will be assumed to be 1, 
			// delta works should be true when the connection has beed sampled by the bsdf (like during the random walk), for deterministic connection, it should be false
			//////////////////////////////////
			template <bool DENSITY_AREA=true>
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
				else
				{
					pdf_solid_angle = hit.geometry->getMaterial()->pdf(hit, to_vertex);
				}
				if constexpr (DENSITY_AREA)
				{
					double res = pdf_solid_angle / dist2;
					
					if (type != Type::Camera)
					{
						res *= std::abs(hit.primitive_normal * to_vertex);
					}
					
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
		__forceinline unsigned int randomWalk(Scene const& scene, Math::Sampler& sampler, VertexStack& res, Ray ray, RGBColor beta, const double pdf, const unsigned int max_depth)const
		{
			Hit hit;
			Vertex* prev = &res.top(); //used for reverse pdf
			double pdf_solid_angle_prev = pdf;
			double cos_prev = std::abs(prev->hit.primitive_normal * ray.direction());
			int nv = 0;
			for (nv = 0; nv < max_depth; ++nv)
			{
				if (scene.full_intersection(ray, hit))
				{
					const double dist2 = hit.z * hit.z;
					Vertex current_vertex(Vertex::Type::Surface, beta, hit, hit.geometry->getMaterial()->delta());
					res.push(current_vertex);
					const double cos_vertex = std::abs(ray.direction() * hit.primitive_normal);
					const double conversion = cos_vertex / dist2;
					if (MODE == TransportMode::Radiance)
						res.end()->pdf_radiance = pdf_solid_angle_prev * conversion;
					else
						res.end()->pdf_importance = pdf_solid_angle_prev * conversion;

					//sample next dir
					DirectionSample next_dir;
					hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, next_dir, sampler);
					prev = &current_vertex;

					cos_prev = std::abs(next_dir.direction * hit.primitive_normal);
					beta = beta * next_dir.bsdf * cos_prev / next_dir.pdf;
					if (beta.isBlack())
						break;
					ray = Ray(hit.point, next_dir.direction);
					pdf_solid_angle_prev *= next_dir.pdf;
				}
				else
					break;
			}
			return nv;
		}

		__forceinline void computeReverseProbabilities(VertexStack& sub_path, TransportMode const& mode)const
		{
			// TODO
		}

		unsigned int traceCameraSubPath(Scene const& scene, Math::Sampler& sampler, VertexStack& res, Ray const& ray)const
		{
			Hit hit;
			hit.normal = hit.primitive_normal = ray.direction();
			hit.point = scene.m_camera.getPosition();

			Vertex cam_vertex(Vertex::Type::Camera, 1, hit, false);
			res.push(cam_vertex);
			cam_vertex.pdf_importance = 1;

			return randomWalk<TransportMode::Importance>(scene, sampler, res, ray, scene.m_camera.We(ray.direction()) / scene.m_camera.pdfWeSolidAngle(ray.direction()), scene.m_camera.pdfWeSolidAngle(ray.direction()), max_camera_depth)+ 1;
		}

		unsigned int traceLightSubPath(Scene const& scene, Math::Sampler& sampler, VertexStack& res)const
		{
			SurfaceLightSample sls;
			sampleLight(scene, sls, sampler);


			Hit hit;
			hit.point = sls.vector;
			hit.normal = hit.primitive_normal = sls.normal;
			hit.tex_uv = sls.uv;
			hit.geometry = sls.geo;

			//connectVertexToCamera(scene, beta, hit, lvs, true);
			res.push(Vertex(Vertex::Type::Light, 1, hit, hit.geometry->getMaterial()->delta())); // init light
			
			res.begin()->pdf_radiance = sls.pdf; // area density

			//generate a direction
			Math::RandomDirection Le_sampler(&sampler, sls.normal, 1);
			DirectionSample dir_from_light = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);

			double cosl = std::abs(sls.normal * dir_from_light.direction);
			double dist2 = hit.z * hit.z;
			//double conversion = cosl / dist2;
			//RGBColor beta = sls.geo->getMaterial()->Le(true, sls.uv) / sls.pdf;

			//beta = beta * cosl / dir_from_light.pdf;

			RGBColor beta = dir_from_light.bsdf * std::abs(dir_from_light.direction * sls.normal) / (sls.pdf * dir_from_light.pdf);


			Ray ray(sls.vector, dir_from_light.direction);
			if (beta.isBlack())
				return 0;

			return randomWalk<TransportMode::Radiance>(scene, sampler, res, ray, beta, sls.pdf * dir_from_light.pdf, max_light_depth) + 1;
		}
		


		void computeSample(Scene const& scene, double u, double v, __in Math::Sampler& sampler, __out RGBColor& pixel_res, LightVertexStack& lt_vertices)const
		{
			VertexStack cameraSubPath, LightSubPath;

			Ray ray = scene.m_camera.getRay(u, v);

			traceCameraSubPath(scene, sampler, cameraSubPath, ray);
			traceLightSubPath(scene, sampler, LightSubPath);
			double Pt = 1;
			for (int t = 1; t <= cameraSubPath.size(); ++t)
			{
				Vertex& camera_top = cameraSubPath[t - 1];
				Pt *= camera_top.pdf_importance;
				double Ps = 1;
				for (int s = 0; s <= LightSubPath.size(); ++s)
				{
					if (s + t > m_max_depth + 2)
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
							double pdf = LightSubPath[0].pdf_radiance;
							double weight = MISWeight(cameraSubPath, LightSubPath, scene.m_camera, s, t, Pt, Ps, scene.m_camera.resolution, pdf);
							pixel_res += L * weight;
						}
					}
					else
					{
						Vertex light_top = LightSubPath[s - 1];
						Ps *= light_top.pdf_radiance;
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
							camera_connection = camera_top.hit.geometry->getMaterial()->BSDF(camera_top.hit, -dir);
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
							light_connection = light_top.hit.geometry->getMaterial()->BSDF(light_top.hit, dir);
						}

						L = camera_top.beta * camera_connection * G * light_connection * light_top.beta;
						L *= MISWeight(cameraSubPath, LightSubPath, scene.m_camera, s, t, Pt, Ps, scene.m_camera.resolution);

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
		//Computes the MIS weights for the bidirectional path tracer
		// - the last parameter if the probability of sampling the last point on the camera sub path if s == 0 (pure path tracing), else it is not necessary 
		////////////////////////////////////////////////////////////////
		double MISWeight(
			VertexStack& cameras, VertexStack& lights,
			const Camera& camera,
			const int main_s, const int main_t,
			const double Pt, const double Ps,
			double resolution, double pdf_sampling_point = -1)const
		{
			// TODO
			return 1.0 / (double)(main_s + main_t);
			return 1;
		}

	public:

		BidirectionalIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			BidirectionalBase(sample_per_pixel, width, height),
			max_camera_depth(1),
			max_light_depth(1)
		{
			
		}

		virtual void setDepth(unsigned int d)override
		{
			Integrator::setDepth(d);
			max_camera_depth = d;
			max_light_depth = d;
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

							computeSample(scene, u, v, sampler, pixel, lvs);

							m_frame_buffer[x][y] += pixel;
							for (LightVertex const& lv : lvs)
							{
								int lx = lv.uv[0] * visu.width();
								int ly = lv.uv[1] * visu.height();
								m_frame_buffer[lx][ly] += lv.light;
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

							computeSample(scene, u, v, sampler, pixel, lvs);

							m_frame_buffer[x][y] += pixel;
							for (LightVertex const& lv : lvs)
							{
								int lx = lv.uv[0] * m_frame_buffer.width();
								int ly = lv.uv[1] * m_frame_buffer.height();
								m_frame_buffer[lx][ly] += lv.light;
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

						computeSample(scene, u, v, sampler, pixel, lvs);

						m_frame_buffer[x][y] += pixel;
						for (LightVertex const& lv : lvs)
						{
							int lx = lv.uv[0] * visu.width();
							int ly = lv.uv[1] * visu.height();
							m_frame_buffer[lx][ly] += lv.light;
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
			//size_t x, y;

			//bool print_zero = false;

			//while (true)
			//{
			//	std::cout << "x= ";
			//	std::cin >> x;
			//	std::cout << "y= ";
			//	std::cin >> y;

			//	std::vector<RGBColor> direct_samples(m_maximum_pass);
			//	std::vector<LightVertexStack> lt_samples(m_maximum_pass);
			//	size_t width = visu.width(), height = visu.height();

			//	///OMP_PARALLEL_FOR
			//	for (size_t passPerPixelCounter = 0; passPerPixelCounter < m_pass_per_pixel; ++passPerPixelCounter)
			//	{
			//		for (size_t sub_x = 0; sub_x < m_sub_pixel_samples; ++sub_x)
			//		{
			//			double xp = subPixel(sub_x);
			//			double u = (x + xp) / double(width);
			//			for (size_t sub_y = 0; sub_y < m_sub_pixel_samples; ++sub_y)
			//			{
			//				double yp = subPixel(sub_y);
			//				double v = (y + yp) / double(height);
			//				size_t pass = passPerPixelCounter * m_sub_pixel_samples * m_sub_pixel_samples + sub_x * m_sub_pixel_samples + sub_y;
			//				size_t seed = pixelSeed(x, y, width, height, pass);
			//				Math::Sampler sampler(seed);

			//				RGBColor & direct_sample = direct_samples[pass];
			//				LightVertexStack & lt_sample = lt_samples[pass];
			//				computeSample(scene, u, v, sampler, direct_sample, lt_sample);
			//				int chunchunmaru = 0;

			//			}
			//		}
			//	}

			//	RGBColor sum = 0;
			//	std::cout << "--------------------------------------------------------------------------" << std::endl;
			//	for (size_t i = 0; i < direct_samples.size(); ++i)
			//	{
			//		RGBColor const& direct_sample = direct_samples[i];
			//		LightVertexStack const& lt_sample = lt_samples[i];
			//		bool printed = false;
			//		if (!(print_zero && direct_sample.isBlack()))
			//		{
			//			std::cout << "direct sample " << i << ":\t" << direct_sample << std::endl;
			//			printed = true;
			//		}

			//		int j = 0;
			//		for (LightVertex const& lv : lt_sample)
			//		{
			//			if (!(print_zero && lv.light.isBlack()))
			//			{
			//				std::cout << "light  sample " << i << ", " << j << " @ " << lv.uv << ":\t" << lv.light << std::endl;
			//				printed = true;
			//			}
			//			++j;
			//		}

			//		if (printed)
			//		{
			//			std::cout << "--------------------------------------------------------------------------" << std::endl;
			//		}
			//	}

			//	//visu.plot(x, y, mean);
			//	visu.update();

			//}
		}
	};
}