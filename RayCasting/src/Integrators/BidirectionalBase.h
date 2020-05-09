#pragma once


#include <Integrators\Integrator.h>
#include <System/ScopedAssignment.h>

namespace Integrator
{
	class BidirectionalBase: public Integrator
	{
	protected:
		Image::Image < RGBColor, Image::IMAGE_ROW_MAJOR> m_frame_buffer;

		void resizeFrameBuffer(size_t w, size_t h)
		{
			m_frame_buffer.resize(w, h);
		}

		struct LightVertex {
			RGBColor light;
			Math::Vector2f uv;
		};

		using LightVertexStack = StackN<LightVertex>;


		__forceinline bool cameraVisibility(Scene const& scene, Math::Vector3f const& point)const
		{
			Math::Vector3f dir = scene.m_camera.m_position - point;
			Ray ray(point, dir);
			return !scene.intersectionCloser(ray, dir.norm());
		}


		__forceinline bool visibility(Scene const& scene, Math::Vector3f const& p, Math::Vector3f const& q)const
		{
			Math::Vector3f dir = p - q;
			Ray ray(q, dir);
			//I really don't like it
			return !scene.intersectionCloser(ray, dir.norm() - 0.0001);
		}

		__forceinline bool cameraSeeSkybox(Scene const& scene, Math::Vector3f const& dir)const
		{
			return scene.noIntersection(Ray(scene.m_camera.m_position, dir));
		}

		static void samplePointDisk(Math::Vector3f const& center, double radius, double radius2, Math::Vector3f const& normal, Math::Sampler& sampler, Math::Vector3f& res, double& pdf)
		{
			Math::Vector3f tg = Math::Vector3f(1, 0, 0);
			tg = tg - normal * (normal * tg);
			if (tg.norm() < std::numeric_limits<double>::epsilon() * 10)
			{
				tg = { 0, 1, 0 };
				tg = tg - normal * (normal * tg);
				if (tg.norm() < std::numeric_limits<double>::epsilon() * 10)
				{
					tg = { 0, 0, 1 };
					tg = tg - normal * (normal * tg);
					assert(tg.norm() > std::numeric_limits<double>::epsilon() * 10);
				}
			}
			tg = tg.normalized();
			Math::Vector3f ctg = normal ^ tg;

			double angle = Math::twoPi * sampler.generateContinuous<double>();
			double rho = sqrt(sampler.generateContinuous<double>()) * (radius);

			res = center + (ctg * cos(angle) + tg * sin(angle)) * rho;

			pdf = 1.0 / (Math::pi * radius2);
		}


		void showFrame(Visualizer::Visualizer& visu, size_t total)const
		{
			if (visu.visible())
			{
				OMP_PARALLEL_FOR
					for (long x = 0; x < m_frame_buffer.width(); ++x)
					{
						for (size_t y = 0; y < m_frame_buffer.height(); ++y)
						{
							visu.plot(x, y, m_frame_buffer(x, y) / double(total));
						}
					}
			}
		}



	public:

		BidirectionalBase(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			m_frame_buffer(width, height)
		{
			m_sample_per_pixel = sample_per_pixel;
		}




		enum TransportMode { Importance, Radiance };
		struct Vertex
		{
			enum Type { Camera, Light, Surface } type;

			RGBColor beta;

			//maybe use something lighter than a full hit?
			Hit hit;

			//wether the bsdf of this has a delta distribution, which would make it unconnectable
			bool delta;

			double rev_pdf, fwd_pdf;


			Vertex(Type t, RGBColor b, Hit const& h, bool d) :
				type(t),
				beta(b),
				hit(h),
				delta(d)
			{}

			Vertex() {}

			Math::Vector3f dir_to_vertex(const Vertex* vert = nullptr)const
			{
				if (vert)
					return (vert->hit.point - hit.point).normalized();
				return 0.0;
			}

			Math::Vector3f const& position()const
			{
				return hit.point;
			}

			Math::Vector3f const& omega_o()const
			{
				return hit.to_view;
			}

			const Geometry::Material* material()const
			{
				return hit.geometry->getMaterial();
			}

			const Geometry::GeometryBase* geometry()const
			{
				return hit.geometry;
			}

			const Geometry::Camera& camera()const
			{
				assert(type == Type::Camera);
				return *hit.camera;
			}

			const Math::Vector3f pNormal()const
			{
				return hit.primitive_normal;
			}

			const Math::Vector3f sNormal()const
			{
				return hit.normal;
			}


			//////////////////////////////////
			// returns the probability of sampling the direction from this to next, knowing this has been sampled from prev
			// the probability returned is in area density
			// the function should handle most cases (all for now)
			// delta_works: if true, the pdf returned by the delta pdf will be assumed to be 1, 
			// delta works should be true when the connection has beed sampled by the bsdf (like during the random walk), for deterministic connection, it should be false
			//////////////////////////////////
			template <TransportMode MODE, bool DENSITY_AREA = true>
			double pdf(Vertex const& next, Math::Vector3f const& wo = Math::Vector3f(0, 0, 0))const
			{
				double pdf_solid_angle;
				Math::Vector3f to_vertex = next.position() - position();
				const double dist2 = to_vertex.norm2();
				const double dist = std::sqrt(dist2);
				to_vertex /= dist;
				if (type == Type::Camera)
				{
					pdf_solid_angle = camera().pdfWeSolidAngle<false>(to_vertex);
				}
				else if (type == Type::Light || (wo.norm2() == 0))
				{
					pdf_solid_angle = hit.geometry->getMaterial()->pdfLight(hit, to_vertex);
				}
				else
				{
					if (material()->delta())
						return 0;
					else
						pdf_solid_angle = material()->pdf(hit, to_vertex, wo, MODE == TransportMode::Radiance);
				}
				if constexpr (DENSITY_AREA)
				{
					double res = pdf_solid_angle / dist2;

					if (next.type != Type::Camera)
					{
						res *= std::abs(next.pNormal() * to_vertex);
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
			double cos_prev = std::abs(prev->pNormal() * ray.direction());
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
					vertex.fwd_pdf = pdf_solid_angle * cos_vertex / dist2;

					//sample next direction
					DirectionSample next_dir;
					double pdf_rev;
					hit.geometry->getMaterial()->sampleBSDF(vertex.hit, next_dir, sampler, MODE == TransportMode::Radiance, &pdf_rev);

					prev->rev_pdf = pdf_rev * cos_prev / dist2;

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
			camera_vertex.fwd_pdf = 1;
			camera_vertex.hit.normal = camera_vertex.hit.primitive_normal = ray.direction();
			camera_vertex.hit.point = scene.m_camera.getPosition();
			camera_vertex.hit.camera = &scene.m_camera;
			double pdf_sa = scene.m_camera.pdfWeSolidAngle<true>(ray.direction());
			RGBColor beta = scene.m_camera.We<true>(ray.direction()) / pdf_sa;
			int nv = randomWalk<TransportMode::Importance>(scene, sampler, res, ray, beta, pdf_sa, m_max_len - 1);

			camera_vertex.rev_pdf = 0;
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

			light_vertex.fwd_pdf = sls.pdf;
			light_vertex.beta = 1.0 / sls.pdf;

			//generate a direction
			Math::RandomDirection Le_sampler(&sampler, sls.normal, 1);
			DirectionSample next_dir = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);

			Ray ray(sls.vector, next_dir.direction);

			RGBColor beta = next_dir.bsdf * std::abs(next_dir.direction * sls.normal) / (sls.pdf * next_dir.pdf);

			int nv = randomWalk<TransportMode::Radiance>(scene, sampler, res, ray, beta, (next_dir.pdf), m_max_len - 2);

			return 1 + nv;
		}


		// Returns the sum of the ratios of the pdf except for technique (main_s, main_t)
		double sumRatioVC(Path& cameras, Path& lights, const int main_s, const int main_t, double s1_pdf)const
		{
			const double resolution = cameras[0].camera().resolution;
			double sum = 0;
			//expand the camera sub path
			{
				double ri = 1.0;
				for (int s = main_s; s >= 1; --s)
				{
					const Vertex& camera_end = lights[s - 1];
					const Vertex* light_end = s == 1 ? nullptr : &lights[s - 2];
					ri *= camera_end.rev_pdf / (camera_end.fwd_pdf);

					const double actual_ri = s != 2 ? ri :
						ri * s1_pdf / light_end->fwd_pdf;

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
					ri *= light_end.rev_pdf / light_end.fwd_pdf;

					double actual_ri = ri;
					if (main_s == 0 && t == main_t)
						actual_ri = ri * s1_pdf / light_end.rev_pdf;

					if (!(light_end.delta || camera_end.delta))
					{
						double ni = (t == 2 ? resolution : 1); // account for the extra samples of the light tracer
						sum += (actual_ri * ni);
					}
				}
			}
			return sum;
		}

		double ComputeSumRatioVC(Path& cameras, Path& lights, const int main_s, const int main_t, double s1_pdf, double * buffer)const
		{
			double sum = 0;
			const double resolution = cameras[0].camera().resolution;
			//expand the camera sub path
			{
				double ri = 1.0;
				for (int s = main_s; s >= 1; --s)
				{
					const Vertex& camera_end = lights[s - 1];
					const Vertex* light_end = s == 1 ? nullptr : &lights[s - 2];
					ri *= camera_end.rev_pdf / (camera_end.fwd_pdf);

					const double actual_ri = s != 2 ? ri :
						ri * s1_pdf / light_end->fwd_pdf;

					if (!(camera_end.delta || (light_end && light_end->delta)))
					{
						buffer[s - 1] = actual_ri;
						sum += actual_ri;
					}
					else
						buffer[s - 1] = 0;
				}
			}

			//expand the light subpath
			{
				double ri = 1.0;
				int s = main_s + 1;
				for (int t = main_t; t >= 2; --t)
				{
					const Vertex& light_end = cameras[t - 1];
					const Vertex& camera_end = cameras[t - 2];
					ri *= light_end.rev_pdf / light_end.fwd_pdf;

					double actual_ri = ri;
					if (main_s == 0 && t == main_t)
						actual_ri = ri * s1_pdf / light_end.rev_pdf;

					if (!(light_end.delta || camera_end.delta))
					{
						double ni = (t == 2 ? resolution : 1); // account for the extra samples of the light tracer
						buffer[s] = (actual_ri * ni);
						sum += (actual_ri * ni);
					}
					else
						buffer[s] = 0;
					++s;
				}
			}
			return sum;
		}


		////////////////////////////////////////////////////////////////
		//Computes te MIS weights for Vertex Connection
		// - the last parameter if the probability of sampling the last point on the camera sub path if s == 0 (pure path tracing), else it is not necessary 
		////////////////////////////////////////////////////////////////
		double VCbalanceWeight(
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

			if (xtm)
			{
				xtm_pdf_rev_sa = { &xtm->rev_pdf, xt->pdf<TransportMode::Importance, true>(*xtm, xt->dir_to_vertex(ys)) };
			}

			if (ysm)
			{
				ysm_pdf_rev_sa = { &ysm->rev_pdf, ys->pdf<TransportMode::Radiance, true>(*ysm, ys->dir_to_vertex(xt)) };
			}

			const double actual_ni = main_t == 1 ? cameras[0].camera().resolution : 1;
			const double actual_main_ri = (main_s == 1 ? s1_pdf / lights[0].fwd_pdf : 1);


			double sum = actual_main_ri * actual_ni;

			sum += sumRatioVC(cameras, lights, main_s, main_t, s1_pdf);

			double weight = (actual_main_ri * actual_ni) / sum;
			return weight;
		}

		

	};
}