#ifndef _Geometry_Scene_H
#define _Geometry_Scene_H

#include <windows.h>
#include <Geometry/Shapes/Geometry.h>
#include <Geometry/PointLight.h>
#include <Visualizer/Visualizer.h>
#include <Geometry/Camera.h>
#include <Geometry/BoundingBox.h>
#include <System/aligned_allocator.h>
#include <Math/Constant.h>
#include <queue>
#include <functional>
#include <random>
#include <Geometry/LightSampler.h>
#include <cmath>
#include <limits>
#include <Geometry/medium.h>
#include <stack>
#include <utils.h>
#include <valarray>
#include <tbb/parallel_for_each.h>
#include <functional>
#include <Geometry/Hit.h>
#include <System/BoundedStack.h>
#include <tbb/atomic.h>
#include "DirectionalLight.h"
#include <Geometry/Shapes/Sphere.h>
#include <Math/Sampler.h>
#include <settings.h>
#include <Geometry/EnvironmentMap.h>
#include <Geometry/BVH.h>
#include <Geometry/Shapes/Disk.h>
#include <System/Parallel.h>
#include <map>
#include <unordered_map>

namespace Geometry
{

	class Scene
	{
	public:
		/// \brief	The scene geometry (basic representation without any optimization).
		//need to use pointers to keep informations about the class
		::std::deque<GeometryCollection *> m_geometries;

		std::vector<Sphere> m_spheres;

		std::vector<Disk> m_disks;
		
		/// \brief	The lights.
		//std::vector<PointLight> m_lights;

		Camera m_camera;

		BoundingBox m_sceneBoundingBox;

		
		std::vector<GeometryBase*> m_surface_lights;

		medium m_medium;

		size_t m_number_triangle = 0;

		EnvironmentMap m_envmap;

#ifdef COUNT_RAYS
		mutable tbb::atomic<size_t> ray_counter=0;
#endif



		BVH<Triangle> m_triangle_bvh;
		BVH<Sphere> m_sphere_bvh;
		BVH<Disk> m_disk_bvh;
		



		Scene(RGBColor const& sky_color=0, medium med=medium()):
			m_medium(med),
			m_envmap(sky_color)
		{}

		/// <summary>
		/// Prints stats about the geometry associated with the scene
		/// </summary>
		void printStats()
		{
			::std::cout << "Scene: " << m_number_triangle << " triangles" << ::std::endl;
		}

		/// <summary>
		/// Computes the scene bounding box.
		/// </summary>
		/// <returns></returns>
		const BoundingBox & getBoundingBox()
		{
			return m_sceneBoundingBox;
		}



		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::add(const GeometryCollection & geometry)
		///
		/// \brief	Adds a geometry to the scene.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		///
		/// \param	geometry The geometry to add.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void add(GeometryCollection * geometry)
		{
			geometry->recomputeBox();
			if (geometry->getVertices().size() == 0) { return; }
			m_geometries.push_back(geometry) ;
			m_geometries.back()->computeVertexNormals(Math::piDiv4/2);
			m_number_triangle += geometry->getTriangles().size();

			if (m_geometries.back()->getMaterial()->is_emissive())
			{
				m_surface_lights.push_back(m_geometries.back());
			}
			m_sceneBoundingBox.update(geometry->box());
		}

		void add(Sphere const& s)
		{
			m_spheres.push_back(s);
			if (m_spheres.back().getMaterial()->is_emissive())
			{
				m_surface_lights.push_back(&m_spheres.back());
			}
			m_sceneBoundingBox.update(s.box());
		}

		void add(Disk const& d)
		{
			m_disks.push_back(d);
			if (m_disks.back().getMaterial()->is_emissive())
			{
				m_surface_lights.push_back(&m_disks.back());
			}
			m_sceneBoundingBox.update(d.box());
		}


		void add(const PointLight & light)
		{
			std::cerr << "Warning, Point lights are not currently supported!" << std::endl;
			//m_lights.push_back(light) ;
		}


		void setCamera(Camera const & cam)
		{
			m_camera = cam ;
		}


		bool full_intersection(Ray const& ray, Hit & res)const
		{
			if (empty())
				return false;

			Intersection<Triangle> rti;
			m_triangle_bvh.intersection(ray, rti);

			Intersection<Sphere> rsi;
			m_sphere_bvh.intersection(ray, rsi);

			Intersection<Disk> rdi;
			m_disk_bvh.intersection(ray, rdi);

			if (rti.valid() || rsi.valid() || rdi.valid())
			{
				if (rti < rsi)
				{
					if (rdi < rti)
					{
						rdi.fill(res, ray);
					}
					else
					{
						rti.fill(res, ray);
					}
				}
				else
				{
					if (rdi < rsi)
					{
						rdi.fill(res, ray);
					}
					else
					{
						rsi.fill(res, ray);
					}
				}

				return true;
			}
			else
			{
				return false;
			}
		}
		

		bool intersectionCloser(Ray const& ray, double d)const
		{
			Intersection<Triangle> rti;
			m_triangle_bvh.intersection(ray, rti);
			if (rti.valid() && rti.t() < d)
			{
				return true;
			}

			Intersection<Sphere> rsi;
			m_sphere_bvh.intersection(ray, rsi);
			if (rsi.valid() && rsi.t() < d)
			{
				return true;
			}
			Intersection<Disk> rdi;
			m_disk_bvh.intersection(ray, rdi);
			if (rdi.valid() && rdi.t() < d)
			{
				return true;
			}
			return false;
		}


		bool noIntersection(Ray const& ray)const
		{
			Intersection<Triangle> rti;
			if (m_triangle_bvh.intersection(ray, rti))
			{
				return false;
			}
			Intersection<Sphere> rsi;
			if (m_sphere_bvh.intersection(ray, rsi))
			{
				return false;
			}
			Intersection<Disk> rdi;
			if (m_disk_bvh.intersection(ray, rdi))
			{
				return false;
			}
			return true;
		}

		bool cameraVisibility(Math::Vector3f const& point)const
		{
			Math::Vector3f dir = m_camera.m_position - point;
			Ray ray(point, dir);
			return !intersectionCloser(ray, dir.norm());
		}

		bool visibility(Math::Vector3f const& p, Math::Vector3f const& q)const
		{
			Math::Vector3f dir = p - q;
			Ray ray(q, dir);
			return !intersectionCloser(ray, dir.norm() - 0.000001); // not stonks
		}

		bool cameraSeesEnvMap(Math::Vector3f const& dir)const
		{
			return noIntersection(Ray(m_camera.m_position, dir));
		}

		bool empty()const
		{
			return m_geometries.empty() && m_spheres.empty() && m_disks.empty();
		}


	public:

		///////////////////////////////
		// Computes the number of intersection between the ray and the acceleration structure
		// This version shows the bounding boxes of the geometrys
		///////////////////////////////
		unsigned int count_bb_intersections(Ray const& ray)const
		{
			//count the bounding boxes of the BVHs
			return m_triangle_bvh.countBB(ray) + m_sphere_bvh.countBB(ray) + m_disk_bvh.countBB(ray);

			//count the bouding boxes of the geometries
			double t_entry, t_exit;
			unsigned int res = 0;
			for (auto const& p : m_geometries)
			{
				if (p->box().intersect(ray, 0, std::numeric_limits<double>::max(), t_entry, t_exit))
				{
					++res;
				}
			}
			for (auto const& s : m_spheres)
			{
				if (s.box().intersect(ray, 0, std::numeric_limits<double>::max(), t_entry, t_exit))
				{
					++res;
				}
			}
			for (auto const& d : m_disks)
			{
				if (d.box().intersect(ray, 0, std::numeric_limits<double>::max(), t_entry, t_exit))
				{
					++res;
				}
			}
			return res;
		}


		//////////////////////////////////////
		// counts the bounding boxes of the acceleration structures for the ray, and convets it into a color
		//////////////////////////////////////
		RGBColor count_box(Ray const& ray)const
		{
			int n = 2*count_bb_intersections(ray);
			double red(0), green(0), blue(0);

			if (n < 256)
			{
				blue = n;
			}
			else if (n < 512)
			{
				blue = 512 - n;
				red = n - 256;
			}
			else
			{
				red = 3 * 256 - n;
				green = n - 512;
			}
			

			RGBColor res(red, green, blue);
			res = res / 255.;
			return res;
		}



	public:
		//get the background color of a ray lost in the skybox
		RGBColor getBackgroundColor(Math::Vector3f direction)const
		{
			return m_envmap.color(direction);
		}

		RGBColor getBackgroundColor()const
		{
			return m_envmap.color();
		}

		void setBackgroundColor(RGBColor const& col)
		{
			m_envmap.setColor(col);
		}

		inline bool intersection_light(Ray const& ray, double light_t)const
		{
			//TODO Add optimisation
			Hit hit;
			if (full_intersection(ray, hit))
			{
				return hit.z < light_t;
			}
			return false;
		}

		//size_t num_shadow = 0;



	protected:


		//////////////////////////////////////////////
		//Computes the visibily between a point and a light
		//////////////////////////////////////////////
		bool is_lighted(Math::Vector3f const& point, PointLight const& l)const
		{
			Math::Vector3f to_light = l.position() - point;

			double dist_to_light_2 = to_light.norm2();
			double dist_to_light = sqrt(dist_to_light_2);


			Ray ray(point, to_light / dist_to_light);

			

			if (intersection_light(ray, dist_to_light))
			{
				return 0;
			}
			else
			{
				// light is reached
				return 1;
			}
		}


	

		
		/*
		/////////////////////////////////////////////////
		//Compute the shadow texture of a triangle,
		//TODO:
		// - Optimisation (allocate if necessary, use more cache efficient representation)
		// - CLean the artefacts
		////////////////////////////////////////////////
		void compute_triangle_shadow(Triangle & triangle, size_t width, size_t height)
		{
			if (triangle.m_light_tex != nullptr)
			{
				delete triangle.m_light_tex;
			}
			triangle.m_light_tex = new TriangleMultipleTexture<bool>(2 * m_lights.size(), width, height);
			TriangleMultipleTexture<bool> & light_texture = *(triangle.m_light_tex);

			for (unsigned char facing = 0; facing < 2; ++facing)
			{
				
				for (unsigned int light_index = 0; light_index < m_lights.size(); ++light_index)
				{
					PointLight const& l = m_lights[light_index];
					bool facing_light = false;
					for (unsigned int i = 0; i < 3; ++i)
					{
						Math::Vector3f vertex = triangle.vertex(i);
						Math::Vector3f normal = triangle.getVertexNormal(i, facing);

						Math::Vector3f to_light = l.position() - vertex;

						if (to_light * normal > 0.0)
						{
							facing_light = true;
							break;
						}
					}

					size_t tex_index = triangle.get_light_tex_index(light_index, facing);

					
					const double mid_pixel = 0.;
					
					
					if (facing_light)
					{
						for (size_t v = 0; v < light_texture.height(); ++v)
						{
							double v_value = light_texture.get_v(v + mid_pixel);
							for (size_t u = 0; u < light_texture.width(tex_index, v); ++u)
							{
								
								
								double u_value = light_texture.get_u(u + mid_pixel);
								Math::Vector3f sampled_point = triangle.samplePoint(u_value, v_value);
								Math::Vector3f sampled_normal = triangle.sampleNormal(u_value, v_value, facing);

								Math::Vector3f to_light = l.position() - sampled_point;

								if (to_light * sampled_normal < 0.0)
								{
									light_texture.set_pixel(tex_index, u, v, false);
								}
								else
								{

									light_texture.set_pixel(tex_index, u, v, is_lighted(sampled_point, l));
								}
							}
						}
						if (light_texture.set_uniform(tex_index))
						{
							
						}
						else
						{
							num_shadow += 1;

						}
					}
					else
					{
						light_texture.clear(tex_index);
						light_texture.create_default_values();
						light_texture.set_default(tex_index, false);
					}
					
				}
			}
		}
		*/
	public:



		void preCompute(unsigned int n_max, double ct, double ci)
		{
			{
				std::vector<const Triangle*> triangles;
				triangles.reserve(m_number_triangle);
				for (unsigned int i = 0; i < m_geometries.size(); ++i)
				{
					m_geometries[i]->recomputeBox();
					for (unsigned int j = 0; j < m_geometries[i]->getTriangles().size(); ++j)
					{
						triangles.push_back(&m_geometries[i]->getTriangles()[j]);
					}
				}
				m_triangle_bvh.compute(n_max, ct, ci, triangles);
			}
			{
				std::vector<const Sphere*> spheres(m_spheres.size());
				OMP_PARALLEL_FOR
					for (long i = 0; i < spheres.size(); ++i)
					{
						spheres[i] = &m_spheres[i];
					}
				m_sphere_bvh.compute(n_max, ct, ci, spheres);
			}
			{
				std::vector<const Disk*> disks(m_disks.size());
				OMP_PARALLEL_FOR
					for (long i = 0; i < disks.size(); ++i)
					{
						disks[i] = &m_disks[i];
					}
				m_disk_bvh.compute(n_max, ct, ci, disks);
			}
			

			m_sceneBoundingBox = BoundingBox();
			if (!m_triangle_bvh.empty())
			{
				m_sceneBoundingBox = m_triangle_bvh.boundingBox();
			}
			if (!m_sphere_bvh.empty())
			{
				m_sceneBoundingBox.update(m_sphere_bvh.boundingBox());
			}
			if (!m_disk_bvh.empty())
			{
				m_sceneBoundingBox.update(m_disk_bvh.boundingBox());
			}
		}


	protected:

		struct RISCandidate
		{
			double w, p_target;
			SurfaceSample sample;
			RGBColor estimate;
		};

		mutable std::vector<std::vector<RISCandidate>> m_candidates_buffers;

		int m_ris_number_of_candidates;

		std::map<const GeometryBase*, int> m_light_to_index;

	public:	
		void preAllocate()
		{
			m_ris_number_of_candidates = m_surface_lights.size();
			m_candidates_buffers = Parallel::preAllocate(std::vector<RISCandidate>(m_ris_number_of_candidates));
			m_disk_bvh.preAllocate();
			m_triangle_bvh.preAllocate();
			m_sphere_bvh.preAllocate();
		}


		
		/*
		/////////////////////////////////////
		//Pre compute the shadows of all the triangles of the scene
		// - This function can be quite expenseive
		// - The size of the shades are at the scales width_factor and height_factor
		// TODO: add a minimal size
		////////////////////////////////////
		void pre_compute_shadows(double width_factor, double height_factor, size_t min_width=1, size_t min_height=1)
		{
			num_shadow = 0;
			std::vector<Triangle *> triangle_list;
			triangle_list.reserve(m_number_triangle);
			for (auto & p : m_geometries)
			{
				for (Triangle & triangle : p.second.getTriangles())
				{
					triangle_list.push_back(&triangle);
				}
			}

			int step = std::max(1ull, triangle_list.size() / 100ull);

#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < triangle_list.size(); ++i)
			{
				Triangle & triangle = *(triangle_list[i]);

				size_t width = std::max(min_width, (size_t)(width_factor * triangle.uAxis().norm()) + 1);
				size_t height = std::max(min_height, (size_t)(height_factor * triangle.vAxis().norm()) + 1);

				compute_triangle_shadow(triangle, width, height);

				//Printing status
				if (i % step == 0)
				{
					int percent = (i * 100) / triangle_list.size();
#pragma omp critical(cout)
					std::cout << "\r" + progession_bar(i, triangle_list.size(), 50) << std::flush;
				}
			}
			std::cout << "\r" + progession_bar(1, 1, 50) << std::endl;
		}
		*/

		//To complete
		void check_capacity()
		{
			/*const size_t pl = m_lights.size();
			const size_t ls = m_lightSamples * m_surface_lights.size();
			const size_t ii = m_diffuseSamples + m_specularSamples;
			std::cout << "Capacity: " << LightStack::capacity << "\n";
			std::cout << "Ray tracing requires: " << pl + 1 << "\n";
			std::cout << "Path tracing requires: " << ls + ii << "\n";*/
		}





		void compute_light_samplers(unsigned int div=1)
		{
			assert(div > 0);
			//standard geometry surface lights of the scene
			m_surface_lights.clear();
			for (auto it = m_geometries.begin(), end = m_geometries.end(); it != end; ++it)
			{
				GeometryBase* geometry = *(it);
				if (geometry->getMaterial()->is_emissive())
				{
					m_surface_lights.push_back(geometry);
					geometry->build_lights();
				}
			}

			for (Sphere& sphere : m_spheres)
			{
				if (sphere.getMaterial()->is_emissive())
				{
					m_surface_lights.push_back(&sphere);
				}
			}
			std::cout << m_surface_lights.size() << " surface light(s) detected!" << std::endl;
			for (GeometryBase* geo : m_surface_lights)
			{
				geo->divide(div);
			}
			
			//skybox of the scene, for bidir methods
			double longest_radius = 0.0;
			if (empty())
			{
				longest_radius = 1.0;
			}
			else
			{
				for (int x = 0; x < 2; ++x)
				{
					for (int y = 0; y < 2; ++y)
					{
						for (int z = 0; z < 2; ++z)
						{
							Math::Vector3f radius = { m_sceneBoundingBox[x][0], m_sceneBoundingBox[y][1], m_sceneBoundingBox[z][2] };
							radius -= m_sceneBoundingBox.center();
							if (radius.norm() > longest_radius)
							{
								longest_radius = radius.norm();
							}
						}
					}
				}
			}

			m_envmap.setSphere(m_sceneBoundingBox.center(), longest_radius * 1.01);

			m_envmap.preCompute(div);

			m_light_to_index.clear();
			for (int i = 0; i < m_surface_lights.size(); ++i)
			{
				m_light_to_index.emplace(m_surface_lights[i], i);
			}
		}


		//A const update (-_-)
		void update_lights_offset(unsigned int lights_samples)const
		{
			for (const GeometryBase* geometry : m_surface_lights)
			{
				geometry->increment_offset(lights_samples);
			}
			m_envmap.incrementOffset(lights_samples);
		}
		

		void reset_surface_lights()const
		{
			for (const GeometryBase* geometry : m_surface_lights)
			{
				geometry->reset();
			}
			m_envmap.reset();
		}



		bool sampleOneLight(Math::Sampler& sampler, double & pdf, const GeometryBase *& geo)const
		{
			//TODO consider the surface of the lights, the posiion of the hit? 
			if (m_surface_lights.empty())
			{
				return false;
			}
			size_t index = sampler.generate(0, m_surface_lights.size() - 1);
			geo = m_surface_lights[index];
			pdf = 1.0 / double(m_surface_lights.size());
			return true;
		}

		double pdfSamplingLight(const GeometryBase* geo)const
		{
			1.0 / double(m_surface_lights.size());
		}

		void sampleLe(Math::Sampler& sampler, SurfaceSample& sample, int offset=0)const
		{
			assert(!m_surface_lights.empty());
			double pdf;
			sampleOneLight(sampler, pdf, sample.geo);
			sample.geo->sampleLight(sample, sampler, offset);
			sample.pdf *= pdf;
		}

		double pdfSampleLe(const GeometryBase* geo)const
		{
			//= probability of sampling this geonetry * probability of sampling a point on this geometry
			double pdf_geo = 1.0 / double(m_surface_lights.size());
			double pdf_point = geo->pdfSamplingPoint();

			return pdf_geo * pdf_point;
		}

		void sampleLi(Math::Sampler& sampler, SurfaceSample& sample, Hit const& ref, int offset = 0)const
		{
			assert(!m_surface_lights.empty());
			double pdf;
			sampleOneLight(sampler, pdf, sample.geo);
			sample.geo->sampleLight(sample, ref, sampler, offset);
			sample.pdf *= pdf;
		}

		double pdfSampleLi(const GeometryBase* geo, Hit const& hit, Math::Vector3f const& point)const
		{
			double pdf_geo = 1.0 / double(m_surface_lights.size());
			double pdf_point = geo->pdfSamplingPoint(hit, point);
			return pdf_geo * pdf_point;
		}

		int m_ris_estimate_N = 1;

		void sampleLiRIS(Math::Sampler& sampler, SurfaceSample& sample, Hit const& ref, RGBColor * Contribution=nullptr)const
		{
			if (m_ris_number_of_candidates == 1) // M = 1 -> classic sampleLi
			{
				sampleLi(sampler, sample, ref);
				if (Contribution)
				{
					Math::Vector3f to_light = sample.vector - ref.point;
					double dist2 = to_light.norm2();
					to_light /= std::sqrt(dist2);
					double G = std::abs((to_light * ref.primitive_normal) * (to_light * sample.normal)) / dist2;
					*Contribution = ref.geometry->getMaterial()->BSDF(ref, to_light) * 
						sample.geo->getMaterial()->Le(sample.normal, sample.uv, -to_light) * G;
					if (Contribution->anythingWrong())	*Contribution = 0;
				}
				return;
			}

			// Draw the set of candidates
			
			std::vector<RISCandidate>& candidates = m_candidates_buffers[Parallel::tid()];
			double sum = 0;
			int index = 0;
			for (int i = 0; i < m_ris_number_of_candidates; ++i)
			{
				RISCandidate& candidate = candidates[index];
				//int light_index = sampler.generate(0, m_surface_lights.size() - 1);
				int light_index = sampler.generateStratified<double>(0, m_surface_lights.size(), i, m_ris_number_of_candidates);
				const GeometryBase* geo = m_surface_lights[light_index];
				geo->sampleLight(candidate.sample, ref, sampler);
				candidate.sample.pdf /= m_surface_lights.size();
				Math::Vector3f to_light = candidate.sample.vector - ref.point;
				const double dist2 = to_light.norm2();
				to_light /= std::sqrt(dist2);
				double G = std::abs((candidate.sample.normal * to_light) * (ref.primitive_normal * to_light)) / dist2;
				const RGBColor Le = candidate.sample.geo->getMaterial()->Le(candidate.sample.normal, candidate.sample.uv, -to_light);
				const RGBColor bsdf = ref.geometry->getMaterial()->BSDF(ref, to_light, false);
				const RGBColor L = Le * bsdf * G;
				if (L.anythingWrong() || L.isBlack() || L.grey() < 1e-100) // I don't like this, but is creates nan else
				{
					// candidate = 0;
					continue;
				}
				else
				{
					candidate.w = L.grey() / candidate.sample.pdf;
					candidate.p_target = L.grey();
					candidate.estimate = L;
					++index;
				}
				sum += candidate.w;
			}

			// Select one of the candidates by using their weight w as probability

			if (sum == 0)
			{
				if (Contribution) *Contribution = 0;
				return;
			}
			double xi = sampler.generateContinuous<double>(0, sum);
			double partial_sum = 0;
			for (int i = 0; i < candidates.size(); ++i)
			{
				double new_sum = partial_sum + candidates[i].w;
				if (xi >= partial_sum && xi < new_sum)
				{
					sample = candidates[i].sample;
					double p_target = candidates[i].p_target;
					sample.pdf = p_target * m_ris_number_of_candidates / sum;
					if(m_ris_estimate_N > 1) 
					{
						Hit sample_hit;
						sample_hit.point = sample.vector;
						sample_hit.geometry = sample.geo;
						sample.primitive = sample.primitive;
						sample_hit.normal = sample_hit.primitive_normal = sample.normal;
						sample_hit.tex_uv = sample.uv;
						double estimate = pdfRISEstimate(ref, sample_hit, sampler, candidates[i].estimate, m_ris_estimate_N-1);
						sample.pdf = (sample.pdf + estimate * (m_ris_estimate_N - 1)) / m_ris_estimate_N;
					}
					assert(p_target > 0);
					if (Contribution)
						*Contribution = candidates[i].estimate;
					return;
				}
				partial_sum = new_sum;
			}
			assert(false);
		}

		double pdfRISEstimate(Hit const& ref, Hit const& sample, Math::Sampler& sampler, RGBColor const& contribution)const
		{
			return pdfRISEstimate(ref, sample, sampler, contribution, m_ris_estimate_N);
		}

		double pdfRISEstimate(Hit const& ref, Hit const& sample, Math::Sampler & sampler, RGBColor const& contribution, int N)const
		{
			if (N == 0)	return 0;
			if (m_ris_number_of_candidates == 1)
				return pdfSampleLi(sample.geometry, ref, sample.point);
			
			double p_target, p_source;
			p_source = pdfSampleLi(sample.geometry, ref, sample.point);
			p_target = contribution.grey();
			if (p_target == 0)
				return 0;
			double w = p_target / p_source;
			double res = 0;
			int sample_light_index = m_light_to_index.at(sample.geometry);
			std::vector<RISCandidate>& candidates = m_candidates_buffers[Parallel::tid()];
			for (int n = 0; n < N; ++n)
			{
				bool has_seen_sample_light_index = false;
				double sum = w;
				int index = 0;
				for (int i = 0; i < m_ris_number_of_candidates; ++i) 
				{
					if (i == sample_light_index && !has_seen_sample_light_index) // 1 less
					{
						// Warning! this works because the number of ris candidates == number of lights!!!
						assert(m_ris_number_of_candidates == m_surface_lights.size());
						has_seen_sample_light_index = true;
						continue; //discard this index 
					}
					RISCandidate& candidate = candidates[index];
					//int light_index = sampler.generate(0, m_surface_lights.size() - 1);
					int light_index = sampler.generateStratified<double>(0, m_surface_lights.size(), i, m_ris_number_of_candidates);
					const GeometryBase* geo = m_surface_lights[light_index];
					geo->sampleLight(candidate.sample, ref, sampler);
					candidate.sample.pdf /= m_surface_lights.size();
					Math::Vector3f to_light = candidate.sample.vector - ref.point;
					const double dist2 = to_light.norm2();
					to_light /= std::sqrt(dist2);
					double G = std::abs((candidate.sample.normal * to_light) * (ref.primitive_normal * to_light)) / dist2;
					const RGBColor Le = candidate.sample.geo->getMaterial()->Le(candidate.sample.normal, candidate.sample.uv, -to_light);
					const RGBColor bsdf = ref.geometry->getMaterial()->BSDF(ref, to_light, false);
					const RGBColor L = Le * bsdf * G;
					if (L.anythingWrong() || L.isBlack() || L.grey() < 1e-100) // I don't like this, but is creates nan else
					{
						// candidate = 0
						continue; 
					}
					else
					{
						candidate.w = L.grey() / candidate.sample.pdf;
						candidate.p_target = L.grey();
						++index;
					}
					sum += candidate.w;
				}
				assert(has_seen_sample_light_index);
				double pdf = m_ris_number_of_candidates * p_target / sum;
				res += pdf;
			}
			return res / N;
		}

		double pdfSkybox()const
		{
			if (empty())
			{
				return 1.0;
			}
			return !m_envmap.color().isBlack() ? 0.5 : 0;
		}

		double radius()const
		{
			return m_envmap.radius();
		}

		double radius2()const
		{
			return m_envmap.radius2();
		}

		void setBackgroundPixels(Image::Image<RGBColor> const& img)
		{
			m_envmap.pixels() = img;
		}

	} ;
}

#endif
