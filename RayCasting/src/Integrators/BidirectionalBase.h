#pragma once


#include <Integrators\Integrator.h>

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
			return !scene.intersectionCloser(ray, dir.norm() - 0.001);
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

	};
}