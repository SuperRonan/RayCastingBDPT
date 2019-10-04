#pragma once

#include <Geometry/RGBColor.h>
#include <Geometry/Sphere.h>
#include <Image/Image.h>
#include <cassert>
#include <vector>
#include <algorithm>

namespace Geometry
{
	class EnvironmentMap
	{
	protected:

		struct ImportancePixel
		{
			unsigned int x, y;
			
			double sum;

			//operator for sorting and finding
			bool operator<(ImportancePixel const& other)const
			{
				return sum < other.sum;
			}

			bool operator>(ImportancePixel const& other)const
			{
				return sum > other.sum;
			}

			bool operator<(double const& _sum)const
			{
				return sum < _sum;
			}

			bool operator>(double const& _sum)const
			{
				return sum > _sum;
			}
		};


		Sphere m_sphere;

		Image::Image<RGBColor> m_texture;
		//RGBColor m_texture_mean;
		std::vector<ImportancePixel>  m_importance_pixels;

		double importance(unsigned int i, unsigned int j)const
		{
			RGBColor col = color(Math::Vector<unsigned int, 2>(i, j));
			double theta = (double(j) + 0.5) / double(m_texture.height()) * Math::pi;
			return importance(col, theta);
		}

		double importance(RGBColor const& col, double theta)const
		{
			return col.brightness() * sin(theta);
		}

	public:

		EnvironmentMap(RGBColor color = 0):
			m_sphere(0.0, 0, new Material(color))
		{

		}

		void setColor(RGBColor const& color)
		{
			m_sphere.getMaterial()->setEmissive(color);
		}

		//returns the "general color" of the env map
		const RGBColor & color()const
		{
			return m_sphere.getMaterial()->getEmissive();
		}

		bool isTextured()const
		{
			return m_texture.valid();
		}

		RGBColor textureMean()const
		{
			if (isTextured())
			{
				return m_texture.mean();
			}
			else
			{
				return 1.0;
			}
		}

		RGBColor mean()const
		{
			return textureMean() * color();
		}
		

		void setSphere(Math::Vector3f const& pos, double radius)
		{
			m_sphere.setCenter(pos);
			m_sphere.setRadius(radius);
		}

		template <class INT>
		RGBColor pixel(Math::Vector<INT, 2> const& uv)const
		{
			if (isTextured())
			{
				assert(uv[0] >= 0);
				assert(uv[1] >= 0);
				assert(uv[0] < m_texture.width());
				assert(uv[1] < m_texture.height());
				return m_texture[uv[0]][uv[1]];
			}
			return 1.0;
		}

		template <class INT>
		RGBColor color(Math::Vector<INT, 2> const& uv)const
		{
			return pixel(uv) * color();
		}

		template <class INT>
		RGBColor colorSafe(Math::Vector<INT, 2> uv)const
		{
			if (uv[0] < 0)
				uv[0] = 0;
			if (uv[1] < 0)
				uv[1] = 0;
			if (uv[0] >= m_texture.width())
				uv[0] = m_texture.width() - 1;
			if (uv[1] >= m_texture.height())
				uv[1] = m_texture.height() - 1;
			return color(uv);
		}

		RGBColor color(Math::Vector2f const& uv)const
		{
			return colorSafe(Math::Vector<int, 2>(m_texture.width() * uv[0], m_texture.height()*uv[1]));
		}


		RGBColor color(Math::Vector3f const& direction)const
		{
			Math::Vector3f spc = spherical_coordinate(direction);
			Math::Vector2f uv = { spc[2] / Math::twoPi, spc[1] / Math::pi};
			return color(uv);
		}




		void preCompute(unsigned int div = 1)
		{
			assert(div != 0);
			m_sphere.divide(div);
			m_sphere.build_lights();
			if (isTextured())
			{
				double sum = 0;
				m_importance_pixels.resize(0);
				m_importance_pixels.reserve(m_texture.size());

				for (unsigned int i = 0; i < m_texture.width(); ++i)
				{
					for (unsigned int j = 0; j < m_texture.height(); ++j)
					{
						const RGBColor& pixel = m_texture[i][j];
						if (!pixel.isBlack())
						{
							double importnce = importance(i, j);
							sum += importnce;
							ImportancePixel imp = { i, j, sum };
							m_importance_pixels.push_back(imp);
						}
					}
				}
				m_importance_pixels.shrink_to_fit();


				
			}
			else
			{
				
				
			}
		}


		bool sampleDirection(SurfaceLightSample & res, Math::Sampler & sampler)const
		{
			if (color().isBlack())
			{
				return false;
			}
			if (isTextured() & true)
			{
				{
					//sample a random pixel among the non 0 pixels with considering the importance
					double sum = m_importance_pixels.back().sum;
					double xi = sampler.generateStratified<double>(m_sphere.offset(), m_sphere.divisions());
					auto imp_ptr = ::std::lower_bound(m_importance_pixels.cbegin(), m_importance_pixels.cend(), xi * sum);
					ImportancePixel const& imp = *imp_ptr;
					
					int x = imp.x;
					int y = imp.y;
					double impo = importance(x, y);
					double pdf_pixel = (double)m_texture.size() * impo / sum;
					
					const double xi1 = sampler.generateContinuous<double>();
					const double xi2 = sampler.generateContinuous<double>();

					//get uv coordinates
					double u = (double(x) + xi1) / (double)m_texture.width();
					double v = (double(y) + xi2) / (double)m_texture.height();

					double inclination = v * Math::pi;
					double azimuth = u * Math::twoPi;

					Math::Vector3f dir = Math::Vector3f::make_sphere(inclination, azimuth + Math::pi);

					res.vector = m_sphere.center() + dir * m_sphere.radius();
					res.uv = { u, v };
					res.normal = dir;
					res.geo = &m_sphere;
					res.pdf = pdf_pixel / (sin(inclination) * Math::twoPi * Math::pi);
					return true;
				}




				const double sum = m_importance_pixels.back().sum;
				//stratify the sample
				const double xi = sampler.generateStratified<double>(0, sum, m_sphere.offset(), m_sphere.divisions());
				//const auto pixel_ptr = ::std::lower_bound(m_importance_pixels.cbegin(), m_importance_pixels.cend(), xi);
				const auto pixel_ptr = &m_importance_pixels[xi / sum * m_importance_pixels.size()];
				//if (pixel_ptr == m_importance_pixels.cend())
				{
					//I think this should never happen
					//return false;
				}
				ImportancePixel const& imp = *pixel_ptr;
				res.geo = &m_sphere;//If that makes any sense
				const double impo = importance(imp.x, imp.y);
				//const double pdf = impo / sum;
				const double pdf = 1.0 / m_importance_pixels.size();

				const double xi1 = sampler.generateContinuous<double>();
				const double xi2 = sampler.generateContinuous<double>();

				const double inclination = (double(imp.y) + xi1) / double(m_texture.height()) * Math::pi;
				const double azimuth = (double(imp.x) + xi2) / double(m_texture.width()) * Math::twoPi + Math::pi;
				const Math::Vector3f direction = Math::Vector3f::make_sphere(inclination, azimuth);
				res.vector = m_sphere.center() + direction * m_sphere.radius();
				res.uv = { azimuth / Math::twoPi, inclination / Math::pi };
				res.normal = direction;
				res.pdf = pdf / (Math::pi * Math::pi * 2 * sin(inclination));
			}
			else
			{
				m_sphere.sampleLight(res, sampler);
				res.pdf = 1.0 / (Math::pi * 4);
			}
			return true;
		}

		void incrementOffset(unsigned int i)const
		{
			m_sphere.increment_offset(i);
		}

		void reset()const
		{
			m_sphere.reset();
		}

		double radius()const
		{
			return m_sphere.radius();
		}

		double radius2()const
		{
			return m_sphere.radius2();
		}
		
		Image::Image<RGBColor> const& pixels()const
		{
			return m_texture;
		}

		Image::Image<RGBColor>& pixels()
		{
			return m_texture;
		}
		
	};
}