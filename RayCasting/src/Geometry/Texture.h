#ifndef _Geometry_Texture_H
#define _Geometry_Texture_H

#include <string>
#include <SOIL.h>
#include <Geometry/RGBColor.h>
#include <Math/Vectorf.h>

namespace Geometry
{
	/// <summary>
	/// A texture class. This class loads a texture from the disk.
	/// </summary>
	class Texture
	{
	private:
		int m_width;
		int m_height;
		unsigned char * m_data;

		Texture(const Texture & other)
		{
			std::cerr << "Error, copying a texture!" << std::endl;
		}


	public:
		Texture(const ::std::string & filename)
		{
			m_data = SOIL_load_image(filename.c_str(), &m_width, &m_height, 0, SOIL_LOAD_RGB);
			if (m_data == NULL)
			{
				::std::cerr << "Invalid texture file: " << filename << ::std::endl;
			}
			if (m_width <= 0 || m_height <= 0)
			{
				::std::cerr << "Error on the size of the texture!" << std::endl;
				::std::cerr << "Error on the size of the texture!" << std::endl;
			}
		}


		~Texture()
		{
			SOIL_free_image_data(m_data);
		}

		bool isValid() const
		{
			return m_data != NULL;
		}

		RGBColor pixel(int x, int y) const
		{
			while (x < 0) { x += m_width; }
			while (y < 0) { y += m_height; }
			x = x%m_width;
			y = y%m_height;
			int offset = y * 3 * m_width + x * 3;
			unsigned char r = m_data[offset];
			unsigned char g = m_data[offset + 1];
			unsigned char b = m_data[offset + 2];
			return RGBColor(double(r) / 255.0f, double(g) / 255.0f, double(b) / 255.0f);
		}



		RGBColor pixel(Math::Vector2f const & v) const
		{
			return pixel(int(v[0] * m_width), int(v[1] * m_height));
		}

		RGBColor safe_pixel(Math::Vector2f const& v)const
		{
			if (isValid())
			{
				return pixel(v);
			}
			return 1;
		}
	};
}

#endif