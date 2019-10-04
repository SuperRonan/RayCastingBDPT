#pragma once

#include <Math/Vector.h>
#include <Math/Vectorf.h>
#include <vector>
#include <cassert>


namespace Geometry
{
	//TODO multi triangle textures 
	template <class T, class precicion=double>
	class TriangleTexture
	{
	protected:

		size_t m_width;
		size_t m_height;
		std::vector<T> m_data;

		unsigned int * m_map;

		
		static unsigned int to_int(precision p) 
		{
			return ((unsigned int)p) + 1;
		}

		void clean()
		{
			if (m_map != nullptr)
			{
				delete[] m_map;
				m_map = nullptr;
			}
		}

		void full_clean()
		{
			clean();
			m_data = std::vector<T>();
			m_width = m_height = 0;
		}

		construct()
		{
			assert(width > 0);
			assert(height > 0);
			m_map[0] = 0;

			precision linear_factor = (precision)m_width / (precicion)m_height;

			for (size_t pixel_v = 1; pixel_v < m_height; ++pixel_v)
			{
				size_t prev = pixel_v - 1;

				precision prev_line_size = m_width - linear_factor * prev;
				unsigned int prev_line_size_i = to_int(prev_line_size);


				m_map[pixel_v] = m_map[prev] + prev_line_size_i;
			}

			size_t total = m_map[pixel_v] + to_int(m_width - linear_factor * (m_height - 1));

			m_data = std::vector<T>(total);
		}

	public:
		T default_value;

		TriangleTexture(T const& p_default = T()) :
			m_width(0),
			m_height(0),
			m_map(nullptr), 
			default_value(p_default)
		{}

		TriangleTexture(size_t width, size_t height, T const& p_default = T()) :
			m_width(width),
			m_height(height),
			m_map(new unsigned int[height]),
			default_value(p_default)
		{
			construct();
		}

		/*
		what if other is not valid???
		*/
		TriangleTexture(TriangleTexture const& other) :
			m_width(other.m_width),
			m_height(other.m_height),
			m_map(new unsigned int[other.m_height]),
			m_data(other.m_data),
			default_value(other.default_value)
		{
			std::memcpy(m_map, other.m_map, m_height);
		}


		~TriangleTexture()
		{
			clean();
		}

		void resize(size_t width, size_t height)
		{
			clean();
			m_width = width;
			m_height = height;

			m_map = new unsigned int[m_height];

			construct();
		}

		T const& get_brut_safe(precision u, precision v)const
		{
			if (!valid())	return default_value;

			return get_brut(u, v);
		}

		T const& get_brut(precision u, precision v)const
		{
			int u_pixel = u * m_width;
			int v_pixel = v * m_height;

			if (u >= m_width)	u = m_width - 1;
			if (v >= m_height)	v = m_height - 1;

			if (u < 0)	u = 0;
			if (v = 0)	v = 0;

			return m_data[m_map[v] + u];
		}

		T const& get_brut_unsafe(precision u, precision v)const
		{
			return m_data[m_map[v] + u];
		}

		void clear()
		{
			full_clean();
		}

		T & access_pixel(size_t u, size_t v)
		{
			assert(u < m_width);
			assert(v < m_height);
			return m_data[m_map[v] + u];
		}

		bool valid()const
		{
			return m_map != nullptr;
		}

		size_t width() const { return m_width; }
		size_t height() const { return m_height; }

		precision get_u(precision pixel_u)const
		{
			return pixel_u / m_width;
		}

		precision get_v(precision pixel_v)const
		{
			return pixel_v / m_height;
		}

	};
}
