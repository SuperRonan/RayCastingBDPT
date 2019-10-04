#pragma once

#include <Math/Vector.h>
#include <Math/Vectorf.h>

namespace Geometry
{
	template <class T, class precision=double>
	class generic_texture
	{
	protected:


		size_t m_height, m_width;
		size_t m_size;
		T * m_data;
		




		inline size_t index(size_t u, size_t v)const
		{
			return (m_height * u) + v;
		}


		//Maybe add (later) other pre computed stuff
		//like the rows

	public:

		T default_value;

		generic_texture(T const& p_default_value=T()):
			m_height(0), m_width(0), m_size(0), m_data(nullptr), default_value(p_default_value)
		{}

		generic_texture(size_t width, size_t height, T const& p_default_value = T()) :
			m_width(width), m_height(height), m_size(width * height), m_data(new T[width * height]), default_value(p_default_value)
		{}


		//TODO, what if the other is not valid?
		generic_texture(generic_texture const& other)
			:m_height(other.m_height), m_width(other.m_width), m_size(other.m_size)
		{
			//std::cout << "copying" << std::endl;
			if (m_size > 0)
			{
				m_data = new T[m_size];
				std::memcpy(m_data, other.m_data, m_size);
			}
			
		}

		~generic_texture()
		{
			if (m_data != nullptr)
			{
				delete[] m_data;
			}
		}


		void resize(size_t width, size_t height)
		{
			if (m_data != nullptr)
			{
				delete[] m_data;
			}

			m_width = width;
			m_height = height;
			m_size = width * height;

			m_data = new T[m_size];
		}

		/*
		T & get_brut(Math::Vector2f const& uv)
		{
			return get_brut(uv[0], uv[1]);
		}

		T & get_brut(precision up, precision vp)
		{
			int u = up * m_width;
			int v = vp * m_height;

			u = std::min(u, m_width - 1);
			v = std::min(v, m_height - 1);

			u = std::max(u, 0);
			v = std::max(v, 0);

			return m_data[u * m_width + v];
		}
		*/

		T const& get_brut(Math::Vector2f const& uv)const
		{
			return get_brut(uv[0], uv[1]);
		}

		T const& get_brut_safe(precision up, precision vp)const
		{
			if (!valid())	return default_value;
			
			return get_brut(up, vp);
		}

		T const& get_brut(precision up, precision vp)const
		{
			int u = up * m_width;
			int v = vp * m_height;

			u = std::min(u, (int)(m_width - 1));
			v = std::min(v, (int)(m_height - 1));

			u = std::max(u, 0);
			v = std::max(v, 0);

			return m_data[index(u, v)];
		}



		void clear()
		{
			if (m_data != nullptr)
			{
				delete[] m_data;
			}

			m_size = m_width = m_height = 0;
		}


		void set_value(size_t u, size_t v, T const& value)
		{
			m_data[index(u, v)] = value;
		}

		T & access_pixel(size_t u, size_t v)
		{
			return m_data[index(u, v)];
		}


		bool valid()const
		{
			return m_data != nullptr;
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