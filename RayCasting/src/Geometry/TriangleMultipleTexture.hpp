#pragma once


#pragma once

#include <Math/Vector.h>
#include <Math/Vectorf.h>
#include <vector>
#include <cassert>


namespace Geometry
{
	////////////////////////////////////////////////////
	//
	// Stores a Triangle texture, optimized for Triangle, shallow textures, and booleans
	//
	////////////////////////////////////////////////////
	template <class T, class precision = double>
	class TriangleMultipleTexture
	{
	protected:

		size_t m_width;
		size_t m_height;

		
		size_t m_size;
		
		std::vector<T> * m_data;

		unsigned int * m_map;

		std::vector<T> m_default_value;

		/*
		static unsigned int to_int(precision p)
		{
			return ((unsigned int)p) + 1;
		}
		*/

		void clean()
		{
			if (m_map != nullptr)
			{
				assert(m_data != nullptr);
				delete[] m_map;
				delete[] m_data;
				m_map = nullptr;
				m_data = nullptr;
			}
			m_default_value = std::vector<T>();
			
		}

		void full_clean()
		{
			clean();
			m_width = m_height = m_size = 0;
		}

		void construct()
		{
			assert(m_width > 0);
			assert(m_height > 0);
			//m_map = new unsigned int[m_height];
			m_map[0] = 0;

			//std::cout << "building texture of size: " << m_width << "x" << m_height << std::endl;

			precision linear_factor = (precision)m_width / (precision)m_height;


			for (size_t pixel_v = 1; pixel_v < m_height; ++pixel_v)
			{
				size_t prev = pixel_v - 1;

				precision prev_line_size = m_width - linear_factor * prev;
				unsigned int prev_line_size_i = (unsigned int)(prev_line_size);


				m_map[pixel_v] = m_map[prev] + prev_line_size_i;
			}

			size_t total = m_map[m_height-1] + (unsigned int)(m_width - linear_factor * (m_height - 1));

			m_data = new std::vector<T>[m_size];
			
			std::fill(m_data, m_data + m_size, std::vector<T>(total));

			
		}


		size_t index(size_t u, size_t v)const
		{
			return m_map[v] + u;
		}

	public:

		TriangleMultipleTexture() :
			m_width(0),
			m_height(0),
			m_size(0),
			m_map(nullptr),
			m_data(nullptr),
			m_default_value()
		{}

		TriangleMultipleTexture(size_t p_size, size_t width, size_t height) :
			m_width(width),
			m_height(height),
			m_size(p_size),
			m_map(new unsigned int[height]),
			m_default_value()
		{
			assert(m_size > 0);
			assert(m_width > 0);
			assert(m_height > 0);
			construct();
		}


		/*
		what if other is not valid???
		*/
		TriangleMultipleTexture(TriangleMultipleTexture const& other) :
			m_width(other.m_width),
			m_height(other.m_height),
			m_map(nullptr),
			m_data(nullptr),
			m_default_value(other.m_default_value)
		{
			if (other.m_map != nullptr)
			{
				assert(other.m_data != nullptr);
				assert(m_size > 1);

				m_map = new unsigned int[m_height];
				std::memcpy(m_map, other.m_map, m_height);
				
				m_data = new std::vector <T>[m_size];
				std::copy(other.m_data, other.m_data + m_size, m_data);
			}
			

		}


		~TriangleMultipleTexture()
		{
			clean();
		}

		void resize(size_t p_size, size_t width, size_t height)
		{
			//TODO Faster if it is the same size...
			clean();
			m_width = width;
			m_height = height;

			assert(m_size > 0);
			assert(m_width > 0);
			assert(m_height > 0);

			m_map = new unsigned int[m_height];

			construct();
		}



		/////////////////////////////////////////
		//TODO no need to interpolate each time
		///////////////////////////////////////////
		precision get_interpolation(const size_t i, precision u, precision v)const
		{
			u *= m_width;
			v *= m_height;
			const int int_u = u;
			const int int_v = v;

			
			
			

			const precision top_left = get_pixel_safe(i, int_u, int_v);
			const precision top_right = get_pixel_safe(i, int_u, int_v + 1);
			const precision down_left = get_pixel_safe(i, int_u+1, int_v);
			const precision down_right = get_pixel_safe(i, int_u+1, int_v + 1);

			

			const precision relativ_u = abs(u - floor(u));
			const precision relativ_v = abs(v - floor(v));

			const precision top = top_right * relativ_v + top_left * (1.0 - relativ_v);
			const precision down = down_right * relativ_v + down_left * (1.0 - relativ_v);

			const precision res = down * relativ_u + top * (1.0 - relativ_u);

			return res;
		}




		T  get_brut_safe(size_t i, precision u, precision v)const
		{
			if (!texture_loaded() || m_data[i].empty())
			{
				assert(i < m_default_value.size());
				return default_value(i);
			}

			return get_brut(i, u, v);
		}



		T  get_brut(size_t i, precision u, precision v)const
		{
			assert(u >= 0 && u <= 1);
			assert(v >= 0 && v <= 1);
			int u_pixel = u * m_width;
			int v_pixel = v * m_height;


			if (v_pixel >= m_height)
			{
				v_pixel = (int)m_height - 1;
			}
			else if (v_pixel < 0)
			{
				v_pixel = 0;
			}

			if (u_pixel >= width(i, v_pixel)) 
			{
				u_pixel = (int)width(i, v_pixel) - 1; 
			}
			else if (u_pixel < 0)
			{
				u_pixel = 0;
			}

			assert(i < m_size);
			assert(m_data != nullptr);


			return get_pixel(i, u_pixel, v_pixel);
		}

		T get_brut_unsafe(size_t i, precision u, precision v)const
		{
			assert(i < m_size);
			assert(m_data != nullptr);
			assert((int)u >= 0);
			//assert((int)(u) < width((int)v);
			assert((int)v >= 0);
			assert((int)v < m_height);
			int u_pixel = u * m_width;
			int v_pixel = v * m_height;
			return get_pixel(i, u_pixel, v_pixel);
		}


		void clear(size_t i)
		{
			assert(i < m_size);
			assert(m_data != nullptr);
			m_data[i] = std::vector<T>();
		}

		void clear()
		{
			full_clean();
		}


		T const& get_pixel_safe(size_t i, int u, int v)const
		{
			int s_u = u, s_v = v;
			if (!texture_loaded() || m_data[i].empty())
			{
				assert(i < m_default_value.size());
				return default_value(i);
			}

			if (v < 0)
			{
				v = 0;
			}
			else if(v >= m_height)
			{
				v = m_height - 1;
			}

			if (u < 0)
			{
				u = 0;
			}
			else if(u >= width(i, v))
			{
				u = width(i, v) - 1;
			}
			

			return get_pixel(i, u, v);
		}

		T const& get_pixel(size_t i, size_t u, size_t v)const
		{
			assert(u < m_width);
			assert(v < m_height);
			assert(i < m_size);
			assert(m_data != nullptr);
			assert(index(u, v) < m_data[i].size());

			std::vector<T>::const_reference pixel = m_data[i][index(u, v)];
			return pixel;
		}

		void set_pixel(size_t i, size_t u, size_t v, T  value)
		{
			assert(u < m_width);
			assert(u < width(i, v));
			assert(v < m_height);
			assert(i < m_size);
			assert(m_data != nullptr);
			size_t id = index(u, v);
			assert(id < m_data[i].size());
			
			std::vector<T>::reference pixel = m_data[i][id];
			pixel = value;
		}

		bool texture_loaded()const
		{
			bool res = m_map != nullptr;
			assert(res == (m_data != nullptr));

			return res;
		}

		bool default_values_loaded()const
		{
			return !m_default_value.empty();
		}

		bool create_default_values()
		{
			if (default_values_loaded())
			{
				return false;
			}
			else
			{
				m_default_value = std::vector<T>(m_size);
				return true;
			}
		}

		T const& default_value(size_t i)const
		{
			assert(default_values_loaded());
			std::vector<T>::const_reference ref = m_default_value[i];
			return ref;
		}

		void set_default(size_t i, T const& value)
		{
			assert(default_values_loaded());
			std::vector<T>::reference ref = m_default_value[i];
			ref = value;
		}





		//I need to add some parameters
		size_t width(size_t i, size_t h) const 
		{ 
			if (h == m_height - 1)
			{
				return m_data[i].size() - m_map[h];
			}
			return m_map[h + 1] - m_map[h]; 
		}
		
		size_t height() const 
		{ 
			return m_height; 
		}

		precision get_u(precision pixel_u)const
		{
			return pixel_u / m_width;
		}

		precision get_v(precision pixel_v)const
		{
			return pixel_v / m_height;
		}





		bool set_uniform(size_t i)
		{
			
			assert(m_data != nullptr);
			assert(i < m_size);
			assert(!m_data[i].empty());

			std::vector<T> & data = m_data[i];
			
			T const& possible_value = data.front();
			
			bool uniform = !std::any_of(data.cbegin(), data.cend(), [&possible_value](T const& t) {return t != possible_value; });
		
			if (uniform)
			{
				m_data[i] = std::vector<T>();
				create_default_values();
				
				std::vector<T>::reference ref = m_default_value[i];
				ref = possible_value;
				return true;
			}
			else
			{
				return false;
			}


			

		}

		void print_infos(std::ostream & out)const
		{
			out << "texture: " << m_width << ", " << m_height << std::endl;
			out << "[";
			for (size_t i = 0; i < m_height; ++i)
			{
				out << m_map[i] << ", ";
			}
			out << "]" << std::endl;
			out << "[";
			for (size_t i = 0; i < m_size; ++i)
			{
				out << m_data[i].size() <<", ";
			}
			out << "]" << std::endl;
		}

		void print_data(std::ostream & out, size_t i)const
		{
			for (unsigned int v = 0; v < m_height; ++v)
			{
				for (unsigned int u = 0; u < width(i, v); ++u)
				{
					T pixel = get_pixel(i, u, v);
					out << pixel;
					out << " ";
				}
				out<<std::endl;
			}
		}


	};
}
