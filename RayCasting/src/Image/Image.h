#pragma once
#include <cassert>
#include <algorithm>
#include <Math/Vector.h>


namespace Image
{
	bool constexpr IMAGE_ROW_MAJOR = true;
	bool constexpr IMAGE_COL_MAJOR = false;
	///////////////////////////////////
	// Class representing a 2D Image (of anything)
	///////////////////////////////////
	template <class T, bool ROW_MAJOR=IMAGE_COL_MAJOR>
	class Image
	{
	protected:

		size_t m_size, m_width, m_height;
		

	public:

		T* m_data;

	public:

		

		Image(size_t width, size_t height):
			m_size(width*height),
			m_width(width),
			m_height(height),
			m_data(new T[m_size])
		{}

		Image(Image const& other) :
			m_size(other.m_size),
			m_width(other.m_width),
			m_height(other.m_height)
		{
			if (other.valid())
			{
				m_data = new T[m_size];
				//or a std::copy?
				std::memcpy(m_data, other.m_data, m_size * sizeof(T));
			}
			else
			{
				m_width = m_height = m_size = 0;
				m_data = nullptr;
			}
		}

		Image():
			m_size(0),
			m_width(0),
			m_height(0),
			m_data(nullptr)
		{}


		~Image()
		{
			if(m_data)
			{
				delete[] m_data;
			}
		}

		Image& operator=(Image const& other)
		{
			if (m_size == other.m_size && valid())
			{
				assert(valid());
				assert(other.valid());
				m_width = other.m_width;
				m_height = other.m_height;
				std::memcpy(m_data, other.m_data, m_size * sizeof(T));
				return *this;
			}
			if (!other.valid())
			{
				if (valid())
				{
					delete[] m_data;
					m_size = other.m_size;
					m_width = other.m_width;
					m_height = other.m_height;
					m_data = other.m_data;
					assert(m_data == nullptr);
				}
			}
			else
			{
				if (valid())
				{
					delete[] m_data;
				}
				m_size = other.m_size;
				m_width = other.m_width;
				m_height = other.m_height;
				m_data = new T[m_size];
				std::memcpy(m_data, other.m_data, m_size * sizeof(T));
			}
			return *this;
		}

		void resize(size_t w, size_t h)
		{
			assert(w != 0);
			assert(h != 0);
			if (!valid())
			{
				m_width = w;
				m_height = h;
				m_size = w * h;
				m_data = new T[m_size];
			}
			else
			{
				assert(m_data != nullptr);
				if (width() != w || height() != h)
				{
					m_width = w;
					m_height = h;
					if (m_size != w * h)
					{
						delete[] m_data;
						m_size = m_width * m_height;
						m_data = new T[m_size];
					}
				}
			}
		}


		bool valid()const
		{
			return m_data != nullptr;
		}


		size_t size()const
		{
			return m_size;
		}

		size_t width()const
		{
			return m_width;
		}

		size_t height()const
		{
			return m_height;
		}
		

		size_t index(size_t i, size_t j)const
		{
			if constexpr (ROW_MAJOR)
				return i * m_height + j;
			else
				return j * m_width + i;
		}

		Math::Vector<size_t, 2> indices(size_t index)const
		{
			if constexpr (ROW_MAJOR)
				return { index / m_height, index % m_height };
			else
				return { index % m_width, index / m_width };
		}
		

		// -i |j
		inline T const& pixel(size_t i, size_t j)const
		{
			assert(i < m_width);
			assert(j < m_height);
			return m_data[index(i, j)];
		}

		// -i |j
		inline T & pixel(size_t i, size_t j)
		{
			assert(i < m_width);
			assert(j < m_height);
			return m_data[index(i, j)];
		}
		
		inline T const& operator()(size_t i, size_t j)const
		{
			return pixel(i, j);
		}

		inline T & operator()(size_t i, size_t j)
		{
			return pixel(i, j);
		}

		inline T const& operator()(Math::Vector<size_t, 2> const& p)const
		{
			return pixel(p[0], p[1]);
		}

		inline T& operator()(Math::Vector<size_t, 2> const& p)
		{
			return pixel(p[0], p[1]);
		}

		inline T const& operator[](size_t index)const
		{
			return m_data[index];
		}

		inline T& operator[](size_t index)
		{
			return m_data[index];
		}

		
		////should be used with an other []
		//inline const T * operator[](size_t i)const
		//{
		//	assert(i < m_width);
		//	return m_data + i * m_height;
		//}

		////should be used with an other []
		//inline  T* operator[](size_t i)
		//{
		//	assert(i < m_width);
		//	return m_data + i * m_height;
		//}


		inline void fill(T const& t = 0)
		{
			assert(m_data);
			std::fill(m_data, m_data + m_size, t);
		}

		T* begin()
		{
			return m_data;
		}

		T* end()
		{
			return m_data + m_size;
		}

		const T* begin()const
		{
			return m_data;
		}

		const T* end()const
		{
			return m_data + m_size;
		}

		const T* cbegin()const
		{
			return m_data;
		}

		const T* cend()const
		{
			return m_data + m_size;
		}

		Image& operator-=(Image const& other)
		{
			assert(m_size == other.m_size);
			for (size_t i = 0; i < m_size; ++i)
			{
				m_data[i] -= other.m_data[i];
			}
			return *this;
		}

		Image operator-(Image const& other)const
		{
			Image<T> res = *this;
			return res -= other;
		}

		Image& operator+=(Image const& other)
		{
			assert(m_size == other.m_size);
			for (size_t i = 0; i < m_size; ++i)
			{
				m_data[i] += other.m_data[i];
			}
			return *this;
		}

		Image operator+(Image const& other)const
		{
			Image<T> res = *this;
			return res += other;
		}

		T mean()const
		{
			T sum = 0;
			for (size_t i=0; i<size(); ++i)
			{
				sum += m_data[i];
			}
			return sum / size();
		}


		static T eqm(Image<T> const&  a, Image<T> const& b)
		{
			assert(a.size() == b.size());
			T sum = 0;
			for (size_t i = 0; i < a.size(); ++i)
			{
				T dif = a.m_data[i] - b.m_data[i];
				sum += dif * dif;
			}
			return sum / a.size();
		}

		static T etm(Image<T> const& a, Image<T> const& b)
		{
			assert(a.size() == b.size());
			T sum = 0;
			for (size_t i = 0; i < a.size(); ++i)
			{
				T dif = a.m_data[i] - b.m_data[i];
				sum += dif * dif * dif;
			}
			return sum / a.size();
		}

		template <class out_t>
		out_t& print(out_t& out)const
		{
			for (size_t i = 0; i < m_width; ++i)
			{
				for (size_t j = 0; j < m_height; ++j)
				{
					out << (*this)[i][j] << " ";
				}
				out << std::endl;
			}
			return out;
		}

		template <class INT>
		__forceinline bool inBounds(Math::Vector<INT, 2> const& px)const noexcept 
		{
			return px[0] >= 0 && px[0] < m_width && px[1] >= 0 && px[1] < m_height;
		}
	};



}