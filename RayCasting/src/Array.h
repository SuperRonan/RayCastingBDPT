#pragma once

#include <cassert>
#include <algorithm>

namespace utils
{
	template<class T>
	class Array
	{
	protected:
		size_t m_size;
		T * m_data;

	public:

		Array():
			m_size(0),
			m_data(nullptr)
		{}

		Array(size_t p_size) :
			m_size(p_size),
			m_data(new T[p_size])
		{
			assert(m_size > 0);
		}

		Array(size_t p_size, T const& default_value):
			m_size(p_size),
			m_data(new T[p_size])
		{
			assert(m_size > 0);
			std::fill(m_data, m_data + m_size, default_value);
		}

		Array(Array const& other) :
			m_size(other.m_size),
			m_data(nullptr)
		{
			if (m_size > 0)
			{
				m_data = new T[m_size];
				std::copy(other.m_data, other.m_data + m_size, m_data);
			}
		}

		~Array()
		{
			if (m_size > 0)
			{
				assert(m_data != nullptr);
				delete[] m_data;
			}
		}


	};
}