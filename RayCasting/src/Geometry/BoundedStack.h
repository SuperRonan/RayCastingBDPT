#pragma once

#include "DirectionalLight.h"

namespace Geometry
{

	class GeometryBase;

	template <size_t N, class T>
	class BoundedStack
	{


	public:

		constexpr static size_t capacity = N;

	protected:

		T m_data[capacity];

		size_t m_size=0;



	public:

		void reset()
		{
			m_size = 0;
		}

		T const& operator[](int i)const
		{
			assert(i >= 0);
			assert(i < size());
			return m_data[i];
		}

		T & operator[](int i)
		{
			assert(i >= 0);
			assert(i < size());
			return m_data[i];
		}

		size_t size()const
		{
			return m_size;
		}

		bool empty()const
		{
			return m_size == 0;
		}

		bool full()const
		{
			return m_size == capacity;
		}

		void push(T const& l)
		{
			assert(!full());
			m_data[m_size] = l;
			++m_size;
		}

		T const& top()const
		{
			assert(!empty());
			return m_data[m_size-1];
		}

		T & top()
		{
			assert(!empty());
			return m_data[m_size - 1];
		}

		void pop()
		{
			assert(!empty());
			--m_size;
		}

		void grow()
		{
			assert(!full());
			++m_size;
		}

		void grow(size_t n)
		{
			assert(!full());
			m_size += n;
			assert(m_size <= capacity);
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

	};

	template <class T>
	using StackN = BoundedStack<12, T>;
}