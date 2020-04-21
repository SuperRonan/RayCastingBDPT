#pragma once

#include <cassert>

template <typename T>
class Merger
{
protected:

	T* m_left = nullptr, * m_right = nullptr;
	size_t m_left_size = 0, m_right_size = 0;
	bool m_left_reversed = false, m_right_reversed = false;

	////////////////////////////////////////////////
	//   ----------------------------------------------------------------------
	//	|left[0]|left[1]|...|left[left_size-1]|right[0]|...|right[right_size-1]|
	//	 ----------------------------------------------------------------------
	// and the sub arrays can be reversed if needed
	////////////////////////////////////////////////

public:

	Merger()
	{}

	Merger(const Merger& other) = default;

	Merger(T* left, size_t ls, T* right, size_t rs, bool left_reverse = false, bool right_reverse = false) :
		m_left(left),
		m_right(right),
		m_left_size(ls),
		m_right_size(rs),
		m_left_reversed(left_reverse),
		m_right_reversed(right_reverse)
	{}

	T& operator[](size_t index)
	{
		assert(index >= 0);
		assert(index < size());
		if (index < m_left_size)
		{
			// Go fetch in the left sub array
			if (m_left_reversed)
			{
				return m_left[m_left_size - index - 1];
			}
			else
			{
				return m_left[index];
			}
		}
		else
		{
			// Go fetch in the right sub array
			index = index - m_left_size;
			assert(index < m_right_size);
			if (m_right_reversed)
			{
				return m_right[m_right_size - index - 1];
			}
			else
			{
				return m_right[index];
			}
		}
	}

	T const& operator[](size_t index)const
	{
		assert(index >= 0);
		assert(index < size());
		if (index < m_left_size)
		{
			// Go fetch in the left sub array
			if (m_left_reversed)
			{
				return m_left[m_left_size - index - 1];
			}
			else
			{
				return m_left[index];
			}
		}
		else
		{
			// Go fetch in the right sub array
			index = index - m_left_size;
			assert(index < m_right_size);
			if (m_right_reversed)
			{
				return m_right[m_right_size - index - 1];
			}
			else
			{
				return m_right[index];
			}
		}
	}

	size_t size()const
	{
		return m_left_size + m_right_size;
	}

};
