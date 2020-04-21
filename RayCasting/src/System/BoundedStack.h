#pragma once

template <size_t N, class T>
class BoundedStack
{
public:

	constexpr static size_t s_capacity = N;

	size_t capacity()const
	{
		return s_capacity;
	}

protected:

	T m_data[s_capacity];

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

	void resize(size_t s)
	{
		assert(s <= capacity());
		m_size = s;
	}

	bool empty()const
	{
		return m_size == 0;
	}

	bool full()const
	{
		return m_size == capacity();
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

	const T* data()const
	{
		return begin();
	}

	T* data()
	{
		return begin();
	}

};

template <class T>
using StackN = BoundedStack<12, T>;

template <class T>
class BoundedStack<-1, T>
{
protected:

	unsigned int m_capacity;

public:

	size_t capacity()const
	{
		return m_capacity;
	}

protected:

	T * m_data;

	size_t m_size = 0;

public:

	BoundedStack() = delete;

	BoundedStack(T * data, unsigned int capacity):
		m_capacity(capacity),
		m_data(data),
		m_size(0)
	{}

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

	T& operator[](int i)
	{
		assert(i >= 0);
		assert(i < size());
		return m_data[i];
	}

	size_t size()const
	{
		return m_size;
	}

	void resize(size_t s)
	{
		assert(s <= capacity());
		m_size = s;
	}

	bool empty()const
	{
		return m_size == 0;
	}

	bool full()const
	{
		return m_size == capacity();
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
		return m_data[m_size - 1];
	}

	T& top()
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

	const T* data()const
	{
		return begin();
	}

	T* data()
	{
		return begin();
	}

};