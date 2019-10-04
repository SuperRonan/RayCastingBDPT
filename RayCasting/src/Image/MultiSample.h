#pragma once


namespace Image
{
	/////////////////////////////////
	//A class representing a "collection" of samples
	// TODO try another implemetation
	/////////////////////////////////
	template <class T, class uint=unsigned int>
	class MultiSample
	{
	protected:

		T m_total;
		uint m_n;

	public:

		MultiSample() :
			m_total(0),
			m_n(0)
		{}

		MultiSample(T const& init) :
			m_total(init),
			m_n(1)
		{}

		uint n()const
		{
			return m_n;
		}

		T const& total()const
		{
			return m_total;
		}

		T mean()const
		{
			return total() / n();
		}

		void add(T const& sample)
		{
//#pragma omp atomic update
			m_total += sample;
//#pragma omp atomic update
			++m_n;
		}

		void reset()
		{
			m_total = 0;
			m_n = 0;
		}

		operator T() const
		{
			return mean();
		}

	};
}