#pragma once

#include <random>
namespace Math
{

#define DISCRET_SAMPLER_UNBIAS

	//////////////////////////////////
	//A sampler that generates numbers uniformaly
	//////////////////////////////////
	class Sampler
	{
	protected:

		std::mt19937_64 m_rng;

		using ull = uint_fast64_t;

	public:

		Sampler(ull seed=0):
			m_rng(seed)
		{}

		ull min()const
		{
			return m_rng.min();
		}

		ull max()const
		{
			return m_rng.max();
		}


		///////////////////////////
		//generates uniformaly an ull in the range [min(), max()]
		///////////////////////////
		ull generate()
		{
			return m_rng();
		}

		///////////////////////////
		//generates uniformaly an ull in the range [min, max]
		///////////////////////////
		ull generate(ull min, ull max)
		{
			if (min == max)	return min;
#ifdef DISCRET_SAMPLER_UNBIAS
			ull len = max - min + 1;
			ull dev = Sampler::max() / len;
			ull prod = dev * len;
			while (1)
			{
				ull n = generate();
				if (n < prod)
				{
					return (n / dev) + min;
				}
			}
#else
			return generate<double>(double(min), double(max));
#endif
		}

		///////////////////////////
		//generates uniformaly an ull in the range [0, 1]
		///////////////////////////
		template <class floot=double>
		floot generateContinuous()
		{
			return floot(generate()) / floot(max());
		}

		///////////////////////////
		//generates uniformaly an ull in the range [min, max]
		///////////////////////////
		template <class floot=double>
		floot generateContinuous(floot min, floot max)
		{
			return generateContinuous<floot>() * (max - min) + min;
		}

		template <class floot=double>
		floot generateStratified(unsigned int i, unsigned int div)
		{
			floot xi = generateContinuous<floot>();
			return (i + xi) / floot(div);
		}

		template <class floot = double>
		floot generateStratified(floot min, floot max, unsigned int i, unsigned int div)
		{
			return generateStratified<floot>(i, div) * (max - min) + min;
		}

		void seed(ull s)
		{
			m_rng.seed(s);
		}
	};
}