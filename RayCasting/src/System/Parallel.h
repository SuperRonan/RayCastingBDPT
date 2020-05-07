#pragma once

#include <settings.h>
#include <omp.h>
#include <vector>

constexpr auto DEFAULT_INCREMENTATION = [](int& i) {++i; };


class Parallel
{
protected:


public:

	static void init()
	{
		
	}

	static void setNumThread(int n)
	{
		omp_set_num_threads(std::max(1, std::min(getCPUThreads(), n)));
	}

	static int getNumThread()
	{
		return omp_get_max_threads();
	}

	static int getCPUThreads()
	{
		return omp_get_num_procs();
	}

	static int getCurrentNumThread()
	{
		return omp_get_num_threads();
	}

	template <class T>
	__forceinline static std::vector<T> preAllocate(T const& t=T())
	{
		return std::vector<T>(getNumThread(), t);
	}

	static int tid()
	{
		return omp_get_thread_num();
	}

	template <class Function>
	__forceinline static void ParallelFor(int min, int max, Function const& function)
	{
		OMP_PARALLEL_FOR
			for (int i = min; i < max; ++i)
			{
				function(i);
			}
	}
};