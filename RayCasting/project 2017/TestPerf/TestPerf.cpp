#include <iostream>
#include <Image/Image.h>
#include <settings.h>
#include <System/Parallel.h>


double function(int i, int j)
{
	double res = 0;
	auto hasher = std::hash<int>();
	for (int n = 0; n < 15; ++n)
	{
		res += hasher(i + j + n);
	}
	return res;
}

void testPerf()
{
	int N = 2048 * 2 * 2 * 2;
	Image::Image<double, true> image(N, N), other(N, N);
	std::cout << "Image: row major" << std::endl;
	std::cout << "One thread" << std::endl;
	std::cout << "Write only" << std::endl;
	std::cout << "top loop: row" << std::endl;
	tic();
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			image(i, j) = function(i, j);
		}
	}
	toc();

	std::cout << "top loop: col" << std::endl;
	tic();
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			image(j, i) = function(i, j);
		}
	}
	toc();
	std::cout << "Read/Write" << std::endl;
	std::cout << "top loop: row" << std::endl;
	tic();
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			other(i, j) = image(i, j);
		}
	}
	toc();

	std::cout << "top loop: col" << std::endl;
	tic();
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			other(j, i) = image(j, i);
		}
	}
	toc();

	std::cout << "Now with multi threading" << std::endl;
	Parallel::init();
	Parallel::setNumThreads(Parallel::getCPUThreads());
	std::cout << "with " << Parallel::getNumThreads() << " threads" << std::endl;

	std::cout << "top mt loop: row" << std::endl;
	tic();
	OMP_PARALLEL_FOR
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			image(i, j) = function(i, j);
		}
	}
	toc();

	std::cout << "top mt loop: col" << std::endl;
	tic();
	OMP_PARALLEL_FOR
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			image(j, i) = function(i, j);
		}
	}
	toc();
	std::cout << "Read/Write" << std::endl;
	std::cout << "top mt loop: row" << std::endl;
	tic();
	OMP_PARALLEL_FOR
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			other(i, j) = image(i, j);
		}
	}
	toc();

	std::cout << "top mt loop: col" << std::endl;
	tic();
	OMP_PARALLEL_FOR
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			other(j, i) = image(j, i);
		}
	}
	toc();
}


int main()
{
	testPerf();
}

