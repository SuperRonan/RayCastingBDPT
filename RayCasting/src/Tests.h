#pragma once

#include <iostream>
#include <Image/Image.h>
#include <settings.h>
#include <System/Parallel.h>
#include <Geometry/RGBColor.h>
#include <Image/ImWrite.h>
#include <Math/Sampler.h>
#include <vector>
#include <atomic>

using RGBColor = Geometry::RGBColor;


template <typename floot>
void testPrecision()
{
	std::cout.precision(20);
	floot f = 1;
	for (unsigned int i = 0; i < 1024; ++i)
	{
		floot tmp = f + 1;
		std::cout << i << " : ";
		if (f != tmp)
		{
			std::cout << "Ok for " << f << std::endl;
		}
		else
		{
			std::cout << "WRONG for " << f << std::endl;
			break;
		}
		f *= 2;
	}
}

void testWaveLength()
{
	Image::Image<RGBColor> img(290, 50);
	Parallel::ParallelFor(0, img.width(), [&img](const int i)
		{
			double l = i + 380;
			RGBColor c = RGBColor::fromWaveLength(l);
			for (int j = 0; j < img.height(); ++j)
			{
				img(i, j) = c;
			}
		});
	Image::ImWrite::write(img);

	img.fill(img.mean());
	Image::ImWrite::write(img);
}

void testRIS(size_t seed = 0)
{
	Math::Sampler sampler(seed);
	struct Sample
	{
		double x, pdf;
		Sample(double x, double pdf) : x(x), pdf(pdf) {};
		Sample() = default;
		operator double() { return x; };
	};
	struct Candidate
	{
		double pdf, target, f;
		double x;
		Candidate(double x, double p, double t, double f) :pdf(p), target(t), f(f), x(x) {};
		Candidate() = default;
		double w()const { return target / pdf; };
	};
	const int M = 2;
	std::vector<Candidate> candidates(M);

	const auto function = [](double x) {return x; };
	const auto sourcePdf = [](double x) {return 1; };
	const auto targetPdf = [](double x) {return x; };
	const auto sampleSource = [&](Math::Sampler& sampler) {double x = sampler.generateContinuous<double>(); return Sample(x, sourcePdf(x)); };
	const auto sampleRIS = [&](Math::Sampler& sampler)
	{
		double sum = 0;
		for (int i = 0; i < M; ++i)
		{
			Sample x = sampleSource(sampler);
			candidates[i] = Candidate(x, x.pdf, targetPdf(x), function(x));
			sum += candidates[i].w();
		}

		double xi = sampler.generateContinuous<double>(0, sum);
		double partial_sum = 0;
		Sample res(0, 0);
		for (int i = 0; i < M; ++i)
		{
			if (xi >= partial_sum && xi < (candidates[i].w() + partial_sum))
			{
				res = Sample(candidates[i].x, candidates[i].target * M / sum);
			}
			partial_sum += candidates[i].w();
		}
		return res;
	};
	const auto RISAnalyticPdf = [&M](double x)
	{
		if (M == 2)
		{
			return 2 * x * (log((x + 1) / x));
		}

		assert(false);
	};
	const auto RISEstimatedPdf = [&](double x, Math::Sampler& sampler)
	{
		double x_w = targetPdf(x) / sourcePdf(x);
		double sum = x_w;
		for (int i = 1; i < M; ++i)
		{
			Sample x = sampleSource(sampler);
			sum += targetPdf(x) / x.pdf;
		}
		return targetPdf(x) * M / sum;
	};
	int K = 128;
	int N = 128;
	double final_estimate = 0;
	for (int k = 0; k < K; ++k)
	{
		double integral = 0;
		for (int i = 0; i < N; ++i)
		{
			Sample x = sampleRIS(sampler);
			double ris_analytic_pdf = RISAnalyticPdf(x);
			//std::cout << "RIS pdf: " << x.pdf << ", RIS Analytic pdf: " << ris_analytic_pdf << std::endl;
			integral += function(x) / x.pdf;
		}
		std::cout << "integral: " << integral / N << std::endl;
		final_estimate += integral;
	}
	std::cout << "Final Integral estimate: " << final_estimate / (N * K) << std::endl;
}

template <class Float, bool second_iteration = false>
Float fastInvSqrt(Float f)
{
	if constexpr (std::is_same<Float, float>::value)
	{
		int32_t i;
		float x2;
		x2 = f * 0.5f;
		i = *(int32_t*)&f;
		i = 0x5f3759df - (i >> 1); // wtf
		f = *(float*)&i;
		f *= 1.5f - (x2 * f * f);
		if constexpr (second_iteration)
			f *= 1.5f - (x2 * f * f);
		return f;
	}
	else //if constexpr (std::is_same<Float, double>::value)
	{
		int64_t i;
		double x2;
		x2 = f * 0.5f;
		i = *(int64_t*)&f;
		i = 0x5FE6EC85E7DE30DA - (i >> 1); // wtf
		f = *(double*)&i;
		f *= 1.5f - (x2 * f * f);
		if constexpr (second_iteration)
			f *= 1.5f - (x2 * f * f);
		return f;
	}
}

template <class Float>
void testInvSqrt()
{
	Math::Sampler sampler;

	for (int i = 0; i < 100; ++i)
	{
		Float f = sampler.generateContinuous<Float>(0, 2);
		Float exact1 = Float(1.0) / (std::sqrt(f));
		Float exact2 = std::sqrt(Float(1.0) / (f));
		Float exact = (exact1 + exact2) * Float(0.5);

		Float fast1 = fastInvSqrt<Float, false>(f);
		Float fast2 = fastInvSqrt<Float, true>(f);

		std::cout << "f: " << f << std::endl;
		std::cout << "exact: " << exact << std::endl;
		std::cout << "fast1: " << fast1 << std::endl;
		std::cout << "fast2: " << fast2 << std::endl;
	}

}

void testParallel()
{
	int nthread = 4 * 2 * 2;
	Parallel::setNumThreads(nthread);
	for (int iter = 0; iter < 100; ++iter)
	{
		int N = 10000;
		std::vector<std::atomic<double>> vec = { 0, 0 };
		Math::Sampler sampler;
		Parallel::ParallelFor(0, N, [&](int i)
			{
				int id = sampler.generate(0, 1);
				vec[id].store(vec[id] + 1);
			});
		double res = vec[0] + vec[1];
		if (res != N)
			std::cout << res << std::endl;
	}
}