#pragma once

/*
#include <Image/Image.h>
#include <Geometry/RGBColor.h>
#include <atomic>
#include <cassert>
#include <settings.h>
#include <Math/Vectorf.h>

namespace Math
{
	
	double sum(arma::colvec const& vec)
	{
		double res = 0;
		for (double d : vec)
		{
			res += d;
		}
		return res;
	}


	//for BDPT only
	//Implementation of the direct estimator (the progressive will come after that) 
	class OptiMISSolver
	{
	protected:

		using uint = unsigned int;
		static constexpr uint COLOR_SAMPLES = 3;

		//maybe keep track of the number of samples for each technique, it might become more releveant with the contributions of light tracing or delta (maybe)

		// size of the std::vector of the tech matrix, we need to store only half of it since it symetrical
		uint techHalfSize(uint N)const
		{
			return ((N + 1) * N) / 2;
		}

		uint techIndex(uint i, uint j)const
		{
			if (i < j) return j * m_num_tech + i;
			else return i * m_num_tech + j;
		}

		uint m_num_tech;

		//std::vector<std::atomic<double>> m_technique_matrix;

		//stores the contribution vector for each color channel
		//std::vector<std::atomic<double>> m_contribution_vector;

		

		//contain the computed estimate
		//+ is pre allocted
		arma::mat m_A;
		arma::colvec m_b[COLOR_SAMPLES];

		mutable arma::colvec m_alpha[COLOR_SAMPLES];

	public:

		OptiMISSolver(uint N) :
			m_num_tech(N),
			//m_technique_matrix(techHalfSize(N)),
			//m_contribution_vector(COLOR_SAMPLES* N),
			m_A(N, N)
		{
			m_A.fill(0);
			for (uint s = 0; s < COLOR_SAMPLES; ++s)
			{
				m_b[s] = arma::colvec(N);
				m_b[s].fill(0);
				m_alpha[s] = arma::colvec(N);
				m_alpha[s].fill(0);
			}
		}

		//////////////////////////////////////////////////////////////////////////\\
		// the raw estimate is just f(x)										   \\
		// wights is the regular balance MIS weights							   //
		// sum_pdf is the sum of the all the pdf the sample (for each technique)  //
		// The estimate is added to the vectors of atomic float					 //
		//////////////////////////////////////////////////////////////////////////
		void addEstimate(const Geometry::RGBColor& raw_estimate, const double * weights, double sum_pdf)
		{
			for (uint i = 0; i < m_num_tech; ++i)
			{
				for (uint j = 0; j <= i; ++j)
				{
					double tmp = weights[i] * weights[j];
					m_A(i, j) += tmp;
					m_A(j, i) += tmp;
				}

				for (uint s = 0; s < COLOR_SAMPLES; ++s)
				{
					m_b[s](i) += raw_estimate[s] * weights[i] / sum_pdf;
				}
			}
		}

		void loop()
		{
			//unused in the direct estimator
			//used in the progressive esimator
		}

		Geometry::RGBColor solve()const
		{
			arma::mat Aplus = arma::pinv(m_A);
			Geometry::RGBColor res = 0;
			for (uint s = 0; s < COLOR_SAMPLES; ++s)
			{
				m_alpha[s] = Aplus * m_b[s];
				double chan = sum(m_alpha[s]);
				res[s] = chan;
			}
			return res;
		}

	};
	

	
	template<class OMISS>
	class ImageSolver
	{
	protected:

		using uint = unsigned int;

		uint m_num_techs;

		std::vector<OMISS> m_solvers;

		mutable Image::Image<Geometry::RGBColor>* m_target;

		uint idx(uint w, uint h)const noexcept 
		{
			return w * m_target->height() + h;
		}

	public:

		ImageSolver(Image::Image<Geometry::RGBColor>* target, uint N):
			m_num_techs(N),
			m_solvers(target->width() * target->height(), N),
			m_target(target)
		{
			assert(m_target != nullptr);
		}

		void solveAll()const
		{
			OMP_PARALLEL_FOR
			for (uint i = 0; i < m_target->width(); ++i)
			{
				for (uint j = 0; j < m_target->height(); ++j)
				{
					m_target->operator[](i)[j] += m_solvers[idx(i, j)].solve();
				}
			}
		}

		void addEstimate(Vector<int, 2> const & px, Geometry::RGBColor const& raw_estimate, double* weights, double sum_pdf)
		{
			if (m_target->inBounds(px))
			{
				m_solvers[idx(px[0], px[1])].addEstimate(raw_estimate, weights, sum_pdf);
			}
			else
			{

			}
		}

		//uv: normalized position of the pixel € [0, 1[²
		void addEstimate(Vector2f const& uv, Geometry::RGBColor const& raw_estimate, const double* weights, double sum_pdf)
		{
			Vector<int, 2> px = { uv[0] * m_target->width(), uv[1] * m_target->height() };
			addEstimate(px, raw_estimate, weights, sum_pdf);
		}

		void loop()
		{
			OMP_PARALLEL_FOR
				for (uint i = 0; i < m_target->width(); ++i)
				{
					for (uint j = 0; j < m_target->height(); ++j)
					{
						m_solvers[idx(i, j)].loop();
					}
				}
		}

	};
	
}

*/
