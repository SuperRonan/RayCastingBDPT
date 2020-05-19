#pragma once

#include <Image/Image.h>
#include <Geometry/RGBColor.h>
#include <omp.h>
#include <cassert>
#include <Math\Vector.h>

#define PRINT(var) std::cout << #var << ": " << var << std::endl;

template <bool MAJOR>
class ImageEstimator
{
protected:

	const int m_numtechs;

	const int m_width, m_height;

	int PixelTo1D(int i, int j)const
	{
		return Image::Image<int, MAJOR>::index(i, j, m_width, m_height);
	}

public:

	ImageEstimator(int N, int width, int height) :
		m_numtechs(N),
		m_width(width),
		m_height(height)
	{}

	virtual void setOverSample(int techIndex, int n)
	{}

	virtual void addEstimate(Geometry::RGBColor const& balanceEstimate, double* balanceWeights, int tech_index, Math::Vector2f const& uv) = 0;

	virtual void addOneTechniqueEstimate(Geometry::RGBColor const& balanceEstimate, int tech_index, Math::Vector2f const& uv) = 0;

	virtual void addZeroEstimate(int tech_index, Math::Vector2f const& uv) = 0;

	virtual void loop() = 0;

	virtual void solve(Image::Image<Geometry::RGBColor, MAJOR>& res, int iterations) = 0;

	virtual void debug(int iterations, bool col_sum, bool matrix, bool vec, bool alpha)const
	{}
};

template <bool MAJOR>
class BalanceEstimatorImage: public ImageEstimator<MAJOR>
{
protected:

	Image::Image<Geometry::RGBColor, MAJOR> m_image;

public:

	BalanceEstimatorImage(int N, int width, int height) :
		ImageEstimator(N, width, height),
		m_image(width, height)
	{
		m_image.fill(0);
	}


	virtual void addEstimate(Geometry::RGBColor const& balanceEstimate, double* balanceWeights, int tech_index, Math::Vector2f const& uv) override
	{
		Math::Vector<int, 2> pixel = { uv[0] * m_width, uv[1] * m_height };
		m_image(pixel) += balanceEstimate;
	}

	virtual void addZeroEstimate(int tech_index, Math::Vector2f const& uv) override
	{}

	virtual void loop() override
	{}

	virtual void addOneTechniqueEstimate(Geometry::RGBColor const& balanceEstimate, int tech_index, Math::Vector2f const& uv) override
	{
		Math::Vector<int, 2> pixel = { uv[0] * m_width, uv[1] * m_height };
		m_image(pixel) += balanceEstimate;
	}

	virtual void solve(Image::Image<Geometry::RGBColor, MAJOR>& res, int iterations) override
	{
		Parallel::ParallelFor(0, m_width * m_height,
			[&](int pixel)
			{
				res[pixel] += m_image[pixel] / iterations;
			});
	}
};

template <bool MAJOR>
class DirectEstimatorImage: public ImageEstimator<MAJOR>
{
protected:

	using MatrixT = arma::mat;
	using VectorT = arma::vec;
	using Float = double;
	using AtomicUInt = unsigned int;// std::atomic<unsigned int>;
	using AtomicFloat = Float;// std::atomic<Float>;

	const int msize;

	size_t matTo1D(int row, int col) const {
		// return row * numTechs + col;
		if (col > row) std::swap(row, col);
		return (row * (row + 1)) / 2 + col;
	}
	
	struct PixelData
	{
		AtomicFloat* techMatrix;
		// The three channels are one after the other [RRRRRRRR GGGGGGGG BBBBBBBBB]
		AtomicFloat* contribVector;
		AtomicUInt* sampleCount;

		PixelData(AtomicFloat* tm, AtomicFloat* cv, AtomicUInt* sc) :
			techMatrix(tm),
			contribVector(cv),
			sampleCount(sc)
		{}
	};

	PixelData getPixelData(int index)
	{
		const size_t offset = index * m_pixel_data_size;
		char* address = (m_data.data() + offset);
		PixelData res((AtomicFloat*)address, (AtomicFloat*)(address + m_vector_ofsset), (AtomicUInt*)+(address + m_counter_offset));
		return res;
	}

	PixelData getPixelData(Math::Vector2f const& uv)
	{
		const Math::Vector<int, 2> pixel(uv[0] * m_width, uv[1] * m_height);
		const int index = PixelTo1D(pixel[0], pixel[1]);
		return getPixelData(index);
	}

	PixelData getPixelData(int index)const 
	{
		const size_t offset = index * m_pixel_data_size;
		const char* address = (m_data.data() + offset);
		PixelData res((AtomicFloat*)address, (AtomicFloat*)(address + m_vector_ofsset), (AtomicUInt*)+(address + m_counter_offset));
		return res;
	}

	PixelData getPixelData(Math::Vector2f const& uv)const 
	{
		const Math::Vector<int, 2> pixel(uv[0] * m_width, uv[1] * m_height);
		const int index = PixelTo1D(pixel[0], pixel[1]);
		return getPixelData(index);
	}

	std::string debugFolder()const
	{
		return RESULT_FOLDER + "debug/";
	}

	// In Byte
	const unsigned int m_pixel_data_size;
	const unsigned int m_vector_ofsset;
	const unsigned int m_counter_offset;

	static_assert(sizeof(char) == 1);
	std::vector<char> m_data;

	//describes how many times more samples each technique has
	std::vector<unsigned int> m_over_samples;


public:

	DirectEstimatorImage(int N, int width, int height):
		ImageEstimator(N, width, height),
		msize(N * (N+1) / 2),
		m_pixel_data_size(msize * sizeof(AtomicFloat) + 3 * m_numtechs * sizeof(AtomicFloat) + m_numtechs * sizeof(AtomicUInt)),
		m_vector_ofsset(msize * sizeof(AtomicFloat)),
		m_counter_offset(msize * sizeof(AtomicFloat) + 3 * m_numtechs * sizeof(AtomicFloat))
	{
		int res = width * height;
		m_data = std::vector<char>(res * m_pixel_data_size, (char)0);

		m_over_samples = std::vector<unsigned int>(m_numtechs, 1);
	}

	virtual void setOverSample(int techIndex, int n) override
	{
		m_over_samples[techIndex] = n;
	}

	void checkSample(Geometry::RGBColor const& balanceEstimate, Float* balanceWeights, int tech_index)
	{
		for (int i = 0; i < m_numtechs; ++i)
		{
			Float weight = balanceWeights[i];
			if (weight < 0 || std::isnan(weight) || std::isinf(weight))
			{
				std::cout << weight << std::endl;
				__debugbreak();
			}
		}
	}

	virtual void addEstimate(Geometry::RGBColor const& balanceEstimate, Float* balanceWeights, int tech_index, Math::Vector2f const& uv) override
	{
		PixelData data = getPixelData(uv);
		++data.sampleCount[tech_index];
		for (int i = 0; i < m_numtechs; ++i)
		{
			for (int j = 0; j <= i; ++j)
			{
				const int mat_index = matTo1D(i, j);
				Float tmp = balanceWeights[i] * balanceWeights[j];
				data.techMatrix[mat_index] = data.techMatrix[mat_index] + tmp;
			}
		}
		if (!balanceEstimate.isBlack())
		{
			for (int k = 0; k < 3; ++k)
			{
				AtomicFloat* vector = data.contribVector + k * m_numtechs;
				for (int i = 0; i < m_numtechs; ++i)
				{
					double tmp = balanceWeights[i] * balanceEstimate[k];
					vector[i] = vector[i] + tmp;
				}
			}
		}
	}

	
	virtual void addZeroEstimate(int tech_index, Math::Vector2f const& uv) override
	{
		PixelData data = getPixelData(uv);
		++data.sampleCount[tech_index];
		int mat_index = matTo1D(tech_index, tech_index);
		data.techMatrix[mat_index] = data.techMatrix[mat_index] + 1.0;
	}

	virtual void addOneTechniqueEstimate(Geometry::RGBColor const& balanceEstimate, int tech_index, Math::Vector2f const& uv) override
	{
		PixelData data = getPixelData(uv);
		++data.sampleCount[tech_index];
		int mat_index = matTo1D(tech_index, tech_index);
		data.techMatrix[mat_index] = data.techMatrix[mat_index] + 1.0;
		for (int k = 0; k < 3; ++k)
		{
			AtomicFloat* vector = data.contribVector + k * m_numtechs;
			vector[tech_index] = vector[tech_index] + balanceEstimate[k];
		}
	}

	virtual void loop() override
	{}

	inline void fillMatrix(MatrixT& matrix, PixelData const& data, int iterations)const
	{
		for (int i = 0; i < m_numtechs; ++i)
		{
			for (int j = 0; j < i; ++j)
			{
				Float elem = data.techMatrix[matTo1D(i, j)];
				if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
				matrix(i, j) = elem;
				matrix(j, i) = elem;
			}
			matrix(i, i) = data.techMatrix[matTo1D(i, i)];
			// Unsampled samples
			size_t expected = m_over_samples[i] * iterations;
			size_t actually = data.sampleCount[i];
			matrix(i, i) += (Float)(expected - actually);
		}
	}

	virtual void saveColSum(int iterations)const
	{
		Image::Image<double, MAJOR> img(m_width * m_numtechs, m_height);
		int resolution = m_width * m_height;
		std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
		Parallel::ParallelFor(0, resolution, [&](int pixel)
			{
				int tid = Parallel::tid();
				Math::Vector<int, 2> indices = Image::Image<double, MAJOR>::indices(pixel, m_width, m_height);
				PixelData data = getPixelData(pixel);
				MatrixT& matrix = matrices[tid];
				fillMatrix(matrix, data, iterations);
				for (int i = 0; i < m_numtechs; ++i)
				{
					size_t expected = m_over_samples[i];
					double sum = 0;
					for (int j = 0; j < m_numtechs; ++j)
					{
						sum += matrix(i, j);
					}
					sum /= iterations;
					double diff = (expected - sum) / expected;
					img(indices[0] * m_numtechs + i, indices[1]) = diff;
				}
			});
		Image::ImWrite::writeEXR(debugFolder() + "OptiMIScolSum" + std::to_string(m_numtechs) + ".exr", img);
	}

	virtual void saveMatrices(int iterations)const
	{
		Image::Image<double, MAJOR> img(m_width * m_numtechs, m_height*m_numtechs);
		int resolution = m_width * m_height;
		std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
		Parallel::ParallelFor(0, resolution, [&](int pixel)
			{
				int tid = Parallel::tid();
				Math::Vector<int, 2> indices = Image::Image<double, MAJOR>::indices(pixel, m_width, m_height);
				PixelData data = getPixelData(pixel);
				MatrixT& matrix = matrices[tid];
				fillMatrix(matrix, data, iterations);
				for (int i = 0; i < m_numtechs; ++i)
				{
					for (int j = 0; j < m_numtechs; ++j)
					{
						img(indices[0] * m_numtechs + i, indices[1] * m_numtechs + j) = matrix(i, j) / iterations;
					}
				}
			});
		Image::ImWrite::writeEXR(debugFolder() + "OptiMISMatrices" + std::to_string(m_numtechs) + ".exr", img);
	}

	virtual void saveVectors(int iterations)const
	{
		Image::Image<Geometry::RGBColor, MAJOR> img(m_width * m_numtechs, m_height);
		int resolution = m_width * m_height;
		Parallel::ParallelFor(0, resolution, [&](int pixel)
			{
				int tid = Parallel::tid();
				Math::Vector<int, 2> indices = Image::Image<double, MAJOR>::indices(pixel, m_width, m_height);
				PixelData data = getPixelData(pixel);
				for (int i = 0; i < m_numtechs; ++i)
				{
					for (int k = 0; k < 3; ++k)
					{
						img(indices[0] * m_numtechs + i, indices[1])[k] = (data.contribVector + k*m_numtechs)[i] / iterations;
					}
				}
			});
		Image::ImWrite::writeEXR(debugFolder() + "OptiMISVectors" + std::to_string(m_numtechs) + ".exr", img);
	}

	virtual void saveAlphas(int iterations)const
	{
		Image::Image<Geometry::RGBColor, MAJOR> img(m_width * m_numtechs, m_height);
		VectorT MVector(m_numtechs);
		std::copy(m_over_samples.begin(), m_over_samples.end(), MVector.begin());
		std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
		std::vector<VectorT> vectors = Parallel::preAllocate(VectorT(m_numtechs));
		int resolution = m_width * m_height;
		Parallel::ParallelFor(0, resolution, [&](int pixel)
			{
				int tid = Parallel::tid();
				Math::Vector<int, 2> indices = Image::Image<double, MAJOR>::indices(pixel, m_width, m_height);
				MatrixT& matrix = matrices[tid];
				VectorT& vector = vectors[tid];

				PixelData data = getPixelData(pixel);

				bool matrix_done = false;

				for (int k = 0; k < 3; ++k)
				{
					bool isZero = true;
					const AtomicFloat* contribVector = data.contribVector + k * m_numtechs;
					for (int i = 0; i < m_numtechs; ++i)
					{
						Float elem = contribVector[i];
						if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
						vector[i] = elem;
						isZero = isZero & (contribVector[i] == 0);
					}

					if (!isZero)
					{
						if (!matrix_done)
						{
							fillMatrix(matrix, data, iterations);
							matrix = pinv(matrix);
							matrix_done = true;
						}
						// Now vector is alpha
						vector = matrix * vector;
						for (int i = 0; i < m_numtechs; ++i)
						{
							vector[i] *= m_over_samples[i];
							img(indices[0] * m_numtechs + i, indices[1])[k] = vector[i];
						}
					}
				}
				
			});
		Image::ImWrite::writeEXR(debugFolder() + "OptiMISAlphas" + std::to_string(m_numtechs) + ".exr", img);
	}

	virtual void debug(int iterations, bool col_sum, bool matrix, bool vec, bool alpha)const override
	{
		if (col_sum)
			saveColSum(iterations);
		if (matrix)
			saveMatrices(iterations);
		if (vec)
			saveVectors(iterations);
		if (alpha)
			saveAlphas(iterations);
	}

	virtual void solve(Image::Image<Geometry::RGBColor, MAJOR>& res, int iterations) override
	{
		VectorT MVector(m_numtechs);
		std::copy(m_over_samples.begin(), m_over_samples.end(), MVector.begin());
		std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
		std::vector<VectorT> vectors = Parallel::preAllocate(VectorT(m_numtechs));
		int resolution = m_width * m_height;
		Parallel::ParallelFor(0, resolution, [&](int pixel)
		{
			int tid = Parallel::tid();
			MatrixT& matrix = matrices[tid];
			VectorT& vector = vectors[tid];

			Geometry::RGBColor color = 0;

			PixelData data = getPixelData(pixel);

			bool matrix_done = false;

			for (int k = 0; k < 3; ++k)
			{
				bool isZero = true;
				const AtomicFloat* contribVector = data.contribVector + k * m_numtechs;
				for (int i = 0; i < m_numtechs; ++i)
				{
					Float elem = contribVector[i];
					if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
					vector[i] = elem;
					isZero = isZero & (contribVector[i] == 0);
				}

				if (!isZero)
				{
					if (!matrix_done)
					{
						fillMatrix(matrix, data, iterations);
						matrix = pinv(matrix);
						matrix_done = true;
					}
					// Now vector is alpha
					vector = matrix * vector;

					double estimate = arma::dot(vector, MVector);
					color[k] = estimate;
				}
			}
			res[pixel] += color;
		});
	}
};