#pragma once

#include <Image/Image.h>
#include <Geometry/RGBColor.h>
#include <omp.h>
#include <cassert>
#include <Math\Vector.h>

#define PRINT(var) std::cout << #var << ": " << var << std::endl;


////////////////////////////////////////////////////////////////////\\
// Base class for Optimal MIS image solvers							//
// Defines the 3 main functions splitting the estimator algorithms //
////////////////////////////////////////////////////////////////////
class OptimalSolverImage {

public:

	using RGBColor = Geometry::RGBColor;

	using Float = double;

	using MatrixT = arma::mat;
	using VectorT = arma::vec;
	using Float = double;
	using AtomicUInt = std::atomic<unsigned int>;
	using AtomicFloat = std::atomic<Float>;

	

protected:
    const int numTechs;
	const int spp;
    int width, height;

	std::vector<AtomicFloat> techMatrices;

	std::vector<AtomicUInt> LTSamples;

	// a vector of all the contribution vectors of all th pixels
	std::vector<AtomicFloat> contribVectors;

    int PixelTo1D(int x, int y) {
        return (x) +
               (y) * width;
    }

	AtomicFloat* getPixelContribVectors(int i, int j) {
		size_t offset = size_t(PixelTo1D(i, j)) * 3 * numTechs;

		return contribVectors.data() + offset;
		// return &contribVectors[ PixelTo1D(i, j) * Spectrum::nSamples *
		// numTechs];
	}

	AtomicFloat* getPixelTechMatrix(int i, int j) {
		int msize = numTechs * (numTechs + 1) / 2;
		size_t offset = size_t(PixelTo1D(i, j)) * msize;
		return techMatrices.data() + offset;
	}

	size_t matTo1D(int row, int col) const {
		// return row * numTechs + col;
		if (col > row) std::swap(row, col);
		return (row * (row + 1)) / 2 + col;
	}

public:

	const bool useDirect = true;
	const bool useLT;


    OptimalSolverImage(int numTechs, int width, int height, int spp, bool useLT)
        : numTechs(numTechs),
          width(width),
          height(height),
		  spp(spp),
		  useLT(useLT)
	{
		int msize = numTechs * (numTechs + 1) / 2;
		int res = width * height;

		techMatrices = std::vector<AtomicFloat>(msize * res);

		if (useLT)
		{
			LTSamples = std::vector<AtomicUInt>(res);
			std::fill(LTSamples.begin(), LTSamples.end(), 0);
		}

		contribVectors = std::vector<AtomicFloat>(res * numTechs * 3);
	}

	bool techIsLT(int techIndex)const
	{
		return useLT && (techIndex == (numTechs - 1));
	}

	int numberOfSamples(int techIndex)const
	{
		return techIsLT(techIndex) ? width * height : 1;
	}

    ///////////////////////////////////////////////////////////////////////////////////////////\\
	// Should be called every time a sample is drawn										    \\
	// Update the matrix and the vector with the drawn sample and the extra info provided	    //
    // Also usualy computes a part of the estimate <F>o when using the
    // progressive estimator   //
    ////////////////////////////////////////////////////////////////////////////////////////////
	void AddEstimate(
		const Math::Vector2f& p,
		const RGBColor& f,
		const Float* pdfs,
		double sumqi,
		int techIndex
	)
	{
		Math::Vector<int, 2> pPixel(p[0]*width, p[1]*height);
		{
			// Update the matrix of p
			AtomicFloat* techMatrix = getPixelTechMatrix(pPixel[0], pPixel[1]);
			for (int i = 0; i < numTechs; ++i) 
			{
				for (int j = 0; j <= i; ++j) 
				{
					const int mat_index = matTo1D(i, j);
					double tmp = pdfs[i] * pdfs[j] / (sumqi * sumqi);
					if (std::isnan(tmp))
						__debugbreak();
					techMatrix[mat_index] = techMatrix[mat_index] + tmp;
				}
			}
		}
		if (useLT && techIndex == (numTechs-1))
			LTSamples[PixelTo1D(pPixel[0], pPixel[1])]++;

		// Update only the contribution vector of the pixel
		// If we want to add filters, this is certainly where it should be done
		if (!f.isBlack()) {
			AtomicFloat* pixelContribVectors =
				getPixelContribVectors(pPixel[0], pPixel[1]);

			for (int k = 0; k < 3; ++k) {
				for (int i = 0; i < numTechs; ++i) {
					double tmp = f[k] * pdfs[i] / (sumqi * sumqi);
					pixelContribVectors[k * numTechs + i] = pixelContribVectors[k * numTechs + i] + tmp;
				}
			}
		}

		// Progressive estimator: update the estimate
		// Again, if we want to add filters, something will have to be done here
		if (!useDirect) {
			// TODO
			__debugbreak();
		}
	}

	void AddZeroEstimate(
		const Math::Vector2f& p,
		int techIndex
	)
	{
		double ni = numberOfSamples(techIndex);
		Math::Vector<int, 2> pPixel(p[0]*width, p[1]*height);
		{
			// Update the matrix of p
			AtomicFloat* techMatrix = getPixelTechMatrix(pPixel[0], pPixel[1]);

			double tmp = 1.0 / (ni * ni);
			techMatrix[matTo1D(techIndex, techIndex)] = techMatrix[matTo1D(techIndex, techIndex)] + tmp;
			
		}
		if (techIsLT(techIndex))
			LTSamples[PixelTo1D(pPixel[0], pPixel[1])]++;


		// Progressive estimator: update the estimate
		// Again, if we want to add filters, something will have to be done here
		if (!useDirect) {
			// TODO
			__debugbreak();
		}
	}


    ////////////////////////////////////////////////////////////////////////////////////\\
	// Should be called at the end, to get the final optimal result						 \\
	// Computes the last line after the loop in the estimators algorithms				 //
    // Solves the system and returns sun alpha for the direct estimator
    // // Returns the avg of the estimates already computed for the progressive
    // estimator //
    ////////////////////////////////////////////////////////////////////////////////////
	template <bool MAJOR>
	void DevelopFilm(Image::Image<RGBColor, MAJOR>* film, int numIterations)
	{
		if (useDirect) {
			// pre allocate the vector and the matrix

			std::vector<MatrixT> mats(omp_get_num_threads()+16, MatrixT(numTechs, numTechs));
			MatrixT LTfiller(numTechs, numTechs);
			std::vector<VectorT> vecs(omp_get_num_threads()+16, VectorT(numTechs));
			
			//std::vector<MatrixT> pinvs(omp_get_num_threads(), MatrixT(numTe;

			double nlt2 = Float(width * height) * Float(width * height);

			LTfiller.fill(0);
			LTfiller(numTechs - 1, numTechs - 1) = 1.0 / nlt2;


OMP_PARALLEL_FOR
			for (int x = 0; x < width; ++x)
			{
				int tid = omp_get_thread_num();
				MatrixT& mat = mats[tid];
				VectorT& vec = vecs[tid];
				for (int y = 0; y < height; ++y)
				{
					const AtomicFloat* pixelTechMatrix = getPixelTechMatrix(x, y);
					for (int i = 0; i < numTechs; ++i) 
					{
						for (int j = 0; j < i; ++j)
						{
							mat(i, j) = pixelTechMatrix[matTo1D(i, j)];
							mat(j, i) = pixelTechMatrix[matTo1D(i, j)];
							
							{
								if (std::isnan(mat(i, j)))
									__debugbreak();
								if (std::isnan(mat(j, i)))
									__debugbreak();
								if (std::isinf(mat(i, j)))
									__debugbreak();
								if (std::isinf(mat(j, i)))
									__debugbreak();
							}
						}
						mat(i, i) = pixelTechMatrix[matTo1D(i, i)];
						{
							if (std::isnan(mat(i, i)))
								__debugbreak();
							if (std::isinf(mat(i, i)))
								__debugbreak();
						}
					}

					if (useLT)
					{
						Float NLTZero = size_t(width) * size_t(height) * size_t(spp) - LTSamples[PixelTo1D(x, y)];

						mat += (LTfiller * NLTZero);

					}
					try
					{
						mat = arma::pinv(mat);
					}
					catch (std::exception const& e)
					{
						std::cout << e.what() << std::endl;
						mat = MatrixT(numTechs, numTechs);
						continue;
					}
					const AtomicFloat* pixelContribVectors = getPixelContribVectors(x, y);
					RGBColor estimate = 0;
					for (int k = 0; k < 3; ++k) {
						// fill the vector
						for (int i = 0; i < numTechs; ++i) {
							vec[i] = pixelContribVectors[k * numTechs + i];
						}
						
						vec = mat * vec;

						if (useLT)
						{
							//vec[numTechs - 1] *= 1 * spp;
						}

						estimate[k] = sum(vec);
						//estimate[k] = vec[vec.size() - 1];
					}
					(*film)(x, y) += estimate;

				}
			}
		}
		else {
			// TODO
		}
	}

    //////////////////////////////////////////////////////////////////////////////////\\
	// Used only for the progressive estimator										   \\
	// This should be called once all the samples have been drawn					   //
    // This is when the optimal estimator <F>o is comnputed
    // // This is also here that the linear system is solved every once in a
    // while		 // and the alpha vector is updated
    // //
    /////////////////////////////////////////////////////////////////////////////////
	void Loop()
	{
		if (!useDirect) {
			// TODO
		}
	}
};

class OptimalSolverImagePBRT {

public:

	using RGBColor = Geometry::RGBColor;

	using Float = double;
	using Film = Image::Image<RGBColor>;

	using MatrixT = arma::mat;
	using VectorT = arma::vec;
	using Float = double;
	using AtomicUInt = std::atomic<unsigned int>;
	using AtomicFloat = std::atomic<Float>;



protected:
	const int numTechs;
	const int spp;
	int width, height;

	std::vector<AtomicFloat> techMatrices;

	std::vector<AtomicUInt> LTSamples;

	// a vector of all the contribution vectors of all th pixels
	std::vector<AtomicFloat> contribVectors;

	int PixelTo1D(int x, int y) {
		return (x)+
			(y)* width;
	}

	AtomicFloat* getPixelContribVectors(int i, int j) {
		size_t offset = size_t(PixelTo1D(i, j)) * 3 * numTechs;

		return contribVectors.data() + offset;
		// return &contribVectors[ PixelTo1D(i, j) * Spectrum::nSamples *
		// numTechs];
	}

	AtomicFloat* getPixelTechMatrix(int i, int j) {
		int msize = numTechs * (numTechs + 1) / 2;
		size_t offset = size_t(PixelTo1D(i, j)) * msize;
		return techMatrices.data() + offset;
	}

	size_t matTo1D(int row, int col) const {
		// return row * numTechs + col;
		if (col > row) std::swap(row, col);
		return (row * (row + 1)) / 2 + col;
	}

public:

	const bool useDirect = true;
	const bool useLT;


	OptimalSolverImagePBRT(int numTechs, int width, int height, int spp, bool useLT)
		: numTechs(numTechs),
		width(width),
		height(height),
		spp(spp),
		useLT(useLT)
	{
		int msize = numTechs * (numTechs + 1) / 2;
		int res = width * height;

		techMatrices = std::vector<AtomicFloat>(msize * res);

		if (useLT)
		{
			LTSamples = std::vector<AtomicUInt>(res);
			std::fill(LTSamples.begin(), LTSamples.end(), 0);
		}

		contribVectors = std::vector<AtomicFloat>(res * numTechs * 3);
	}

	bool techIsLT(int techIndex)const
	{
		return useLT && (techIndex == (numTechs - 1));
	}

	int numberOfSamples(int techIndex)const
	{
		return techIsLT(techIndex) ? width * height : 1;
	}

	///////////////////////////////////////////////////////////////////////////////////////////\\
	// Should be called every time a sample is drawn										    \\
	// Update the matrix and the vector with the drawn sample and the extra info provided	    //
	// Also usualy computes a part of the estimate <F>o when using the
	// progressive estimator   //
	////////////////////////////////////////////////////////////////////////////////////////////
	void AddEstimate(
		const Math::Vector2f& p,
		const RGBColor& f,
		const Float* pdfs,
		double sumqi,
		int techIndex
	)
	{
		Math::Vector<int, 2> pPixel(p[0] * width, p[1] * height);
		{
			// Update the matrix of p
			AtomicFloat* techMatrix = getPixelTechMatrix(pPixel[0], pPixel[1]);
			for (int i = 0; i < numTechs; ++i)
			{
				for (int j = 0; j <= i; ++j)
				{
					const int mat_index = matTo1D(i, j);
					double tmp = (pdfs[i] * numberOfSamples(i) * pdfs[j] * numberOfSamples(j)) / (sumqi * sumqi);
					if (std::isnan(tmp))
						__debugbreak();
					techMatrix[mat_index] = techMatrix[mat_index] + tmp;
				}
			}
		}
		if (useLT && techIndex == (numTechs - 1))
			LTSamples[PixelTo1D(pPixel[0], pPixel[1])]++;

		// Update only the contribution vector of the pixel
		// If we want to add filters, this is certainly where it should be done
		if (!f.isBlack()) {
			AtomicFloat* pixelContribVectors =
				getPixelContribVectors(pPixel[0], pPixel[1]);

			for (int k = 0; k < 3; ++k) {
				for (int i = 0; i < numTechs; ++i) {
					double tmp = (f[k] * pdfs[i] * numberOfSamples(i)) / (sumqi * sumqi);
					pixelContribVectors[k * numTechs + i] = pixelContribVectors[k * numTechs + i] + tmp;
				}
			}
		}

		// Progressive estimator: update the estimate
		// Again, if we want to add filters, something will have to be done here
		if (!useDirect) {
			// TODO
			__debugbreak();
		}
	}

	void AddZeroEstimate(
		const Math::Vector2f& p,
		int techIndex
	)
	{
		double ni = numberOfSamples(techIndex);
		Math::Vector<int, 2> pPixel(p[0] * width, p[1] * height);
		{
			// Update the matrix of p
			AtomicFloat* techMatrix = getPixelTechMatrix(pPixel[0], pPixel[1]);

			double tmp = 1.0 / (ni * ni);
			techMatrix[matTo1D(techIndex, techIndex)] = techMatrix[matTo1D(techIndex, techIndex)] + tmp;

		}
		if (techIsLT(techIndex))
			LTSamples[PixelTo1D(pPixel[0], pPixel[1])]++;


		// Progressive estimator: update the estimate
		// Again, if we want to add filters, something will have to be done here
		if (!useDirect) {
			// TODO
			__debugbreak();
		}
	}


	////////////////////////////////////////////////////////////////////////////////////\\
	// Should be called at the end, to get the final optimal result						 \\
	// Computes the last line after the loop in the estimators algorithms				 //
	// Solves the system and returns sun alpha for the direct estimator
	// // Returns the avg of the estimates already computed for the progressive
	// estimator //
	////////////////////////////////////////////////////////////////////////////////////
	void DevelopFilm(Film* film, int numIterations)
	{
		if (useDirect) {
			// pre allocate the vector and the matrix

			std::vector<MatrixT> mats(omp_get_num_threads() + 16, MatrixT(numTechs, numTechs));
			MatrixT LTfiller(numTechs, numTechs);
			std::vector<VectorT> vecs(omp_get_num_threads() + 16, VectorT(numTechs));

			double nlt2 = Float(width * height) * Float(width * height);

			LTfiller.fill(0);
			LTfiller(numTechs - 1, numTechs - 1) = 1.0;


			OMP_PARALLEL_FOR
				for (int x = 0; x < width; ++x)
				{
					int tid = omp_get_thread_num();
					MatrixT& mat = mats[tid];
					VectorT& vec = vecs[tid];
					for (int y = 0; y < height; ++y)
					{
						const AtomicFloat* pixelTechMatrix = getPixelTechMatrix(x, y);
						for (int i = 0; i < numTechs; ++i)
						{
							for (int j = 0; j < i; ++j)
							{
								mat(i, j) = pixelTechMatrix[matTo1D(i, j)];
								mat(j, i) = pixelTechMatrix[matTo1D(i, j)];

								{
									if (std::isnan(mat(i, j)))
										__debugbreak();
									if (std::isnan(mat(j, i)))
										__debugbreak();
									if (std::isinf(mat(i, j)))
										__debugbreak();
									if (std::isinf(mat(j, i)))
										__debugbreak();
								}
							}
							mat(i, i) = pixelTechMatrix[matTo1D(i, i)];
							{
								if (std::isnan(mat(i, i)))
									__debugbreak();
								if (std::isinf(mat(i, i)))
									__debugbreak();
							}
						}

						if (useLT)
						{
							Float NLTZero = size_t(width) * size_t(height) * size_t(spp) - LTSamples[PixelTo1D(x, y)];

							mat += (LTfiller * NLTZero);

						}
						try
						{
							mat = arma::pinv(mat);
						}
						catch (std::exception const& e)
						{
							std::cout << e.what() << std::endl;
							mat = MatrixT(numTechs, numTechs);
							continue;
						}
						const AtomicFloat* pixelContribVectors = getPixelContribVectors(x, y);
						RGBColor estimate = 0;
						for (int k = 0; k < 3; ++k) {
							// fill the vector
							for (int i = 0; i < numTechs; ++i) {
								vec[i] = pixelContribVectors[k * numTechs + i];
							}

							vec = mat * vec;

							if (useLT)
							{
								vec[numTechs - 1] *= (width*height);
							}

							estimate[k] = sum(vec);
							//estimate[k] = vec[vec.size() - 1];
						}
						(*film)(x, y) += estimate;

					}
				}
		}
		else {
			// TODO
		}
	}

	//////////////////////////////////////////////////////////////////////////////////\\
	// Used only for the progressive estimator										   \\
	// This should be called once all the samples have been drawn					   //
	// This is when the optimal estimator <F>o is comnputed
	// // This is also here that the linear system is solved every once in a
	// while		 // and the alpha vector is updated
	// //
	/////////////////////////////////////////////////////////////////////////////////
	void Loop()
	{
		if (!useDirect) {
			// TODO
		}
	}
};


class BalanceSolverImage {

public:

	using RGBColor = Geometry::RGBColor;

	using Float = double;
	using Film = Image::Image<RGBColor>;

	using MatrixT = arma::mat;
	using VectorT = arma::vec;
	using Float = double;
	using AtomicUInt = std::atomic<unsigned int>;
	using AtomicFloat = std::atomic<Float>;



protected:
	const int numTechs;
	int width, height;

	
	std::vector<RGBColor> m_result;
	

	int PixelTo1D(int x, int y) {
		return (x) + (y)* width;
	}

	size_t matTo1D(int row, int col) const {
		// return row * numTechs + col;
		if (col > row) std::swap(row, col);
		return (row * (row + 1)) / 2 + col;
	}

public:

	const bool useDirect = true;
	const bool useLT;


	BalanceSolverImage(int numTechs, int width, int height, int spp, bool useLT)
		: numTechs(numTechs),
		width(width),
		height(height),
		useLT(useLT)
	{
		int msize = numTechs * (numTechs + 1) / 2;
		int res = width * height;

		m_result.resize(res);
		std::fill(m_result.begin(), m_result.end(), 0);
	}

	///////////////////////////////////////////////////////////////////////////////////////////\\
	// Should be called every time a sample is drawn										    \\
	// Update the matrix and the vector with the drawn sample and the extra info provided	    //
	// Also usualy computes a part of the estimate <F>o when using the
	// progressive estimator   //
	////////////////////////////////////////////////////////////////////////////////////////////
	void AddEstimate(
		const Math::Vector2f& p,
		const RGBColor& f,
		const Float* pdfs,
		double sumqi,
		int techIndex
	)
	{
		double ni = (useLT && techIndex == numTechs - 1) ? width * height : 1;
		double qi = pdfs[techIndex] * ni;
		int x = p[0] * width, y = p[1]* height;
		m_result[PixelTo1D(x, y)] += f / sumqi;
		//m_result[PixelTo1D(x, y)] += f / (qi * numTechs);
		//if (techIndex == 1)
		//	m_result[PixelTo1D(x, y)] += f / qi;
	}

	void AddZeroEstimate(
		const Math::Vector2f& p,
		int techIndex
	)
	{
		
	}


	////////////////////////////////////////////////////////////////////////////////////\\
	// Should be called at the end, to get the final optimal result						 \\
	// Computes the last line after the loop in the estimators algorithms				 //
	// Solves the system and returns sun alpha for the direct estimator
	// // Returns the avg of the estimates already computed for the progressive
	// estimator //
	////////////////////////////////////////////////////////////////////////////////////
	template <bool MAJOR>
	void DevelopFilm(Image::Image<RGBColor, MAJOR>* film, int numIterations)
	{
		if (useDirect) {


#pragma omp parallel for schedule(dynamic)
			for (int x = 0; x < width; ++x)
			{
				for (int y = 0; y < height; ++y)
				{
					(*film)(x, y) += m_result[PixelTo1D(x, y)] / numIterations;
				}
			}
		}
		else {
			// TODO
		}
	}

	//////////////////////////////////////////////////////////////////////////////////\\
	// Used only for the progressive estimator										   \\
	// This should be called once all the samples have been drawn					   //
	// This is when the optimal estimator <F>o is comnputed
	// // This is also here that the linear system is solved every once in a
	// while		 // and the alpha vector is updated
	// //
	/////////////////////////////////////////////////////////////////////////////////
	void Loop()
	{
		if (!useDirect) {
			// TODO
		}
	}
};

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

	virtual void addZeroEstimate(int tech_index, Math::Vector2f const& uv) = 0;

	virtual void loop() = 0;

	virtual void solve(Image::Image<Geometry::RGBColor, MAJOR>& res, int iterations) = 0;
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


	virtual void addEstimate(Geometry::RGBColor const& balanceEstimate, double* balanceWeights, int tech_index, Math::Vector2f const& uv)
	{
		Math::Vector<int, 2> pixel = { uv[0] * m_width, uv[1] * m_height };
		m_image(pixel) += balanceEstimate;
	}

	virtual void addZeroEstimate(int tech_index, Math::Vector2f const& uv)
	{}

	virtual void loop()
	{}

	virtual void solve(Image::Image<Geometry::RGBColor, MAJOR>& res, int iterations)
	{
		Parallel::ParallelFor(0, m_width * m_height,
			[&](int pixel)
			{
				res[pixel] += m_image[pixel] / iterations;
			});
	}
};

template <bool MAJOR>
class OptimalEstimatorImage: public ImageEstimator<MAJOR>
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

	// In Byte
	const unsigned int m_pixel_data_size;
	const unsigned int m_vector_ofsset;
	const unsigned int m_counter_offset;

	static_assert(sizeof(char) == 1);
	std::vector<char> m_data;

	//describes how many times more samples each technique has
	std::vector<unsigned int> m_over_samples;


public:

	OptimalEstimatorImage(int N, int width, int height):
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

	void setOverSample(int techIndex, int n)
	{
		m_over_samples[techIndex] = n;
	}

	void addEstimate(Geometry::RGBColor const& balanceEstimate, Float* balanceWeights, int tech_index, Math::Vector2f const& uv)
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

	
	void addZeroEstimate(int tech_index, Math::Vector2f const& uv)
	{
		PixelData data = getPixelData(uv);
		++data.sampleCount[tech_index];
		int mat_index = matTo1D(tech_index, tech_index);
		data.techMatrix[mat_index] = data.techMatrix[mat_index] + 1.0;
	}

	void loop()
	{}

	inline void fillMatrix(MatrixT& matrix, PixelData const& data, int iterations)const
	{
		for (int i = 0; i < m_numtechs; ++i)
		{
			for (int j = 0; j < i; ++j)
			{
				Float elem = data.techMatrix[matTo1D(i, j)];
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

	void solve(Image::Image<Geometry::RGBColor, MAJOR>& res, int iterations)
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
					vector[i] = contribVector[i];
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