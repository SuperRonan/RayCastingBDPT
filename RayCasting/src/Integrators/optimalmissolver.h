#pragma once

#include <Image/Image.h>
#include <Geometry/RGBColor.h>
#include <omp.h>
#include <cassert>

#define PRINT(var) std::cout << #var << ": " << var << std::endl;


////////////////////////////////////////////////////////////////////\\
// Base class for Optimal MIS image solvers							//
// Defines the 3 main functions splitting the estimator algorithms //
////////////////////////////////////////////////////////////////////
class OptimalSolverImage {

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
	Film* film;

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

	size_t To1D(int row, int col) const {
		// return row * numTechs + col;
		if (col > row) std::swap(row, col);
		return (row * (row + 1)) / 2 + col;
	}

public:

	const bool useDirect = true;
	const bool useLT = true;


    OptimalSolverImage(int numTechs, Film* img)
        : numTechs(numTechs),
          width(img->width()),
          height(img->height()),
          film(img)
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

    ///////////////////////////////////////////////////////////////////////////////////////////\\
	// Should be called every time a sample is drawn										    \\
	// Update the matrix and the vector with the drawn sample and the extra info provided	    //
    // Also usualy computes a part of the estimate <F>o when using the
    // progressive estimator   //
    ////////////////////////////////////////////////////////////////////////////////////////////
	void AddEstimate(
		const Math::Vector2f& p,
		const RGBColor& balanceEstimate,
		const Float* balanceWeights,
		const Float* Pdfs,
		int techIndex
	)
	{
		
		Math::Vector<int, 2> pPixel(p[0], p[1]);
		{
			// Update the matrix of p
			AtomicFloat* techMatrix = getPixelTechMatrix(pPixel[0], pPixel[1]);
			for (int i = 0; i < numTechs; ++i) {
				assert(!(std::isnan(balanceWeights[i]) ||
					std::isinf(balanceWeights[i])));
				for (int j = 0; j <= i; ++j) {
					double tmp = (balanceWeights[i] * balanceWeights[j]);
					techMatrix[To1D(i, j)] = techMatrix[To1D(i, j)] + tmp;
				}
			}
		}
		if (useLT)
			LTSamples[PixelTo1D(pPixel[0], pPixel[1])]++;

		// Update only the contribution vector of the pixel
		// If we want to add filters, this is certainly where it should be done
		if (!balanceEstimate.isBlack()) {
			AtomicFloat* pixelContribVectors =
				getPixelContribVectors(pPixel[0], pPixel[1]);

			for (int k = 0; k < 3; ++k) {
				for (int i = 0; i < numTechs; ++i) {
					double tmp = balanceEstimate[k] * balanceWeights[i];
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
	/*
	void AddEstimate(
		const Math::Vector2f& p,
		const RGBColor& balanceEstimate,
		const Float* balanceWeights,
		const Float* Pdfs,
		int techIndex
	)
	{
		Math::Vector<int, 2> pPixel(p[0], p[1]);
		{
			// Update the matrix of p
			AtomicFloat* techMatrix = getPixelTechMatrix(pPixel[0], pPixel[1]);
			for (int i = 0; i < numTechs; ++i) {
				assert(!(std::isnan(balanceWeights[i]) ||
					std::isinf(balanceWeights[i])));
				for (int j = 0; j <= i; ++j) {
					double tmp = (balanceWeights[i] * balanceWeights[j]);
					techMatrix[To1D(i, j)] = techMatrix[To1D(i, j)] + tmp;
				}
			}
		}
		if (useLT)
			LTSamples[PixelTo1D(pPixel[0], pPixel[1])]++;

		// Update only the contribution vector of the pixel
		// If we want to add filters, this is certainly where it should be done
		if (!balanceEstimate.isBlack()) {
			AtomicFloat* pixelContribVectors =
				getPixelContribVectors(pPixel[0], pPixel[1]);

			for (int k = 0; k < 3; ++k) {
				for (int i = 0; i < numTechs; ++i) {
					double tmp = balanceEstimate[k] * balanceWeights[i];
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
	*/

    ////////////////////////////////////////////////////////////////////////////////////\\
	// Should be called at the end, to get the final optimal result						 \\
	// Computes the last line after the loop in the estimators algorithms				 //
    // Solves the system and returns sun alpha for the direct estimator
    // // Returns the avg of the estimates already computed for the progressive
    // estimator //
    ////////////////////////////////////////////////////////////////////////////////////
	void DevelopFilm(int numIterations, bool balanceMis = false)
	{
		if (useDirect) {
			// pre allocate the vector and the matrix

			std::vector<MatrixT> mats(omp_get_num_threads(), MatrixT(numTechs, numTechs));
			MatrixT LTfiller(numTechs, numTechs);
			std::vector<VectorT> vecs(omp_get_num_threads(), VectorT(numTechs));
			
			//std::vector<MatrixT> pinvs(omp_get_num_threads(), MatrixT(numTe;


			LTfiller.fill(0);
			LTfiller(numTechs - 1, numTechs - 1) = 1;

#pragma omp parallel for schedule(dynamic)
			for (int x = 0; x < width; ++x)
			{
				int tid = omp_get_thread_num();
				MatrixT& mat = mats[tid];
				VectorT& vec = vecs[tid];
				for (int y = 0; y < height; ++y)
				{
					const AtomicFloat* pixelTechMatrix =
						getPixelTechMatrix(x, y);
					for (int i = 0; i < numTechs; ++i) {
						for (int j = 0; j <= i; ++j) {
							mat(i, j) = pixelTechMatrix[To1D(i, j)];
							mat(j, i) = pixelTechMatrix[To1D(i, j)];
						}
						mat(i, i) = pixelTechMatrix[To1D(i, i)];
					}

					if (useLT)
					{
						mat += LTfiller * Float(width * height -
							LTSamples[PixelTo1D(x, y)]); /*
	/
		   Float(width * height);*/
					}

					const auto pinv = arma::pinv(mat);

					const AtomicFloat* pixelContribVectors = getPixelContribVectors(x, y);
					RGBColor estimate = 0;
					for (int k = 0; k < 3; ++k) {
						// fill the vector
						for (int i = 0; i < numTechs; ++i) {
							vec[i] = pixelContribVectors[k * numTechs + i];
						}
						
						vec = pinv * vec;
						
						if (useLT)
						{
							vec[vec.size() - 1] *= Float(width * height) / numIterations;
							//vec[vec.size() - 1] = 0;
						}


						estimate[k] = sum(vec);
						//estimate[k] = vec[vec.size() - 1];
					}
					(*film)[x][y] = estimate;

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
