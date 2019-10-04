#ifndef _Math_RandomDirection_H
#define _Math_RandomDirection_H

#include <math.h>
#include <stdlib.h>
#include <Math/sobol.h>
#include <Math/Constant.h>
#include <Math/Quaternion.h>
#include <Math/Sampler.h>

namespace Math
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	RandomDirection
	///
	/// \brief	Random direction sampling. the sampling is biased by a cosine distribution, useful for
	/// 		respecting a BRDF distribution (diffuse or specular).
	/// 		
	/// \author	F. Lamarche, University of Rennes 1.
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class RandomDirection
	{
	protected:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	static ::std::pair<double,double> RandomDirection::randomPolar(double n=1.0)
		///
		/// \brief	Random sampling of spherical coordinates.
		///
		/// \author	F. Lamarche, University of Rennes 1.
		/// \date	04/12/2013
		///
		/// \param	n	(optional) The shininess (1.0 if diffuse).
		///
		/// \return	Random spherical coordinates repecting a cos^n distribution.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		::std::pair<double,double> randomPolar(double n=1.0) const
		{
			//double rand1 = random() ;
			double rand1 = sampler->generateContinuous<double>();
			double p = pow(rand1, 1/(n+1)) ;
			double theta = acos(p) ;
			//double rand2 = random() ;
			double rand2 = sampler->generateContinuous<double>();
			double phy = (double)(2*pi*rand2) ;
			return ::std::make_pair(theta, phy) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	static Math::Vector3 RandomDirection::getVector(double theta, double phy)
		///
		/// \brief	Computes the normalized vector corresponding to provided spherical coordinates.
		/// 		
		/// \author	F. Lamarche, University of Rennes 1.
		/// \date	04/12/2013
		///
		/// \param	theta	
		/// \param	phy  	
		///
		/// \return	The vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		static Math::Vector3f getVector(double theta, double phy)
		{
			return Math::makeVector(sin(theta)*cos(phy), sin(theta)*sin(phy), cos(theta)) ;
		}

	public:


	protected:

		/// \brief	The main direction for sampling.
		Math::Vector3f m_direction ;

		/// \brief	A direction normal to the direction vector.
		Math::Vector3f m_directionNormal ;

		/// \brief	The specular coefficient.
		double m_n ;

		Math::Sampler* sampler;

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RandomDirection::RandomDirection(Math::Vector3 const & direction, double n=1.0)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, University of Rennes 1.
		/// \date	04/12/2013
		///
		/// \param	direction	The main direction of the random sampling.
		/// \param	n		 	n The shininess of the surface (1.0 is diffuse component, the 
		/// 					shininess otherwise)
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RandomDirection(Math::Sampler* p_sampler, Math::Vector3f const & direction, double n=1.0)
			: m_direction(direction.normalized()), m_n(n), sampler(p_sampler)
		{
			// We compute a vector normal to the main direction
			m_directionNormal = Math::makeVector(1.0f,0.0f,0.0f) ;
			m_directionNormal = m_directionNormal - m_direction*(m_direction*m_directionNormal) ;
			if(m_directionNormal.norm()<std::numeric_limits<double>::epsilon()*10.0)
			{
				m_directionNormal = Math::makeVector(0.0f,1.0f,0.0f) ;
				m_directionNormal =  m_directionNormal - m_direction*(m_direction*m_directionNormal) ;
				if(m_directionNormal.norm()<std::numeric_limits<double>::epsilon()*10.0)
				{
					m_directionNormal = Math::makeVector(0.0f,0.0f,1.0f) ;
					m_directionNormal =  m_directionNormal - m_direction*(m_direction*m_directionNormal) ;

				}	
			}
			m_directionNormal = m_directionNormal.normalized();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector3 RandomDirection::generate() const
		///
		/// \brief	Generate a random direction respecting a cosine^n distribution.
		///
		/// \author	F. Lamarche, University of Rennes 1.
		/// \date	04/12/2013
		///
		/// \return	The random direction.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector3f generate() const
		{
			::std::pair<double,double> perturbation = randomPolar(m_n) ;
			Math::Quaternion<double> q1(m_directionNormal, perturbation.first) ;
			Math::Quaternion<double> q2(m_direction, perturbation.second) ;
			Math::Quaternion<double> result = q2.rotate(q1.rotate(m_direction)) ;
			return result.v() ;
		}

	};
}

#endif