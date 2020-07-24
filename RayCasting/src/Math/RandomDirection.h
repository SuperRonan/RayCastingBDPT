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
	class OldRandomDirection
	{
	protected:


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


		OldRandomDirection(Math::Sampler* p_sampler, Math::Vector3f const & direction, double n=1.0)
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


		Vector3f generate() const
		{
			::std::pair<double,double> perturbation = randomPolar(m_n) ;
			Math::Quaternion<double> q1(m_directionNormal, perturbation.first) ;
			Math::Quaternion<double> q2(m_direction, perturbation.second) ;
			Math::Quaternion<double> result = q2.rotate(q1.rotate(m_direction)) ;
			return result.v() ;
		}

	};

	class RandomDirection
	{
	protected:


		::std::pair<double, double> polar(double rand1, double rand2, double n = 1.0) const
		{
			double p = pow(rand1, 1 / (n + 1));
			double theta = acos(p);
			//double rand2 = random() ;
			double phi = (double)(2 * pi * rand2);
			return ::std::make_pair(theta, phi);
		}


		static Math::Vector3f getVector(double theta, double phy)
		{
			return Math::makeVector(sin(theta) * cos(phy), sin(theta) * sin(phy), cos(theta));
		}

	public:


	protected:

		/// \brief	The main direction for sampling.
		Math::Vector3f m_direction;

		/// \brief	A direction normal to the direction vector.
		Math::Vector3f m_directionNormal;

		/// \brief	The specular coefficient.
		double m_n;


	public:


		RandomDirection(Math::Vector3f const& direction, double n = 1.0)
			: m_direction(direction.normalized()), m_n(n)
		{
			// We compute a vector normal to the main direction
			m_directionNormal = Math::makeVector(1.0f, 0.0f, 0.0f);
			m_directionNormal = m_directionNormal - m_direction * (m_direction * m_directionNormal);
			if (m_directionNormal.norm() < std::numeric_limits<double>::epsilon() * 10.0)
			{
				m_directionNormal = Math::makeVector(0.0f, 1.0f, 0.0f);
				m_directionNormal = m_directionNormal - m_direction * (m_direction * m_directionNormal);
				if (m_directionNormal.norm() < std::numeric_limits<double>::epsilon() * 10.0)
				{
					m_directionNormal = Math::makeVector(0.0f, 0.0f, 1.0f);
					m_directionNormal = m_directionNormal - m_direction * (m_direction * m_directionNormal);

				}
			}
			m_directionNormal = m_directionNormal.normalized();
		}


		Vector3f generate(double rand1, double rand2) const
		{
			::std::pair<double, double> perturbation = polar(rand1, rand2);
			Math::Quaternion<double> q1(m_directionNormal, perturbation.first);
			Math::Quaternion<double> q2(m_direction, perturbation.second);
			Math::Quaternion<double> result = q2.rotate(q1.rotate(m_direction));
			return result.v();
		}

	};
}

#endif