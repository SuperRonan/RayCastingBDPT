#pragma once

#include <Math/Sampler.h>
#include <Math/Quaternion.h>
#include <Math/Vectorf.h>
#include <Math/Constant.h>

namespace Math
{
	class SolidAngleSampler
	{
	protected:

		Vector3f m_normal;
		Vector3f m_tangeant;

		double m_max_theta;
		double m_cos_max_theta;

	public:

		SolidAngleSampler(Vector3f const& direction, double max_angle, double cos_max_angle):
			m_normal(direction),
			m_max_theta(max_angle),
			m_cos_max_theta(cos_max_angle)
		{
			assert(m_max_theta <= piDiv2);
			m_tangeant = Math::makeVector(1.0f, 0.0f, 0.0f);
			m_tangeant = m_tangeant - m_normal * (m_normal * m_tangeant);
			if (m_tangeant.norm() < std::numeric_limits<double>::epsilon() * 10.0)
			{
				m_tangeant = Math::makeVector(0.0f, 1.0f, 0.0f);
				m_tangeant = m_tangeant - m_normal * (m_normal * m_tangeant);
				if (m_tangeant.norm() < std::numeric_limits<double>::epsilon() * 10.0)
				{
					m_tangeant = Math::makeVector(0.0f, 0.0f, 1.0f);
					m_tangeant = m_tangeant - m_normal * (m_normal * m_tangeant);

				}
			}
			m_tangeant = m_tangeant.normalized();
		}

		std::pair<double, double> polar(double xi1, double xi2)const
		{
			double azimuth = twoPi * xi2; 
			xi1 = 1 -(xi1 * (1 - m_cos_max_theta));
			double inclination = acos(xi1);
			//double inclination = xi1 * m_max_theta;
			return { inclination, azimuth };
		}

		//instead of taking a sampler, it takes the result of what a sampler could give, this is done for segmetation of spheres
		Vector3f generate(double xi1, double xi2)const
		{
			::std::pair<double, double> perturbation = polar(xi1, xi2);
			Math::Quaternion<double> q1(m_tangeant, perturbation.first);
			Math::Quaternion<double> q2(m_normal, perturbation.second);
			Math::Quaternion<double> result = q2.rotate(q1.rotate(m_normal));
			return result.v();
		}

		double surface()const
		{
			return twoPi * (1 - m_cos_max_theta);
		}
	};
}