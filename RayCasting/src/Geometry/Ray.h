#ifndef _Geometry_Ray
#define _Geometry_Ray

#include <Math/Vectorf.h>
#include <iostream>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Ray
	///
	/// \brief	A ray with a source and a direction.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	09/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Ray
	{
	protected:
		/// \brief	Source for the ray.
		Math::Vector3f m_source ;
		/// \brief	The direction of the ray (normalized).
		Math::Vector3f m_direction ;
		Math::Vector3f m_invDirection ;
		int m_sign[3] ;
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Ray::Ray(Math::Vector3 const & source, Math::Vector3 const & direction)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	09/12/2013
		///
		/// \param	source   	The source of the ray.
		/// \param	direction	The direction of the ray.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Ray(Math::Vector3f const & source, Math::Vector3f const & direction)
			: m_source(source), m_direction(direction/direction.norm())
		{
			//m_invDirection = m_direction.simdInv() ;
			m_invDirection = Math::makeVector(1.0f / m_direction[0], 1.0f / m_direction[1], 1.0f / m_direction[2]);
			m_sign[0] = m_direction[0]<0.0 ;
			m_sign[1] = m_direction[1]<0.0 ;
			m_sign[2] = m_direction[2]<0.0 ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Math::Vector3 & Ray::source() const
		///
		/// \brief	Gets the source of the ray.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	09/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Math::Vector3f & source() const
		{ return m_source ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Math::Vector3 & Ray::direction() const
		///
		/// \brief	Gets the direction of the ray.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	09/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Math::Vector3f & direction() const
		{ return m_direction ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector3 Ray::invDirection() const
		///
		/// \brief	Gets the inverse direction of the ray (1.0/dx, 1.0/dy, 1.0/dz). Useful for ray box 
		/// 		intersection computation.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	09/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Math::Vector3f & invDirection() const
		{
			return m_invDirection ;
			//return m_direction.simdInv() ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Ray::project(Math::Vector3 const & point, double & t, Math::Vector3 & delta)
		///
		/// \brief	Projects a point on the ray.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	09/12/2013
		///
		/// \param	point		 	The point.
		/// \param [in,out]	t	 	The parameter such as the projection is source()+direction()*t.
		/// \param [in,out]	delta	The vector such as the projection is source()+delta.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void project(Math::Vector3f const & point, double & t, Math::Vector3f & delta)
		{
			t=(point-source())*direction() ;
			delta = (point-source())-direction()*t ;
		}





		Math::Vector3f sample_point(double t)const noexcept
		{
			return m_source + m_direction * t;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	int * Ray::getSign() const
		///
		/// \brief	Gets a three dimensional array containing the sign of each coordinate (1 if negative, 
		/// 		0 if positive).
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	10/12/2013
		///
		/// \return	null if it fails, else the sign.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const int * getSign() const
		{
			return m_sign ;
		}
	} ;

	inline std::ostream & operator<< (std::ostream & out, Ray const & ray)
	{
		out<<"Ray ("<<ray.source()<<","<<ray.direction()<<")" ;
		return out ;
	}
}

#endif
