#pragma once
#include <Math/Vectorf.h>
#include "RGBColor.h"
namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	PointLight
	///
	/// \brief	A point light.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class DirectionalLight
	{
	protected:
		/// \brief	The light direction (norm).
		Math::Vector3f m_direction;

		/// \brief	The light color.
		RGBColor m_color;

	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	PointLight::PointLight(Math::Vector3 const & position, RGBColor const & color)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	position	The position.
		/// \param	color   	The color.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		DirectionalLight(Math::Vector3f const & direction = Math::makeVector(0.0, 0.0, 1.0), RGBColor const & color = RGBColor())
			: m_direction(direction), m_color(color)
		{}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const RGBColor & PointLight::color() const
		///
		/// \brief	Gets the light color.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const RGBColor & color() const
		{
			return m_color;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Math::Vector3 & PointLight::position() const
		///
		/// \brief	Gets the position of the light.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Math::Vector3f & direction() const
		{
			return m_direction;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	virtual RGBColor PointLight::color(RayTriangleIntersection const & intersection)
		///
		/// \brief	Compute the color of an intersection lightened by this light.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	intersection A valid intersection between a ray and a triangle.
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/*
		virtual RGBColor color(RayTriangleIntersection const & intersection)
		{
			std::cerr << "PointLight::color Not implemented [TODO]" << std::endl;
			return RGBColor(0, 1, 0);
		}
		*/
	};
}