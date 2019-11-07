#ifndef _Geometry_camera_H
#define _Geometry_camera_H

#include <Math/Vectorf.h>
#include <Math/Quaternion.h>
#include <Geometry/Ray.h>
#include <math.h>
#include <Math/Constant.h>
#include <utils.h>


#define USE_CAMERA_PLANE


namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Camera
	///
	/// \brief	A pinhole camera.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Camera
	{
	public:
		/// \brief	The camera position.
		Math::Vector3f m_position ;
		/// \brief	The aim of the camera.
		Math::Vector3f m_target ;
		/// \brief	Distance of the focal plane.
		static double constexpr m_planeDistance = 1.0;
		/// \brief	Width of the projection rectangle.
		double	      m_planeWidth ;
		/// \brief	Height of the projection rectangle.
		double		  m_planeHeight ;


		/// \brief	The front vector of the camera.
		Math::Vector3f m_front ;
		/// \brief	The right vector.
		Math::Vector3f m_right ;
		/// \brief	The down vector.
		Math::Vector3f m_down ;

		double resolution;


	
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Camera::computeParameters()
		///
		/// \brief	Calculates the camera parameters.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void computeParameters()
		{
			m_front = m_target-m_position ;
			m_front = m_front/m_front.norm() ;
			
			m_right = Math::Quaternion<double>(Math::makeVector(0.0f, 0.0f, 1.0f), -3.14159265f/2.0f).rotate(m_front) ;
			m_right[2] = 0.0f;
			m_right = m_right/m_right.norm() ;
			
			m_down  = m_front ^ m_right ;
			m_down  = m_down/m_down.norm() ;
		}

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Camera::Camera(Math::Vector3 const & position = Math::Vector3(0.0, 0.0, 0.0),
		/// 	Math::Vector3 const & target = Math::Vector3(0.0, 1.0, 0.0), double planeDistance=1.0,
		/// 	double planeWidth=1.0, double planeHeight=1.0)
		///
		/// \brief	Constructeur de la caméra.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	position	 	The camera position
		/// \param	target		 	The target of the camera
		/// \param	planeDistance	La distance of the focal plane.
		/// \param	planeWidth   	Width of the projection rectangle.
		/// \param	planeHeight  	Height of the projection rectangle.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Camera(Math::Vector3f const & position = Math::makeVector(0.0f, 0.0f, 0.0f), 
			   Math::Vector3f const & target = Math::makeVector(0.0f, 1.0f, 0.0f),
			   double planeDist=1,
			   double planeWidth=1.0f, double planeHeight=1.0f): 
			m_position(position), 
			m_target(target), 
		    m_planeWidth(planeWidth / planeDist), 
			m_planeHeight(planeHeight / planeDist)
		{
			computeParameters() ;
		}

		/// <summary>
		/// Translates the camera in local coordinates (X = right, Y = front, Z=up).
		/// </summary>
		/// <param name="translation">The translation vector.</param>
		void translateLocal(Math::Vector3f const & translation)
		{
			Math::Vector3f trans =m_right*translation[0] + m_front*translation[1] - m_down*translation[2];
			m_position = m_position + trans;
			m_target = m_target + trans;
		}


		void update_both(Math::Vector3f const & position, Math::Vector3f const & target)
		{
			m_position = position;
			m_target = target;
			computeParameters();
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Camera::setPosition(Math::Vector3 const & position)
		///
		/// \brief	Sets the camera position position.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	position	The new camera position.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void setPosition(Math::Vector3f const & position)
		{
			m_position = position ;
		}

		Math::Vector3f const& getPosition()const
		{
			return m_position;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Camera::setTarget(Math::Vector3 const & target)
		///
		/// \brief	Sets the target.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	target	The new target.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void setTarget(Math::Vector3f const & target)
		{
			m_target = target ;
			computeParameters() ;
		}

		void scale_screen(double scale)
		{
			m_planeWidth *= scale;
			m_planeHeight *= scale;

#ifndef USE_CAMERA_PLANE
			if (m_planeHeight > Math::pi)
			{
				m_planeHeight = Math::pi;
			}
			if (m_planeWidth > Math::twoPi)
			{
				m_planeWidth = Math::twoPi;
			}
#endif

		}

		//returns the position on the screen of the point p (in world space)
		__forceinline Math::Vector2f screen_position(Math::Vector3f p)const
		{
			p -= m_position;
			return screen_direction(p.normalized());

		}

		__forceinline bool validRaster(Math::Vector2f const& raster)const
		{
			return raster[0] >= 0 && raster[0] <= 1 && raster[1] >= 0 && raster[1] <= 1;
		}

		__forceinline Math::Vector2f raster(Math::Vector3f const& dir)const
		{
			return screen_direction(dir);
		}



#ifdef USE_CAMERA_PLANE


		Ray getRay(double coordX, double coordY) const
		{
			return Ray(m_position, m_front - m_right * (0.5-coordX) * m_planeWidth - m_down * (0.5 - coordY) * m_planeHeight) ;
		}


		__forceinline Math::Vector2f screen_direction(Math::Vector3f const& direction)const
		{
			double w = (direction * m_front);
			if (w < 0)
				return { -1, -1 };
			Math::Vector3f pw = m_front * w;

			double u = ((direction - pw) * m_right) / m_planeWidth * m_planeDistance / (w);
			double v = ((direction - pw) * (-m_down)) / m_planeHeight * m_planeDistance / (w);
			return { u + 0.5, 1 - v - 0.5 };
		}

		__forceinline double We(Math::Vector3f const& dir)const
		{
			return pdfWeSolidAngle(dir);
		}


		//returns the pdf of sampling a direction
		//assumes dir is normalized
		//return a solid angle pdf
		double pdfWeSolidAngle(Math::Vector3f const& dir)const
		{
			if (validRaster(raster(dir)))
			{
				double cost = dir * m_front;
				return pdfWeArea() / (cost * cost * cost);
			}
			return 0;
		}


		__forceinline double pdfWeArea()const
		{
			return resolution / (m_planeWidth * m_planeHeight);
		}
#else

		Ray getRay(double coordX, double coordY) const
		{
			double inclination = Math::piDiv2 + (0.5 - coordY) * m_planeHeight;
			double azimuth = (0.5 - coordX) * m_planeWidth;
			Math::Vector3f sphere_dir = Math::Vector3f::make_sphere(inclination, azimuth);
			Math::Vector3f dir = m_front * sphere_dir[0] - m_right * sphere_dir[1] + m_down * sphere_dir[2];
			return Ray(m_position, dir.normalized());
		}


		__forceinline Math::Vector2f screen_direction(Math::Vector3f const& direction)const
		{
			//first, we need the direction in the camera basis
			Math::Vector3f camera_direction = {(direction * m_front), -direction * m_right, -(direction * m_down)};

			Math::Vector2f inc_azi = inclinationAzimuthOfNormalized(camera_direction);

			return Math::Vector2f(0.5 - (inc_azi[1]) / m_planeWidth, 0.5+ (inc_azi[0] - Math::piDiv2) / m_planeHeight);
		}



		__forceinline double We(Math::Vector3f const& dir)const
		{
			return pdfWeSolidAngle(dir);
		}


		//returns the pdf of sampling a direction
		//assumes dir is normalized
		//return a solid angle pdf
		double pdfWeSolidAngle(Math::Vector3f const& dir)const
		{
			if (validRaster(raster(dir)))
			{
				double cos_inclination = dir * m_down;
				double sin_inclination = std::sqrt(1 - cos_inclination * cos_inclination);
				return pdfWeArea() / sin_inclination;
			}
			return 0;
		}


		__forceinline double pdfWeArea()const
		{
			return resolution / (m_planeWidth * m_planeHeight);
		}
#endif

	} ;
}

#endif
