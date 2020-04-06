#ifndef _Geometry_Cube_H
#define _Geometry_Cube_H

#include <Geometry/Geometry.h>
#include <Geometry/Square.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Cube
	///
	/// \brief	A cube centered in (0,0,0), 1.0 height/width....
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Cube : public GeometryCollection
	{
	protected:
		Square * m_square[6] ;

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Cube::Cube(Material * material)
		///
		/// \brief	Constructeur
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	material	The material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Cube(Material* material, Math::Vector3f const& origin = { 0, 0, 0 }, Math::Vector3f const& u = { 1, 0, 0 }, Math::Vector3f const& v = { 0, 1, 0 }, Math::Vector3f const& w = { 0, 0, 1 })
			: GeometryCollection(material)
		{
			{//down
				Square sq0(material, origin, u, v);
				sq0.translate(-w/2);
				merge(sq0);
			}

			{//up
				Square sq1(material, origin, v, u);
				sq1.translate(w/2);
				//sq1.rotate(Math::Quaternion<double>(Math::makeVector(1.0f, 0.0f, 0.0f), (double)M_PI));
				merge(sq1);
			}

			{//right
				Square sq2(material, origin, w, u);
				sq2.translate(-v/2);
				//sq2.rotate(Math::Quaternion<double>(Math::makeVector(1.0f, 0.0f, 0.0f), (double)M_PI));
				merge(sq2);
			}

			{//left
				Square sq3(material, origin, u, w);
				sq3.translate(v/2);
				//sq3.rotate(Math::Quaternion<double>(Math::makeVector(1.0f, 0.0f, 0.0f), (double)-M_PI / 2.0f));
				merge(sq3);
			}

			{//front
				Square sq4(material, origin, v, w);
				sq4.translate(-u/2);
				//sq4.rotate(Math::Quaternion<double>(Math::makeVector(0.0f, 1.0f, 0.0f), (double)M_PI / 2.0f));
				merge(sq4);
			}

			{//back
				Square sq5(material, origin, w, v);
				sq5.translate(u/2);
				//sq5.rotate(Math::Quaternion<double>(Math::makeVector(0.0f, 1.0f, 0.0f), (double)M_PI));
				merge(sq5);
			}
			//ca marche peu
			// Déclaration des carrés constituant les faces
			//m_square[0] = new Square(material) ;
			//m_square[0]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[1] = new Square(material) ;
			//m_square[1]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[1]->rotate(Math::Quaternion(Math::Vector3(1.0f,0.0f,0.0f), (double)M_PI/2.0f)) ;
			//m_square[2] = new Square(material) ;
			//m_square[2]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[2]->rotate(Math::Quaternion(Math::Vector3(1.0f,0.0f,0.0f), (double)M_PI)) ;
			//m_square[3] = new Square(material) ;
			//m_square[3]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[3]->rotate(Math::Quaternion(Math::Vector3(1.0f,0.0f,0.0f), (double)-M_PI/2.0f)) ;
			//m_square[4] = new Square(material) ;
			//m_square[4]->translate(Math::Vector3(0.0,0.0,0.5)) ;
			//m_square[4]->rotate(Math::Quaternion(Math::Vector3(0.0f,1.0f,0.0f), (double)M_PI/2.0f)) ;
			//m_square[5] = new Square(material) ;
			//m_square[5]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[5]->rotate(Math::Quaternion(Math::Vector3(0.0f,1.0f,0.0f), (double)-M_PI/2.0f)) ;
			// Ajout des géométries dans cette géométrie
			//merge(m_square[0]) ;
			//merge(m_square[1]) ;
			//merge(m_square[2]) ;
			//merge(m_square[3]) ;
			//merge(m_square[4]) ;
			//merge(m_square[5]) ;
		}
		
	} ;
}

#endif
