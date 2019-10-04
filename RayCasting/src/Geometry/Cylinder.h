#ifndef _Geometry_Cylinder_H
#define _Geometry_Cylinder_H

#include <Geometry/Geometry.h>
#include <Geometry/Material.h>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	CylinderFabrice
	///
	/// \brief	A cylinder centered in (0,0,0), 1.0 height, 1.0 diameter.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class CylinderFabrice : public GeometryCollection
	{
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	CylinderFabrice::CylinderFabrice(int nbDiv, double scaleDown, double scaleUp, Material * material)
		///
		/// \brief	Constructor
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	nbDiv	 	Number of base circle subdivisions.
		/// \param	scaleDown	Radius of the bottom circle.
		/// \param	scaleUp  	Radius of the top circle.
		/// \param	material 	The material.
		/// \param d1			Add the top circle
		/// \param d2			Add the bottom circle
		////////////////////////////////////////////////////////////////////////////////////////////////////
		CylinderFabrice(int nbDiv, double scaleDown, double scaleUp, Material * material, bool d1=true, bool d2=true):
			GeometryCollection(material)
		{

			DiskFabrice disk1(nbDiv, material);
			disk1.scale(scaleUp);
			disk1.translate(Math::makeVector(0.0f, 0.0f, 0.5f));
			if (d1)
			{	
				merge(disk1);
			}

			DiskFabrice disk2(nbDiv, material);
			disk2.scale(scaleDown);
			disk2.translate(Math::makeVector(0.0f, 0.0f, -0.5f));
			if (d2)
			{
				merge(disk2);
			}

			
			

			for(int cpt=1 ; cpt<=nbDiv ; cpt++)
			{
				addTriangle(disk1.getVertices()[cpt], disk1.getVertices()[(cpt)%nbDiv+1], disk2.getVertices()[cpt]) ;
				addTriangle(disk1.getVertices()[(cpt)%nbDiv+1], disk2.getVertices()[cpt], disk2.getVertices()[(cpt)%nbDiv+1]) ;
				//addTriangle(new Triangle(disk1.getVertices()[cpt], disk1.getVertices()[(cpt+1)%nbDiv], 
				//			disk2.getVertices()[cpt], material)) ;
				//addTriangle(new Triangle(disk1.getVertices()[(cpt+1)%nbDiv], disk2.getVertices()[cpt], 
				//			disk2.getVertices()[(cpt+1)%nbDiv], material)) ;
			}
		}
	};
}

#endif
