#ifndef _Rennes1_Math_Matrix4x4_H
#define _Rennes1_Math_Matrix4x4_H

#include <Math/Vector.h>
#include <cassert>
#include <math.h>
#include <stdexcept>

namespace Math
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Matrix4x4
	///
	/// \brief	Matrix 4x4 stored in row major format.
	///
	/// \author	Fabrice Lamarche, University of Rennes 1
	/// \date	22/12/2011
	////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class Float>
	class Matrix4x4
	{
	protected:
		//! Matrix data stored in row major format.
		Math::Vector<Float,4> m_data[4] ;

	public:
		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Constructor. Creates a 0 matrix.
		/// 
		/// \return 
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		Matrix4x4()
		{
			Math::Vector<Float,4> zero = makeVector(Float(), Float(), Float(), Float()) ;
			m_data[0] = zero ;
			m_data[1] = zero ;
			m_data[2] = zero ;
			m_data[3] = zero ;
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Constructor intializing matrix with fiven rows
		/// 
		/// \param row0 The first row
		/// \param row1 The second row
		/// \param row2 The third row
		/// \param row3 The fourth row
		/// 	   
		/// \return 
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		Matrix4x4(Math::Vector<Float,4> const & row0, Math::Vector<Float,4> const & row1, Math::Vector<Float,4> const & row2, Math::Vector<Float,4> const & row3)
		{
			m_data[0] = row0 ;
			m_data[1] = row1 ;
			m_data[2] = row2 ;
			m_data[3] = row3 ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float,4> const & Matrix4x4::getRow(int index) const
		///
		/// \brief	Gets a row.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \param	index	Zero-based index of the row.
		///
		/// \return	The row.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float,4> const & getRow(int index) const
		{
			assert(index>=0 && index<4) ;
			return m_data[index] ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Matrix4x4::setRow(const Math::Vector<Float,4> & row, int index)
		///
		/// \brief	Sets a row of the matrix.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/03/2016
		///
		/// \param	row  	The row.
		/// \param	index	Zero-based index of the row.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void setRow(const Math::Vector<Float,4> & row, int index)
		{
			assert(index>=0 && index<4) ;
			m_data[index] = row ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float, 4> Matrix4x4::getColumn(int index) const
		///
		/// \brief	Gets a column.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \param	index	Zero-based index of the column.
		///
		/// \return	The column.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float, 4> getColumn(int index) const
		{
			assert(index>=0 && index<4) ;
			return makeVector(m_data[0][index], m_data[1][index], m_data[2][index], m_data[3][index]) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Matrix4x4::setColumn(Math::Vector<Float, 4> const & column, int index)
		///
		/// \brief	Sets a column of the matrix.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/03/2016
		///
		/// \param	column	The column.
		/// \param	index 	Zero-based index of the column.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void setColumn(Math::Vector<Float, 4> const & column, int index)
		{
			assert(index>=0 && index<4) ;
			m_data[0][index]=column[0] ;
			m_data[1][index]=column[1] ;
			m_data[2][index]=column[2] ;
			m_data[3][index]=column[3] ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Float & Matrix4x4::operator() (int row, int column)
		///
		/// \brief	 Gets a value in the matrix.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Float & operator() (int row, int column)
		{
			assert(row>=0 && row<4) ;
			assert(column>=0 && column<4) ;
			return m_data[row][column] ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Float & Matrix4x4::operator() (int row, int column) const
		///
		/// \brief	Gets a value in the matrix.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Float & operator() (int row, int column) const
		{
			assert(row>=0 && row<4) ;
			assert(column>=0 && column<4) ;
			return m_data[row][column] ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Matrix4x4 Matrix4x4::operator+(Matrix4x4 const & matrix) const
		///
		/// \brief	Addition operator.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \param	matrix	The matrix.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Matrix4x4 operator+(Matrix4x4 const & matrix) const
		{
			return Matrix4x4(m_data[0]+matrix.m_data[0], m_data[1]+matrix.m_data[1], m_data[2]+matrix.m_data[2], m_data[3]+matrix.m_data[3]) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Matrix4x4 Matrix4x4::operator-(Matrix4x4 const & matrix) const
		///
		/// \brief	Soustraction operator.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \param	matrix	The matrix.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Matrix4x4 operator-(Matrix4x4 const & matrix) const
		{
			return Matrix4x4(m_data[0]-matrix.m_data[0], m_data[1]-matrix.m_data[1], m_data[2]-matrix.m_data[2], m_data[3]-matrix.m_data[3]) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Matrix4x4 Matrix4x4::operator*(Matrix4x4 const & matrix) const
		///
		/// \brief	Multiplication operator.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \param	matrix	The matrix.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Matrix4x4 operator*(Matrix4x4 const & matrix) const
		{
			Matrix4x4 result ;
			for(int otherColumn=0 ; otherColumn<4 ; ++otherColumn)
			{
				Math::Vector<Float,4> column = matrix.getColumn(otherColumn) ;
				for(int thisRow=0 ; thisRow<4 ; ++thisRow)
				{
					result(thisRow, otherColumn) = column*getRow(thisRow) ;
				}
			}
			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Matrix4x4 Matrix4x4::operator*(Float const & value) const
		///
		/// \brief	Multiplication operator.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \param	value	The value.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Matrix4x4 operator*(Float const & value) const
		{
			return Matrix4x4(getRow(0)*value, getRow(1)*value, getRow(2)*value, getRow(3)*value) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float,4> Matrix4x4::operator* (Math::Vector<Float,4> const & value) const
		///
		/// \brief	Multiplication operator.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \param	value	The value.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float,4> operator* (Math::Vector<Float,4> const & value) const
		{
			return makeVector(getRow(0)*value, getRow(1)*value, getRow(2)*value, getRow(3)*value) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float,3> Matrix4x4::operator* (Math::Vector<Float, 3> const & value) const
		///
		/// \brief	Multiplication operator.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \param	value	The value.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float,3> operator* (Math::Vector<Float, 3> const & value) const
		{
			Math::Vector<Float,4> result = (*this)*makeVector(value[0], value[1], value[2], Float(1.0)) ;
			return makeVector(result[0], result[1], result[2])/result[3] ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Matrix4x4::toBuffer(Float * buffer)
		///
		/// \brief	Write this matrix into a buffer compatible with OpenGL (i.e. in a column major matrix). 
		/// 		The buffer must be able to contain 16 Float values.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	02/03/2016
		///
		/// \param [in,out]	buffer	The buffer that will be filled.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void toBuffer(Float * buffer) const
		{
			for(int cpt=0 ; cpt<4 ; ++cpt)
			{
				Vector<Float,4> column = getColumn(cpt) ;
				::std::copy(column.begin(), column.end(), buffer) ;
				buffer += 4 ;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Matrix4x4::loadFromBuffer(Float * buffer)
		///
		/// \brief	Loads this matrix from data in the buffer. Buffer must contain the matrix in column major
		/// 		representation, this ensures a compatibility with OpenGL
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	02/03/2016
		///
		/// \param [in,out]	buffer	If non-null, the buffer.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void loadFromBuffer(Float * buffer)
		{
			m_data[0] = makeVector(buffer[0], buffer[4], buffer[8], buffer[12]) ;
			m_data[1] = makeVector(buffer[1], buffer[5], buffer[9], buffer[13]) ;
			m_data[2] = makeVector(buffer[2], buffer[6], buffer[10], buffer[14]) ;
			m_data[3] = makeVector(buffer[3], buffer[7], buffer[11], buffer[15]) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Matrix4x4 Matrix4x4::inverse() const
		///
		/// \brief	Gets the inverse of this matrix.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	02/03/2016
		///
		/// \return	The inverse of this matrix.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Matrix4x4 inverse() const
		{
			// Uses the code of MESA (gluInvertMatrix)
			Float m[16] ; 
			Float inv[16], det;
			int i;

			toBuffer(m) ;

			// Trust the compiler to optimize this thing ;)
			inv[0] = m[5]  * m[10] * m[15] - m[5]  * m[11] * m[14] - m[9]  * m[6]  * m[15] + m[9]  * m[7]  * m[14] + m[13] * m[6]  * m[11] - m[13] * m[7]  * m[10];  //
			inv[4] = -m[4]  * m[10] * m[15] + m[4]  * m[11] * m[14] + m[8]  * m[6]  * m[15] - m[8]  * m[7]  * m[14] - m[12] * m[6]  * m[11] + m[12] * m[7]  * m[10]; //
			inv[8] = m[4]  * m[9] * m[15] - m[4]  * m[11] * m[13] - m[8]  * m[5] * m[15] + m[8]  * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
			inv[12] = -m[4]  * m[9] * m[14] + m[4]  * m[10] * m[13] + m[8]  * m[5] * m[14] - m[8]  * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
			inv[1] = -m[1]  * m[10] * m[15] + m[1]  * m[11] * m[14] + m[9]  * m[2] * m[15] - m[9]  * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
			inv[5] = m[0]  * m[10] * m[15] - m[0]  * m[11] * m[14] - m[8]  * m[2] * m[15] + m[8]  * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
			inv[9] = -m[0]  * m[9] * m[15] + m[0]  * m[11] * m[13] + m[8]  * m[1] * m[15] - m[8]  * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
			inv[13] = m[0]  * m[9] * m[14] - m[0]  * m[10] * m[13] - m[8]  * m[1] * m[14] + m[8]  * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
			inv[2] = m[1]  * m[6] * m[15] - m[1]  * m[7] * m[14] - m[5]  * m[2] * m[15] + m[5]  * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
			inv[6] = -m[0]  * m[6] * m[15] + m[0]  * m[7] * m[14] + m[4]  * m[2] * m[15] - m[4]  * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
			inv[10] = m[0]  * m[5] * m[15] - m[0]  * m[7] * m[13] - m[4]  * m[1] * m[15] + m[4]  * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
			inv[14] = -m[0]  * m[5] * m[14] + m[0]  * m[6] * m[13] + m[4]  * m[1] * m[14] - m[4]  * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
			inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
			inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
			inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
			inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

			det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

			if (det == 0)
				throw ::std::runtime_error("Zero determinant!");

			det = Float(1.0) / det;

			for (i = 0; i < 16; i++)
				inv[i] = inv[i] * det;

			Matrix4x4 result ;
			result.loadFromBuffer(inv) ;

			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Matrix4x4 Matrix4x4::transpose() const
		///
		/// \brief	Gets the transpose of the current matrix.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/03/2016
		///
		/// \return	The transpose of the current matrix.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Matrix4x4 transpose() const
		{
			Matrix4x4 result ;
			for(int cpt=0 ; cpt<4 ; ++cpt)
			{
				result.setColumn(getRow(cpt), cpt) ;
			}
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Creates an identity matrix.
		/// 
		/// \return An identity matrix.
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		static Matrix4x4 getIdentity()
		{
			return Matrix4x4(makeVector(Float(1.0),Float(0.0),Float(0.0),Float(0.0)), 
								makeVector(Float(0.0),Float(1.0),Float(0.0),Float(0.0)),
								makeVector(Float(0.0),Float(0.0),Float(1.0),Float(0.0)),
								makeVector(Float(0.0),Float(0.0),Float(0.0),Float(1.0))
				);
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Creates a translation matrix in homogeneous coordinates.
		/// 
		/// \param The translation vector (3D)
		/// \return A translation matrix in homogeneous coordinates.
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		static Matrix4x4 getTranslation(Math::Vector<Float,3> const & trans)
		{
			return Matrix4x4(makeVector(Float(1.0),Float(0.0),Float(0.0),Float(trans[0])), 
								makeVector(Float(0.0),Float(1.0),Float(0.0),Float(trans[1])),
								makeVector(Float(0.0),Float(0.0),Float(1.0),Float(trans[2])),
								makeVector(Float(0.0),Float(0.0),Float(0.0),Float(1.0))
				);
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Creates a rotation around X matrix.
		/// 
		/// \param angle Rotation angle
		/// \return 
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		static Matrix4x4 getRotationX(Float const & angle)
		{
			Float cosine = Float(cos(angle)) ;
			Float sine = Float(sin(angle)) ;
			return Matrix4x4(makeVector(Float(1.0),Float(0.0),Float(0.0),Float(0.0)), 
								makeVector(Float(0.0),cosine   ,-sine    ,Float(0.0)),
								makeVector(Float(0.0),sine     ,cosine   ,Float(0.0)),
								makeVector(Float(0.0),Float(0.0),Float(0.0),Float(1.0))
				);
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Creates a rotation around Y matrix.
		/// 
		/// \param angle The rotation angle
		/// \return 
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		static Matrix4x4 getRotationY(Float const & angle)
		{
			Float cosine = Float(cos(angle)) ;
			Float sine = Float(sin(angle)) ;
			return Matrix4x4(makeVector(cosine   ,Float(0.0),sine     ,Float(0.0)), 
								makeVector(Float(0.0),Float(1.0),Float(0.0),Float(0.0)),
								makeVector(-sine    ,Float(0.0),cosine   ,Float(0.0)),
								makeVector(Float(0.0),Float(0.0),Float(0.0),Float(1.0))
				);
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Creates a rotation around Z matrix.
		/// 
		/// \param angle The rotation angle.
		/// \return 
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		static Matrix4x4 getRotationZ(Float const & angle)
		{
			Float cosine = Float(cos(angle)) ;
			Float sine = Float(sin(angle)) ;
			return Matrix4x4(makeVector(cosine   ,-sine    ,Float(0.0),Float(0.0)), 
								makeVector(sine     ,cosine   ,Float(0.0),Float(0.0)),
								makeVector(Float(0.0),Float(0.0),Float(1.0),Float(0.0)),
								makeVector(Float(0.0),Float(0.0),Float(0.0),Float(1.0))
				);
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Creates a scale matrix.
		/// 
		/// \param sx X scale
		/// \param sy Y scale
		/// \param sz Z scale
		/// \return A scale matrix.
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		static Matrix4x4 getScale(Float const & sx, Float const & sy, Float const & sz)
		{
			return Matrix4x4(makeVector(sx		  ,Float(0.0),Float(0.0),Float(0.0)), 
							 makeVector(Float(0.0),sy	     ,Float(0.0),Float(0.0)),
							 makeVector(Float(0.0),Float(0.0),sz		,Float(0.0)),
							 makeVector(Float(0.0),Float(0.0),Float(0.0),Float(1.0))
				);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	static Matrix4x4 Matrix4x4::getRotation(Vector<Float,3> const & axis, Float angle)
		///
		/// \brief	Gets the rotation matrix corresponding axis angle rotation.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	02/03/2016
		///
		/// \param	axis 	The axis.
		/// \param	angle	The angle.
		///
		/// \return	The matrix corresponding to the axis angle rotation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		static Matrix4x4 getRotation(Vector<Float,3> const & axis, Float angle)
		{
			Math::Vector<Float, 3> normalized = axis.normalized() ;
			Float c = Float(cos(angle)), s = Float(sin(angle)), t = Float(1)-c, x = normalized[0], y = normalized[1], z=normalized[2] ;
			return Matrix4x4(makeVector(t*x*x + c  , t*x*y - z*s, t*x*z + y*s, Float(0.0)),
							 makeVector(t*x*y + z*s, t*y*y + c  , t*y*z - x*s, Float(0.0)),
							 makeVector(t*x*z - y*s, t*y*z + x*s, t*z*z + c  , Float(0.0)),
							 makeVector(Float(0.0) , Float(0.0) , Float(0.0) , Float(1.0))
				) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	static Matrix4x4 Matrix4x4::getNormalTransform(const Matrix4x4 & transformVertex)
		///
		/// \brief	If transformVertex is a vertex transformation matrix, the returned matrix is the
		/// 		associated normal transformation matrix.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/03/2016
		///
		/// \param	transformVertex	The transformVertex.
		///
		/// \return	The normal transform.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		static Matrix4x4 getNormalTransform(const Matrix4x4 & transformVertex)
		{
			return transformVertex.inverse().transpose() ;
		}
	};
}


#endif
