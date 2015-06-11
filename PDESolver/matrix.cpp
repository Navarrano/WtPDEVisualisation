
#include "matrix.h"


template<class T>
Matrix<T> operator*(const Matrix<T>& mat1, double d)
{
	Matrix<T> mat(mat1);

	for (int i=0; i<mat.nrows; i++)
		for (int j=0; j<mat.ncols; j++) mat[i][j] = mat1[i][j]*d;
	return mat;
}

template<class T>
Matrix<T> operator*(double d, const Matrix<T>& mat1)
{
	Matrix<T> mat(mat1);

	for (int i=0; i<mat.nrows; i++)
		for (int j=0; j<=mat.ncols; j++) mat[i][j] = mat1[i][j]*d;
	return mat;
}


Matrix<double> operator*(const Matrix<double>& mat1, double d)
{
	Matrix<double> mat(mat1);

	for (int i=0; i<mat.nrows; i++)
		for (int j=0; j<mat.ncols; j++) mat[i][j] = mat1[i][j]*d;
	return mat;
}


Matrix<double> operator*(double d, const Matrix<double>& mat1)
{
	Matrix<double> mat(mat1);

	for (int i=0; i<mat.nrows; i++)
		for (int j=0; j<mat.ncols; j++) mat[i][j] = mat1[i][j]*d;
	return mat;
}


template<>
Matrix<double> Matrix<Point1D>::GetX() const
{
	Matrix<double> m(nrows,ncols);

	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			m[i][j] = (*this)[i][j].GetX();

	return m;
}

template<>
Matrix<double> Matrix<Point2D>::GetX() const
{
	Matrix<double> m(nrows,ncols);

	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			m[i][j] = (*this)[i][j].GetX();

	return m;
}

template<>
Matrix<double> Matrix<Point2D>::GetY() const
{
	Matrix<double> m(nrows,ncols);

	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			m[i][j] = (*this)[i][j].GetY();

	return m;
}

template<>
Matrix<double> Matrix<Point3D>::GetX() const
{
	Matrix<double> m(nrows,ncols);

	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			m[i][j] = ((*this)[i][j]).GetX();

	return m;
}

template<>
Matrix<double> Matrix<Point3D>::GetY() const
{
	Matrix<double> m(nrows,ncols);

	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			m[i][j] = ((*this)[i][j]).GetY();

	return m;
}


template<>
Matrix<double> Matrix<Point3D>::GetZ() const
{
	Matrix<double> m(nrows,ncols);

	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			m[i][j] = ((*this)[i][j]).GetZ();

	return m;
}


template<>
Matrix<Point2D>::Matrix(const Matrix<double>& m1, const Matrix<double>& m2)
{
	Matrix<Point2D> m(m1.GetNumRows(),m1.GetNumCols());	

	for (int i=0; i<m1.GetNumRows(); i++) 
		for (int j=0; j<m1.GetNumCols(); j++) m[i][j] = Point2D(m1[i][j],m2[i][j]);

	*this=m;
}

template<>
Matrix<Point3D>::Matrix(const Matrix<double>& m1, const Matrix<double>& m2, const Matrix<double>& m3)
{
	Matrix<Point3D> m(m1.GetNumRows(),m1.GetNumCols());	
	for (int i=0; i<m1.GetNumRows(); i++) 
		for (int j=0; j<m1.GetNumCols(); j++) m[i][j] = Point3D(m1[i][j],m2[i][j],m3[i][j]);

	*this=m;
}
