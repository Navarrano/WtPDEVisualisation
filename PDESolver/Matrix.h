
#ifndef MATRIX
#define MATRIX

#include "point.h"
template<class T>
class FMatrix;


template<class T>
class Matrix : public Vector<Vector<T> > {
private:
	int nrows;	// start row index
	int ncols;	 // end row index
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("class double")) 
			return std::string("class Matrix<double>");
		else {
			std::string s(typeid(T).name()), s1("class Matrix<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	Matrix();
	Matrix(int R, int C);	
	Matrix(int R, int C, const T& val);
	Matrix(const Vector<T>& v);
	Matrix(const FMatrix<T>& mat);
	Vector<T> GetRow(int index) const;
	Vector<T> GetCol(int index) const;
	Vector<T> GetCol(int indrow, int indcol) const;
	Vector<T> GetRow(int indrow, int indcol) const;
	Matrix<double> GetX() const;
	Matrix<double> GetY() const;
	Matrix<double> GetZ() const;
	Matrix<T>& operator=(double d);
	int GetNumRows() const { return nrows;}
	int GetNumCols() const { return ncols;}
	Matrix(const Matrix<double>& m1, const Matrix<double>& m2, const Matrix<double>& m3);
	Matrix(const Matrix<double>& m1, const Matrix<double>& m2);
	template<class T1>
	Matrix(const Matrix<T1>& m)
	{
		Matrix<T> m1(m.GetNumRows(),m.GetNumCols());	
		for (int i=0; i<m.GetNumRows(); i++) 
			for (int j=0; j<m.GetNumCols(); j++) m1[i][j] = m[i][j];
		*this=m1;
	}
	Vector<T> CreateVector() const;
	Matrix<T>& InsertRow(const Vector<T>& v, int index);
	Matrix<T>& InsertCol(const Vector<T>& v, int index);
	Matrix<T>& InsertRow(const Vector<T>& v, int indrow, int indcol);
	Matrix<T>& InsertCol(const Vector<T>& v, int indrow, int indcol);
	static bool CheckIndicesMult(const Matrix<T>& t1, const Matrix<T>& t2);
	static bool CheckIndicesAdd(const Matrix<T>& t1, const Matrix<T>& t2);
	Vector<T> CreateKroneckerVector() const;
	friend Matrix<T> operator*(const Matrix<T>& mat1, double d);
	friend Matrix<T> operator*(double d, const Matrix<T>& mat1);
	friend Matrix<double> operator*(const Matrix<double>& mat1, double d);
	friend Matrix<double> operator*(double d, const Matrix<double>& mat1);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};


template<class T>
bool Matrix<T>::CheckIndicesMult(const Matrix<T>& t1, const Matrix<T>& t2) 
{
	return (t1.ncols == t2.nrows);
}

template<class T>
bool Matrix<T>::CheckIndicesAdd(const Matrix<T>& t1, const Matrix<T>& t2) 
{
	return (t1.nrows == t2.nrows && t1.ncols == t2.ncols);
}


template<class T>
Matrix<T>::Matrix() : Vector<Vector<T> >(), nrows(0), ncols(0) {}


template<class T>
Matrix<T>::Matrix(int R, int C) : Vector<Vector<T> >(R, Vector<T>(C)), nrows(R), ncols(C) 
{
}


template<class T>
Matrix<T>::Matrix(const Vector<T>& v) : Vector<Vector<T> >(1, v), nrows(1), ncols(v.GetNum()) {
}


template<class T>
Matrix<T>::Matrix(int R, int C, const T& val) : Vector<Vector<T> >(R, Vector<T>(C, val)), nrows(R), ncols(C) {
}

template<class T>
Matrix<T>::Matrix(const FMatrix<T>& mat) : Vector<Vector<T> >(mat.GetNumRows(), Vector<T>(mat.GetNumCols())), nrows(mat.GetNumRows()),
										ncols(mat.GetNumCols())
{
	for (int i=mat.GetMinRow(); i<=mat.GetMaxRow(); i++)
		for (int j=mat.GetMinCol(); j<=mat.GetMaxCol(); j++) (*this)[i-mat.GetMinRow()][j-mat.GetMinCol()]=mat[i][j];
}


template<class T>
Matrix<T>& Matrix<T>::operator=(double d) 
{
	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++) (*this)[i][j]=d;

	return *this;
}


template<class T>
Vector<T> Matrix<T>::CreateVector() const
{
	Vector<T> v(nrows*ncols+1,0.0);

	int count=1;
	for (int j=1; j<ncols; j++) 
		for (int i=1; i<nrows; i++) {
			v[count]= (*this)[i][j];
			count++;
		}

	return v;
}

template<class T>
Vector<T> Matrix<T>::CreateKroneckerVector() const
{
	Vector<T> v(nrows*ncols,0.0);

	int count=0;
	for (int j=0; j<ncols; j++) 
		for (int i=0; i<nrows; i++) {
			v[count]= (*this)[i][j];
			count++;
		}

	return v;
}

template<class T>
Vector<T> Matrix<T>::GetRow(int index) const
{
	Vector<T> v(ncols);
	for (int j=0; j<ncols; j++) v[j] = (*this)[index][j];
	return v;
}

template<class T>
Vector<T> Matrix<T>::GetCol(int index) const
{
	Vector<T> v(nrows);
	for (int i=0; i<nrows; i++) v[i] = (*this)[i][index];
	return v;
}


template<class T>
Vector<T> Matrix<T>::GetRow(int indrow, int indcol) const
{
	Vector<T> v(ncols-indcol);
	for (int j=indcol; j<ncols; j++) v[j] = (*this)[indrow][j];
	return v;
}

template<class T>
Vector<T> Matrix<T>::GetCol(int indrow, int indcol) const
{
	Vector<T> v(nrows-indrow);
	for (int i=indrow; i<nrows; i++) v[i] = (*this)[i][indcol];
	return v;
}



template<class T>
Matrix<T>& Matrix<T>::InsertRow(const Vector<T>& v, int index)
{
	for (int j=0; j<ncols; j++) (*this)[index][j]=v[j];
	return *this;
}


template<class T>
Matrix<T>& Matrix<T>::InsertCol(const Vector<T>& v, int index)
{
	for (int i=0; i<nrows; i++) (*this)[i][index] = v[i];
	return *this;
}



template<class T>
Matrix<T>& Matrix<T>::InsertRow(const Vector<T>& v, int indrow, int indcol)
{
	for (int j=indcol; j<ncols; j++) (*this)[indrow][j]=v[j];
	return *this;
}


template<class T>
Matrix<T>& Matrix<T>::InsertCol(const Vector<T>& v, int indrow, int indcol)
{
	for (int i=indrow; i<nrows; i++) (*this)[i][indcol] = v[i];
	return *this;
}



template <class T>
void Matrix<T>::read(std::istream& is)
{
	Vector<Vector<T> >::read(is);
}



template <class T>
void Matrix<T>::readfile(std::ifstream& ifs)
{
	Vector<Vector<T> >::readfile(ifs);
}

template <class T>
void Matrix<T>::writefile(std::ofstream& ofs)
{
	Vector<Vector<T> >::writefile(ofs);
}


template <class T>
void Matrix<T>::write(std::ostream& os) 
{
	Vector<Vector<T> >::write(os);
}



#endif