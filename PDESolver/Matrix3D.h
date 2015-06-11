#ifndef MAT3D
#define MAT3D

#include "matrix.h"


template<class T>
class Matrix3D : public Vector<Matrix<T> >
{
private:
	int nrows, ncols, npols;
	virtual ObjectID Identity() const 
	{ 
			if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class Matrix3D<double>");
		else {
			std::string s(typeid(T).name()), s1("class Matrix3D<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	Matrix3D();
	Matrix3D(int R, int C, int H);
	Matrix3D(int R, int C, int H, const T& val);
	Matrix3D(const Vector<Matrix<T> >& v, int num);
	template<class T1>
	Matrix3D(const Matrix3D<T1>& m)
	{
		Matrix3D<T> m1(m.GetNumRows(),m.GetNumCols(),m.GetNumPols());	
		for (int i=0; i<m.GetNumRows(); i++) 
			for (int j=0; j<m.GetNumCols(); j++) 
				for (int k=0; k<m.GetNumPols(); k++) m1[k][i][j] = m[k][i][j];
		*this=m1;
	}
	Matrix<T> GetUV(int index) const;
	Matrix<T> GetUW(int index) const;
	Matrix<T> GetVW(int index) const;
	Vector<T> GetU(int indv, int indw) const;
	Vector<T> GetV(int indu, int indw) const;
	Vector<T> GetW(int indu, int indv) const;
	Matrix3D<T>& InsertVW(const Matrix<T>& v, int index);	
	Matrix3D<T>& InsertUV(const Matrix<T>& v, int index);
	Matrix3D<T>& InsertUW(const Matrix<T>& v, int index);
	Matrix3D<T>& InsertU(const Vector<T>& v, int ind1, int ind2);
	Matrix3D<T>& InsertV(const Vector<T>& v, int ind1, int ind2);
	Matrix3D<T>& InsertW(const Vector<T>& v, int ind1, int ind2);
	int GetNumRows() const { return nrows; }
	int GetNumCols() const { return ncols; }
	int GetNumPols() const { return npols; }
	static bool CheckIndicesMult(const Matrix3D<T>& t1, const Matrix<double>& t2);
	static bool CheckIndicesMult(const Matrix<double>& t1, const Matrix3D<T>& t2);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

template<class T>
bool Matrix3D<T>::CheckIndicesMult(const Matrix3D<T>& t1, const Matrix<double>& t2) 
{
	return (t1.GetNumCols() == t2.GetNumRows() || t1.GetNumPols() == t2.GetNumRows());
}

template<class T>
bool Matrix3D<T>::CheckIndicesMult(const Matrix<double>& t1, const Matrix3D<T>& t2) 
{
	return (t1.GetNumCols() == t2.GetNumRows() || t1.GetNumCols() == t2.GetNumPols());
}


template<class T>
Matrix3D<T>::Matrix3D(const Vector<Matrix<T> >& v, int num) : Vector<Matrix<T> >(num,Matrix<T>(v[0].GetNumRows(),v[0].GetNumCols())),
nrows(v[0].GetNumRows()), ncols(v[0].GetNumCols()), npols(num) 

{
	for (int i=0; i<num; i++) (*this)[i] = v[i];
}

template<class T>
Matrix3D<T>::Matrix3D() :
Vector<Matrix<T> >(), nrows(0), ncols(0), npols(0) 
{
}



template<class T>
Matrix3D<T>::Matrix3D(int R, int C, int H, const T& val) :
Vector<Matrix<T> >(H, Matrix<T>(R,C,val)), nrows(R), ncols(C), npols(H) {
}


template<class T>
Matrix3D<T>::Matrix3D(int R, int C, int H) :
Vector<Matrix<T> >(H, Matrix<T>(R,C)), nrows(R), ncols(C), npols(H) {
}


template<class T>
Matrix<T> Matrix3D<T>::GetUV(int index) const
{
	Matrix<T> uv(nrows, ncols);
	for (int i=0; i<nrows; i++)
	 	for (int j=0; j<ncols; j++) uv[i][j] = (*this)[index][i][j];
	return uv;
}

template<class T>
Matrix<T> Matrix3D<T>::GetUW(int index) const
{
	Matrix<T> uw(nrows, npols);
	for (int i=0; i<nrows; i++)
	 	for (int k=0; k<npols; k++) uw[i][k] = (*this)[k][i][index];
	return uw;
}

template<class T>
Matrix<T> Matrix3D<T>::GetVW(int index) const
{
	Matrix<T> vw(ncols,npols);
	for (int j=0; j<ncols; j++)
	 	for (int k=0; k<npols; k++) vw[j][k] = (*this)[k][index][j];
	return vw;
}

template<class T>
Vector<T> Matrix3D<T>::GetU(int indv, int indw) const
{
	Vector<T> u(nrows);
	for (int i=0; i<nrows; i++) u[i] = (*this)[indw][i][indv];
	return u;
}


template<class T>
Vector<T> Matrix3D<T>::GetV(int indu, int indw) const
{
	Vector<T> v(ncols);
	for (int j=0; j<ncols; j++) v[j] = (*this)[indw][indu][j];
	return v;
}

template<class T>
Vector<T> Matrix3D<T>::GetW(int indu, int indv) const
{
	Vector<T> w(npols);
	for (int k=0; k<npols; k++) w[k] = (*this)[k][indu][indv];
	return w;
}


template<class T>
Matrix3D<T>& Matrix3D<T>::InsertUV(const Matrix<T>& v, int index)
{
	for (int i=0; i<nrows; i++) 
		for (int j=0; j<npols; j++) (*this)[index][i][j]=v[i][j];
	return *this;
}


template<class T>
Matrix3D<T>& Matrix3D<T>::InsertVW(const Matrix<T>& v, int index)
{
	for (int j=0; j<ncols; j++) 
		for (int k=0; k<npols; k++) (*this)[k][index][j]=v[j][k];
	return *this;
}

template<class T>
Matrix3D<T>& Matrix3D<T>::InsertUW(const Matrix<T>& v, int index)
{
	for (int i=0; i<nrows; i++) 
		for (int k=0; k<npols; k++) (*this)[k][i][index]=v[i][k];
	return *this;
}

template<class T>
Matrix3D<T>& Matrix3D<T>::InsertU(const Vector<T>& v, int ind1, int ind2)
{
	for (int i=0; i<nrows; i++) (*this)[ind2][i][ind1]=v[i];
	return *this;
}

template<class T>
Matrix3D<T>& Matrix3D<T>::InsertV(const Vector<T>& v, int ind1, int ind2)
{
	for (int j=0; j<ncols; j++) (*this)[ind2][ind1][j]=v[j];
	return *this;
}


template<class T>
Matrix3D<T>& Matrix3D<T>::InsertW(const Vector<T>& v, int ind1, int ind2)
{
	for (int k=0; k<npols; k++) (*this)[k][ind1][ind2]=v[k];
	return *this;
}


template <class T>
void Matrix3D<T>::read(std::istream& is)
{
	Vector<Matrix<T> >::read(is);
}


template <class T>
void Matrix3D<T>::readfile(std::ifstream& ifs)
{
	Vector<Matrix<T> >::readfile(ifs);
}


template <class T>
void Matrix3D<T>::write(std::ostream& os) 
{
	Vector<Matrix<T> >::write(os);
}


template <class T>
void Matrix3D<T>::writefile(std::ofstream& ofs) 
{
	Vector<Matrix<T> >::writefile(ofs);
}



#endif

