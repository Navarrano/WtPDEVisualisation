
#ifndef FMATRIX
#define FMATRIX

#include "point.h"

template<class T>
class FMatrix : public FVector<FVector<T> > {
private:
	int r1;	// start row index
	int r2;	 // end row index
	int c1;  // start col index
	int c2;  // end col index
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("class double")) 
			return std::string("class FMatrix<double>");
		else {
			std::string s(typeid(T).name()), s1("class FMatrix<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	FMatrix();
	FMatrix(int R1, int R2, int C1, int C2);
	FMatrix(int R1, int R2, int C1, int C2, const T& val);
	FMatrix(int R, int C, const T& val);
	FVector<T> GetRow(int index) const;
	FVector<T> GetCol(int index) const;
	int GetNumRows() const { return r2-r1+1;}
	int GetNumCols() const { return c2-c1+1;}
	int GetMinRow() const { return r1; }
	int GetMaxRow() const { return r2; }
	int GetMinCol() const { return c1; }
	int GetMaxCol() const { return c2; }
	static bool CheckIndices(const FMatrix<T>& t1, const FMatrix<T>& t2);
	friend FMatrix<T> operator*(const FMatrix<T>& mat1, double d);
	friend FMatrix<T> operator*(double d, const FMatrix<T>& mat1);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

template<class T>
bool FMatrix<T>::CheckIndices(const FMatrix<T>& t1, const FMatrix<T>& t2) 
{
	return (t1.c1 == t2.r1 && t1.c2 == t2.r2);
}


template<class T>
FMatrix<T>::FMatrix() : FVector<FVector<T> >(), r1(0), r2(0), c1(0), c2(0) {}

template<class T>
FMatrix<T>::FMatrix(int R1, int R2, int C1, int C2): FVector<FVector<T> >(R1, R2, FVector<T>(C1,C2)), r1(R1), r2(R2), c1(C1), c2(C2) {
}

/*template<class T>
FMatrix<T>::FMatrix(int R, int C) : FVector<FVector<T> >(R, FVector<T>(C)), r1(0), r2(R-1), c1(0), c2(C-1) {
}*/

template<class T>
FMatrix<T>::FMatrix(int R1, int R2, int C1, int C2, const T& val): FVector<FVector<T> >(R1, R2, FVector<T>(C1, C2, val)), r1(R1), r2(R2), c1(C1), c2(C2) {
}

template<class T>
FMatrix<T>::FMatrix(int R, int C, const T& val) : FVector<FVector<T> >(0, R-1, FVector<T>(0, C-1, val)), r1(0), r2(R-1), c1(0), c2(C-1) {
}

template<class T>
FVector<T> FMatrix<T>::GetRow(int index) const
{
	FVector<T> v(c1, c2);
	for (int j=c1; j<=c2; j++) v[j] = (*this)[index][j];
	return v;
}

template<class T>
FVector<T> FMatrix<T>::GetCol(int index) const
{
	FVector<T> v(r1, r2);
	for (int i=r1; i<=r2; i++) v[i] = (*this)[i][index];
	return v;
}


template<class T>
FMatrix<T> operator*(const FMatrix<T>& mat1, double d)
{
	FMatrix<T> mat(mat1);

	for (int i=mat1.r1; i<=mat.r2; i++)
		for (int j=mat.c1; j<=mat.c2; j++) mat[i][j] = mat1[i][j]*d;
	return mat;
}

template<class T>
FMatrix<T> operator*(double d, const FMatrix<T>& mat1)
{
	FMatrix<T> mat(mat1);

	for (int i=mat1.r1; i<=mat.r2; i++)
		for (int j=mat.c1; j<=mat.c2; j++) mat[i][j] = mat1[i][j]*d;
	return mat;
}


template <class T>
void FMatrix<T>::read(std::istream& is)
{
	if (r2 ==  0 && c2 == 0) {
		int ncol, nrow;

		std::cout << "How many rows in the FMatrix?\n";
		is >> nrow;
		std::cout << "How many columns in the FMatrix?\n";
		is >> ncol;
		FMatrix<T> m(0,nrow-1,0,ncol-1);
		std::cout << "input your FMatrix row by row\n";
		for (int i=m.r1; i<=m.r2; i++) {
			std::cout << "row" << i << "\n";
			for (int j=m.c1; j<=m.c2; j++)
				is >> m[i][j];
		}
		*this = m;
	} else {
		for (int i=r1; i<=r2; i++) {
			std::cout << "row" << i << "\n";
			for (int j=c1; j<=c2; j++)
				is >> (*this)[i][j];
		}
	}
}



template <class T>
void FMatrix<T>::readfile(std::ifstream& ifs)
{
	if (r2 == 0 && c2 == 0) {
		int nrow, ncol;

		ifs >> nrow;
		ifs >> ncol;
	
		FMatrix<T> m(0,nrow-1,0,ncol-1);
		for (int i=m.r1; i<=m.r2; i++) {
			for (int j=m.c1; j<=m.c2; j++)
			ifs >> m[i][j];
		}
		*this = m;
	} else {
		for (int i=r1; i<=r2; i++) 
			for (int j=c1; j<=c2; j++)
				ifs >> (*this)[i][j];
	}	
}


template <class T>
void FMatrix<T>::write(std::ostream& os) 
{
	os << "Your FMatrix is" << "\n";
	for (int i=r1; i<=r2; i++) {
		Vector<T> v(GetRow(i));
		v.write(os);
	}
}


template <class T>
void FMatrix<T>::writefile(std::ofstream& ofs) 
{
	ofs << "Your FMatrix is" << "\n";
	for (int i=r1; i<=r2; i++) {
		Vector<T> v(GetRow(i));
		v.writefile(ofs);
	}
}


#endif