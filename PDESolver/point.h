#ifndef POINT
#define POINT

#include "vector.h"

class Point1D : public TextObject, public FTextObject {
	double x;
	virtual ObjectID Identity() const  { return std::string("class Point1D"); }
public:
	Point1D();
	Point1D(double X);
	Point1D(const Point2D& pt);
	Point1D(const Point3D& pt);
	double GetX() const;
	operator double() const;
	Point1D& operator=(const Point1D& p);
	Point1D& operator=(double d);
	operator FVector<double>() const;
	Point1D(const FVector<double>& v);
	Point1D(const Vector<double>& v);
	friend Point1D operator*(double d, const Point1D& p);
	friend Point1D operator*(const Point1D& p, double d);
	friend Point1D operator*(const Point1D& p1, const Point1D& p2);
	friend Point1D operator*(const Point1D& p, double d);
	friend Point1D operator-(const Point1D& p, double d);
	friend Point1D operator-(double d, const Point1D& p);
	friend Point1D operator-(const Point1D& p1, const Point1D& p2);
	friend Point1D operator+(double d, const Point1D& p);
	friend Point1D operator+(const Point1D& p, double d);
	friend Point1D operator+(const Point1D& p1, const Point1D& p2);
	friend Point1D operator-(const Point1D& p, double d);
	friend Point1D operator-(double d, const Point1D& p);
	friend Point1D operator-(const Point1D& p1, const Point1D& p2);
	friend Point1D operator/(const Point1D& p, double d);
	friend bool operator>(const Point1D& p, double d);
	friend bool operator>(const Point1D& p1, const Point1D& p2);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};


class Point2D : public TextObject, public FTextObject {
	double x;
	double y;
	virtual ObjectID Identity() const { return std::string("class Point2D"); }
public:
	Point2D();
	Point2D(double X, double Y);
	Point2D(double d);
	Point2D(const Point1D& pt);
	Point2D(const Point3D& pt);
	double GetX() const;
	double GetY() const;
	operator double() const;
	Point2D& operator=(const Point2D& p);
	Point2D& operator=(const Point3D& p);
	Point2D& operator=(double d);
	operator FVector<double>() const;
	Point2D(const FVector<double>& v);
	Point2D(const Vector<double>& v);
	friend Point2D operator*(double d, const Point2D& p);
	friend Point2D operator*(const Point2D& p, double d);
	friend Point2D operator*(const Point2D& p1, const Point2D& p2);
	friend Point2D operator*(const Point2D& p, double d);
	friend Point2D operator*(const Point2D& p1, const Point1D& p2);
	friend Point2D operator+(const Point2D& p, double d);
	friend Point2D operator+(double d, const Point2D& p);
	friend Point2D operator+(const Point2D& p1, const Point2D& p2);
	friend Point2D operator-(const Point2D& p, double d);
	friend Point2D operator-(double d, const Point2D& p);
	friend Point2D operator-(const Point2D& p1, const Point2D& p2);
	friend Point2D operator/(const Point2D& p, double d);
	friend bool operator>(const Point2D& p, double d);
	friend bool operator>(const Point2D& p1, const Point2D& p2);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

class Point3D :  public TextObject, public FTextObject {
	double x;
	double y;
	double z;
	virtual ObjectID Identity() const { return std::string("class Point3D"); }
public:
	Point3D();
	Point3D(double X, double Y, double Z);
	Point3D(double d);
	Point3D(const FVector<double>& v);
	Point3D(const Vector<double>& v);
	Point3D(const Point2D& pt);
	Point3D(const Point1D& pt);
	Point3D(const Point4D& pt);
	double GetX() const;
	double GetY() const;
	double GetZ() const;
	operator FVector<double>() const;
	operator double() const;
	Point3D& operator=(const Point3D& p);
	Point3D& operator=(const FVector<double>& vec);
	Point3D& operator=(const Point4D& p);
	Point3D& operator=(double d);
	double Product(const Point3D& p) const;
	Point3D operator-();
	double Length() const;
	friend Point3D operator*(double d, const Point3D& p);
	friend Point3D operator*(const Point3D& p, double d);
	friend Point3D operator*(const Point3D& p1, const Point3D& p2);
	friend Point3D operator*(const Point3D& p1, const Point2D& p2);
	friend Point3D operator*(const Point3D& p1, const Point1D& p2);
	friend Point3D operator*(const Point3D& p, double d);
	friend Point3D operator+(const Point3D& p, double d);
	friend Point3D operator+(double d, const Point3D& p);
	friend Point3D operator+(const Point3D& p1, const Point3D& p2);
	friend Point3D operator+=(const Point3D& p1, const Point3D& p2);
	friend Point3D operator-(const Point3D& p, double d);
	friend Point3D operator-(double d, const Point3D& p);
	friend Point3D operator-(const Point3D& p1, const Point3D& p2);
	friend Point3D operator-=(const Point3D& p1, const Point3D& p2);
	friend Point3D operator/(const Point3D& p, double d);
	friend bool operator>(const Point3D& p, double d);
	friend bool operator>(const Point3D& p1, const Point3D& p2);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

class Point4D :  public TextObject, public FTextObject {
	double x;
	double y;
	double z;
	double w;
	virtual ObjectID Identity() const { return std::string("class Point4D"); }
public:
	Point4D();
	Point4D(double X, double Y, double Z, double W);
	Point4D(double d);
	Point4D(const FVector<double>& v);
	Point4D(const Vector<double>& v);
	Point4D(const Point3D& pt);
	Point4D(const Point2D& pt);
	Point4D(const Point1D& pt);
	double GetX() const;
	double GetY() const;
	double GetZ() const;
	double GetW() const;
	operator FVector<double>() const;
	operator double() const;
	Point4D& operator=(const Point4D& p);
	Point4D& operator=(const Point3D& p);
	Point4D& operator=(double d);
	double Product(const Point4D& p) const;
	friend Point4D operator*(double d, const Point4D& p);
	friend Point4D operator*(const Point4D& p, double d);
	friend Point4D operator*(const Point4D& p1, const Point4D& p2);
	friend Point4D operator*(const Point4D& p1, const Point3D& p2);
	friend Point4D operator*(const Point4D& p1, const Point2D& p2);
	friend Point4D operator*(const Point4D& p1, const Point1D& p2);
	friend Point4D operator*(const Point4D& p, double d);
	friend Point4D operator+(const Point4D& p, double d);
	friend Point4D operator+(double d, const Point4D& p);
	friend Point4D operator+(const Point4D& p1, const Point4D& p2);
	friend Point4D operator+=(const Point4D& p1, const Point4D& p2);
	friend Point4D operator-(const Point4D& p, double d);
	friend Point4D operator-(double d, const Point4D& p);
	friend Point4D operator-(const Point4D& p1, const Point4D& p2);
	friend Point4D operator-=(const Point4D& p1, const Point4D& p2);
	friend Point4D operator/(const Point4D& p, double d);
	friend bool operator>(const Point4D& p, double d);
	friend bool operator>(const Point4D& p1, const Point4D& p2);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};


template<class T>
class FVector : public TextObject, public FTextObject {
private:
	int min;            // Minimum (base) index
	int max;            // Maximum index
	T *pdata;           // Pointer to first element
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("class double"))
			return std::string("class FVector<double>");
		else {
			std::string s(typeid(T).name()), s1("class FVector<");
			return s1.append(s,0,s.size())+">";
		}
	}
protected:
	void Init(int Min, int Max, const T& val);
	void Init(int Min, int Max);
public:
	FVector();
	FVector(int Min, int Max);
	FVector(int Min, int Max, const T& val);
	FVector(int Num);
	FVector(const FVector<T>& copy);
	FVector(const Vector<T>& copy);
	FVector(double d);
	FVector(const Point1D& p);
	FVector(const Point2D& p);
	FVector(const Point3D& p);
	operator double() const;
	operator Point3D() const;
	operator Point2D() const;
	operator Point1D() const;
	int GetMin() const;
	int GetMax() const;
	int GetNum() const;
	~FVector();
	operator T*();
	FVector<T>& operator=(const FVector<T>& v);
	T& operator[](int i);
	const T& operator[](int i) const;
	T& operator()(int i) const;
	friend FVector<T> operator*(const FVector<T>& vec, double d);
	friend FVector<T> operator*(double d, const FVector<T>& vec);
	friend FVector<T> operator+(const FVector<T>& vec1, const FVector<T>& vec2);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};


template<class T>
FVector<T>::operator Point3D() const
{
	return Point3D((*this)[min],(*this)[min+1],(*this)[min+2]);
}

template<class T>
FVector<T>::operator Point2D() const
{
	return Point2D((*this)[min],(*this)[min+1]);
}

template<class T>
FVector<T>::operator Point1D() const
{
	return Point1D((*this)[min]);
}

template<class T>
FVector<T>::operator double() const
{
	return (*this)[min];
}

template<class T>
FVector<T>::FVector(double d)
{
	Init(1,1,d);
}

template<class T>
FVector<T>::FVector(const Point1D& p)
{
	Init(1,1,p.GetX());
}


template<class T>
FVector<T>::FVector(const Point2D& p)
{
	Init(1,2);
	pdata[0] = p.GetX();
	pdata[1] = p.GetY();
}

template<class T>
FVector<T>::FVector(const Point3D& p)
{
	Init(1,3);
	pdata[0] = p.GetX();
	pdata[1] = p.GetY();
	pdata[2] = p.GetZ();
}


template<class T>
FVector<T>::FVector() : min(-1), max(0), pdata(0) {}

template<class T>
void FVector<T>::Init(int Min, int Max, const T& val)
{
	Init(Min, Max);
	for (int i=0; i<=max-min; i++) pdata[i] = val;
}

template<class T>
void FVector<T>::Init(int Min, int Max)
{
	min = Min;
	max = Max;
	if (min > max) {
	    min=max=0;
		pdata=0;  // Object construction failed! default created
	} else pdata = new T[max-min+1];  // Allocate memory to FVector	
}


template <class T>
void FVector<T>::write(std::ostream& os) 
{
	for (int i=min; i<=max; i++) {
		os << (*this)[i];
		os << " ";
	}
	os << "\n";
//	std::copy(v.begin(),v.end(), std::ostream_iterator<T>(os, " "));
}


template <class T>
void FVector<T>::read(std::istream& is)
{
	if (min == -1) {
		int num;
		std::cout << "How many elements in your FVector" << "\n";
		is >> num;
		Init(0, num-1);
		std::cout << "Input your FVector" << "\n";
		for (int i=0; i<num; i++) is >> pdata[i];
	} else for (int i=0; i<max-min+1; i++) is >> pdata[i];
} 


template <class T>
void FVector<T>::writefile(std::ofstream& ofs) 
{
	for (int i=min; i<=max; i++) {
		ofs << (*this)[i]; 
		ofs << " ";
	}
	ofs << "\n";
}

template <class T>
void FVector<T>::readfile(std::ifstream& ifs)
{
	if (min == -1) {
		int num;
		ifs >> num;

		Init(0, num-1);
		for (int i=0; i<max-min+1; i++) ifs >> pdata[i];
	} else for (int i=0; i<max-min+1; i++) ifs >> pdata[i];
} 


template<class T>
FVector<T>::FVector(int Min, int Max, const T& val)	{
	Init(Min, Max, val);
}

template<class T>
FVector<T>::FVector(int Min, int Max)	{
	Init(Min, Max);
}

template<class T>
FVector<T>::FVector(int Num)	{
	Init(0,Num-1);
}



template<class T>
FVector<T>::FVector(const FVector<T> &copy)	{
	Init(copy.GetMin(), copy.GetMax());
	for(int i=min; i<=max; i++)
		(*this)[i] = copy[i]; 
}

template<class T>
FVector<T>::FVector(const Vector<T> &copy)	{
	Init(0, copy.GetNum()-1);
	for(int i=min; i<=max; i++)
		(*this)[i] = copy[i]; 
}


template<class T>	
FVector<T>::~FVector()	{
	delete[] pdata;  // Note special form of delete[]
}

template<class T>
inline int FVector<T>::GetMin() const { return min;}

template<class T>
inline int FVector<T>::GetMax() const { return max; }

template<class T>
inline int FVector<T>::GetNum() const { return max-min+1;}


template<class T>
FVector<T>& FVector<T>::operator=(const FVector<T> &copy)
{
	if (this != &copy) {   // Can't copy self to self
		if (pdata) delete[] pdata;
		Init(copy.GetMin(), copy.GetMax());
		if (pdata) {
			for(int i=min; i<=max; i++) {
				(*this)[i] = copy[i];
			}
		}
	}
	return *this;
}

template<class T>
FVector<T>::operator T*()
{
	return pdata;
}


template<class T>
T& FVector<T>::operator() (int i) const
{
	if ((i < min) || (i > max)) 
		i = min;  // Range error causes index to equal base
	return pdata[i - min];
}



template<class T>
T& FVector<T>::operator[] (int i) 
{
	if ((i < min) || (i > max)) 
		i = min;  // Range error causes index to equal base
	return pdata[i - min];
}


template<class T>
const T& FVector<T>::operator[] (int i) const
{
	if ((i < min) || (i > max)) 
		i = min;  // Range error causes index to equal base
	return pdata[i - min];
}

template<class T>
FVector<T> operator+(const FVector<T>& vec1, const FVector<T>& vec2)
{
	FVector<T> res(vec1.min,vec1.max);
	for (int i=vec1.min ; i<=vec1.max; i++) res[i] = vec1[i]+vec2[i];
	return res;
}


template<class T>
FVector<T> operator*(double d, const FVector<T>& vec)
{
	FVector<T> res(vec.min, vec.max);

	for (int i=vec.min; i<=vec.max; i++) res[i] = d*vec[i];
	return res;
}

template<class T>
FVector<T> operator*(const FVector<T>& vec,double d)
{
	FVector<T> res(vec.min, vec.max);

	for (int i=vec.min; i<=vec.max; i++) res[i] = d*vec[i];
	return res;
}



#endif
