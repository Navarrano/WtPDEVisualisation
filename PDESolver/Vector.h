#ifndef VECTOR
#define VECTOR


#include "object.h"

template<class T>
class FVector;



class Point1D;
class Point2D;
class Point3D;
class Point4D;

template<class T>
class Vector : public TextObject, public FTextObject {
public:
	Vector();
	Vector(int n);
	Vector(int n, const T& t);
	Vector(const Vector<T>& v);
	Vector(const std::vector<T>&);
	Vector(const FVector<T>& v);
	std::vector<T> convert() const;
	~Vector();
	int GetNum() const;
	int size() const;
	Vector(const Point1D& p);
	Vector(const Point2D& p);
	Vector(const Point3D& p);
	Vector(const Point4D& p);

	Vector<T>& operator=(const Vector<T>& v);
	T& operator[](int i);
	const T& operator[](int i) const;
	T* begin() { return pdata; }
	T* end() { return pdata+num; }
	const T* begin() const { return pdata; }
	const T* end() const { return pdata+num; }
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
private:
	int num;			// num of elements
	T* pdata;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class Vector<double>");
		else {
			std::string s(typeid(T).name()), s1("class Vector<");
			return s1.append(s,0,s.size())+">";
		}
	}
};




template<class T>
Vector<T>::Vector() : num(0), pdata(0)
{
}


template<class T>
Vector<T>::Vector(int n) : num(n), pdata(new T[num])
{
}


template<class T>
Vector<T>::~Vector()
{
	delete [] pdata;
}

template<class T>
T& Vector<T>::operator[](int i)
{
	return pdata[i];
}

template<class T>
const T& Vector<T>::operator[](int i) const
{
	return pdata[i];
}


template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& v)
{
	if (this == &v) return *this;
	if (v.num) {
		num = v.num;
		delete [] pdata;
		pdata = new T[num];
		for (int i=0; i<num; i++) pdata[i] = v.pdata[i];
	}
	return *this;
}

template<class T>
Vector<T>::Vector(int n, const T& t) : num(n), pdata(new T[num])
{
	for (int i=0; i<num; i++) pdata[i] = t;
}

template<class T>
Vector<T>::Vector(const Point1D& p) : num(1), pdata(new T[1])
{
	(*this)[0] = p.GetX();
}

template<class T>
Vector<T>::Vector(const Point2D& p) : num(2), , pdata(new T[2])
{
	(*this)[0] = p.GetX();
	(*this)[1] = p.GetY();
}

template<class T>
Vector<T>::Vector(const Point3D& p) : num(3), pdata(new T[3])
{
	(*this)[0] = p.GetX();
	(*this)[1] = p.GetY();
	(*this)[2] = p.GetZ();
}


template<class T>
Vector<T>::Vector(const Point4D& p) : num(4), pdata(new T[4])
{
	(*this)[0] = p.GetX();
	(*this)[1] = p.GetY();
	(*this)[2] = p.GetZ();
	(*this)[3] = p.GetW();
}

template<class T>
Vector<T>::Vector(const Vector<T>& v) : num(v.num), pdata(new T[num])
{
	for (int i=0; i<num; i++) pdata[i]=v.pdata[i];
}

template<class T>
Vector<T>::Vector(const std::vector<T>& v) : num(v.size()), pdata(new T[num])
{
	for (int i=0; i<num; i++) pdata[i]=v[i];
}



template<class T>
Vector<T>::Vector(const FVector<T>& v)
{
	Vector<T> w(v.GetNum());

	for (int i=v.GetMin(); i<=v.GetMax(); i++) w[i-v.GetMin()]=v[i];

	*this = w;
}

template<class T>
inline int Vector<T>::size() const
{
	return num;
}


template<class T>
inline int Vector<T>::GetNum() const
{
	return num;
}


template<class T>
std::vector<T> Vector<T>::convert() const
{
	std::vector<T> v(num);

	for (int i = 0; i < num; i++) v[i] = (*this)[i];
	return v;
}


template <class T>
void Vector<T>::write(std::ostream& os)
{	
	for (int i=0; i<num; i++) os << (*this)[i] << " ";
		os << "\n";
	//std::copy(begin(),end(),  std::ostream_iterator<T>(std::cout, " "));
}


template <class T>
void Vector<T>::read(std::istream& is)
{
		if (num == 0) {
		int n;
		std::cout << "How many elements in your vector" << "\n";
		is >> n;
		Vector<T> v(n);
		num = n;
		std::cout << "Input your vector" << "\n";
		for (int i=0; i<num; i++) is >> v[i];
		*this = v;
	} else {
		std::string s(typeid(T).name());
		std::cerr << "input vector of " << s  << "\n";
		Vector<T> v(num);
		for (int i=0; i<num; i++) is >> v[i];
		*this = v;
	}
} 


template <class T>
void Vector<T>::writefile(std::ofstream& ofs)
{
//	ofs << num;
	for (int i=0; i<num; i++) {
		ofs << (*this)[i] << " ";
	}
	ofs << "\n";
}


template <class T>
void Vector<T>::readfile(std::ifstream& ifs)
{
	if (num == 0) {
		int n;
		ifs >> n;
		Vector<T> v(n);
		num = n;
		for (int i=0; i<num; i++) ifs >> v[i];
//		reserve(n);
		*this = v;
	} else {
		//ifs >> num;
		for (int i=0; i<num; i++) ifs >> (*this)[i];
	}
} 




#endif