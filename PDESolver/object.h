#ifndef __OBJECT
#define __OBJECT
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <ctime>
#include <algorithm>
#include "ptr.h"


#define M1 259200l
#define IA1 7141
#define IC1 54773
#define RM1 1.0/M1
#define M2 134456l
#define IA2 8121
#define IC2 28411
#define RM2 1.0/M2
#define M3 243000l
#define IA3 4561
#define IC3 51349


template<class T>
class Vector;

class TextObject;

class FTextObject;

template<class T>
class CPoints;


template<class T>
class Matrix;


template<class T>
class FVector;

template<class T>
class FMatrix;


template <class T>
class CPoints;

class Orthonormal;

template<class T>
class Curve;

template<class T>
class Surf;

template<class T>
class Vol;

class Point1D;

class Point2D;

class Point3D;

class Point4D;


class Math;


template<class T>
class FnCurvObject;

template<class T>
class FnSurfObject;


// CURVES

class BezCurvBasisFunc;

class BezCurvBasisFuncSet;

class BspCurvBasisFunc;

class BspCurvBasisFuncSet;

class PolyCurvBasisFunc;

class PolyCurvBasisFuncSet;


template<class T>
class BezCurv;

template<class T>
class PolyCurv;


template<class T>
class BspCurv;

template<class T>
class CompBezCurv;

template<class T>
class CompPolyCurv;


class BspCurvDouble;

class FBspCurv;

class PBspCurv2D;

class PBspCurv3D;


// SURFACES

class BezSurfBasisFunc;

class BezSurfBasisFuncSet;

class BspSurfBasisFunc;

class BspSurfBasisFuncSet;

class PolySurfBasisFunc;

class PolySurfBasisFuncSet;


template<class T>
class PolySurf;

template<class T>
class BezSurf;

template<class T>
class BspSurf;

template<class T>
class CompBezSurf;

template<class T>
class CompPolySurf;


class KnotSet;

class BspSurfDouble;


class FBspSurf;

class PBspSurf3D;



#define polyTol 0.0000001
const double Pi = 4.0*atan(1.0);




const double pi = 4.0*atan(1.0);
const double cpntTol = 0.0001;
const double knotTol=0.000001;

class TextObject;
class FTextObject;

typedef std::string ObjectID;

typedef TextObject * PTextObject;
typedef TextObject & RTextObject;

typedef FTextObject * PFTextObject;
typedef FTextObject & RFTextObject;

class FTextObject {
	friend std::ifstream & operator>> (std::ifstream& ifs, RFTextObject r);
	friend std::ifstream & operator>> (std::ifstream& ifs, PFTextObject p);
	friend std::ofstream & operator<< (std::ofstream& ofs, FTextObject& r);
	friend std::ofstream & operator<< (std::ofstream& ofs, PFTextObject p);
public:
	virtual void readfile(std::ifstream& ifs) = 0;
	virtual void writefile(std::ofstream& ofs)  = 0;
};

class TextObject {
	friend std::ostream & operator<< (std::ostream& os, RTextObject r);
	friend std::ostream & operator<< (std::ostream& os, PTextObject p);
	friend std::istream & operator>> (std::istream& is, RTextObject r);
	//friend std::istream & operator>> (std::istream& is, PTextObject p);
public:
	virtual void write(std::ostream& os)= 0;
	virtual void read(std::istream &is) = 0;
};

template<class T>
class Curve : public TextObject, public FTextObject {
	virtual ObjectID Identity() const = 0;
public:
	virtual T operator()(double) const = 0;
	virtual T operator()(int, double) const = 0;
	virtual T Derive(int, double) const = 0;
	virtual double GetLeftLimit() const = 0;
	virtual double GetRightLimit() const = 0;
	virtual void readfile(std::ifstream& ifs) = 0;
	virtual void writefile(std::ofstream& ofs)  = 0;
	virtual void write(std::ostream& os)  = 0;
	virtual void read(std::istream &is) = 0;
};

template<class T>
class Surf : public TextObject, public FTextObject {
	virtual ObjectID Identity() const = 0;
public:
	virtual T operator()(double, double) const = 0;
	virtual T Derive(int, int, double, double) const =0;
	virtual T operator()(int, int, double, double) const = 0;
	virtual double GetLeftLimitU() const = 0;
	virtual double GetRightLimitU() const = 0;
	virtual double GetLeftLimitV() const = 0;
	virtual double GetRightLimitV() const = 0;
	virtual void readfile(std::ifstream& ifs) = 0;
	virtual void writefile(std::ofstream& ofs)   = 0;
	virtual void write(std::ostream& os)  = 0;
	virtual void read(std::istream &is) = 0;
};

template<class T>
class Vol : public TextObject, public FTextObject {
public:
	virtual ObjectID Identity() const = 0;
	virtual T operator()(double, double, double) const = 0;
	virtual T operator()(int, int, int, double, double, double) const = 0;
	virtual T Derive(int, int, int, double, double, double) const = 0;
	virtual double GetLeftLimitU() const = 0;
	virtual double GetRightLimitU() const = 0;
	virtual double GetLeftLimitV() const = 0;
	virtual double GetRightLimitV() const = 0;
	virtual double GetLeftLimitW() const = 0;
	virtual double GetRightLimitW() const = 0;
	virtual void readfile(std::ifstream& ifs) = 0;
	virtual void writefile(std::ofstream& ofs) = 0;
	virtual void write(std::ostream& os) = 0;
	virtual void read(std::istream &is) = 0;
};



#endif