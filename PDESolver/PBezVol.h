
#ifndef PBEZVOL
#define PBEZVOL


#include "bezvol.h"


class BezVolDouble : public BezVol<double>
{
	virtual ObjectID Identity() const { return std::string("class BezVol"); }
public:
	BezVolDouble();
	BezVolDouble(const BezVol<double>& v);
	BezVolDouble(const Matrix3D<double>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw);
	BezVolDouble(const Matrix3D<double>& Cpts, int Ordu, int Ordv, int Ordw);
	BezVolDouble(const Matrix3D<double>& Cpts, int Ordu, int Ordv, int Ordw, double Lu, double Ru, double Lv, double Rv, double Lw, double Rw);
	BezVolDouble(int Ordu, int Ordv, int Ordw);
};

class FBezVol : public BezVol<Point1D>
{
	virtual ObjectID Identity() const { return std::string("class BezVol"); }
public:
	FBezVol();
	FBezVol(const BezVol<Point1D>& v);
	FBezVol(const Matrix3D<Point1D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw);
	FBezVol(const Matrix3D<Point1D>& Cpts, int Ordu, int Ordv, int Ordw);
	FBezVol(const Matrix3D<Point1D>& Cpts, int Ordu, int Ordv, int Ordw, double Lu, double Ru, double Lv, double Rv, double Lw, double Rw);
	FBezVol(int Ordu, int Ordv, int Ordw);
};

class PBezVol3D : public BezVol<Point3D>
{
	virtual ObjectID Identity() const { return std::string("class PBezVol3D"); }
public:
	PBezVol3D();
	PBezVol3D(const BezVol<Point3D>& v);
	PBezVol3D(const Matrix3D<Point3D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw);
	PBezVol3D(const Matrix3D<Point3D>& Cpts, int Ordu, int Ordv, int Ordw);
	PBezVol3D(const Matrix3D<Point3D>& Cpts, int Ordu, int Ordv, int Ordw, double Lu, double Ru, double Lv, double Rv, double Lw, double Rw);
	PBezVol3D(int Ordu, int Ordv, int Ordw);
	PBezVol3D(const FBezVol& v);
};

#endif