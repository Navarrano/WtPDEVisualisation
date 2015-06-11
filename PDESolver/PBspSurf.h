#ifndef PBSPSURF
#define PBSPSURF


#include "bspsurf.h"


class BspSurfDouble : public BspSurf<double>
{
	virtual ObjectID Identity() const { return std::string("class BspSurfDouble"); }
public:
	BspSurfDouble();
	BspSurfDouble(const BspSurf<double>& b);
	BspSurfDouble(const Matrix<double>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv, int Numu, int Numv);
	BspSurfDouble(const Matrix<double>& Cpts, int Ordu, int Ordv, int Numu, int Numv);
};


class FBspSurf : public BspSurf<Point1D>
{
	virtual ObjectID Identity() const { return std::string("class FBspSurf"); }
public:
	FBspSurf();
	FBspSurf(const BspSurf<Point1D>& b);
	FBspSurf(const BspSurf<double>& b);
	FBspSurf(const Matrix<Point1D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv, int Numu, int Numv);
	FBspSurf(const Matrix<Point1D>& Cpts, int Ordu, int Ordv, int Numu, int Numv);
	
	FBspSurf ComputePDESolution(double L, double W, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, Vector<Matrix<double> >& surface, Vector<int>& natgeom, Matrix<double>& bound) const;
};


class PBspSurf3D : public BspSurf<Point3D>
{
	virtual ObjectID Identity() const { return std::string("class PBspSurf3D"); }
public:
	PBspSurf3D();
	PBspSurf3D(const BspSurf<Point3D>& s);
	PBspSurf3D(const BspSurf<Point4D>& s);
	PBspSurf3D(const Matrix<Point3D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv, int Numu, int Numv);
	PBspSurf3D(const Matrix<Point3D>& Cpts, int Ordu, int Ordv, int Numu, int Numv);
	PBspSurf3D(const FBspSurf& b);
	PBspSurf3D(const FBspSurf& s1, const FBspSurf& s2, const FBspSurf& s3);
};



#endif
