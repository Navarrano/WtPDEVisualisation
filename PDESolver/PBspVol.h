#ifndef PBSPVOL
#define PBSPVOL

#include "BspVol.h"
#include "PBspSurf.h"
#include "PBspCurv.h"

class BspVolDouble : public BspVol<double>
{
	virtual ObjectID Identity() const { return std::string("class BspVolDouble"); }
public:
	BspVolDouble();

	BspVolDouble(const BspVol<double>& v);
	BspVolDouble(const Matrix3D<double>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw);
	BspVolDouble(const Matrix3D<double>& Cpts, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw);

};


class FBspVol : public BspVol<Point1D>
{
	virtual ObjectID Identity() const { return std::string("class FBspVol"); }
public:
	FBspVol();
	
	FBspVol(const BspVol<Point1D>& v);
	FBspVol(const Matrix3D<Point1D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw);
	FBspVol(const Matrix3D<Point1D>& Cpts, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw);

	double Poisson3DError(int numx, int numy, int numz);
	double Laplace3D1Error(int numx, int numy, int numz);
	double Laplace3D2Error(int numx, int numy, int numz);
	double Laplace3D1ErrorHeatFlowZ0(int numx, int numy, int numz);
	double Laplace3D2ErrorHeatFlowZ0(int numx, int numy, int numz);

	FBspVol ComputeSolution(double L, double W, double H, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, Vector<Matrix3D<double> >& volume, Vector<int>& natgeom,
										Matrix3D<double>& bound) const;
	FBspVol ComputeSolution1(double L, double W, double H, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, Vector<Matrix3D<double> >& volume, Vector<int>& natgeom,
										Matrix3D<double>& bound) const;
	FBspVol ComputeFiniteElementNew(double L, double W, double H, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, BspVol<double>& volume, Vector<int>& natgeom,
										Matrix3D<double>& bound) const;
};



class PBspVol3D : public BspVol<Point3D>
{
	virtual ObjectID Identity() const { return std::string("class PBspVol3D"); }
public:
	PBspVol3D();
	
	PBspVol3D(const BspVol<Point3D>& v);
	PBspVol3D(const Matrix3D<Point3D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw);
	PBspVol3D(const Matrix3D<Point3D>& Cpts, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw);
	PBspVol3D(const FBspVol& v);
};


#endif