
#include "pbezvol.h"

BezVolDouble::BezVolDouble() : BezVol<double>() {}
BezVolDouble::BezVolDouble(const BezVol<double>& v) : BezVol<double>(v) {}
BezVolDouble::BezVolDouble(const Matrix3D<double>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw) :
BezVol<double>(Cpts, Ktsu, Ktsv, Ktsw, Ordu, Ordv, Ordw) {}
BezVolDouble::BezVolDouble(const Matrix3D<double>& Cpts, int Ordu, int Ordv, int Ordw) :
BezVol<double>(Cpts,Ordu,Ordv,Ordw) {}
BezVolDouble::BezVolDouble(const Matrix3D<double>& Cpts, int Ordu, int Ordv, int Ordw, double Lu, double Ru, double Lv, double Rv, double Lw, double Rw) :
BezVol<double>(Cpts,Ordu,Ordv,Ordw,Lu,Ru,Lv,Rv,Lw,Rw) {}
BezVolDouble::BezVolDouble(int Ordu, int Ordv, int Ordw) : BezVol<double>(Ordu,Ordv,Ordw) {}




FBezVol::FBezVol() : BezVol<Point1D>() {}
FBezVol::FBezVol(const BezVol<Point1D>& v) : BezVol<Point1D>(v) {}
FBezVol::FBezVol(const Matrix3D<Point1D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw) :
BezVol<Point1D>(Cpts, Ktsu, Ktsv, Ktsw, Ordu, Ordv, Ordw) {}
FBezVol::FBezVol(const Matrix3D<Point1D>& Cpts, int Ordu, int Ordv, int Ordw) :
BezVol<Point1D>(Cpts,Ordu,Ordv,Ordw) {}
FBezVol::FBezVol(const Matrix3D<Point1D>& Cpts, int Ordu, int Ordv, int Ordw, double Lu, double Ru, double Lv, double Rv, double Lw, double Rw) :
BezVol<Point1D>(Cpts,Ordu,Ordv,Ordw,Lu,Ru,Lv,Rv,Lw,Rw) {}
FBezVol::FBezVol(int Ordu, int Ordv, int Ordw) : BezVol<Point1D>(Ordu,Ordv,Ordw) {}





PBezVol3D::PBezVol3D() : BezVol<Point3D>() {}
PBezVol3D::PBezVol3D(const BezVol<Point3D>& v) : BezVol<Point3D>(v) {}
PBezVol3D::PBezVol3D(const Matrix3D<Point3D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw) :
BezVol<Point3D>(Cpts, Ktsu, Ktsv, Ktsw, Ordu, Ordv, Ordw) {}
PBezVol3D::PBezVol3D(const Matrix3D<Point3D>& Cpts, int Ordu, int Ordv, int Ordw) :
BezVol<Point3D>(Cpts,Ordu,Ordv,Ordw) {}
PBezVol3D::PBezVol3D(const Matrix3D<Point3D>& Cpts, int Ordu, int Ordv, int Ordw, double Lu, double Ru, double Lv, double Rv, double Lw, double Rw) :
BezVol<Point3D>(Cpts,Ordu,Ordv,Ordw,Lu,Ru,Lv,Rv,Lw,Rw) {}
PBezVol3D::PBezVol3D(int Ordu, int Ordv, int Ordw) : BezVol<Point3D>(Ordu,Ordv,Ordw) {}





PBezVol3D::PBezVol3D(const FBezVol& v)
{
	Matrix3D<Point4D> m1(v.GetNumU(),v.GetNumV(),v.GetNumW());
	Vector<double> v1(v.GetKnotSetU().ComputeKnotSetAver());
	Vector<double> v2(v.GetKnotSetV().ComputeKnotSetAver());
	Vector<double> v3(v.GetKnotSetV().ComputeKnotSetAver());

	for (int k=0; k<v.GetNumW(); k++) 
		for (int i=0; i<v.GetNumU(); i++) 
			for (int j=0; j<v.GetNumV(); j++) 
				m1[k][i][j] = Point4D(v1[i],v2[j],v3[k],v.GetCPoints()[k][i][j].GetX());

	*this = PBezVol3D(m1,v.GetKnotsU(),v.GetKnotsV(),v.GetKnotsW(),v.GetOrdU(),v.GetOrdV(),v.GetOrdW());
}

