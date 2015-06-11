#ifndef PBSPCURV
#define PBSPCURV

#include "bspcurv.h"


class BspCurvDouble : public BspCurv<double>
{
	virtual ObjectID Identity() const { return std::string("class BspCurvDouble"); }
public:
	BspCurvDouble();
	BspCurvDouble(const BspCurv<double>& v);
	BspCurvDouble(const Vector<double>& Cpts, const Vector<double>& Kts, int Ord, int Num);
	BspCurvDouble(const Vector<double>& Cpts, int Ord, int Num);
};


class FBspCurv : public BspCurv<Point1D>
{
	virtual ObjectID Identity() const { return std::string("class FBspCurv"); }
public:
	FBspCurv();
	FBspCurv(const BspCurv<Point1D>& v);
	FBspCurv(const BspCurv<double>& b);
	FBspCurv(const Vector<Point1D>& Cpts, const Vector<double>& Kts, int Ord, int Num);
	FBspCurv(const Vector<Point1D>& Cpts, int Ord, int Num);
	FBspCurv ComputeLeastSquares(const Vector<Point1D>& y, const Vector<double>& tau, const Vector<double>& W, int m, int k, int s, double rho, double& err);
	FBspCurv ComputeSmoothCurv(double u1, double u2, Vector<int>& v, int num) const;
	FBspCurv ComputeFiniteElement(double C, double L, double Q, double QL, double F, double FL) const;
	FBspCurv ComputeFiniteElementNew(double L, double C, Vector<int>& nums, Vector<double>& supports, Vector<double>& clamped, Vector<double>& ptforce, Vector<double>& disforce) const;
	FBspCurv ComputeFiniteElementNew(double L, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, Vector<Vector<double> >& curve, Vector<double>& bound) const;
	FBspCurv ComputeFiniteElementNew(double L, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, BspCurv<double>& curve, Vector<double>& bound) const;
	FBspCurv ComputeDiffEquation1(double L, double C, Vector<int>& nums, Vector<double>& supports, Vector<double>& clamped, Vector<double>& ptforce, Vector<double>& disforce) const;
	FBspCurv ComputeDiffEquation2(double L, double C, Vector<int>& nums, Vector<double>& supports, Vector<double>& clamped, Vector<double>& ptforce, Vector<double>& disforce) const;
	FBspCurv ComputeDiffEquation3(double L, double C, Vector<int>& nums, Vector<double>& supports, Vector<double>& clamped, Vector<double>& ptforce, Vector<double>& disforce) const;
	FBspCurv ComputeLeastSquaresPer(const Vector<Point1D>& y, const Vector<double>& tau, const Vector<double>& w, int m, int ord, int nu, int s, double& err);
	FBspCurv ComputeWeightedSpline1(const Vector<double>& tau, const Vector<double>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	FBspCurv ComputeWeightedSpline2(const Vector<double>& tau, const Vector<double>& y, int l, const Vector<double>& nu, const Vector<double>& w);		
	FBspCurv ComputeWeightedSpline3(const Vector<double>& tau, const Vector<double>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	FBspCurv ComputeWeightedSpline4(const Vector<double>& tau, const Vector<double>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	FBspCurv ComputeLeastSquaresSpline(std::ifstream& ifs, std::ofstream& ofs);
	FBspCurv ComputeInterpolatingSpline(std::ifstream& ifs, std::ofstream& ofs);
};

class PBspCurv2D : public BspCurv<Point2D>
{
	virtual ObjectID Identity() const { return std::string("class PBspCurv2D"); }
public:
	PBspCurv2D();
	PBspCurv2D(const BspCurv<Point2D>& v);
	PBspCurv2D(const Vector<Point2D>& Cpts, const Vector<double>& Kts, int Ord, int Num);
	PBspCurv2D(const Vector<Point2D>& Cpts, int Ord, int Num);
	PBspCurv2D(const FBspCurv& b);
	PBspCurv2D ComputeSmoothCurv(double u1, double u2, Vector<int>& v, int num) const;
	PBspCurv2D ComputeLeastSquares(const Vector<Point2D>& y, const Vector<double>& tau, const Vector<double>& w, int m, int ord, int s, double rho, double& err);
//	PBspCurv2D ComputeFiniteElementNew(double L, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, const PBspCurv2D& bcurv, Vector<double>& bound) const;
	PBspCurv2D ComputeLeastSquaresPer(const Vector<Point2D>& y, const Vector<double>& tau, const Vector<double>& w, int m, int ord, int nu, int s, double& err);
	PBspCurv2D ComputeWeightedSpline1(const Vector<double>& tau, const Vector<Point2D>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	PBspCurv2D ComputeWeightedSpline2(const Vector<double>& tau, const Vector<Point2D>& y, int l, const Vector<double>& nu, const Vector<double>& w);		
	PBspCurv2D ComputeWeightedSpline3(const Vector<double>& tau, const Vector<Point2D>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	PBspCurv2D ComputeWeightedSpline4(const Vector<double>& tau, const Vector<Point2D>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	PBspCurv2D ComputeLeastSquaresSpline(std::ifstream& ifs, std::ofstream& ofs);
	PBspCurv2D ComputeInterpolatingSpline(std::ifstream& ifs, std::ofstream& ofs);
};

class PBspCurv3D : public BspCurv<Point3D>
{
	virtual ObjectID Identity() const { return std::string("class PBspCurv3D"); }
public:
	PBspCurv3D();
	PBspCurv3D(const BspCurv<Point3D>& v);
	PBspCurv3D(const BspCurv<Point4D>& v);
	PBspCurv3D(const Vector<Point3D>& Cpts, const Vector<double>& Kts, int Ord, int Num);
	PBspCurv3D(const Vector<Point3D>& Cpts, int Ord, int Num);
	PBspCurv3D(const FBspCurv& b);
	PBspCurv3D(const PBspCurv2D& b);
	PBspCurv3D PerturbControlPoints(int u1, int u2) const;
	PBspCurv3D ComputeSmoothCurv(double u1, double u2, Vector<int>& v, int num) const;
	PBspCurv3D ComputeLeastSquares(const Vector<Point3D>& y, const Vector<double>& tau, const Vector<double>& w, int m, int ord, int s, double rho, double& err);
	PBspCurv3D ComputeLeastSquaresPer(const Vector<Point3D>& y, const Vector<double>& tau, const Vector<double>& w, int m, int ord, int nu, int s, double& err);
	PBspCurv3D ComputeWeightedSpline1(const Vector<double>& tau, const Vector<Point3D>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	PBspCurv3D ComputeWeightedSpline2(const Vector<double>& tau, const Vector<Point3D>& y, int l, const Vector<double>& nu, const Vector<double>& w);		
	PBspCurv3D ComputeWeightedSpline3(const Vector<double>& tau, const Vector<Point3D>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	PBspCurv3D ComputeWeightedSpline4(const Vector<double>& tau, const Vector<Point3D>& y, int l, const Vector<double>& nu, const Vector<double>& w);
	PBspCurv3D ComputeLeastSquaresSpline(std::ifstream& ifs, std::ofstream& ofs);
	PBspCurv3D ComputeInterpolatingSpline(std::ifstream& ifs, std::ofstream& ofs);
	Vector<Point2D> ComputeDifferenceArray(const PBspCurv3D& c, int m) const;
	Point3D ComputeUnitTangent(double u) const;
	PBspCurv3D ComputeTangent() const;
	Point3D ComputeTangentAt(double u) const;
	Point3D ComputeCurvatureVec(double u) const;
	double ComputeCurvature(double u) const;
	Vector<Point2D> ComputeCurvatureArray(int m) const;
	Point3D ComputePrincipalNormal(double u) const;
	double MaxDeviation(const PBspCurv3D& curv, int m) const;
	double MaxDeviation1(const PBspCurv3D& curv, int m) const;
	PBspCurv3D ComputeSmoothCurv1(double u1, double u2, double sm) const;
	PBspCurv3D ComputeSmoothCurv2(double u1, double u2, double sm) const;
	PBspCurv3D ComputeSmoothCurv3(double u1, double u2, double sm) const;
	PBspCurv3D ComputeSmoothCurv4(double u1, double u2, double sm) const;
	PBspCurv3D ComputeSmoothCurv5(double u1, double u2, double sm) const;
};


#endif