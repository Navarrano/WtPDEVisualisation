#ifndef BSPVOL
#define BSPVOL

#include "matrix3D.h"
#include "knotset.h"
#include "compbezvol.h"
#include "bspsurf.h"

template<class T>
class BspVol;

template<class T>
class CompBezVol;

template<class T>
class CompPolyVol;

template<class T>
class PolyVol;

template<class T>
class BezVol;

class BspVolBasisFunc : public Vol<double> {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int ordw;	// order in w
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<BspVol<double> > b;	// BspVol representation of basis function

	// private functions
	BspVol<double> CreateBspVol() const;	// creates the BspVol
	virtual ObjectID Identity() const { return std::string("class BspVolBasisFunc"); }
public:
	// constructors
	BspVolBasisFunc();
	BspVolBasisFunc(int Ordu, int Ordv, int Ordw, const Vector<double>& Ktsu, 
	const Vector<double>& Ktsv, const Vector<double>& Ktsw);

	// access functions
	int ComputeDimU() const;
	int ComputeDimV() const;
	int ComputeDimW() const;
	BspVol<double> GetBspVol() const;		// get BspVol representation

	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	Vector<double> GetKnotsW() const;

	// evaluators
	double Eval(double u, double v, double w) const;	// evaluate basis function
	virtual double operator()(double u, double v, double w) const;
	virtual double operator()(int, int, int, double u, double v, double w) const;
	
	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	
	virtual double Derive(int, int , int, double, double, double) const;
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
	virtual double GetLeftLimitW() const;
	virtual double GetRightLimitW() const;
};


class BspVolBasisFuncSet : public TextObject, public FTextObject {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int ordw;
	int numu;
	int numv;
	int numw;
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<Matrix3D<BspVolBasisFunc> > b;
	virtual ObjectID Identity() const { return std::string("class BspVolBasisFuncSet"); }
public:
	// constructors
	BspVolBasisFuncSet();
	BspVolBasisFuncSet(int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw);

	// evaluators
	Matrix3D<double> Eval(double u, double v, double w) const;	// evaluate basis function
	Matrix3D<double> operator()(double u, double v, double w) const;
	Matrix3D<double> operator()(int, int, int, double u, double v, double w) const;
	Matrix3D<double> Derive(int, int, int, double u, double v, double w) const;
	
	Matrix<double> ComputeUBasisMatrix() const;
	Matrix<double> ComputeVBasisMatrix() const;
	Matrix<double> ComputeWBasisMatrix() const;
	Matrix3D<double> CreateIntegral(double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateIntegralNewU(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateIntegralNewV(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateIntegralNewW(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<Matrix3D<double> > CreateMatrixIntegral(int levu, int levv, int levw) const;
	Matrix3D<Matrix3D<double> > CreateMatrixIntegral(int levu, int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisation1(int levu, int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisation1(int levu, int levv, int levw) const;
	Matrix3D<double> CreateMatrixMinimisation(int levu, int levv, int levw) const;
	Matrix3D<double> CreateMatrixMinimisation(int levu, int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationUV(int levu, int levv, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationUW(int levu, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationVW(int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationU(int levu, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationV(int levv, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationV1(int levv, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationW(int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix<double> CreateMatrixKnotAveragesU() const;
	Matrix<double> CreateMatrixKnotAveragesV() const;
	Matrix<double> CreateMatrixKnotAveragesW() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	Vector<double> GetKnotsW() const;

	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};


template<class T>
class BspVol : public Vol<T> {
private:
	// data
	int ordu;
	int ordv;
	int ordw;
	int numu;
	int numv;
	int numw;
	Ptr<Matrix3D<T> > cpts;
	Ptr<Vector<double> > ktsu;
	Ptr<Vector<double> > ktsv;
	Ptr<Vector<double> > ktsw;
	Ptr<KnotSet> ksetu;
	Ptr<KnotSet> ksetv;
	Ptr<KnotSet> ksetw;
	double lu1, lu2, lv1, lv2, lw1, lw2;

	// private functions

	Matrix3D<T> DeriveCPointsU(int level) const;
	Matrix3D<T> DeriveCPointsV(int level) const;
	Matrix3D<T> DeriveCPointsW(int level) const;
	BspVol<T> DeriveUV(int levu, int levv) const;
	BspVol<T> DeriveUW(int levu, int levw) const;
	BspVol<T> DeriveVW(int levv, int levw) const;
	T DeriveU(int level,double u, double v, double w) const;
	T DeriveV(int level,double u, double v, double w) const;
	T DeriveW(int level,double u, double v, double w) const;
	T DeriveUV(int levu,int levv, double u, double v) const;
	T DeriveUW(int levu,int levw, double u, double v, double w) const;
	T DeriveVW(int levv,int levw, double u, double v, double w) const;
	BspVol<T> RemovePossibleKnotsUV(double tol=knotTol) const;
	BspVol<T> RemovePossibleKnotsUW(double tol=knotTol) const;
	BspVol<T> RemovePossibleKnotsVW(double tol=knotTol) const;
	Matrix3D<T> ElevateCPointsU(int level) const;
	Matrix3D<T> ElevateCPointsV(int level) const;
	Matrix3D<T> ElevateCPointsW(int level) const;
	BspVol<T> ElevateUV(int levu, int levv) const;
	BspVol<T> ElevateUW(int levu, int levw) const;
	BspVol<T> ElevateVW(int levv, int levw) const;
	BspVol<T> IntegrateU() const;
	BspVol<T> IntegrateV() const;
	BspVol<T> IntegrateW() const;
	BspVol<T> IntegrateUV() const;
	BspVol<T> IntegrateUV1() const;
	BspVol<T> IntegrateUV2() const;
	BspVol<T> IntegrateUW() const;
	BspVol<T> IntegrateVW() const;
	BspVol<T> IntegrateUW1() const;
	BspVol<T> IntegrateVW1() const;
	T Integrate2(double u1, double u2, double v1, double v2, double w1, double w2) const;
	T Integrate1(double u1, double u2, double v1, double v2, double w1, double w2) const;
	T IntegrateU(double u1, double u2) const;
	T IntegrateV(double v1, double v2) const;
	T IntegrateW(double w1, double w2) const;
	T IntegrateUV(double u1, double u2, double v1, double v2) const;
	T IntegrateUW(double u1, double u2, double w1, double w2) const;
	T IntegrateVW(double v1, double v2, double w1, double w2) const;
	Matrix3D<T> IntegrateCPointsU() const;
	Matrix3D<T> IntegrateCPointsV() const;
	Matrix3D<T> IntegrateCPointsW() const;
	Matrix3D<T> IntegrateCPointsUV() const;
	Matrix3D<T> IntegrateCPointsUW() const;
	Matrix3D<T> IntegrateCPointsVW() const;
	BspVol<T> SubdivideU(double u1, double u2) const;
	BspVol<T> SubdivideV(double v1, double v2) const;
	BspVol<T> SubdivideW(double w1, double w2) const;
	Matrix3D<T> SubdivideCPointsU(double u1, double u2) const;
	Matrix3D<T> SubdivideCPointsV(double v1, double v2) const;
	Matrix3D<T> SubdivideCPointsW(double w1, double w2) const;
	Matrix3D<T> SubdivideCPointsU(int level) const;
	Matrix3D<T> SubdivideCPointsV(int level) const;
	Matrix3D<T> SubdivideCPointsW(int level) const;
	Matrix3D<T> SubdivideCPointsVW(int levv, int levw) const;
	Matrix3D<T> SubdivideCPointsUW(int levu, int levw) const;
	Matrix3D<T> SubdivideCPointsUV(int levu, int levv) const;
	BspVol<T> SubdivideUV(double u1, double u2, double v1, double v2) const;
	BspVol<T> SubdivideUW(double u1, double u2, double w1, double w2) const;
	BspVol<T> SubdivideVW(double u1, double u2, double v1, double v2) const;
	BspVol<T> SubdivideU(int level) const;
	BspVol<T> SubdivideV(int level) const;
	BspVol<T> SubdivideW(int level) const;
	BspVol<T> SubdivideUV(int levu, int levv) const;
	BspVol<T> SubdivideUW(int levu, int levw) const;
	BspVol<T> SubdivideVW(int levv, int levw) const;
	BspVol<T> InsertKnotU(double u) const;
	BspVol<T> InsertKnotU(double u, int level) const;
	BspVol<T> InsertKnotU(const Vector<double>& Kts, int n) const;
	BspVol<T> InsertKnotV(double v) const;
	BspVol<T> InsertKnotV(double v, int level) const;
	BspVol<T> InsertKnotV(const Vector<double>& Kts, int n) const;
	BspVol<T> InsertKnotW(double w) const;
	BspVol<T> InsertKnotW(double w, int level) const;
	BspVol<T> InsertKnotW(const Vector<double>& Kts, int n) const;
	BspVol<T> InsertKnotUV(double u, double v) const;
	BspVol<T> InsertKnotUV(double u, double v, int levu, int levv) const;
	BspVol<T> InsertKnotUV(const Vector<double>& Ktsu, int n, const Vector<int>& Ktsv, int m) const;
	BspVol<T> InsertKnotUV(const Vector<double>& Ktsu, const Vector<int>& multu, int N, const Vector<double>& Ktsv, const Vector<int>& multv, int m) const;
	BspVol<T> InsertKnotUW(double u, double w) const;
	BspVol<T> InsertKnotUW(double u, double w, int levu, int levw) const;
	BspVol<T> InsertKnotUW(const Vector<double>& Ktsu, int n, const Vector<int>& Ktsw, int m) const;
	BspVol<T> InsertKnotUW(const Vector<double>& Ktsu, const Vector<int>& multu, int N, const Vector<double>& Ktsw, const Vector<int>& multw, int m) const;
	BspVol<T> InsertKnotVW(double v, double w) const;
	BspVol<T> InsertKnotVW(double v, double w, int levv, int levw) const;
	BspVol<T> InsertKnotVW(const Vector<double>& Ktsv, int n, const Vector<int>& Ktsw, int m) const;
	BspVol<T> InsertKnotVW(const Vector<double>& Ktsv, const Vector<int>& multv, int N, const Vector<double>& Ktsw, const Vector<int>& multw, int m) const;
	BspVol<T> InsertKnotUVW(double u, double v, double w) const;
	BspVol<T> InsertKnotUVW(double u, double v, double w, int levu, int levv, int levw) const;
	BspVol<T> InsertKnotUVW(const Vector<double>& Ktsu, int n, const Vector<int>& Ktsv, int m, const Vector<double>& Ktsw, int p) const;
	BspVol<T> MakeCompatableU(const BspVol<T>& b) const;
	BspVol<T> MakeCompatableV(const BspVol<T>& b) const;
	BspVol<T> MakeCompatableW(const BspVol<T>& b) const;
	BspVol<T> MakeCompatableUV(const BspVol<T>& b) const;
	BspVol<T> MakeCompatableUW(const BspVol<T>& b) const;
	BspVol<T> MakeCompatableVW(const BspVol<T>& b) const;
	CompBezVol<T> ConvertCompBezVol2() const;
	BspVol<T> KnotRemovalU(double knot) const;
	BspVol<T> KnotRemovalV(double knot) const;
	BspVol<T> KnotRemovalW(double knot) const;
	BspVol<T> RemovePossibleKnotsU(double tol=knotTol) const;
	BspVol<T> RemovePossibleKnotsV(double tol=knotTol) const;
	BspVol<T> RemovePossibleKnotsW(double tol=knotTol) const;
	bool IsKnotRemovableU(double knot, double tol=cpntTol) const;
	bool IsKnotRemovableV(double knot, double tol=cpntTol) const;
	bool IsKnotRemovableW(double knot, double tol=cpntTol) const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class BspVol<double>");
		else {
			std::string s(typeid(T).name()), s1("class BspVol<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	BspVol();
	BspVol(const Matrix3D<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, 
		const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw,  int Numu, int Numv, int Numw);
	BspVol(const Matrix3D<T>& Cpts, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw);
	BspVol(const PolyVol<T>& p, const KnotSet& KsetU, const KnotSet& KsetV, const KnotSet& KsetW);
	BspVol(const PolyVol<T>& p, int Ordu, int Ordv, int Ordw, const KnotSet& KsetU, const KnotSet& KsetV, const KnotSet& KsetW);

	// access functions
	int GetOrdU() const;
	int GetNumU() const;
	int GetOrdV() const;
	int GetNumV() const;
	int GetOrdW() const;
	int GetNumW() const;
	Matrix3D<T> GetCPoints() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	Vector<double> GetKnotsW() const;
	KnotSet GetKnotSetU() const;
	KnotSet GetKnotSetV() const;
	KnotSet GetKnotSetW() const;
	BezVol<T> GetVol(int i, int j, int k) const; 
	virtual double GetLeftLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetLeftLimitW() const;

	virtual double GetRightLimitU() const;
	virtual double GetRightLimitV() const;
	virtual double GetRightLimitW() const;

	// evaluators
	virtual T operator()(double u, double v, double w) const;
	virtual T operator()(int, int, int, double u, double v, double w) const;
	T Eval(double u, double v, double w) const;
	T Eval1(double u, double v, double w) const;
	Matrix3D<T> Eval2(double u, double v, double w) const;
	Matrix3D<T> EvalDeriv(int levu, int levv, int levw, double u, double v, double w) const;
	Matrix3D<T> ComputeKnotSetAverValues() const;
	Matrix3D<T> ComputePoints(int n, int m, int p) const;
	
	// derivatives
	Matrix3D<T> DeriveCPoints(int levu, int levv, int levw) const;
	BspVol<T> Derive(int levu, int levv, int levw) const;
	virtual T Derive(int levu,int levv, int levw, double u, double v, double w) const;
	BspVol<T> DeriveU(int level) const; 
	BspVol<T> DeriveV(int level) const;
	BspVol<T> DeriveW(int level) const;
	
	// knot removal
	BspVol<T> RemovePossibleKnots(double tol=knotTol) const;
	
	// degree elevation
	Matrix3D<T> ElevateCPoints(int levu, int levv, int levw) const;
	BspVol<T> Elevate(int levu, int levv, int levw) const;
	BspVol<T> ElevateU(int level) const;
	BspVol<T> ElevateV(int level) const;
	BspVol<T> ElevateW(int level) const;

	// integration
	BspVol<T> Integrate() const;
	T Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<T> IntegrateCPoints() const;
	BspVol<T> IntegrateU(int level) const;
	BspVol<T> IntegrateV(int level) const;
	BspVol<T> IntegrateW(int level) const;


	// subdivision
	Matrix3D<T> SubdivideCPoints(int levu, int levv, int levw) const;
	BspVol<T> Subdivide(double u1, double u2, double v1, double v2, double w1, double w2) const;
	BspVol<T> Subdivide(int levu, int levv, int levw) const;

	// addition and subtraction
	BspVol<T> Add(const BspVol<T>& b) const;
	BspVol<T> Subtract(const BspVol<T>& b) const;

	// knot insertion
	BspVol<T> InsertKnotU(const Vector<double>& Kts, const Vector<int>& mult, int N) const;
	BspVol<T> InsertKnotV(const Vector<double>& Kts, const Vector<int>& mult, int N) const;
	BspVol<T> InsertKnotW(const Vector<double>& Kts, const Vector<int>& mult, int N) const;
	BspVol<T> InsertKnot(const Vector<double>& Ktsu, const Vector<int>& Multu, int N, const Vector<double>& Ktsv, const Vector<int>& Multv, int m, const Vector<double>& Ktsw, const Vector<int>& Multw, int p) const;

	BspVol<T> MakeKnotSetCompatableW(const KnotSet& KsetW) const;
	BspVol<T> MakeKnotSetCompatableV(const KnotSet& KsetV) const;
	BspVol<T> MakeKnotSetCompatableU(const KnotSet& KsetU) const;
	BspVol<T> MakeKnotSetCompatable(const KnotSet& KsetU, const KnotSet& KsetV, const KnotSet& KsetW) const;


	// make compatable with another volume


	template<class T1>
	BspVol<T> MakeCompatable(const BspVol<T1>& b) const
	{
		return MakeCompatableW(MakeCompatableUV(b));
	}
	
	template<class T1>
	BspVol<T> MakeBreakCompatable(const BspVol<T1>& b) const
	{
		return MakeBreakCompatableW(MakeBreakCompatableV(MakeBreakCompatableU(b)));
	}

	template<class T1>
	BspVol<T> MakeBreakCompatableU(const BspVol<T1>& b) const
	{
		// assume they are both of the same degree
		// create knot set for current object
	
		// create normalised knot set for current object relative to b
		KnotSet norm = (*ksetu).Normalise(b.GetKnotSetU());
	
		// add knots to current object in knot set for b not present 
		// in knot set for Current object
		KnotSet add1 = norm.KnotUnion(b.GetKnotSetU());
	
		// find differences in knots
		KnotSet kdiff1 = add1.KnotDifference(norm);
	
		if (kdiff1.GetNum() > 0) {
			// insert knots into current object to get compatibility with b
			return 	BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
		} else return BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw);
	}

	template<class T1>
	BspVol<T> MakeBreakCompatableV(const BspVol<T1>& b) const
	{
		// assume they are both of the same degree
		// create knot set for current object
	
		// create normalised knot set for current object relative to b
		KnotSet norm = (*ksetv).Normalise(b.GetKnotSetV());
	
		// add knots to current object in knot set for b not present 
		// in knot set for Current object
		KnotSet add1 = norm.KnotUnion(b.GetKnotSetV());
	
		// find differences in knots
		KnotSet kdiff1 = add1.KnotDifference(norm);
	
		if (kdiff1.GetNum() > 0) {
			// insert knots into current object to get compatibility with b
			return 	BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
		} else return BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw);
	}
	
	template<class T1>
	BspVol<T> MakeBreakCompatableW(const BspVol<T1>& b) const
	{
		// assume they are both of the same degree
		// create knot set for current object
	
		// create normalised knot set for current object relative to b
		KnotSet norm = (*ksetw).Normalise(b.GetKnotSetW());
	
		// add knots to current object in knot set for b not present 
		// in knot set for Current object
		KnotSet add1 = norm.KnotUnion(b.GetKnotSetW());
	
		// find differences in knots
		KnotSet kdiff1 = add1.KnotDifference(norm);
	
		if (kdiff1.GetNum() > 0) {
			// insert knots into current object to get compatibility with b
			return 	BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw).InsertKnotW(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
		} else return BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw);
	}


	BspVol<T> MakeBreakCompatable(const KnotSet& KsetU, const KnotSet& KsetV, const KnotSet& KsetW) const;
	BspVol<T> MakeBreakCompatableU(const KnotSet& KsetU) const;
	BspVol<T> MakeBreakCompatableV(const KnotSet& KsetV) const;
	BspVol<T> MakeBreakCompatableW(const KnotSet& KsetW) const;

	// get isoparametric curve and surfaces
	BspSurf<T> GetIsoparametricU(double u) const;
	BspSurf<T> GetIsoparametricV(double v) const;
	BspSurf<T> GetIsoparametricW(double w) const;
	BspCurv<T> GetIsoparametricUV(double u, double v) const;
	BspCurv<T> GetIsoparametricUW(double u, double w) const;
	BspCurv<T> GetIsoparametricVW(double v, double w) const;


	// conversions	
	CompBezVol<T> ConvertCompBezVol() const;
	CompPolyVol<T> ConvertCompPolyVol() const;

	// products
	template<class T1>
	BspVol<T> Product(const BspVol<T1>& b) const
	{
		// make local copies of the two surfaces
		BspVol<T> d(*this);
		BspVol<T1> e(b);

		// make the surfaces compatable in u and v

		if (!(*ksetu).IsSameDistinctKnotSet(b.GetKnotSetU())) {
			// needs to add distinct knots from b not in current object
			// into current object (independent of order)
			d = MakeBreakCompatable(b);
			e = b.MakeBreakCompatable(d);
		}
	
		// convert to CompBezVol
		CompBezVol<T> e1 = d.ConvertCompBezVol();
		CompBezVol<T1> e2 = e.ConvertCompBezVol();

		// multiply
		CompBezVol<T> prod = e1.Product(e2);
  
		// remove all possible knots in u and v
		return prod.ConvertBspVol();
		//return (prod.ConvertBspVol1()).RemovePossibleKnotsUVW();
	}

	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};

// CONSTRUCTORS

// constructor builds a BspVol from a Matrix of control points, knots in u and v
// orders in u and v and number of control points in u and v
template<class T>
BspVol<T>::BspVol(const Matrix3D<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, 
	int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw) :
			ordu(Ordu), ordv(Ordv), ordw(Ordw), numu(Numu), numv(Numv), numw(Numw), 
				cpts(new Matrix3D<T>(Cpts)), ktsu(new Vector<double>(Ktsu)), 
				ktsv(new Vector<double>(Ktsv)), ktsw(new Vector<double>(Ktsw)),
				ksetu(new KnotSet(*ktsu,ordu,ordu+numu)),
				ksetv(new KnotSet(*ktsv,ordv,ordv+numv)),
				ksetw(new KnotSet(*ktsw,ordw,ordw+numw)), lu1((*ktsu)[ordu-1]), lu2((*ktsu)[numu]),lv1((*ktsv)[ordv-1]), lv2((*ktsv)[numv]),lw1((*ktsw)[ordw-1]), lw2((*ktsw)[numw])
{
}


// default constructor
template<class T>
BspVol<T>::BspVol() : ordu(0), ordv(0), ordw(0), numu(0), numv(0), numw(0), 
cpts(), ktsu(), ktsv(), ktsw() { }


// constructor builds a BspVol from a Matrix of control points, 
// orders in u and v and number of control points in u and v
// creates the knots using a simple algorithm
template<class T>
BspVol<T>::BspVol(const Matrix3D<T>& Cpts, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw) :
			ordu(Ordu), ordv(Ordv), ordw(Ordw), numu(Numu), numv(Numv), numw(Numw), 
			cpts(new Matrix3D<T>(Cpts))
{
	ktsu = new Vector<double>(Math::CreateKnots(Ordu,Numu));
	ktsv = new Vector<double>(Math::CreateKnots(Ordv,Numv));
	ktsw = new Vector<double>(Math::CreateKnots(Ordw,Numw));
	ksetu = new KnotSet(*ktsu,ordu,ordu+numu);
	ksetv = new KnotSet(*ktsv,ordv,ordv+numv);
	ksetw = new KnotSet(*ktsw,ordw,ordw+numw);
	lu1= (*ktsu)[ordu-1];
	lu2= (*ktsu)[numu];
	lv1= (*ktsv)[ordv-1];
	lv2 =(*ktsv)[numv];
	lw1 = (*ktsw)[ordw-1]; 
	lw2= (*ktsw)[numw];
}

template<class T>
BspVol<T>::BspVol(const PolyVol<T>& p, const KnotSet& KsetU, const KnotSet& KsetV, const KnotSet& KsetW) : ordu(KsetU.GetOrd()), ordv(KsetV.GetOrd()), ordw(KsetW.GetOrd()),
				numu(KsetU.GetNum()-KsetU.GetOrd()), numv(KsetV.GetNum()-KsetV.GetOrd()), numw(KsetW.GetNum()-KsetW.GetOrd()),
					cpts(new Matrix3D<T>(p.ElevateUVW(KsetU.GetOrd()-p.GetOrdU(),KsetV.GetOrd()-p.GetOrdV(),KsetW.GetOrd()-p.GetOrdW()).ConvertBspVol().MakeKnotSetCompatable(KsetU,KsetV,KsetW).GetCPoints())), ktsu(new Vector<double>(KsetU.GetKnots())), ktsv(new Vector<double>(KsetV.GetKnots())), ktsw(new Vector<double>(KsetW.GetKnots())), 
					ksetu(new KnotSet(KsetU)), ksetv(new KnotSet(KsetV)), ksetw(new KnotSet(KsetW)),lu1((*ktsu)[ordu-1]), lu2((*ktsu)[numu]),lv1((*ktsv)[ordv-1]), lv2((*ktsv)[numv]),lw1((*ktsw)[ordw-1]), lw2((*ktsw)[numw])
{

}


template<class T>
BspVol<T>::BspVol(const PolyVol<T>& p, int Ordu, int Ordv, int Ordw, const KnotSet& KsetU, const KnotSet& KsetV, const KnotSet& KsetW) : ordu(Ordu), ordv(Ordv), ordw(Ordw)
{
	*this = p.ConvertBspVol().MakeBreakCompatable(KsetU,KsetV,KsetW);
}



// ACCESS FUNCTIONS

// get the order of the BspVol in u
template<class T>
inline int BspVol<T>::GetOrdU() const { return ordu; }

// get the order of the BspVol in v
template<class T>
inline int BspVol<T>::GetOrdV() const { return ordv; }


// get the order of the BspVol in w
template<class T>
inline int BspVol<T>::GetOrdW() const { return ordw; }

// get the number of control points in u
template<class T>
inline int BspVol<T>::GetNumU() const { return numu; }


// get the number of control points in v
template<class T>
inline int BspVol<T>::GetNumV() const { return numv; }

// get the number of control points in w
template<class T>
inline int BspVol<T>::GetNumW() const { return numw; }


// get the control points
template<class T>
inline Matrix3D<T> BspVol<T>::GetCPoints() const { return *cpts; }

// get the knot vector in u
template<class T>
inline Vector<double> BspVol<T>::GetKnotsU() const { return *ktsu; }


// get the knot vector in v
template<class T>
inline Vector<double> BspVol<T>::GetKnotsV() const { return *ktsv; }

// get the knot vector in w
template<class T>
inline Vector<double> BspVol<T>::GetKnotsW() const { return *ktsw; }

// get the knot vector in w
template<class T>
inline KnotSet BspVol<T>::GetKnotSetU() const { return *ksetu; }

// get the knot vector in w
template<class T>
inline KnotSet BspVol<T>::GetKnotSetV() const { return *ksetv; }

// get the knot vector in w
template<class T>
inline KnotSet BspVol<T>::GetKnotSetW() const { return *ksetw; }


// get the knot vector
template<class T>
double BspVol<T>::GetLeftLimitU() const { return (*ktsu)[ordu-1]; }


// get the knot vector
template<class T>
double BspVol<T>::GetLeftLimitV() const { return (*ktsv)[ordv-1]; }

// get the knot vector
template<class T>
double BspVol<T>::GetLeftLimitW() const { return (*ktsw)[ordw-1]; }


// get the knot vector
template<class T>
double BspVol<T>::GetRightLimitU() const { return (*ktsu)[numu]; }


// get the knot vector
template<class T>
double BspVol<T>::GetRightLimitV() const { return (*ktsv)[numv]; }

// get the knot vector
template<class T>
double BspVol<T>::GetRightLimitW() const { return (*ktsw)[numw]; }

// Get the i,j pacth of the Volace and return as a BspVol
// should be as a BezVol?
template<class T>       
BezVol<T> BspVol<T>::GetVol(int i, int j, int k) const
{
	// subdivide at appropriate knots
	BspVol<T> b = Subdivide((*ksetu).GetDistinctKnots()[i-1],(*ksetu).GetDistinctKnots()[i],
				(*ksetv).GetDistinctKnots()[j-1],(*ksetv).GetDistinctKnots()[j],
				(*ksetw).GetDistinctKnots()[k-1],(*ksetw).GetDistinctKnots()[k]);

	return BezVol<T>(b.GetCPoints(),b.GetKnotsU(),b.GetKnotsV(),b.GetKnotsW(),ordu,ordv,ordw);
}



// EVALUATORS

// evaluate the BspVol at the point x using de Boor algorithm
template<class T>
T BspVol<T>::operator()(double u, double v, double w) const
{
//	if (u < (*ktsu)[ordu-1] || u > (*ktsu)[numu]) return 0.0;
//	if (v < (*ktsv)[ordv-1] || v > (*ktsv)[numv]) return 0.0;
//	if (w < (*ktsw)[ordw-1] || w > (*ktsw)[numw]) return 0.0;

	Vector<T> v2(numu);

	// use tensor product evaluation
	// evaluate in v and create points in u
	for (int i=0; i<numu; i++) v2[i]= BspSurf<T>((*cpts).GetVW(i),*ktsv,*ktsw,ordv,ordw,numv,numw)(v,w);

	// evaluate the resulting curve
	return BspCurv<T>(v2,*ktsu,ordu,numu)(u);
}


// evaluate the BspVol using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T BspVol<T>::Eval(double u, double v, double w) const
{
	return (*this)(u,v,w);//ConvertCompBezVol().Eval(u,v,w);
}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
T BspVol<T>::Eval1(double u, double v, double w) const
{
	Vector<double> uvec = ksetu->CreateVectorInterp(u);
	Vector<double> vvec = ksetv->CreateVectorInterp(v);
	Vector<double> wvec = ksetw->CreateVectorInterp(w);

	return Math::mult0(Math::mult3(Math::mult7(uvec,*cpts),vvec),wvec);
}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Matrix3D<T> BspVol<T>::Eval2(double u, double v, double w) const
{
	Matrix<double> mu(ksetu->CreateVectorInterp(u));
	Matrix<double> mv(ksetv->CreateVectorInterp(v));
	Matrix<double> mw(ksetw->CreateVectorInterp(w));

	return Math::mult9(mw,Math::mult5(Math::mult6(mu,*cpts),Math::transpose(mv)));
}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Matrix3D<T> BspVol<T>::EvalDeriv(int levu, int levv, int levw, double u, double v, double w) const
{
	Matrix<double> mu(ksetu->CreateVectorInterpDeriv(levu,u));
	Matrix<double> mv(ksetv->CreateVectorInterpDeriv(levv,v));
	Matrix<double> mw(ksetw->CreateVectorInterpDeriv(levw,w));

	Matrix3D<T> m1 = Math::mult6(Math::mult1(mu,ksetu->CreateMatrixDeriv(levu)),*cpts);
	Matrix3D<T> m2 = Math::mult5(m1,Math::transpose(Math::mult1(mv,ksetv->CreateMatrixDeriv(levv))));
	return Math::mult9(Math::mult1(mw,ksetw->CreateMatrixDeriv(levw)),m2);
}



template<class T>
Matrix3D<T> BspVol<T>::ComputeKnotSetAverValues() const
{
	Matrix3D<T> m(numu,numv,numw);

	Vector<double> uval((*ksetu).ComputeKnotSetAver());
	Vector<double> vval((*ksetv).ComputeKnotSetAver());
	Vector<double> wval((*ksetw).ComputeKnotSetAver());

	for (int i=0; i<numu; i++)
		for (int j=0; j<numv; j++) 
			for (int k=0; k<numw; k++) m[k][i][j] = (*this)(uval[i],vval[j],wval[k]);

	return m;
}


// KNOT REMOVAL


// remove all possible knots from the Volace
// remove possible u knots then v knots
template<class T>
BspVol<T> BspVol<T>::RemovePossibleKnots(double tol) const
{
	return (RemovePossibleKnotsU(tol).RemovePossibleKnotsV(tol)).RemovePossibleKnotsW(tol);
}

/*
template<class T>
BspVol<T> BspVol<T>::MakeBreakCompatableU(const BspVol<T>& b) const
{	
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetu).Normalise(b.GetKnotSetU());
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(b.GetKnotSetU());
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw);
}


template<class T>
BspVol<T> BspVol<T>::MakeBreakCompatableV(const BspVol<T>& b) const
{	
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetv).Normalise(b.GetKnotSetV());
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(b.GetKnotSetV());
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw);
}


template<class T>
BspVol<T> BspVol<T>::MakeBreakCompatableW(const BspVol<T>& b) const
{	
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetw).Normalise(b.GetKnotSetW());
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(b.GetKnotSetW());
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw).InsertKnotW(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw);
}


template<class T> 
BspVol<T> BspVol<T>::MakeBreakCompatable(const BspVol<T>& b) const
{
	return MakeBreakCompatableW(MakeBreakCompatableV(MakeBreakCompatableU(b)));
}



// MAKE COMPATABLE

template<class T> 
BspVol<T> BspVol<T>::MakeCompatable(const BspVol<T>& b) const
{
	return MakeCompatableW(MakeCompatableUV(b));
}
*/

// ADD and SUBTRACT

// add two Volaces together, assuming they are not compatable
// different knot sets and different orders
template<class T>
BspVol<T> BspVol<T>::Add(const BspVol<T>& b) const
{
	// make local copies of the two Volaces
	BspVol<T> a(*this), c(b);

	// make them the same order
	if (ordu > b.GetOrdU()) c = c.ElevateU(ordu-b.GetOrdU());
	if (b.GetOrdU() > ordu) a = a.ElevateU(b.GetOrdU()-ordu);
	if (ordv > b.GetOrdV()) c = c.ElevateV(ordv-b.GetOrdV());
	if (b.GetOrdV() > ordv) a = a.ElevateV(b.GetOrdV()-ordv);
	if (ordw > b.GetOrdW()) c = c.ElevateW(ordw-b.GetOrdW());
	if (b.GetOrdW() > ordw) a = a.ElevateW(b.GetOrdW()-ordw);


	BspVol<T> d1 = a.MakeCompatable(c);
	BspVol<T> d2 = c.MakeCompatable(d1);


	// now add them together
	// add the control points
	Matrix3D<T> temp(d1.GetNumU(),d1.GetNumV(),d1.GetNumW());

	for (int i=0; i<d1.GetNumU(); i++)
		for (int j=0; j<d1.GetNumV(); j++)
			for (int k=0; k<d1.GetNumW(); k++)
				temp[k][i][j]=d1.GetCPoints()[k][i][j]+d2.GetCPoints()[k][i][j];

	return BspVol<T>(temp,d1.GetKnotsU(),d1.GetKnotsV(),d1.GetKnotsW(),d1.GetOrdU(),d1.GetOrdV(),d1.GetOrdW(),d1.GetNumU(),d1.GetNumV(),d1.GetNumW());
}

// subtract two Volaces together, assuming they are not compatable
// different orders and different knots
template<class T>
BspVol<T> BspVol<T>::Subtract(const BspVol<T>& b) const
{
	// make local copies of the two Volaces
	BspVol<T> a(*this), c(b);

	// make them the same order
	if (ordu > b.GetOrdU()) c = c.ElevateU(ordu-b.GetOrdU());
	if (b.GetOrdU() > ordu) a = a.ElevateU(b.GetOrdU()-ordu);
	if (ordv > b.GetOrdV()) c = c.ElevateV(ordv-b.GetOrdV());
	if (b.GetOrdV() > ordv) a = a.ElevateV(b.GetOrdV()-ordv);
	if (ordw > b.GetOrdW()) c = c.ElevateW(ordw-b.GetOrdW());
	if (b.GetOrdW() > ordw) a = a.ElevateW(b.GetOrdW()-ordw);


	BspVol<T> d1 = a.MakeCompatable(c);
	BspVol<T> d2 = c.MakeCompatable(d1);


	// now add them together
	// add the control points
	Matrix3D<T> temp(d1.GetNumU(),d1.GetNumV(),d1.GetNumW());

	for (int i=0; i<d1.GetNumU(); i++)
		for (int j=0; j<d1.GetNumV(); j++)
			for (int k=0; k<d1.GetNumW(); k++)
				temp[k][i][j]=d1.GetCPoints()[k][i][j]-d2.GetCPoints()[k][i][j];
	
	return BspVol<T>(temp,d1.GetKnotsU(),d1.GetKnotsV(),d1.GetKnotsW(),d1.GetOrdU(),d1.GetOrdV(),d1.GetOrdW(),d1.GetNumU(),d1.GetNumV(),d1.GetNumW());
}



template<class T>	
BspVol<T> BspVol<T>::MakeBreakCompatable(const KnotSet& KsetU, const KnotSet& KsetV, const KnotSet& KsetW) const
{
	return MakeBreakCompatableW(KsetW).MakeBreakCompatableV(KsetV).MakeBreakCompatableU(KsetU);
}

template<class T>		
BspVol<T> BspVol<T>::MakeBreakCompatableU(const KnotSet& KsetU) const
{
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetu).Normalise(KsetU);
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(KsetU);
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw);
}

	
template<class T>	
BspVol<T> BspVol<T>::MakeBreakCompatableV(const KnotSet& KsetV) const
{
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetv).Normalise(KsetV);
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(KsetV);
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw);
}


template<class T>	
BspVol<T> BspVol<T>::MakeBreakCompatableW(const KnotSet& KsetW) const
{
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetw).Normalise(KsetW);
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(KsetW);
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw).InsertKnotW(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw);
}

// DEGREE ELEVATION

// elevate the degree of the BspVol by levu in u and levv in v
template<class T>
Matrix3D<T> BspVol<T>::ElevateCPoints(int levu, int levv, int levw) const
{
	Matrix3D<T> mat = ElevateCPointsUV(levu,levv);
	// create the knot set
	KnotSet kset1 = (*ksetu).CreateKnotSetElevate(levu);
	// create the knot set
	KnotSet kset2 = (*ksetv).CreateKnotSetElevate(levv);
	
	// elevate in u then v the w
	return BspVol<T>(mat,kset1.GetKnots(),kset2.GetKnots(),*ktsw,ordu+levu,ordv+levv,ordw,kset1.GetNum()-(ordu+levu),kset2.GetNum()-(ordv+levv),numw).ElevateCPointsW(levw);
}



// elevate the degree of the BspVol by levu in u and levv in v
template<class T>
BspVol<T> BspVol<T>::Elevate(int levu, int levv, int levw) const
{
	// elevate in u then v the w
	return (ElevateU(levu).ElevateV(levv)).ElevateW(levw);
}


// DERIVATIVES

template<class T>
T BspVol<T>::operator() (int valu, int valv, int valw, double u, double v, double w) const
{
	return Derive(valu,valv,valw,u,v,w);
}



// compute the derivative of the BspVol of order deriv and
// represent the result as another BspVol
template<class T>
BspVol<T> BspVol<T>::Derive(int levu, int levv, int levw) const
{   
	// v derivative then w derivative
	return (DeriveU(levu).DeriveV(levv)).DeriveW(levw);
}


// evaluate the derivative of the BspVol of order deriv at a point
// x. Computes the derivative as a BspVol and then evaluates this at x
template<class T>
T BspVol<T>::Derive(int levu, int levv, int levw, double u, double v, double w) const
{
	return Derive(levu,levv,levw)(u,v,w);
}


// INTEGRATION

// integrate the BVol
// computes the indefinite integral as a BspVol and then evaluates 
// the integral by adding and subtracting the areas
template<class T>
T BspVol<T>::Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	BspVol<T> temp = Integrate();

/*	if (u2 < temp.GetKnotsU()[temp.GetOrdU()-1] || u1 > temp.GetKnotsU()[temp.GetNumU()]) return 0;
	if (v2 < temp.GetKnotsV()[temp.GetOrdV()-1] || v1 > temp.GetKnotsV()[temp.GetNumV()]) return 0;
	if (w2 < temp.GetKnotsW()[temp.GetOrdW()-1] || w1 > temp.GetKnotsW()[temp.GetNumW()]) return 0;


	if (u1 < temp.GetKnotsU()[temp.GetOrdU()-1]) u1 = temp.GetKnotsU()[temp.GetOrdU()-1];
	if (v1 < temp.GetKnotsV()[temp.GetOrdV()-1]) v1 = temp.GetKnotsV()[temp.GetOrdV()-1];
	if (w1 < temp.GetKnotsW()[temp.GetOrdW()-1]) w1 = temp.GetKnotsW()[temp.GetOrdW()-1];
	
	
	if (u2 > temp.GetKnotsU()[temp.GetNumU()]) u2 = temp.GetKnotsU()[temp.GetNumU()];
	if (v2 > temp.GetKnotsV()[temp.GetNumV()]) v2 = temp.GetKnotsV()[temp.GetNumV()];
	if (w2 > temp.GetKnotsW()[temp.GetNumW()]) w2 = temp.GetKnotsW()[temp.GetNumW()];*/

	if (u2 < lu1 || u1 > lu2) return 0;
	if (v2 < lv1 || v1 > lv2) return 0;
	if (w2 < lw1 || w1 > lw2) return 0;

	if (u1 < lu1) u1 = lu1;
	if (v1 < lv1) v1 = lv1;
	if (w1 < lw1) w1 = lw1;
	if (u2 > lu2) u2 = lu2;
	if (v2 > lv2) v2 = lv2;
	if (w2 > lw2) w2 = lw2;
	
	return temp(u2,v2,w2)-temp(u2,v2,w1)-temp(u2,v1,w2)-temp(u1,v2,w2)+
		temp(u2,v1,w1)+temp(u1,v2,w1)+temp(u1,v1,w2)-temp(u1,v1,w1);
}


// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::Integrate() const
{
	// create knot sets
	KnotSet kset1 = (*ksetu).CreateKnotSetIntegrate();
	KnotSet kset2 = (*ksetv).CreateKnotSetIntegrate();
	KnotSet kset3 = (*ksetw).CreateKnotSetIntegrate();
	// create and return the new BspVol
	return BspVol<T>(IntegrateCPoints(),kset1.GetKnots(),kset2.GetKnots(),kset3.GetKnots(),ordu+1,ordv+1,ordw+1,numu+1,numv+1,numw+1);
}       


// compute the indefinite integral of the BspVol as a BspVol
// and return just the control points
template<class T>
Matrix3D<T> BspVol<T>::IntegrateCPoints() const
{
	Matrix3D<T> mat = IntegrateCPointsVW();
	KnotSet kset1 = (*ksetv).CreateKnotSetIntegrate();
	KnotSet kset2 = (*ksetw).CreateKnotSetIntegrate();

	return BspVol<T>(mat,*ktsu,kset1.GetKnots(),kset2.GetKnots(),ordu,ordv+1,ordw+1,numu,numv+1,numw+1).IntegrateCPointsU();
}

/*
// PRODUCT

// product using conversion to piecewise Bezier and multiplication
template<class T>
BspVol<T> BspVol<T>::Product(const BspVol<T>& b) const
{
		// make local copies of the two surfaces
	BspVol<T> d(*this), e(b);

	// make the surfaces compatable in u and v

	if (!(*ksetu).IsSameDistinctKnotSet(b.GetKnotSetU())) {
		// needs to add distinct knots from b not in current object
		// into current object (independent of order)
		d = MakeBreakCompatable(b);
		e = b.MakeBreakCompatable(d);
	}
	

	// convert to CompBezVol
	CompBezVol<T> e1 = d.ConvertCompBezVol();
	CompBezVol<T> e2 = e.ConvertCompBezVol();

	// multiply
	CompBezVol<T> prod = e1.Product(e2);
   
	// remove all possible knots in u and v
	return prod.ConvertBspVol();
	//return (prod.ConvertBspVol1()).RemovePossibleKnotsUVW();
}
*/

// SUBDIVISION


// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPoints(int levu, int levv, int levw) const
{
	Matrix3D<T> mat = SubdivideCPointsVW(levv,levw);
	KnotSet kset1 = (*ksetv).CreateKnotSetSubdivide(levv);
	KnotSet kset2 = (*ksetw).CreateKnotSetSubdivide(levw);

	// subdivide 
	return BspVol<T>(mat,*ktsu,kset1.GetKnots(),kset2.GetKnots(),ordu,ordv,ordw,numu,kset1.GetNum()-ordv,kset2.GetNum()-ordw).SubdivideCPointsU(levu);
}


// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
BspVol<T> BspVol<T>::Subdivide(int levu, int levv, int levw) const
{
	// subdivide u then v
	return (SubdivideU(levu).SubdivideV(levv)).SubdivideW(levw);
}

// subdivide the BspVol between the limits x1 and x2
// and return the new BspVol
template<class T>
BspVol<T> BspVol<T>::Subdivide(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	// subdivideu then v then w
	return (SubdivideU(u1,u2).SubdivideV(v1,v2)).SubdivideW(w1,w2);
}



// KNOT INSERTION
 
// inserts a vector of knots into the BspVol according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspVol<T> BspVol<T>::InsertKnotU(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	if (n <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetu).CreateKnotSetInsert(t,mult,n);
	Matrix3D<T> ncpts(kset.GetNum()-ordu,numv, numw);
	// insert into v curves
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).InsertKnotCPoints(t,mult,n);
			for (int i=0; i<kset.GetNum()-ordu; i++) ncpts[k][i][j]=temp[i];
		}
	return BspVol<T>(ncpts,kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu, numv, numw);
}
	     
 
// inserts a vector of knots into the BspVol according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspVol<T> BspVol<T>::InsertKnotV(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	if (n <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(t,mult,n);
	Matrix3D<T> ncpts(numu,kset.GetNum()-ordv,numw);

	// insert into v curves
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).InsertKnotCPoints(t,mult,n);
			for (int j=0; j<kset.GetNum()-ordv; j++) ncpts[k][i][j]=temp[j];
		}
	return BspVol<T>(ncpts,*ktsu, kset.GetKnots(), *ktsw, ordu, ordv, ordw, numu,kset.GetNum()-ordv,numw);  
}

// inserts a vector of knots into the BspVol according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspVol<T> BspVol<T>::InsertKnotW(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// create knot set
	if (n <= 0) return *this;
	KnotSet kset = (*ksetw).CreateKnotSetInsert(t,mult,n);
	Matrix3D<T> ncpts(numu,numv,kset.GetNum()-ordw);

	// insert into v curves
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).InsertKnotCPoints(t,mult,n);
			for (int k=0; k<kset.GetNum()-ordw; k++) ncpts[k][i][j]=temp[k];
		}
	return BspVol<T>(ncpts,*ktsu, *ktsv, kset.GetKnots(), ordu, ordv,ordw,numu,numv,kset.GetNum()-ordw);    
}

template<class T>
BspVol<T> BspVol<T>::InsertKnot(const Vector<double>& Ktsu, const Vector<int>& multu, int N, const Vector<double>& Ktsv, const Vector<int>& multv, int M, const Vector<double>& Ktsw, const Vector<int>& multw, int P) const
{
	// insert into u then w
	return (InsertKnotU(Ktsu,multu,n).InsertKnotV(Ktsv,multv,m)).InsertKnot(Ktsw,multw,p);
}



// CONVERSION

// convert the BspVol to a composite Bezier Vol
// by converting the u and v curves
template<class T>
CompBezVol<T> BspVol<T>::ConvertCompBezVol() const
{
	// create knotset object

	KnotSet kset1 = (*ksetu).CreateKnotSetCompBezCurv();
	KnotSet kset2 = (*ksetv).CreateKnotSetCompBezCurv();
	KnotSet kset3 = (*ksetw).CreateKnotSetCompBezCurv();

	
	Matrix3D<T> ncpts(kset1.GetNum()-ordu,kset2.GetNum()-ordv,kset3.GetNum()-ordw);


	// convert v curves
	for (int i=0; i<numu; i++) 
		for(int k=0; k<numw; k++) {
			Vector<T> v =BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).ConvertCompBezCurvCPoints();
			for (int j=0; j<kset2.GetNum()-ordv; j++) ncpts[k][i][j]=v[j];
		}

	// convert resulting u curves
	for (int j=0; j<kset2.GetNum()-ordv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> v = BspCurv<T>(ncpts.GetU(j,k), *ktsu, ordu,numu).ConvertCompBezCurvCPoints();
			for (int i=0; i<kset1.GetNum()-ordu; i++) ncpts[k][i][j]=v[i];
		}

	// convert resulting u curves
	for (int i=0; i<kset1.GetNum()-ordu; i++) 
		for (int j=0; j<kset2.GetNum()-ordv; j++) {
			Vector<T> v = BspCurv<T>(ncpts.GetW(i,j), *ktsw, ordw,numw).ConvertCompBezCurvCPoints();
			for (int k=0; k<kset3.GetNum()-ordw; k++) ncpts[k][i][j]=v[k];
		}

	// create and return the CompBezVol
	return CompBezVol<T>(kset1.GetNumDistinct()-1, kset2.GetNumDistinct()-1, kset3.GetNumDistinct()-1, ordu, ordv, ordw, ncpts, kset1.GetKnots(), kset2.GetKnots(), kset3.GetKnots()); 
}


// convert to composite poly Vol
template<class T>
CompPolyVol<T> BspVol<T>::ConvertCompPolyVol() const
{
	// convert to composite Bezier form then to poly form
	return (ConvertCompBezVol().ConvertCompPolyVol());
}


// ISOPARAMETRIC CURVES & SURFACES

// return the isoparametric curve in u as a BspCurv in v
template<class T>
BspSurf<T> BspVol<T>::GetIsoparametricU(double u) const
{
	Matrix<T> ncpts(numv,numw);
	// evaluate the u curves at point u
	for (int j=0; j<numv; j++)
		for (int k=0; k<numw; k++)
			ncpts[j][k] = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).Eval(u);
	// return the v curve
	return BspSurf<T>(ncpts,*ktsv,*ktsw,ordv,ordw,numv,numw);
}
	

// return the isoparametric curve in v as a BspCurv in u
template<class T>
BspSurf<T> BspVol<T>::GetIsoparametricV(double v) const
{
	Matrix<T> ncpts(numu,numw);
	// evaluate the v curves at v
	for (int i=0; i<numu; i++)
		for (int k=0; k<numw; k++)
			ncpts[i][k] = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).Eval(v);
	// return the u curve
	return BspSurf<T>(ncpts,*ktsu,*ktsw,ordu,ordw,numu,numw);
}


// return the isoparametric curve in v as a BspCurv in u
template<class T>
BspSurf<T> BspVol<T>::GetIsoparametricW(double w) const
{
	Matrix<T> ncpts(numu,numv);
	// evaluate the v curves at v
	for (int i=0; i<numu; i++)
		for (int j=0; j<numv; j++)
			ncpts[i][j] = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).Eval(w);
	// return the u curve
	return BspSurf<T>(ncpts,*ktsu,*ktsv,ordu,ordv,numu,numv);
}


// return the isoparametric curve in v as a BspCurv in u
template<class T>
BspCurv<T> BspVol<T>::GetIsoparametricUV(double u, double v) const
{
	Vector<T> ncpts(numw);
	// evaluate the v curves at v
	for (int k=0; i<numw; k++)
		ncpts[k] = BspSurf<T>((*cpts).GetUV(k),*ktsu,*ktsv,ordu,ordv,numu,numv).Eval(u,v);
	// return the u curve
	return BspCurv<T>(ncpts,*ktsw,ordw,numw);
}


// return the isoparametric curve in v as a BspCurv in u
template<class T>
BspCurv<T> BspVol<T>::GetIsoparametricUW(double u, double w) const
{
	Vector<T> ncpts(numv);
	// evaluate the v curves at v
	for (int j=0; j<numv; j++)
		ncpts[j] = BspSurf<T>((*cpts).GetUW(j),*ktsu,*ktsw,ordu,ordw,numu,numw).Eval(u,w);
	// return the u curve
	return BspCurv<T>(ncpts,*ktsv,ordv,numv);
}


// return the isoparametric curve in v as a BspCurv in u
template<class T>
BspCurv<T> BspVol<T>::GetIsoparametricVW(double v, double w) const
{
	Vector<T> ncpts(numu);
	// evaluate the v curves at v
	for (int i=0; i<numu; i++)
		ncpts[i] = BspSurf<T>((*cpts).GetVW(i),*ktsv,*ktsw,ordv,ordw,numv,numw).Eval(v,w);
	// return the u curve
	return BspCurv<T>(ncpts,*ktsu,ordu,numu);
}

// READ and WRITE

template <class T>
void BspVol<T>::write(std::ostream& os) 
{
	os << "RBspline Volume\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "order in w is " << ordw << "\n";
	os << "number of control points in u is " << numu << "\n";
	os << "number of control points in v is " << numv << "\n";
	os << "number of control points in w is " << numw << "\n";
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are \n";
	os << *ktsv;
	
	os << "\nknots in w are \n";
	os << *ktsw;
	
	os << "\ncontrol points are\n";
	os << *cpts;
}

template <class T>
void BspVol<T>::read(std::istream& is)
{
	int Ordu, Ordv, Ordw;
	std::cout << "order of RBspline Volume in u and v and w\n";
	is >> Ordu >> Ordv >> Ordw;
	std::cout << "number of control points in u and v and w\n";
	int Numu, Numv, Numw;
	is >> Numu >> Numv >> Numw;
	Matrix3D<T> Cpts(Numu,Numv,Numw);
	Vector<double> Ktsu(Ordu+Numu), Ktsv(Ordv+Numv), Ktsw(Ordw+Numw);
	std::cout << "knots in u\n";
	is >> Ktsu;
	std::cout << "knots in v\n";
	is >> Ktsv;
	std::cout << "knots in w\n";
	is >> Ktsw;
	
	std::cout << "\ninput control points\n";
	is >> Cpts;
	*this = BspVol<T>(Cpts,Ktsu,Ktsv,Ktsw,Ordu,Ordv,Ordw,Numu,Numv,Numw);
} 


template <class T>
void BspVol<T>::writefile(std::ofstream& ofs)
{

	ofs << ordu << " " << ordv << " " << ordw << "\n";
	ofs << numu << " " << numv << " " <<  numw << "\n";
	ofs << *ktsu;

	ofs << *ktsv;
	

	ofs << *ktsw;

	ofs << *cpts;
}


template<class T>
Matrix3D<T> BspVol<T>::ComputePoints(int n, int m, int p) const
{
	Matrix3D<T> v(n+1,m+1, p+1);

	double stepw = ((*ktsw)[numw]-(*ktsw)[ordw-1])/p;

		// evaluate
	double valw;
	for (int k=0; k<=p; k++) {
		valw = (*ktsw)[ordw-1]+k*stepw;
		BspSurf<T> surf = GetIsoparametricW(valw);
		v.InsertUV(surf.ComputePoints(n,m),k);
	}
	return v;
}


template <class T>
void BspVol<T>::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw;
	int Numu, Numv, Numw;
	ifs >> Ordu >> Ordv >> Ordw;
	ifs >> Numu >> Numv >> Numw;

	Matrix3D<T> Cpts(Numu,Numv,Numw);

	Vector<double> Ktsu(Ordu+Numu), Ktsv(Ordv+Numv), Ktsw(Ordw+Numw);
	ifs >> Ktsu;
	ifs >> Ktsv;
	ifs >> Ktsw;
	ifs >> Cpts;
	
	*this = BspVol<T>(Cpts,Ktsu,Ktsv,Ktsw,Ordu,Ordv,Ordw,Numu,Numv,Numw);
} 


// PRIVATE FUNCTIONS

// remove all possible knots from the Volace
// remove possible u knots then v knots
template<class T>
BspVol<T> BspVol<T>::RemovePossibleKnotsUV(double tol) const
{
	return RemovePossibleKnotsU(tol).RemovePossibleKnotsV(tol);
}

// remove all possible knots from the Volace
// remove possible u knots then v knots
template<class T>
BspVol<T> BspVol<T>::RemovePossibleKnotsUW(double tol) const
{
	return RemovePossibleKnotsU(tol).RemovePossibleKnotsW(tol);
}


// remove all possible knots from the Volace
// remove possible u knots then v knots
template<class T>
BspVol<T> BspVol<T>::RemovePossibleKnotsVW(double tol) const
{
	return RemovePossibleKnotsV(tol).RemovePossibleKnotsW(tol);
}

template<class T> 
BspVol<T> BspVol<T>::MakeCompatableU(const BspVol<T>& b) const
{
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetu).Normalise(b.GetKnotSetU());

	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(b.GetKnotSetU());

	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);

	if (kdiff1.GetNum() > 0) 
		return BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetMult(),kdiff1.GetNumDistinct());
	else return BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw);
}

template<class T> 
BspVol<T> BspVol<T>::MakeCompatableV(const BspVol<T>& b) const
{
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetv).Normalise(b.GetKnotSetV());

	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(b.GetKnotSetV());

	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);

	if (kdiff1.GetNum() > 0) 
		return BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetMult(),kdiff1.GetNumDistinct());
	else return BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw);
}


template<class T> 
BspVol<T> BspVol<T>::MakeCompatableW(const BspVol<T>& b) const
{
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetw).Normalise(b.GetKnotSetW());

	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(b.GetKnotSetW());

	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);

	if (kdiff1.GetNum() > 0) 
		return BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw).InsertKnotW(kdiff1.GetDistinctKnots(),kdiff1.GetMult(),kdiff1.GetNumDistinct());
	else return BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw);
}


template<class T> 
BspVol<T> BspVol<T>::MakeCompatableUV(const BspVol<T>& b) const
{
	return MakeCompatableV(MakeCompatableU(b));
}


template<class T> 
BspVol<T> BspVol<T>::MakeCompatableUW(const BspVol<T>& b) const
{
	return MakeComptableW(MakeCompatableU(b));
}

	
template<class T> 
BspVol<T> BspVol<T>::MakeCompatableVW(const BspVol<T>& b) const
{
	return MakeComptableW(MakeComptableV(b));
}

template<class T> 
BspVol<T> BspVol<T>::MakeKnotSetCompatableU(const KnotSet& KsetU) const
{
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetu).Normalise(KsetU);
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(KsetU);
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
	// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,norm.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,numu,numv,numw);
}

template<class T> 
BspVol<T> BspVol<T>::MakeKnotSetCompatableV(const KnotSet& KsetV) const
{
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetv).Normalise(KsetV);
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(KsetV);
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,*ktsu,norm.GetKnots(),*ktsw,ordu,ordv,ordw,numu,numv,numw);
}

template<class T> 
BspVol<T> BspVol<T>::MakeKnotSetCompatableW(const KnotSet& KsetW) const
{
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetv).Normalise(KsetW);
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(KsetW);
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspVol<T>(*cpts,*ktsu,*ktsv,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw).InsertKnotW(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspVol<T>(*cpts,*ktsu,*ktsw,norm.GetKnots(),ordu,ordv,ordw,numu,numv,numw);
}

template<class T>
BspVol<T> BspVol<T>::MakeKnotSetCompatable(const KnotSet& KsetU, const KnotSet& KsetV, const KnotSet& KsetW) const
{
	return MakeKnotSetCompatableW(KsetW).MakeKnotSetCompatableV(KsetV).MakeKnotSetCompatableU(KsetU);
}


// elevate the degree of the BspVol by level
template<class T>
BspVol<T> BspVol<T>::ElevateU(int level) const
{
	// create knot set
	KnotSet kset = (*ksetu).CreateKnotSetElevate(level);
	// return the new BspVol
	return BspVol<T>(ElevateCPointsU(level), kset.GetKnots(), *ktsv, *ktsw, ordu+level, ordv, ordw, kset.GetNum()-(ordu+level), numv, numw);
}


// elevate the degree of the BspVol by level
template<class T>
Matrix3D<T> BspVol<T>::ElevateCPointsU(int level) const
{
	KnotSet kset = (*ksetu).CreateKnotSetElevate(level);
	// create matrix of control points
	Matrix3D<T> ncpts(kset.GetNum()-(ordu+level),numv,numw);

	// elevate each u curve
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).ElevateCPoints(level);
			// extract control points
			for (int i=0; i<kset.GetNum()-(ordu+level); i++) ncpts[k][i][j]=temp[i];
		}
	// return the new points
	return ncpts;
}


// elevate the degree of the BspVol by level
template<class T>
BspVol<T> BspVol<T>::ElevateV(int level) const
{
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetElevate(level);
	// return the new BspVol
	return BspVol<T>(ElevateCPointsV(level), *ktsu, kset.GetKnots(),*ktsw, ordu, ordv+level, ordw, numu, kset.GetNum()-(ordv+level),  numw);
}


// elevate the degree of the BspVol by level
template<class T>
Matrix3D<T> BspVol<T>::ElevateCPointsV(int level) const
{
	KnotSet kset = (*ksetv).CreateKnotSetElevate(level);
	// create matrix of control points
	Matrix3D<T> ncpts(numu,kset.GetNum()-(ordv+level),numw);

	// elevate each v curve
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).ElevateCPoints(level);
			// extract control points
			for (int j=0; j<kset.GetNum()-(ordv+level); j++) ncpts[k][i][j]=temp[j];
		}
	// return the new BspVol
	return ncpts;
}


// elevate the degree of the BspVol by level
template<class T>
Matrix3D<T> BspVol<T>::ElevateCPointsW(int level) const
{	
	KnotSet kset = (*ksetw).CreateKnotSetElevate(level);
	// create matrix of control points
	Matrix3D<T> ncpts(numu,numv,kset.GetNum()-(ordw+level));

	// elevate each v curve
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).ElevateCPoints(level);
			// extract control points
			for (int k=0; k<kset.GetNum()-(ordw+level); k++) ncpts[k][i][j]=temp[k];
		}
	// return the new BspVol
	return ncpts;
}


// elevate the degree of the BspVol by level
template<class T>
BspVol<T> BspVol<T>::ElevateW(int level) const
{
	// create the knot set
	KnotSet kset = (*ksetw).CreateKnotSetElevate(level);
	// return the new BspVol
	return BspVol<T>(ElevateCPointsW(level), *ktsu, *ktsv, kset.GetKnots(), ordu, ordv, ordw+level, numu, numv, kset.GetNum()-(ordw+level));
}

// compute the u derivative of the BspVol of order deriv and
// represent the result as another BspVol
template<class T>
Matrix3D<T> BspVol<T>::DeriveCPointsU(int level) const
{   
	// compute knot set to find size of control point array
	KnotSet kset=(*ksetu).CreateKnotSetDeriv(level);
	int numcpts = kset.GetNum()-kset.GetOrd(); 

	// create matrix of control points
	Matrix3D<T> ncpts(numcpts,numv,numw);

	// only need to compute knot vector once, use version that
	// returns control points only
	// for derivative of each u curve
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).DeriveCPoints(level);
			// extract control points
			for (int i=0; i<numcpts; i++) ncpts[k][i][j]=temp[i];
		}
	// return the new BspVol
	return ncpts;
}

// compute the u derivative of the BspVol of order deriv and
// represent the result as another BspVol
template<class T>
BspVol<T> BspVol<T>::DeriveU(int level) const
{   
	// compute knot set to find size of control point array
	KnotSet kset=(*ksetu).CreateKnotSetDeriv(level);
	int numcpts = kset.GetNum()-kset.GetOrd(); 
	// return the new BspVol
	return BspVol<T>(DeriveCPointsU(level), kset.GetKnots(), *ktsv, *ktsw, ordu-level, ordv, ordw, numcpts, numv, numw);
}


// compute the derivative of the BspVol of order deriv and
// represent the result as another BspVol
template<class T>
BspVol<T> BspVol<T>::DeriveV(int level) const
{   
	// compute knot set to find size of control point array
	KnotSet kset= (*ksetv).CreateKnotSetDeriv(level);
	
	int numcpts = kset.GetNum()-kset.GetOrd(); 
	
	return BspVol<T>(DeriveCPointsV(level), *ktsu, kset.GetKnots(), *ktsw, ordu, ordv-level, ordw, numu, numcpts, numw);
}

// compute the derivative of the BspVol of order deriv and
// represent the result as another BspVol
template<class T>
Matrix3D<T> BspVol<T>::DeriveCPointsV(int level) const
{   
	// compute knot set to find size of control point array
	KnotSet kset= (*ksetv).CreateKnotSetDeriv(level);
	int numcpts = kset.GetNum()-kset.GetOrd(); 
	
	// create matrix of control points
	Matrix3D<T> ncpts(numu,numcpts,numw);
	
	// for derivative of each v curve
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).DeriveCPoints(level);
			// extract control points
			for (int j=0; j<numcpts; j++) ncpts[k][i][j]=temp[j];
		}
	// return the new BspVol
	return ncpts;
}


// compute the derivative of the BspVol of order deriv and
// represent the result as another BspVol
template<class T>
BspVol<T> BspVol<T>::DeriveW(int level) const
{
	// compute knot set to find size of control point array
	KnotSet kset= (*ksetw).CreateKnotSetDeriv(level);
	int numcpts = kset.GetNum()-kset.GetOrd(); 
	// return the new BspVol
	return BspVol<T>(DeriveCPointsW(level), *ktsu, *ktsv, kset.GetKnots(), ordu, ordv, ordw-level, numu, numv, numcpts);
}



// compute the derivative of the BspVol of order deriv and
// represent the result as another BspVol
template<class T>
Matrix3D<T> BspVol<T>::DeriveCPointsW(int level) const
{
	// compute knot set to find size of control point array
	KnotSet kset= (*ksetw).CreateKnotSetDeriv(level);
	int numcpts = kset.GetNum()-kset.GetOrd(); 
	
	// create matrix of control points
	Matrix3D<T> ncpts(numu,numv,numcpts);
	
	// for derivative of each v curve
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).DeriveCPoints(level);
			// extract control points
			for (int k=0; k<numcpts; k++) ncpts[k][i][j]=temp[k];
		}
	// return the new BspVol
	return ncpts;
}


// evaluate the derivative of the BspVol of order deriv at a point
// x. Computes the derivative as a BspVol and then evaluates this at x
template<class T>
T BspVol<T>::DeriveU(int level, double u, double v, double w) const
{
	return DeriveU(level)(u,v,w);
}

// evaluate the derivative of the BspVol of order deriv at a point
// x. Computes the derivative as a BspVol and then evaluates this at x
template<class T>
T BspVol<T>::DeriveV(int level, double u, double v, double w) const
{
	return DeriveV(level)(u,v,w);
}

// evaluate the derivative of the BspVol of order deriv at a point
// x. Computes the derivative as a BspVol and then evaluates this at x
template<class T>
T BspVol<T>::DeriveW(int level, double u, double v, double w) const
{
	return DeriveW(level)(u,v,w);
}


// compute the indefinite integral of the BspCurv and represent
// it as a BspCurv of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateU(int level) const
{
	BspVol<T> b = IntegrateU();

	if (level >= 0) {
		for (int i=1; i<level; i++) b = b.IntegrateU();
		return b;
	}
	else return *this;
}	

// compute the indefinite integral of the BspCurv and represent
// it as a BspCurv of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateW(int level) const
{
	BspVol<T> b = IntegrateW();

	if (level >= 0) {
		for (int i=1; i<level; i++) b = b.IntegrateW();
		return b;
	}
	else return *this;
}	

// compute the indefinite integral of the BspCurv and represent
// it as a BspCurv of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateV(int level) const
{
	BspVol<T> b = IntegrateV();

	if (level >= 0) {
		for (int i=1; i<level; i++) b = b.IntegrateV();
		return b;
	}
	else return *this;
}	



// integrate the BspVol between the limits u1,u2 and v1,v2. 
template<class T>
T BspVol<T>::Integrate1(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Vector<T> w(numw);

	// Integrate each curve in u
	for (int k=0; k<numw; k++)
		w[k] = BspSurf<T>(GetUV(k),*ktsu,*ktsv,ordu,ordv,numu,numv).Integrate(u1,u2,v1,v2);
	// integrate resulting curve in v
	return BspCurv<T>(w,*ktsw,ordw,numw).Integrate(w1,w2);
}
// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateUV1() const
{
	// create knot sets
	KnotSet ksetu = (*ksetu).CreateKnotSetIntegrate();
	KnotSet ksetv = (*ksetv).CreateKnotSetIntegrate();
	
	// create the array of new control points
	Matrix3D<T> ncpts(numu+1,numv+1,numw);

	// integrate the u curves
	for (int i=0; i<numu; i++) {
		Vector<T> v1 = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).IntegrateCPoints();                      
		for (int k=0; k<numv+1; k++) ncpts[i][k]=v1[k];
	}

	// integrate the v curves
	for (int j=0; j<numv+1; j++) {
		Vector<T> v2 = BspCurv<T>(ncpts.GetCol(j),*ktsu,ordu,numu).IntegrateCPoints();
		for (int k=0; k<numu+1; k++) ncpts[k][j]=v2[k];
	}
	// create and return the BspVol
	return BspVol<T>(ncpts,ksetu.GetKnots(),ksetv.GetKnots(),ordu+1,ordv+1,numu+1,numv+1);
}       

// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateUV2() const
{
	// integrateU then V
	return IntegrateU().IntegrateV();
}       



// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateUW1() const
{
	// integrateU then V
	return IntegrateU().IntegrateW();
}       


// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateVW1() const
{
	// integrateU then V
	return IntegrateV().IntegrateW();
}   


// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
Matrix3D<T> BspVol<T>::IntegrateCPointsU() const
{
	// create the array of control points
	Matrix3D<T> ncpts(numu+1,numv,numw);
	
	// integrate each u curve
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> v = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).IntegrateCPoints();                       
			for (int i=0; i<numu+1; i++) ncpts[k][i][j]=v[i];
		}
	
	// create and return the BspVol
	return ncpts;
}       


// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateU() const
{

	// create knot sets
	KnotSet kset = (*ksetu).CreateKnotSetIntegrate();
	// create and return the BspVol
	return BspVol<T>(IntegrateCPointsU(),kset.GetKnots(),*ktsv,*ktsw,ordu+1,ordv,ordw,numu+1,numv, numw);
}       


// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateV() const
{
	// create knot sets
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();

	// create and return the new BspVol
	return BspVol<T>(IntegrateCPointsV(),*ktsu,kset.GetKnots(),*ktsw,ordu,ordv+1,ordw, numu,numv+1,numw);
}       



// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
Matrix3D<T> BspVol<T>::IntegrateCPointsV() const
{
	// create the matrix of control points
	Matrix3D<T> ncpts(numu,numv+1,numw);

	// integrate each v curve
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> v = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).IntegrateCPoints();                       
			for (int j=0; j<numv+1; j++) ncpts[k][i][j]=v[j];
		}
	// create and return the new BspVol
	return ncpts;
}       


// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateW() const
{
	// create knot sets
	KnotSet kset = (*ksetw).CreateKnotSetIntegrate();
	// create and return the new BspVol
	return BspVol<T>(IntegrateCPointsW(),*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw+1,numu,numv,numw+1);
}       


// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateUV() const
{
	// create knot sets
	KnotSet kset1 = (*ksetu).CreateKnotSetIntegrate();
	KnotSet kset2 = (*ksetv).CreateKnotSetIntegrate();
	
	// create and return the new BspVol
	return BspVol<T>(IntegrateCPointsUV(),kset1.GetKnots(),kset2.GetKnots(),*ktsw,ordu+1,ordv+1,ordw,numu+1,numv+1,numw);
}   

// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateUW() const
{
	// create knot sets
	KnotSet kset1 = (*ksetu).CreateKnotSetIntegrate();
	KnotSet kset2 = (*ksetw).CreateKnotSetIntegrate();
	
	// create and return the new BspVol
	return BspVol<T>(IntegrateCPointsUV(),kset1.GetKnots(),*ktsv,kset2.GetKnots(),ordu+1,ordv,ordw+1,numu+1,numv,numw+1);
}   

// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
BspVol<T> BspVol<T>::IntegrateVW() const
{
	// create knot sets
	KnotSet kset1 = (*ksetv).CreateKnotSetIntegrate();
	KnotSet kset2 = (*ksetw).CreateKnotSetIntegrate();
	
	// create and return the new BspVol
	return BspVol<T>(IntegrateCPointsVW(),*ktsu,kset1.GetKnots(),kset2.GetKnots(),ordu,ordv+1,ordw+1,numu,numv+1,numw+1);
}   

// compute the indefinite integral of the BspVol and represent
// it as a BspVol of one higher degree
template<class T>
Matrix3D<T> BspVol<T>::IntegrateCPointsW() const
{
	// create the matrix of control points
	Matrix3D<T> ncpts(numu,numv,numw+1);
	
	// integrate each v curve
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> v = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).IntegrateCPoints();                       
			for (int k=0; k<numw+1; k++) ncpts[k][i][j]=v[k];
		}
	// create and return the new BspVol
	return ncpts;
}       




// compute the indefinite integral of the BspVol as a BspVol
// and return just the control points
template<class T>
Matrix3D<T> BspVol<T>::IntegrateCPointsUV() const
{
	Matrix3D<T> mat = IntegrateCPointsV();
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();

	return BspVol<T>(mat,*ktsu,kset.GetKnots(),*ktsw,ordu,ordv+1,ordw,numu,numv+1,numw).IntegrateCPointsU();

}


// compute the indefinite integral of the BspVol as a BspVol
// and return just the control points
template<class T>
Matrix3D<T> BspVol<T>::IntegrateCPointsUW() const
{
	Matrix3D<T> mat = IntegrateCPointsW();
	KnotSet kset = (*ksetw).CreateKnotSetIntegrate();

	return BspVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw+1,numu,numv,numw+1).IntegrateCPointsU();

}


// compute the indefinite integral of the BspVol as a BspVol
// and return just the control points
template<class T>
Matrix3D<T> BspVol<T>::IntegrateCPointsVW() const
{
	Matrix3D<T> mat = IntegrateCPointsV();
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();

	return BspVol<T>(mat,*ktsu,kset.GetKnots(),*ktsw,ordu,ordv+1,ordw,numu,numv+1,numw).IntegrateCPointsW();
}
// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
BspVol<T> BspVol<T>::SubdivideU(int level) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(level);
	// return the new BspVol
	return BspVol<T>(SubdivideCPointsU(level), kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu, numv, numw);
}

// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsU(int level) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(level);
	int count = kset.GetNum()-ordu;
	Matrix3D<T> ncpts(count,numv,numw);

	// subdivide each u curve
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) { 
			Vector<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).SubdivideCPoints(level);
		// extract control points
			for (int i=0; i<count; i++) ncpts[k][i][j]=temp[i];
		}
	// return the new BspVol
	return ncpts;
}


// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
BspVol<T> BspVol<T>::SubdivideV(int level) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(level);
	// return the new BspVol
	return BspVol<T>(SubdivideCPointsV(level), *ktsu, kset.GetKnots(), *ktsw, ordu, ordv, ordw, numu, kset.GetNum()-ordv,numw);
}



// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsV(int level) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(level);
	int count = kset.GetNum()-ordv;
	Matrix3D<T> ncpts(numu,count,numw);

	// subdivide each v curve
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).SubdivideCPoints(level);
		// extract control points
			for (int j=0; j<count; j++) ncpts[k][i][j]=temp[j];
		}
	// return the new BspVol
	return ncpts;
}


// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
BspVol<T> BspVol<T>::SubdivideW(int level) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(level);
	// return the new BspVol
	return BspVol<T>(SubdivideCPointsW(level), *ktsu, *ktsv, kset.GetKnots(), ordu, ordv, ordw, numu, numv, kset.GetNum()-ordw);
}



// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsW(int level) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(level);
	int count = kset.GetNum()-ordw;
	Matrix3D<T> ncpts(numu,numv,count);

	// subdivide each v curve
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).SubdivideCPoints(level);
		// extract control points
			for (int k=0; k<count; k++) ncpts[k][i][j]=temp[k];
		}
	// return the new BspVol
	return ncpts;
}


// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsVW(int levv, int levw) const
{
	Matrix3D<T> mat = SudivideCPointsW(levw);
	
	KnotSet  kset = (*ksetw).CreateKnotSetSubdivide(levw);	

	return BspVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw,numu,numv,kset.GetNum()-ordw).SubdivideCPointsV(levv);
}

// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsUV(int levu, int levv) const
{
	Matrix3D<T> mat = SudivideCPointsV(levv);
	
	KnotSet  kset = (*ksetv).CreateKnotSetSubdivide(levv);	

	return BspVol<T>(mat,*ktsu,kset.GetKnots(),*ktsw,ordu,ordv,ordw,numu,kset.GetNum()-ordv,numw).SubdivideCPointsU(levu);
}

// subdivide the BspVol upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspVol after the knot refinement
template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsUW(int levu, int levw) const
{
	Matrix3D<T> mat = SudivideCPointsW(levw);
	
	KnotSet  kset = (*ksetw).CreateKnotSetSubdivide(levw);	

	return BspVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw,numu,numv,kset.GetNum()-ordw).SubdivideCPointsU(levu);
}


// subdivide the BspVol between the limits x1 and x2
// and return the new BspVol
template<class T>
BspVol<T> BspVol<T>::SubdivideU(double u1, double u2) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(u1,u2);
	// return the new BspVol
	return BspVol<T>(SubdivideCPointsU(u1,u2), kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu, numv, numw);
}

template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsU(double u1,double u2) const
{
	KnotSet kset = KnotSet(*ktsu,ordu,ordu+numu).CreateKnotSetSubdivide(u1,u2);
	int count = kset.GetNum()-ordu;
	Matrix3D<T> ncpts(count,numv,numw);

	// subdivide each u curve
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).SubdivideCPoints(u1,u2);
			// extract control points
			for (int i=0; i<count; i++) ncpts[k][i][j]=temp[i];
		}	
	return ncpts;
}


// subdivide the BspVol between the limits x1 and x2
// and return the new BspVol
template<class T>
BspVol<T> BspVol<T>::SubdivideV(double v1, double v2) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);

	// return the new BspVol
	return BspVol<T>(SubdivideCPointsV(v1,v2), *ktsu, kset.GetKnots(), *ktsw, ordu, ordv, ordw, numu, kset.GetNum()-ordv,numw);
}


template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsV(double v1,double v2) const
{
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);
	int count = kset.GetNum()-ordv;
	Matrix3D<T> ncpts(numu,count,numw);

	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).SubdivideCPoints(v1,v2);
			// extract control points
			for (int j=0; j<count; j++) ncpts[k][i][j]=temp[j];
		}
	return ncpts;
}
	

// subdivide the BspVol between the limits x1 and x2
// and return the new BspVol
template<class T>
BspVol<T> BspVol<T>::SubdivideW(double w1, double w2) const
{
	// create subdivision knot set
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(w1,w2);
	// return the new BspVol
	return BspVol<T>(SubdivideCPointsW(w1,w2), *ktsu, *ktsv, kset.GetKnots(), ordu, ordv, ordw, numu, numv, kset.GetNum()-ordw);
}


template<class T>
Matrix3D<T> BspVol<T>::SubdivideCPointsW(double w1,double w2) const
{
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(w1,w2);
	int count = kset.GetNum()-ordw;
	Matrix3D<T> ncpts(numu,numv,count);

	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).SubdivideCPoints(w1,w2);
			// extract control points
			for (int k=0; k<kset.GetNum()-ordw; k++) ncpts[k][i][j]=temp[k];
		}
	return ncpts;
}


// subdivide the BspVol between the limits x1 and x2
// and return the new BspVol
template<class T>
BspVol<T> BspVol<T>::SubdivideUV(double u1, double u2, double v1, double v2) const
{
	// subdivide u then v
	return SubdivideU(u1,u2).SubdivideV(v1,v2);
}


// subdivide the BspVol between the limits x1 and x2
// and return the new BspVol
template<class T>
BspVol<T> BspVol<T>::SubdivideUW(double u1, double u2, double w1, double w2) const
{
	// subdivide u then w
	return SubdivideU(u1,u2).SubdivideW(w1,w2);
}


// subdivide the BspVol between the limits x1 and x2
// and return the new BspVol
template<class T>
BspVol<T> BspVol<T>::SubdivideVW(double v1, double v2, double w1, double w2) const
{
	// subdivide v then w
	return SubdivideV(v1,v2).SubdivideW(w1,w2);
}

// inserts a vector of knots into the BspVol and returns the new
// BspVol
template<class T>
BspVol<T> BspVol<T>::InsertKnotU(const Vector<double>& Kts, int N) const
{
	// find size of array
	if (N <=0) return *this;
	KnotSet kset = (*ksetu).CreateKnotSetInsert(Kts,N);
	Matrix3D<T> ncpts(kset.GetNum()-ordu,numv,numw);

	// insert knots into u curves
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).InsertKnotCPoints(Kts,N);
			for (int i=0; i<kset.GetNum()-ordu; i++) ncpts[k][i][j]=temp[i];
		}
	return BspVol<T>(ncpts,kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu, numv, numw);
}
    
	 
// inserts a knot x into the BspVol  and returns the new BspVol
template<class T>
BspVol<T> BspVol<T>::InsertKnotU(double x) const
{       
	// find size of array
	KnotSet kset = (*ksetu).CreateKnotSetInsert(x);
	Matrix3D<T> ncpts(kset.GetNum()-ordu,numv,numw);
	// insert knots into v curves
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).InsertKnotCPoints(x);
			for (int i=0; i<kset.GetNum()-ordu; i++) ncpts[k][i][j]=temp[i];
		}
	return BspVol<T>(ncpts,kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu, numv, numw);
}                              
      

// inserts a knot x 'level' times into the BspVol and returns
// the new BspVol
template<class T>
BspVol<T> BspVol<T>::InsertKnotU(double x, int level) const
{
	if (level <=0) return *this;
	// create knot set 
	KnotSet kset = (*ksetu).CreateKnotSetInsert(x,level);
	Matrix3D<T> ncpts(kset.GetNum()-ordu,numv,numw);
	// insert into u curves
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).InsertKnotCPoints(x,level);
			for (int i=0; i<kset.GetNum()-ordu; i++) ncpts[k][i][j]=temp[i];
		}
	return BspVol<T>(ncpts,kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu, numv, numw);
}


// inserts a vector of knots into the BspVol and returns the new
// BspVol
template<class T>
BspVol<T> BspVol<T>::InsertKnotV(const Vector<double>& Kts, int N) const
{
	if (N <=0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(Kts,N);
	Matrix3D<T> ncpts(numu,kset.GetNum()-ordv,numw);
	// insert into v curves
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).InsertKnotCPoints(Kts,N);
			for (int j=0; j<kset.GetNum()-ordv; j++) ncpts[k][i][j]=temp[j];
		}
	return BspVol<T>(ncpts,*ktsu, kset.GetKnots(), *ktsw, ordu, ordv,ordw, numu,kset.GetNum()-ordv,numw);
}
    
	 
// inserts a knot x into the BspVol  and returns the new BspVol
template<class T>
BspVol<T> BspVol<T>::InsertKnotV(double x) const
{       
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(x);
	Matrix3D<T> ncpts(numu,kset.GetNum()-ordv,numw);

	// insert into v curves
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).InsertKnotCPoints(x);
			for (int j=0; j<kset.GetNum()-ordv; j++) ncpts[k][i][j]=temp[j];
		}
	return BspVol<T>(ncpts,*ktsu, kset.GetKnots(), *ktsw,ordu, ordv,ordw,numu,kset.GetNum()-ordv,numw);
}                              
      

// inserts a knot x 'level' times into the BspVol and returns
// the new BspVol
template<class T>
BspVol<T> BspVol<T>::InsertKnotV(double x, int level) const
{
	if (level <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(x,level);
	Matrix3D<T> ncpts(numu,kset.GetNum()-ordv,numw);
	// insert into v curves
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).InsertKnotCPoints(x,level);
			for (int j=0; j<kset.GetNum()-ordv; j++) ncpts[k][i][j]=temp[j];
		}
	return BspVol<T>(ncpts,*ktsu, kset.GetKnots(), *ktsw, ordu, ordv,ordw,numu,kset.GetNum()-ordv,numw);
}

// inserts a knot x 'level' times into the BspVol and returns
// the new BspVol
template<class T>
BspVol<T> BspVol<T>::InsertKnotW(double x, int level) const
{
	if (level <=0) return *this;
	// create knot set
	KnotSet kset = (*ksetw).CreateKnotSetInsert(x,level);
	Matrix3D<T> ncpts(numu,numv,kset.GetNum()-ordw);

	// insert into v curves
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).InsertKnotCPoints(x,level);
			for (int k=0; k<kset.GetNum()-ordw; k++) ncpts[k][i][j]=temp[k];
		}
	return BspVol<T>(ncpts,*ktsu, *ktsv, kset.GetKnots(), ordu, ordv,ordw,numu,numv,kset.GetNum()-ordw);    
}

// inserts a knot x 'level' times into the BspVol and returns
// the new BspVol
template<class T>
BspVol<T> BspVol<T>::InsertKnotW(double x) const
{
	// create knot set
	KnotSet kset = (*ksetw).CreateKnotSetInsert(x);
	Matrix3D<T> ncpts(numu,numv,kset.GetNum()-ordw);

	// insert into v curves
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).InsertKnotCPoints(x);
			for (int k=0; k<kset.GetNum()-ordw; k++) ncpts[k][i][j]=temp[k];
		}
	return BspVol<T>(ncpts,*ktsu, *ktsv, kset.GetKnots(), ordu, ordv,ordw,numu,numv,kset.GetNum()-ordw);    
}

template<class T>
BspVol<T> BspVol<T>::InsertKnotW(const Vector<double>& Kts, int N) const
{
	if (N <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetw).CreateKnotSetInsert(Kts,N);
	Matrix3D<T> ncpts(numu,numv,kset.GetNum()-ordw);

	// insert into v curves
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).InsertKnotCPoints(Kts,N);
			for (int k=0; k<kset.GetNum()-ordw; k++) ncpts[k][i][j]=temp[k];
		}		
	return BspVol<T>(ncpts,*ktsu, *ktsv, kset.GetKnots(), ordu, ordv,ordw,numu,numv,kset.GetNum()-ordw);    
}
	     

template<class T>
BspVol<T> BspVol<T>::InsertKnotUV(double u, double v) const
{
	// insert into u then v
	return InsertKnotU(u).InsertKnotV(v);
}
	
template<class T>
BspVol<T> BspVol<T>::InsertKnotUV(double u, double v, int levu, int levv) const
{
	// insert into u then v
	return InsertKnotU(u,levu).InsertKnotV(v,levv);
}

template<class T>
BspVol<T> BspVol<T>::InsertKnotUV(const Vector<double>& Ktsu, int n, const Vector<int>& Ktsv, int m) const
{
	// insert into u then v
	return InsertKnotU(Ktsu,n).InsertKnotV(Ktsv,m);
}


template<class T>
BspVol<T> BspVol<T>::InsertKnotUW(double u, double w) const
{
	// insert into u then v
	return InsertKnotU(u).InsertKnotW(w);
}
	
template<class T>
BspVol<T> BspVol<T>::InsertKnotUW(double u, double w, int levu, int levw) const
{
	// insert into u then v
	return InsertKnotU(u,levu).InsertKnotW(w,levw);
}

template<class T>
BspVol<T> BspVol<T>::InsertKnotUW(const Vector<double>& Ktsu, int n, const Vector<int>& Ktsw, int m) const
{
	// insert into u then v
	return InsertKnotU(Ktsu,n).InsertKnotV(Ktsw,m);
}

	
template<class T>
BspVol<T> BspVol<T>::InsertKnotUW(const Vector<double>& Ktsu, const Vector<int>& multu, int N, const Vector<double>& Ktsw, const Vector<int>& multw, int m) const
{
	// insert into u then w
	return InsertKnotU(Ktsu,multu,n).InsertKnotV(Ktsw,multw,m);
}


template<class T>
BspVol<T> BspVol<T>::InsertKnotUVW(double u, double v, double w) const
{
	// insert into u then v
	return (InsertKnotU(u).InsertKnotV(v)).InsertKnotW(w);
}
	
template<class T>
BspVol<T> BspVol<T>::InsertKnotUVW(double u, double v, double w, int levu, int levv, int levw) const
{
	// insert into u then v
	return (InsertKnotU(u,levu).InsertKnotV(v,levv)).InsertW(w,levw);
}

template<class T>
BspVol<T> BspVol<T>::InsertKnotUVW(const Vector<double>& Ktsu, int n, const Vector<int>& Ktsv, int m, const Vector<double>& Ktsw, int p) const
{
	// insert into u then v
	return (InsertKnotU(Ktsu,n).InsertKnotV(Ktsv,m)).InsertKnotW(Ktsw,p);
}



// conevrt to CompBezVol form using explicit knot insertion in both
// directions
template<class T>
CompBezVol<T> BspVol<T>::ConvertCompBezVol2() const
{
	// create distinct knot and multiplicity arrays
	Vector<double> ndtsu((*ksetu).GetNumDistinct());
	Vector<int> nmultu((*ksetu).GetNumDistinct());

	Vector<double> ndtsv((*ksetv).GetNumDistinct());
	Vector<int> nmultv((*ksetv).GetNumDistinct());

	Vector<double> ndtsw((*ksetw).GetNumDistinct());
	Vector<int> nmultw((*ksetw).GetNumDistinct());


	// create distinct and multiplicity vectors
	// required for insertion
	int count1=0;
	for (int i=1; i<(*ksetu).GetNumDistinct()-1; i++) 
		if ((*ksetu).GetMult()[i] < ordu) {
			ndtsu[count1]=(*ksetu).GetDistinctKnots()[i];
			nmultu[count1]=ordu-(*ksetu).GetMult()[i];
			count1++;
		}

	int count2=0;
	for (int j=1; j<(*ksetv).GetNumDistinct()-1; j++) 
		if ((*ksetv).GetMult()[j] < ordv) {
			ndtsv[count2]=(*ksetv).GetDistinctKnots()[j];
			nmultv[count2]=ordv-(*ksetv).GetMult()[j];
			count2++;
		}

	int count3=0;
	for (int k=1; k<(*ksetw).GetNumDistinct()-1; k++) 
		if ((*ksetw).GetMult()[i] < ordw) {
			ndtsw[count3]=(*ksetw).GetDistinctKnots()[k];
			nmultw[count3]=ordw-(*ksetw).GetMult()[k];
			count3++;
		}


	// insert knots in u and v
	BspVol<T> b = InsertKnotUVW(ndtsu,multu,count1,ndtsv,nmultv,count2,ndtsw,multw,count3);

	// create the CompBezCurv Curve
	return CompBezVol<T>((*ksetu).GetNumDistinct()-1,(*ksetv).GetNumDistinct()-1,(*ksetw).GetNumDistinct()-1,ordu,ordv,ordw,b.GetCPoints(),b.GetKnotsU(), b.GetKnotsV(),b.GetKnotsW());
}
// need to check on the range for the knot
// remove a single knot from the u knot set
template<class T>
BspVol<T> BspVol<T>::KnotRemovalU(double knot) const
{
	KnotSet kset= (*ksetu).CreateKnotSetRemoval(knot);
	// remove knot in u direction
	// create Curves and remove knots 
	Matrix3D<T> ncpts(numu-1,numv, numw);

	// remove knot from each u curve
	for (int j=0; j<numv; j++) 
		for (int k=0; k<numw; k++) {
			BspCurv<T> temp = BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).KnotRemoval(knot);
			// extract control points
			for (int i=0; i<numu-1; i++) ncpts[k][i][j]=temp.GetCPoints()[i];
		}
	// return the new BspVol
	return BspVol<T>(ncpts, kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, numu-1, numv, numw);
}
   
// need to check on the range for the knot
// remove a single knot from the v knot set
template<class T>
BspVol<T> BspVol<T>::KnotRemovalV(double knot) const
{
	KnotSet kset=(*ksetv).CreateKnotSetRemoval(knot);
	// remove knot in v direction
	// create curves and remove knots 
	Matrix3D<T> ncpts(numu,numv-1, numw);
	
	// remove knot from each v curve
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++) {
			BspCurv<T> temp = BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).KnotRemoval(knot);
			// extract control points
			for (int j=0; j<numv-1; j++) ncpts[k][i][j]=temp.GetCPoints()[j];
		}
	// return the new BspVol
	return BspVol<T>(ncpts, *ktsu, kset.GetKnots(), *ktsw, ordu, ordv, ordw, numu, numv-1, numw);
}
   
// need to check on the range for the knot
// remove a single knot from the v knot set
template<class T>
BspVol<T> BspVol<T>::KnotRemovalW(double knot) const
{
	KnotSet kset= (*ksetw).CreateKnotSetRemoval(knot);
	// remove knot in w direction
	// create curves and remove knots 
	Matrix3D<T> ncpts(numu,numv, numw-1);
	
	// remove knot from each v curve
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++) {
			BspCurv<T> temp = BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).KnotRemoval(knot);
			// extract control points
			for (int k=0; k<numw-1; k++) ncpts[k][i][j]=temp.GetCPoints()[k];
		}
	// return the new BspVol
	return BspVol<T>(ncpts, *ktsu, *ktsv, kset.GetKnots(), ordu, ordv, ordw, numu, numv, numw-1);
}


// check to see whether a knot can be removed from the u direction
template<class T>
bool BspVol<T>::IsKnotRemovableU(double knot, double tol) const
{
	// create curves and test knot removal
	// if knot can be removed from all u curves 
	for (int j=0; j<numv; j++)
		for (int k=0; k<numw; k++)
			if (!BspCurv<T>((*cpts).GetU(j,k),*ktsu,ordu,numu).IsKnotRemovable(knot,tol)) return false;
				
	return true;
}
   
// check to see whether a knot can be removed from the v knot set
template<class T>
bool BspVol<T>::IsKnotRemovableV(double knot, double tol) const
{
	// create curves and test knot removal 
	// if knot can be removed from all v curves 
	for (int i=0; i<numu; i++) 
		for (int k=0; k<numw; k++)
			if (!BspCurv<T>((*cpts).GetV(i,k),*ktsv,ordv,numv).IsKnotRemovable(knot,tol)) return false;
				
	return true;
}

// check to see whether a knot can be removed from the v knot set
template<class T>
bool BspVol<T>::IsKnotRemovableW(double knot, double tol) const
{
	// create curves and test knot removal 

	// if knot can be removed from all v curves 
	for (int i=0; i<numu; i++) 
		for (int j=0; j<numv; j++)
			if (!BspCurv<T>((*cpts).GetW(i,j),*ktsw,ordw,numw).IsKnotRemovable(knot,tol)) return false;
				
	return true;
}

   
// remove all possible u knots from the Volace
// within the tolerance tol
template<class T>
BspVol<T> BspVol<T>::RemovePossibleKnotsU(double tol) const
{
	// create local copy
	BspVol<T> temp(*this);

	// for each distinct knot, try to remove it
	for (int i=1; i<(*ksetu).GetNumDistinct()-1; i++)
		for (int j=0; j<(*ksetu).GetMult()[i]; j++)
			if (temp.IsKnotRemovableU((*ksetu).GetDistinctKnots()[i],tol))
				temp = temp.KnotRemovalU((*ksetu).GetDistinctKnots()[i]);
	// return new BezVol
	return temp;
}


// remove all possible v knots from the Vol
template<class T>
BspVol<T> BspVol<T>::RemovePossibleKnotsV(double tol) const
{
	BspVol<T> temp(*this);

	// check for each distinct v knot
	for (int i=1; i<(*ksetv).GetNumDistinct()-1; i++)
		for (int j=0; j<(*ksetv).GetMult()[i]; j++)
			if (temp.IsKnotRemovableV((*ksetv).GetDistinctKnots()[i],tol))
				temp = temp.KnotRemovalV((*ksetv).GetDistinctKnots()[i]);
	// return new BezVol
	return temp;
}

// remove all possible v knots from the Vol
template<class T>
BspVol<T> BspVol<T>::RemovePossibleKnotsW(double tol) const
{
	BspVol<T> temp(*this);

	
	// check for each distinct v knot
	for (int i=1; i<(*ksetw).GetNumDistinct()-1; i++)
		for (int j=0; j<(*ksetw).GetMult()[i]; j++)
			if (temp.IsKnotRemovableW((*ksetw).GetDistinctKnots()[i],tol))
				temp = temp.KnotRemovalW((*ksetw).GetDistinctKnots()[i]);
	// return new BezVol
	return temp;
}




#endif
 
