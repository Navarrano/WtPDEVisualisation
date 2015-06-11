#ifndef BSPSURF
#define BSPSURF

#include "matrix.h"
#include "knotset.h"
#include "compbezsurf.h"

class BspSurfBasisFunc : public Surf<double> {

	// data
	int ordu;	 // order of basis function
	int ordv;
	Ptr<Vector<double> > ktsu;		// knots 
	Ptr<Vector<double> > ktsv;
	Ptr<BspSurf<double> > b;	// BspSurf representation of basis function

	// private functions
	BspSurf<double> CreateBspSurf() const;	// creates the BspSurf
	virtual ObjectID Identity() const { return std::string("class BspSurfBasisFunc"); }
public:

	// constructors
	BspSurfBasisFunc();
	BspSurfBasisFunc(int Ordu, int Ordv, const Vector<double>& Ktsu, const Vector<double>& Ktsv);

	// access functions
	int ComputeDimU() const;
	int ComputeDimV() const;
	BspSurf<double> GetBspSurf() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	double CreateIntegral(double u1, double u2, double v1, double v2) const;


	// evaluators
	double Eval(double u, double v) const;
	virtual double operator()(double u, double v) const;
	virtual double operator()(int, int, double u, double v) const;
	virtual double Derive(int, int, double u, double v) const;
	
	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);

	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
};


class BspSurfBasisFuncSet : public TextObject, public FTextObject {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int numu;
	int numv;
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Matrix<BspSurfBasisFunc> > b;
	virtual ObjectID Identity() const { return std::string("class BspSurfBasisFuncSet"); }
public:
	// constructors
	BspSurfBasisFuncSet();
	BspSurfBasisFuncSet(int Ordu, int Ordv, int Numu, int Numv, const Vector<double>& Ktsu, const Vector<double>& Ktsv);


	// evaluators
	Matrix<double> Eval(double u, double v) const;	// evaluate basis function
	Matrix<double> operator()(double u, double v) const;
	Matrix<double> operator()(int, int, double u, double v) const;
	Matrix<double> Derive(int, int, double u, double v) const;
	Matrix<double> ComputeUBasisMatrix() const;
	Matrix<double> ComputeVBasisMatrix() const;
	Matrix<double> ComputeNMatrix(double u, double v) const;
	Matrix<double> CreateIntegral(double u1, double u2, double v1, double v2) const;
	Matrix<double> CreateIntegralNewU(const BspSurf<double>& s, const BspSurf<double>& norm, double u1, double u2, double v1, double v2) const;
	Matrix<double> CreateIntegralNewV(const BspSurf<double>& s, const BspSurf<double>& norm, double u1, double u2, double v1, double v2) const;
	Matrix<Matrix<double> > CreateMatrixIntegral(int levu, int levv) const;
	Matrix<Matrix<double> > CreateMatrixIntegral(int levu, int levv, double u1, double u2, double v1, double v2) const;
	Matrix<double> CreateMatrixMinimisation(int levu, int levv) const;
	Matrix<double> CreateMatrixMinimisation(int levu, int levv,double u1, double u2, double v1, double v2) const;
	Matrix<double> CreateMatrixMinimisationNewU(int levu, double u1, double u2, double v1, double v2) const;
	Matrix<double> CreateMatrixMinimisationNewV(int levv, double u1, double u2, double v1, double v2) const;
	Matrix<double> CreateMatrixKnotAveragesU() const;
	Matrix<double> CreateMatrixKnotAveragesV() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	Vector<Matrix<double> > GetSmoothingMatrices() const;


	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};

template<class T>
class BspSurf : public Surf<T> {
private:
	// data
	int ordu;
	int ordv;
	int numu;
	int numv;

	// private functions
	Ptr<Matrix<T> > cpts;
	Ptr<Vector<double> > ktsu;
	Ptr<Vector<double> > ktsv;
	Ptr<KnotSet> ksetu;
	Ptr<KnotSet> ksetv;
	double lu1, lu2, lv1, lv2;

	BspSurf<T> MakeCompatableU(const BspSurf<T>& b) const;
	BspSurf<T> MakeCompatableV(const BspSurf<T>& b) const;
	
	
	Matrix<T> ElevateCPointsU(int level) const;
	Matrix<T> ElevateCPointsV(int level) const;


	BspSurf<T> SubdivideU(int level) const;
	BspSurf<T> SubdivideV(int level) const;
	Matrix<T> SubdivideCPointsU(int level) const;
	Matrix<T> SubdivideCPointsV(int level) const;
	BspSurf<T> SubdivideU(double u1, double u2) const;
	BspSurf<T> SubdivideV(double v1, double v2) const;
	Matrix<T> SubdivideCPointsU(double u1, double u2) const;
	Matrix<T> SubdivideCPointsV(double v1, double v2) const;


	Matrix<T> InsertKnotCPointsU(double u) const;
	Matrix<T> InsertKnotCPointsU(double u, int level) const;
	Matrix<T> InsertKnotCPointsU(const Vector<double>& Kts, int n) const;
	Matrix<T> InsertKnotCPointsV(double v) const;
	Matrix<T> InsertKnotCPointsV(double v, int level) const;
	Matrix<T> InsertKnotCPointsV(const Vector<double>& Kts, int n) const;	
	Matrix<T> InsertKnotCPointsUV(double u, double v) const;
	Matrix<T> InsertKnotCPointsUV(double u, double v, int levu, int levv) const;
	Matrix<T> InsertKnotCPointsUV(const Vector<double>& Ktsu, int n, const Vector<int>& Ktsv, int m) const;
	Matrix<T> DeriveCPointsU(int level) const;
	Matrix<T> DeriveCPointsV(int level) const;
	T DeriveU(int level, double u, double v) const;
	T DeriveV(int level, double u, double v) const;
	BspSurf<T> IntegrateU() const;
	BspSurf<T> IntegrateV() const;
	Matrix<T> IntegrateCPointsU() const;
	Matrix<T> IntegrateCPointsV() const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class BspSurf<double>");
		else {
			std::string s(typeid(T).name()), s1("class BspSurf<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	BspSurf();
	BspSurf(const Matrix<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv, int Numu, int Numv);
	BspSurf(const Matrix<T>& Cpts, int Ordu, int Ordv, int Numu, int Numv);
	BspSurf(const PolySurf<T>& p, const KnotSet& KsetU, const KnotSet& KsetV);
	BspSurf(const PolySurf<T>& p, int Ordu, int Ordv, const KnotSet& KsetU, const KnotSet& KsetV);


	// access functions
	int GetOrdU() const;
	int GetNumU() const;
	int GetOrdV() const;
	int GetNumV() const;
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
	Matrix<T> GetCPoints() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	KnotSet GetKnotSetU() const;
	KnotSet GetKnotSetV() const;
	BezSurf<T> GetPatch(int i, int j) const;	
	BspSurfBasisFuncSet GetBasisFuncSet() const;

	// evaluators
	virtual T operator()(double u, double v) const;
	virtual T operator()(int, int, double u, double v) const;
	T Eval(double u, double v) const;
	T Eval1(double u, double v) const;
	Matrix<T> Eval2(double u, double v) const;
	Matrix<T> EvalDeriv(int levu, int levv, double u, double v) const;
	Matrix<T> ComputePoints(int n, int m) const;
	Matrix<T> ComputeKnotSetAverValues() const;

	// knot removal
	BspSurf<T> KnotRemovalU(double knot) const;
	BspSurf<T> KnotRemovalV(double knot) const;
	BspSurf<T> RemovePossibleKnotsU(double tol=knotTol) const;
	BspSurf<T> RemovePossibleKnotsV(double tol=knotTol) const;
	BspSurf<T> RemovePossibleKnots(double tol=knotTol) const;
	bool IsKnotRemovableU(double knot, double tol=cpntTol) const;
	bool IsKnotRemovableV(double knot, double tol=cpntTol) const;

	
	BspSurf<T> MakeKnotSetCompatableU(const KnotSet& KsetU) const;
	BspSurf<T> MakeKnotSetCompatableV(const KnotSet& KsetV) const;
	BspSurf<T> MakeKnotSetCompatable(const KnotSet& KsetU, const KnotSet& KsetV) const;

	BspSurf<T> MakeBreakCompatable(const KnotSet& KsetU, const KnotSet& KsetV) const;
	BspSurf<T> MakeBreakCompatableU(const KnotSet& KsetU) const;
	BspSurf<T> MakeBreakCompatableV(const KnotSet& KsetV) const;
	

	// make surface compatable with another
	template<class T1>
	BspSurf<T> MakeCompatable(const BspSurf<T1>& b) const
	{
		return MakeCompatableU(b).MakeCompatableV(b);
		//return MakeCompatableV(MakeCompatableU(b));
	}
	
	
	template<class T1>
	BspSurf<T> MakeBreakCompatableU(const BspSurf<T1>& b) const
	{
		// assume they are both of the same degree
		// create knot set for current object
	
		// create normalised knot set for current object relative to b
		KnotSet norm(GetKnotSetU());
	
		// add knots to current object in knot set for b not present 
		// in knot set for Current object
		KnotSet add1 = norm.KnotUnion(b.GetKnotSetU());
		// find differences in knots
		KnotSet kdiff1 = add1.KnotDifference(norm);

		if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
			return 	BspSurf<T>(*cpts,norm.GetKnots(),*ktsv,ordu,ordv,numu,numv).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
		} else return BspSurf<T>(*cpts,norm.GetKnots(),*ktsv,ordu,ordv,numu,numv);
	}



	template<class T1>
	BspSurf<T> MakeBreakCompatableV(const BspSurf<T1>& b) const
	{
		// assume they are both of the same degree
		// create knot set for current object
	
		// create normalised knot set for current object relative to b
		KnotSet norm(GetKnotSetV());
	
		// add knots to current object in knot set for b not present 
		// in knot set for Current object
		KnotSet add1 = norm.KnotUnion(b.GetKnotSetV());
	
		// find differences in knots
		KnotSet kdiff1 = add1.KnotDifference(norm);
	
		if (kdiff1.GetNum() > 0) {
			// insert knots into current object to get compatibility with b
			return 	BspSurf<T>(*cpts,*ktsu,norm.GetKnots(),ordu,ordv,numu,numv).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
		} else return BspSurf<T>(*cpts,*ktsu,norm.GetKnots(),ordu,ordv,numu,numv);
	}

	template<class T1>
	BspSurf<T> MakeBreakCompatable(const BspSurf<T1>& b) const
	{
		return MakeBreakCompatableV(b).MakeBreakCompatableU(b);
		//	return MakeBreakCompatableV(MakeBreakCompatableU(b));
	}


	// degree elevation
	Matrix<T> ElevateCPoints(int levu, int levv) const;
	BspSurf<T> Elevate(int levu, int levv) const;
	BspSurf<T> ElevateU(int level) const;
	BspSurf<T> ElevateV(int level) const;

	// product
	template<class T1>
	BspSurf<T> Product(const BspSurf<T1>& b) const
	{
		// make local copies of the two surfaces
		BspSurf<T> d(*this);
		BspSurf<T1> e(b);
		// make the surfaces compatable in u and v
		if (!(*ksetu).IsSameDistinctKnotSet(b.GetKnotSetU())) {
			// needs to add distinct knots from b not in current object
			// into current object (independent of order)
			d = d.MakeBreakCompatableU(e);
			e = e.MakeBreakCompatableU(d);
		}
		if (!(*ksetv).IsSameDistinctKnotSet(b.GetKnotSetV())) {
			// needs to add distinct knots from b not in current object
			// into current object (independent of order)
			d = d.MakeBreakCompatableV(e);
			e = e.MakeBreakCompatableV(d);
		}


		// convert to CompBezSurf
		CompBezSurf<T> e1 = d.ConvertCompBezSurf();
		CompBezSurf<T1> e2 = e.ConvertCompBezSurf();

		// multiply
		CompBezSurf<T> prod = e1.Product(e2);

		// remove all possible knots in u and v
		return prod.ConvertBspSurf().RemovePossibleKnots();
	}
	template<class T1>
	Matrix<T> ProductCPoints(const BspSurf<T1>& b) const
	{
		// make local copies of the two surfaces
		BspSurf<T> a(*this);
		BspSurf<T1> c(b);

		// make the surfaces compatable in u and v
		BspSurf<T> d1 = a.MakeBreakCompatable(c);
		BspSurf<T1> d2 = c.MakeBreakCompatable(d1);

		// convert to CompBezSurf
		CompBezSurf<T> e1 = d1.ConvertCompBezSurf();
		CompBezSurf<T1> e2 = d2.ConvertCompBezSurf();

		// multiply
		CompBezSurf<T> prod = e1.Product(e2);

		// remove all possible knots in u and v
		return ((prod.ConvertBspSurf()).RemovePossibleKnotsUV()).GetCPoints();
	}

	// addition and subtraction
	BspSurf<T> Add(const BspSurf<T>& b) const;
	BspSurf<T> Subtract(const BspSurf<T>& b) const;

	// subdivision
	BspSurf<T> Subdivide(int levu, int levv) const;
	BspSurf<T> SubdivideCPoints(int levu, int levv) const;
	BspSurf<T> Subdivide(double u1, double u2, double v1, double v2) const;
	Matrix<T> SubdivideCPoints(double u1, double u2, double v1, double v2) const;

	// knot insertion	
	BspSurf<T> InsertKnotU(double u) const;
	BspSurf<T> InsertKnotU(double u, int level) const;
	BspSurf<T> InsertKnotU(const Vector<double>& Kts, int n) const;
	BspSurf<T> InsertKnotUV(const Vector<double>& Ktsu, int n, const Vector<double>& Ktsv, int m) const;
	BspSurf<T> InsertKnotUV(const Vector<double>& Ktsu, const Vector<int>& multu, int N, const Vector<double>& Ktsv, const Vector<int>& multv, int m) const;
	BspSurf<T> InsertKnotV(double v) const;
	BspSurf<T> InsertKnotV(double v, int level) const;
	BspSurf<T> InsertKnotV(const Vector<double>& Kts, int n) const;
	BspSurf<T> InsertKnotUV(double u, double v) const;
	BspSurf<T> InsertKnotUV(double u, double v, int levu, int levv) const;
	BspSurf<T> InsertKnotU(const Vector<double>& Kts, const Vector<int>& mult, int N) const;
	BspSurf<T> InsertKnotV(const Vector<double>& Kts, const Vector<int>& mult, int N) const;
	BspSurf<T> InsertKnot(const Vector<double>& Ktsu, const Vector<int>& multu, int N, const Vector<double>& Ktsv, const Vector<int>& multv, int m) const;
	Matrix<T> InsertKnotCPointsU(const Vector<double>& Kts, const Vector<int>& mult, int N) const;
	Matrix<T> InsertKnotCPointsV(const Vector<double>& Kts, const Vector<int>& mult, int N) const;
	Matrix<T> InsertKnotCPointsUV(const Vector<double>& Ktsu, const Vector<int>& multu, int N, const Vector<double>& Ktsv, const Vector<int>& multv, int m) const;

	// derivatives
	BspSurf<T> DeriveU(int level) const;
	BspSurf<T> DeriveV(int level) const;
	BspSurf<T> Derive(int levu, int levv) const;
	Matrix<T> DeriveCPoints(int levu, int levv) const;
	virtual T Derive(int levu, int levv, double u, double v) const;

	// integration
	T IntegrateUV1(double u1, double u2, double v1, double v2) const;
	BspSurf<T> Integrate() const;
	BspSurf<T> IntegrateU(int level) const;
	BspSurf<T> IntegrateV(int level) const;
	Matrix<T> IntegrateCPoints() const;
	T Integrate(double u1, double u2, double v1, double v2) const;
	
	// isoparametric curves
	BspCurv<T> GetIsoparametricU(double u) const;
	BspCurv<T> GetIsoparametricV(double v) const;

	// conversion
	CompBezSurf<T> ConvertCompBezSurf() const;
	CompPolySurf<T> ConvertCompPolySurf() const;

	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};


// CONSTRUCTORS

// constructor builds a BspSurf from a Matrix of control points, knots in u and v
// orders in u and v and number of control points in u and v
template<class T>
BspSurf<T>::BspSurf(const Matrix<T>& Cpts, const Vector<double>& Ktsu, 
			const Vector<double>& Ktsv, int Ordu, int Ordv, int Numu, int Numv) :
			ordu(Ordu), ordv(Ordv), numu(Numu), numv(Numv), 
			cpts(new Matrix<T>(Cpts)), ktsu(new Vector<double>(Ktsu)), 
			ktsv(new Vector<double>(Ktsv)),ksetu(new KnotSet(*ktsu,ordu,ordu+numu)),
			ksetv(new KnotSet(*ktsv,ordv,ordv+numv)), lu1(Ktsu[Ordu-1]), lu2(Ktsu[Numu]), lv1(Ktsv[Ordv-1]), lv2(Ktsv[Numv])
{
}


// default constructor
template<class T>
BspSurf<T>::BspSurf() : ordu(0), ordv(0), numu(0), numv(0), cpts(), ktsu(), ktsv() { }


// constructor builds a BspSurf from a Matrix of control points, 
// orders in u and v and number of control points in u and v
// creates the knots using a simple algorithm
template<class T>
BspSurf<T>::BspSurf(const Matrix<T>& Cpts, int Ordu, int Ordv, int Numu, int Numv) :
			ordu(Ordu), ordv(Ordv), numu(Numu), numv(Numv), cpts(new Matrix<T>(Cpts))
{
	ktsu = new Vector<double>(Math::CreateKnots(numu,ordu));
	ktsv = new Vector<double>(Math::CreateKnots(numv,ordv));
	ksetu = new KnotSet(*ktsu,ordu,ordu+numu);
	ksetv = new KnotSet(*ktsv,ordv,ordv+numv);
	lu1 = (*ktsu)[ordu-1];
	lu2 = (*ktsu)[numu];
	lv1 = (*ktsv)[ordv-1];
	lv2 = (*ktsv)[numv];
}


template<class T>
BspSurf<T>::BspSurf(const PolySurf<T>& p, const KnotSet& KsetU, const KnotSet& KsetV) : ordu(KsetU.GetOrd()), ordv(KsetV.GetOrd()), 
				numu(KsetU.GetNum()-KsetU.GetOrd()), numv(KsetV.GetNum()-KsetV.GetOrd()),
					cpts(new Matrix<T>(p.Elevate(KsetU.GetOrd()-p.GetOrdU(),KsetV.GetOrd()-p.GetOrdV()).ConvertBspSurf().MakeKnotSetCompatable(KsetU,KsetV).GetCPoints())), 
					ktsu(new Vector<double>(KsetU.GetKnots())), ktsv(new Vector<double>(KsetV.GetKnots())), ksetu(new KnotSet(KsetU)), ksetv(new KnotSet(KsetV)),
					lu1((*ktsu)[ordu-1]), lu2((*ktsu)[numu]), lv1((*ktsv)[ordv-1]), lv2((*ktsv)[numv])					
{
}


template<class T>
BspSurf<T>::BspSurf(const PolySurf<T>& p, int Ordu, int Ordv, const KnotSet& KsetU, const KnotSet& KsetV) 
{
	*this = p.ConvertBspSurf().MakeBreakCompatable(KsetU,KsetV);

}



// ACCESS FUNCTIONS

// get the order of the BspSurf in u
template<class T>
int BspSurf<T>::GetOrdU() const { return ordu; }

// get the order of the BspSurf in v
template<class T>
int BspSurf<T>::GetOrdV() const { return ordv; }


// get the number of control points in u
template<class T>
int BspSurf<T>::GetNumU() const { return numu; }

// get the number of control points in v
template<class T>
int BspSurf<T>::GetNumV() const { return numv; }

template<class T>
double BspSurf<T>::GetLeftLimitU() const
{
	return (*ktsu)[ordu-1];
}

template<class T>
double BspSurf<T>::GetRightLimitU() const
{
	return (*ktsu)[numu];
}

template<class T>
double BspSurf<T>::GetLeftLimitV() const
{
	return (*ktsv)[ordv-1];
}

template<class T>
double BspSurf<T>::GetRightLimitV() const
{
	return (*ktsv)[numv];
}

// get the control points
template<class T>
Matrix<T> BspSurf<T>::GetCPoints() const { return *cpts; }

// get the knot vector in u
template<class T>
Vector<double> BspSurf<T>::GetKnotsU() const { return *ktsu; }


// get the knot vector in v
template<class T>
Vector<double> BspSurf<T>::GetKnotsV() const { return *ktsv; }



// get the knot vector in v
template<class T>
KnotSet BspSurf<T>::GetKnotSetU() const { return *ksetu; }


// get the knot vector in v
template<class T>
KnotSet BspSurf<T>::GetKnotSetV() const { return *ksetv; }

template<class T>
BspSurfBasisFuncSet BspSurf<T>::GetBasisFuncSet() const
{
	return BspSurfBasisFuncSet(ordu,ordv,ordu+numu,ordv+numv,*ktsu, *ktsv);
}


// Get the i,j pacth of the surface and return as a BspSurf
template<class T>	
BezSurf<T> BspSurf<T>::GetPatch(int i, int j) const
{
	// create knotset objects
	
	// subdivide at appropriate knots
	BspSurf<T> b = Subdivide((*ksetu).GetDistinctKnots()[i-1],(*ksetu).GetDistinctKnots()[i],
				(*ksetv).GetDistinctKnots()[j-1],(*ksetv).GetDistinctKnots()[j]);
	return BezSurf<T>(b.GetCPoints(),b.GetKnotsU(),b.GetKnotsV(),ordu,ordv);
}


// EVALUATORS

// evaluate the BspSurf at the point x using de Boor algorithm
template<class T>
T BspSurf<T>::operator()(double u, double v) const
{
//	if (u < (*ktsu)[ordu-1] || u > (*ktsu)[numu]) return 0.0;
//	if (v < (*ktsv)[ordv-1] || v > (*ktsv)[numv]) return 0.0;

	Vector<T> v2(numu);

	// use tensor product evaluation
	// evaluate in v and create points in u
	for (int i=0; i<numu; i++) v2[i]= BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv)(v);

	// evaluate the resulting curve
    return BspCurv<T>(v2,*ktsu,ordu,numu)(u);
}

// evaluate the BspSurf at the point x using de Boor algorithm
template<class T>
T BspSurf<T>::Eval(double u, double v) const
{
	if (u < (*ktsu)[ordu-1] || u > (*ktsu)[numu]) return 0.0;
	if (v < (*ktsv)[ordv-1] || v > (*ktsv)[numv]) return 0.0;
	Vector<T> v2(numu);

	// use tensor product evaluation
	// evaluate in v and create points in u
	for (int i=0; i<numu; i++) v2[i]= BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).Eval(v);

	// evaluate the resulting curve
    return BspCurv<T>(v2,*ktsu,ordu,numu).Eval(u);
}

// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
T BspSurf<T>::Eval1(double u, double v) const
{
	Vector<double> uvec = ksetu->CreateVectorInterp(u);
	Vector<double> vvec = ksetv->CreateVectorInterp(v);

	return Math::mult0(Math::mult4(uvec,*cpts),vvec);
}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Matrix<T> BspSurf<T>::Eval2(double u, double v) const
{
	Matrix<double> mu(ksetu->CreateVectorInterp(u));
	Matrix<double> mv(ksetv->CreateVectorInterp(v));

	return Math::mult1(Math::mult1(mu,*cpts),Math::transpose(mv));
}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Matrix<T> BspSurf<T>::EvalDeriv(int levu, int levv, double u, double v) const
{
	Matrix<double> mu(ksetu->CreateVectorInterpDeriv(levu,u));
	Matrix<double> mv(ksetv->CreateVectorInterpDeriv(levv,v));

	Matrix<T> m1 = Math::mult1(Math::mult1(mu,ksetu->CreateMatrixDeriv(levu)),*cpts);
	return Math::mult1(m1,Math::transpose(Math::mult1(mv,ksetv->CreateMatrixDeriv(levv))));
}


template<class T>
Matrix<T> BspSurf<T>::ComputeKnotSetAverValues() const
{
	Matrix<T> m(numu,numv);
	Vector<double> uval((*ksetu).ComputeKnotSetAver());
	Vector<double> vval((*ksetv).ComputeKnotSetAver());


	for (int i=0; i<numu; i++)
		for (int j=0; j<numv; j++) m[i][j] = (*this)(uval[i],vval[j]);

	return m;
}


template<class T>
Matrix<T> BspSurf<T>::ComputePoints(int n, int m) const
{
	Matrix<T> v(n+1,m+1);

	double stepu = ((*ktsu)[numu]-(*ktsu)[ordu-1])/n;

		// evaluate
	double valu;
	for (int i=0; i<=n; i++) {
		valu = (*ktsu)[ordu-1]+i*stepu;
		BspCurv<T> curv = GetIsoparametricU(valu);
		Vector<T> vec = curv.ComputePoints(m);
		v.InsertRow(vec,i);
	}

	return v;
}



// KNOT REMOVAL

// need to check on the range for the knot
// remove a single knot from the u knot set
template<class T>
BspSurf<T> BspSurf<T>::KnotRemovalU(double knot) const
{
	KnotSet kset = (*ksetu).CreateKnotSetRemoval(knot);
	// remove knot in u direction
	// create Curves and remove knots 
	Matrix<T> ncpts(numu-1,numv);

	// remove knot from each u curve
	for (int j=0; j<numv; j++) {
		BspCurv<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).KnotRemoval(knot);
		// extract control points
		for (int i=0; i<numu-1; i++) ncpts[i][j]=temp.GetCPoints()[i];
	}
	// return the new BspSurf
	return BspSurf<T>(ncpts, kset.GetKnots(), *ktsv, ordu, ordv, numu-1, numv);
}

   
// need to check on the range for the knot
// remove a single knot from the v knot set
template<class T>
BspSurf<T> BspSurf<T>::KnotRemovalV(double knot) const
{
	KnotSet kset = (*ksetv).CreateKnotSetRemoval(knot);
	// remove knot in v direction
	// create curves and remove knots 
	Matrix<T> ncpts(numu,numv-1);
	
	// remove knot from each v curve
	for (int i=0; i<numu; i++) {
		BspCurv<T> temp = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).KnotRemoval(knot);
		// extract control points
		for (int j=0; j<numv-1; j++) ncpts[i][j]=temp.GetCPoints()[j];
	}
	// return the new BspSurf
	return BspSurf<T>(ncpts, *ktsu, kset.GetKnots(), ordu, ordv, numu, numv-1);
}
   

// check to see whether a knot can be removed from the u direction
template<class T>
bool BspSurf<T>::IsKnotRemovableU(double knot, double tol) const
{
	// create curves and test knot removal

	// if knot can be removed from all u curves 
	for (int j=0; j<numv; j++) 
			if (!BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).IsKnotRemovable(knot,tol)) return false;
				
	return true;
}
   
// check to see whether a knot can be removed from the v knot set
template<class T>
bool BspSurf<T>::IsKnotRemovableV(double knot, double tol) const
{
	// create curves and test knot removal 

	// if knot can be removed from all v curves 
	for (int i=0; i<numu; i++) 
			if (!BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).IsKnotRemovable(knot,tol)) return false;
				
	return true;
}
   
// remove all possible u knots from the surface
// within the tolerance tol
template<class T>
BspSurf<T> BspSurf<T>::RemovePossibleKnotsU(double tol) const
{
	// create local copy
	BspSurf<T> temp(*this);

	// for each distinct knot, try to remove it
	for (int i=1; i<(*ksetu).GetNumDistinct()-1; i++)
		for (int j=0; j<(*ksetu).GetMult()[i]; j++)
			if (temp.IsKnotRemovableU((*ksetu).GetDistinctKnots()[i],tol))
				temp = temp.KnotRemovalU((*ksetu).GetDistinctKnots()[i]);
	// return new surface
	return temp;
}


// remove all possible v knots from the surface
template<class T>
BspSurf<T> BspSurf<T>::RemovePossibleKnotsV(double tol) const
{
	BspSurf<T> temp(*this);


	// check for each distinct v knot
	for (int i=1; i<(*ksetv).GetNumDistinct()-1; i++)
		for (int j=0; j<(*ksetv).GetMult()[i]; j++)
			if (temp.IsKnotRemovableV((*ksetv).GetDistinctKnots()[i],tol))
				temp = temp.KnotRemovalV((*ksetv).GetDistinctKnots()[i]);
	// return new surface
	return temp;
}

// remove all possible knots from the surface
// remove possible u knots then v knots
template<class T>
BspSurf<T> BspSurf<T>::RemovePossibleKnots(double tol) const
{
	return RemovePossibleKnotsU(tol).RemovePossibleKnotsV(tol);
}



// MAKE COMPATABLE

template<class T>
BspSurf<T> BspSurf<T>::MakeBreakCompatableU(const KnotSet& KsetU) const
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
		return 	BspSurf<T>(*cpts,norm.GetKnots(),*ktsv,ordu,ordv,numu,numv).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspSurf<T>(*cpts,norm.GetKnots(),*ktsv,ordu,ordv,numu,numv);
}

template<class T>
BspSurf<T> BspSurf<T>::MakeBreakCompatableV(const KnotSet& KsetV) const
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
		return 	BspSurf<T>(*cpts,*ktsu,norm.GetKnots(),ordu,ordv,numu,numv).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspSurf<T>(*cpts,*ktsu,norm.GetKnots(),ordu,ordv,numu,numv);
}

template<class T> 
BspSurf<T> BspSurf<T>::MakeBreakCompatable(const KnotSet& KsetU, const KnotSet& KsetV) const
{
	return MakeBreakCompatableV(KsetV).MakeBreakCompatableU(KsetU);
}



// ADDITION and SUBTRACTION	

// add two surfaces together, assuming they are not compatable
// different knot sets and different orders
template<class T>
BspSurf<T> BspSurf<T>::Add(const BspSurf<T>& b) const
{
	// make local copies of the two surfaces
	BspSurf<T> a(*this), c(b);

	// make them the same order
	if (ordu > b.GetOrdU()) c = c.ElevateU(ordu-b.GetOrdU());
	if (b.GetOrdU() > ordu) a = a.ElevateU(b.GetOrdU()-ordu);
	if (ordv > b.GetOrdV()) c = c.ElevateV(ordv-b.GetOrdV());
	if (b.GetOrdV() > ordv) a = a.ElevateV(b.GetOrdV()-ordv);

	BspSurf<T> d1 = a.MakeCompatable(c);
	BspSurf<T> d2 = c.MakeCompatable(d1);

	// now add them together
	// add the control points
	Matrix<T> temp(d1.GetNumU(),d1.GetNumV());

	for (int i=0; i<d1.GetNumU(); i++)
		for (int j=0; j<d1.GetNumV(); j++)
			temp[i][j]=d1.GetCPoints()[i][j]+d2.GetCPoints()[i][j];

	return BspSurf<T>(temp,d1.GetKnotsU(),d1.GetKnotsV(),d1.GetOrdU(),d1.GetOrdV(),d1.GetNumU(),d1.GetNumV());
}

// subtract two surfaces together, assuming they are not compatable
// different orders and different knots
template<class T>
BspSurf<T> BspSurf<T>::Subtract(const BspSurf<T>& b) const
{
	// make local copies of the two surfaces
	BspSurf<T> a(*this), c(b);

	// make them the same order
	if (ordu > b.GetOrdU()) c = c.ElevateU(ordu-b.GetOrdU());
	if (b.GetOrdU() > ordu) a = a.ElevateU(b.GetOrdU()-ordu);
	if (ordv > b.GetOrdV()) c = c.ElevateV(ordv-b.GetOrdV());
	if (b.GetOrdV() > ordv) a = a.ElevateV(b.GetOrdV()-ordv);


	BspSurf<T> d1 = a.MakeCompatable(c);
	BspSurf<T> d2 = c.MakeCompatable(d1);

	// now subtract them
	// subtract the control points
	Matrix<T> temp(d1.GetNumU(),d1.GetNumV());

	for (int i=0; i<d1.GetNumU(); i++)
		for (int j=0; j<d1.GetNumV(); j++)
			temp[i][j]=d1.GetCPoints()[i][j]-d2.GetCPoints()[i][j];

	return BspSurf<T>(temp,d1.GetKnotsU(),d1.GetKnotsV(),d1.GetOrdU(),d1.GetOrdV(),d1.GetNumU(),d1.GetNumV());

}


// DEGREE ELEVATION

// elevate the degree of the BspSurf by levu in u and levv in v
template<class T>
BspSurf<T> BspSurf<T>::Elevate(int levu, int levv) const
{
	// elevate in u then v
	return ElevateU(levu).ElevateV(levv);
}


// elevate the degree of the BspSurf by levu in u and levv in v
template<class T>
Matrix<T> BspSurf<T>::ElevateCPoints(int levu, int levv) const
{
	// elevate in u then v
	Matrix<T> mat = ElevateCPointsU(levu);

	KnotSet kset = (*ksetu).CreateKnotSetElevate(levu);

	return BspSurf<T>(mat,kset.GetKnots(),*ktsv,ordu+levu,ordv,kset.GetNum()-(ordu+levu),numv).ElevateCPointsV(levv);
}



// DERIVATIVES

template<class T>
T BspSurf<T>::operator() (int valu, int valv, double u, double v) const
{
	return Derive(valu,valv,u,v);
}


// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
BspSurf<T> BspSurf<T>::Derive(int levu, int levv) const
{
	return DeriveU(levu).DeriveV(levv);
}

// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
Matrix<T> BspSurf<T>::DeriveCPoints(int levu, int levv) const
{	
	// derive in v
	Matrix<T> ncpts = DeriveCPointsV(levv);
	KnotSet kset = (*ksetv).CreateKnotSetDeriv(levv);

	return BspSurf<T>(ncpts,*ktsu,kset.GetKnots(),ordu,ordv-levv,numu,kset.GetNum()-(ordv-levv)).DeriveCPointsU(levu);
}

// evaluate the derivative of the BezSurf of order deriv at a point
// x. Computes the derivative as a BezSurf and then evaluates this at x
template<class T>
T BspSurf<T>::Derive(int levu, int levv, double u, double v) const
{
	return (DeriveU(levu).DeriveV(levv))(u,v);
}


// INTEGRATION


// compute the indefinite integral of the BezSurf and represent
// it as a BezSurf of one higher degree
template<class T>
BspSurf<T> BspSurf<T>::Integrate() const
{
	return IntegrateU().IntegrateV();
}	

// compute the indefinite integral of the BspCurv and represent
// it as a BspCurv of one higher degree
template<class T>
BspSurf<T> BspSurf<T>::IntegrateU(int level) const
{
	BspSurf<T> b = IntegrateU();

	if (level >= 0) {
		for (int i=1; i<level; i++) b = b.IntegrateU();
		return b;
	}
	else return *this;
}	

// compute the indefinite integral of the BspCurv and represent
// it as a BspCurv of one higher degree
template<class T>
BspSurf<T> BspSurf<T>::IntegrateV(int level) const
{
	BspSurf<T> b = IntegrateV();

	if (level >= 0) {
		for (int i=1; i<level; i++) b = b.IntegrateV();
		return b;
	}
	else return *this;
}	


// integrate the BezSurf between the limits x1 and x2. Computes
// the indefinite integral as a BezSurf and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T BspSurf<T>::Integrate(double u1, double u2, double v1, double v2) const
{
	BspSurf<T> temp = Integrate();

	if (u2 < lu1 || u1 > lu2) return 0;
	if (v2 < lv1 || v1 > lv2) return 0;

	if (u1 < lu1) u1 = lu1;
	if (v1 < lv1) v1 = lv1;
	if (u2 > lu2) u2 = lu2;
	if (v2 > lv2) v2 = lv2;


	return temp(u2,v2)-temp(u1,v2)-temp(u2,v1)+temp(u1,v1);
}


// compute the indefinite integral of the BezSurf as a BezSurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> BspSurf<T>::IntegrateCPoints() const
{
	Matrix<T> ncpts = IntegrateCPointsV();
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();
	return BspSurf<T>(ncpts,*ktsu,kset.GetKnots(),ordu,ordv+1,numu,kset.GetNum()-(ordv+1)).IntegrateCPointsU();
}	



// KNOT INSERTION

// inserts a vector of knots into the BspSurf according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsU(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	if (n <= 0) return *cpts;

	int count=0;
	for (int i=0; i<n; i++) count += mult[i];

	Matrix<T> ncpts(count+numu,numv);
	// insert into v curves
	for (int j=0; j<numv; j++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).InsertKnotCPoints(t,mult,n);
		for (int i=0; i<numu+count; i++) ncpts[i][j]=temp[i];
	}
	return ncpts;
}
             

// inserts a vector of knots into the BspSurf according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotU(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// create knot set
	if (n <= 0) return *this;
	KnotSet kset = (*ksetu).CreateKnotSetInsert(t,mult,n);

	return BspSurf<T>(InsertKnotCPointsU(t,mult,n),kset.GetKnots(),*ktsv,ordu,ordv,kset.GetNum()-ordu,numv);
}


 
// inserts a vector of knots into the BspSurf according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsV(const Vector<double>& t,const Vector<int>& mult,int m) const
{
	if (m <= 0) return *cpts;
	int count=0;
	for (int i=0; i<m; i++) count+= mult[i];

	Matrix<T> ncpts(numu,numv+count);

	// insert into v curves
	for (int i=0; i<numu; i++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).InsertKnotCPoints(t,mult,m);
		for (int j=0; j<numv+count; j++) ncpts[i][j]=temp[j];
	}
	return ncpts;
}

 
// inserts a vector of knots into the BspSurf according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotV(const Vector<double>& t,const Vector<int>& mult,int m) const
{
	if (m <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(t,mult,m);

	return BspSurf<T>(InsertKnotCPointsV(t,mult,m),*ktsu, kset.GetKnots(), ordu, ordv,numu,kset.GetNum()-ordv);	
}
             

template<class T>
BspSurf<T> BspSurf<T>::InsertKnotUV(const Vector<double>& Ktsu, const Vector<int>& multu, int N, const Vector<double>& Ktsv, const Vector<int>& multv, int m) const
{
	// insert into u then v
	return InsertKnotU(Ktsu,multu,n).InsertKnotV(Ktsv,multv,m);
}

	
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsUV(const Vector<double>& Ktsu, const Vector<int>& multu, int n, const Vector<double>& Ktsv, const Vector<int>& multv, int m) const
{
	// insert into u then v
	Matrix<T> mat = InsertKnotCPointsV(Ktsv,multv,m);

	KnotSet kset = (*ksetv).CreateKnotSetInsert(Ktsv,multv,m);

	return BspSurf<T>(mat,*ktsu,kset.GetKnots(),ordu,ordv,numu,kset.GetNum()-ordv).InsertKnotCPointsU(Ktsu,multu,n);
}


// SUBDIVISION

template<class T>
BspSurf<T> BspSurf<T>::Subdivide(int levu, int levv) const
{
	return SubdivideU(levu).SubdivideV(levv);
}

// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
BspSurf<T> BspSurf<T>::Subdivide(double u1, double u2, double v1, double v2) const
{
	return SubdivideU(u1,u2).SubdivideV(v1,v2);
}

// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
Matrix<T> BspSurf<T>::SubdivideCPoints(double u1, double u2, double v1, double v2) const
{
	Matrix<T> ncpts = SubdivideV(v1,v2);
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);
	
	return BspSurf<T>(ncpts,*ktsu,kset.GetKnots(),ordu,ordv,numu,kset.GetNum()-ordv).SubdivideCPointsU(u1,u2);
}


// CONVERSION

// convert the BspSurf to a composite Bezier surf
// by converting the u and v curves
template<class T>
CompBezSurf<T> BspSurf<T>::ConvertCompBezSurf() const
{
	// create knotset object
	KnotSet kset1 = (*ksetu).CreateKnotSetCompBezCurv();
	KnotSet kset2 = (*ksetv).CreateKnotSetCompBezCurv();

	Matrix<T> ncpts(kset1.GetNum()-ordu,kset2.GetNum()-ordv);

	// convert v curves
	for (int i=0; i<numu; i++) {
		Vector<T> v =BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).ConvertCompBezCurvCPoints();
		for (int j=0; j<kset2.GetNum()-ordv; j++) ncpts[i][j]=v[j];
	}

	// convert resulting u curves
	for (int j=0; j<kset2.GetNum()-ordv; j++) {
		Vector<T> v = BspCurv<T>(ncpts.GetCol(j), *ktsu, ordu,numu).ConvertCompBezCurvCPoints();
		for (int i=0; i<kset1.GetNum()-ordu; i++) ncpts[i][j]=v[i];
	}
	// create and return the CompBezSurf
	return CompBezSurf<T>(kset1.GetNumDistinct()-1, kset2.GetNumDistinct()-1, ordu, ordv, ncpts, kset1.GetKnots(), kset2.GetKnots()); 
}


template<class T>
CompPolySurf<T> BspSurf<T>::ConvertCompPolySurf() const
{
	return ConvertCompBezSurf().ConvertCompPolySurf();
}


// ISOPARAMETRIC CURVES

// return the isoparametric curve in v as a BspCurv in u
template<class T>
BspCurv<T> BspSurf<T>::GetIsoparametricV(double v) const
{
	Vector<T> ncpts(numu);
	// evaluate the v curves at v
	for (int i=0; i<numu; i++)
		ncpts[i] = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).Eval(v);
	// return the u curve
	return BspCurv<T>(ncpts,*ktsu,ordu,numu);
}


// return the isoparametric curve in u as a BspCurv in v
template<class T>
BspCurv<T> BspSurf<T>::GetIsoparametricU(double u) const
{
	Vector<T> ncpts(numv);
	// evaluate the v curves at v
	for (int j=0; j<numv; j++)
		ncpts[j] = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).Eval(u);
	// return the u curve
	return BspCurv<T>(ncpts,*ktsv,ordv,numv);
}


// READ and WRITE

template <class T>
void BspSurf<T>::write(std::ostream& os)
{
	os << "Bspline Surface\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "number of control points in u is " << numu << "\n";
	os << "number of control points in v is " << numv << "\n";
	
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are \n";
	os << *ktsv;
	
	os << "\ncontrol points are\n";
	os << *cpts;
}

template <class T>
void BspSurf<T>::read(std::istream& is)
{
	int Ordu, Ordv;
	std::cout << "order of Bspline Surface in u and v\n";
	is >> Ordu >> Ordv;
	std::cout << "number of control points in u and v\n";
	int Numu, Numv;
	is >> Numu >> Numv;
	Matrix<T> Cpts(Numu,Numv);
	Vector<double> Ktsu(Ordu+Numu), Ktsv(Ordv+Numv);
	std::cout << "knots in u\n";
	is >> Ktsu;
	std::cout << "knots in v\n";
	is >> Ktsv;
	std::cout << "input control points\n";
	is >> Cpts;
	*this = BspSurf<T>(Cpts,Ktsu,Ktsv,Ordu,Ordv,Numu,Numv);
} 


template <class T>
void BspSurf<T>::writefile(std::ofstream& ofs)
{
	ofs << ordu << " " << ordv << " " << numu << " " << numv << "\n";
	ofs << *ktsu;
	ofs << *ktsv;
	ofs << *cpts;
}


template <class T>
void BspSurf<T>::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv;
	int Numu, Numv;
	ifs >> Ordu >> Ordv;
	ifs >> Numu >> Numv;
	Matrix<T> Cpts(Numu,Numv);
	Vector<double> Ktsu(Ordu+Numu), Ktsv(Ordv+Numv);
	ifs >> Ktsu;
	ifs >> Ktsv;
	ifs >> Cpts;
	*this = BspSurf<T>(Cpts,Ktsu,Ktsv,Ordu,Ordv,Numu,Numv);
} 



// PRIVATE FUNCTIONS

template<class T> 
BspSurf<T> BspSurf<T>::MakeCompatableU(const BspSurf<T>& b) const
{
	// assume they are both of the same degree
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetu).Normalise(b.GetKnotSetU());

	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(b.GetKnotSetU());

	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	if (kdiff1.GetNum() > 0 ) {
	}
	if (kdiff1.GetNum() > 0) 
		return BspSurf<T>(*cpts,norm.GetKnots(),*ktsv,ordu,ordv,numu,numv).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetMult(),kdiff1.GetNumDistinct());
	else return BspSurf<T>(*cpts,norm.GetKnots(),*ktsv,ordu,ordv,numu,numv);
}

template<class T> 
BspSurf<T> BspSurf<T>::MakeCompatableV(const BspSurf<T>& b) const
{
	// assume they are both of the same degree
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*ksetv).Normalise(b.GetKnotSetV());

	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(b.GetKnotSetV());

	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);

	if (kdiff1.GetNum() > 0) 
		return BspSurf<T>(*cpts,*ktsu,norm.GetKnots(),ordu,ordv,numu,numv).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetMult(),kdiff1.GetNumDistinct());
	else return BspSurf<T>(*cpts,*ktsu,norm.GetKnots(),ordu,ordv,numu,numv);
}


template<class T>
BspSurf<T> BspSurf<T>::MakeKnotSetCompatableU(const KnotSet& KsetU) const
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
		return 	BspSurf<T>(*cpts,norm.GetKnots(),*ktsv,ordu,ordv,numu,numv).InsertKnotU(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspSurf<T>(*cpts,norm.GetKnots(),*ktsv,ordu,ordv,numu,numv);
}


template<class T>
BspSurf<T> BspSurf<T>::MakeKnotSetCompatableV(const KnotSet& KsetV) const
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
		return 	BspSurf<T>(*cpts,*ktsu,norm.GetKnots(),ordu,ordv,numu,numv).InsertKnotV(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspSurf<T>(*cpts,*ktsu,norm.GetKnots(),ordu,ordv,numu,numv);
}

template<class T>
BspSurf<T> BspSurf<T>::MakeKnotSetCompatable(const KnotSet& KsetU, const KnotSet& KsetV) const
{
	return MakeKnotSetCompatableV(KsetV).MakeKnotSetCompatableU(KsetU);
}




// elevate the degree of the BspSurf by level
template<class T>
Matrix<T> BspSurf<T>::ElevateCPointsU(int level) const
{
	// not equal to num + level
	// find num of distinct knots	
	KnotSet kset = (*ksetu).CreateKnotSetElevate(level);
	int count = kset.GetNum()-(ordu+level);
	// create matrix of control points
	Matrix<T> ncpts(count,numv);

	// elevate each u curve
	for (int j=0; j<numv; j++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).ElevateCPoints(level);
		// extract control points
		for (int i=0; i<kset.GetNum()-(ordu+level); i++) ncpts[i][j]=temp[i];
	}
	// return the new BspSurf
	return ncpts;
}


// elevate the degree of the BspSurf by level
template<class T>
BspSurf<T> BspSurf<T>::ElevateU(int level) const
{
	// not equal to num + level
	// find num of distinct knots	
	KnotSet kset = (*ksetu).CreateKnotSetElevate(level);

	// return the new BspSurf
	return BspSurf<T>(ElevateCPointsU(level), kset.GetKnots(), *ktsv, ordu+level, ordv, kset.GetNum()-(ordu+level), numv);
}

// elevate the degree of the BspSurf by level
template<class T>
BspSurf<T> BspSurf<T>::ElevateV(int level) const
{
	// not equal to num + level
	// find num of distinct knots
	KnotSet kset = (*ksetv).CreateKnotSetElevate(level);
	// return the new BspSurf
	return BspSurf<T>(ElevateCPointsV(level), *ktsu, kset.GetKnots(), ordu, ordv+level, numu, kset.GetNum()-(ordv+level));
}



// elevate the degree of the BspSurf by level
template<class T>
Matrix<T> BspSurf<T>::ElevateCPointsV(int level) const
{
	// find num of distinct knots	
	KnotSet kset = (*ksetv).CreateKnotSetElevate(level);
	int count = kset.GetNum()-(ordv+level);
	
	// create the matrix of control points
	Matrix<T> ncpts(numu,count);

	// elevate each v curve
	for (int i=0; i<numu; i++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).ElevateCPoints(level);
		// extract control points
		for (int j=0; j<kset.GetNum()-(ordv+level); j++) ncpts[i][j]=temp[j];
	}
	// return the new BspSurf
	return ncpts;
}

// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
BspSurf<T> BspSurf<T>::DeriveU(int level) const
{   
	if (level <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetu).CreateKnotSetDeriv(level);

	// return the new Bsurf
	return BspSurf<T>(DeriveCPointsU(level), kset.GetKnots(), *ktsv, ordu-level, ordv,kset.GetNum()-(ordu-level),numv);
}


// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
Matrix<T> BspSurf<T>::DeriveCPointsU(int level) const
{   
	if (level <= 0) return *cpts;
	KnotSet kset = KnotSet(*ktsu,ordu,ordu+numu).CreateKnotSetDeriv(level);
	Matrix<T> ncpts(kset.GetNum()-(ordu-level),numv);

	// returns control points only
	for (int j=0; j<numv; j++) {
			Vector<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).DeriveCPoints(level);
			// extract control points
			for (int i=0; i<kset.GetNum()-(ordu-level); i++) ncpts[i][j]=temp[i];
	}
	// return the new control points
	return ncpts;
}



// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
BspSurf<T> BspSurf<T>::DeriveV(int level) const
{
	if (level <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetDeriv(level);

	// return the new Bsurf
	return BspSurf<T>(DeriveCPointsV(level), *ktsu, kset.GetKnots(), ordu, ordv-level,numu, kset.GetNum()-(ordv-level));
}


// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
Matrix<T> BspSurf<T>::DeriveCPointsV(int level) const
{
	if (level <= 0) return *cpts;
	KnotSet kset = (*ksetv).CreateKnotSetDeriv(level);
	Matrix<T> ncpts(numu,kset.GetNum()-(ordv-level));

	// derive in v
	for (int i=0; i<numu; i++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).DeriveCPoints(level);
		// extract control points
		for (int j=0; j<kset.GetNum()-(ordv-level); j++) ncpts[i][j]=temp[j];
	}

	// return the new control points
	return ncpts;
}


// evaluate the derivative of the BezSurf of order deriv at a point
// x. Computes the derivative as a BezSurf and then evaluates this at x
template<class T>
T BspSurf<T>::DeriveU(int level, double u, double v) const
{
	return DeriveU(level)(u,v);
}

// evaluate the derivative of the BezSurf of order deriv at a point
// x. Computes the derivative as a BezSurf and then evaluates this at x
template<class T>
T BspSurf<T>::DeriveV(int level, double u, double v) const
{
	return DeriveV(level)(u,v);
}

// compute the indefinite integral of the BezSurf and represent
// it as a BezSurf of one higher degree
template<class T>
BspSurf<T> BspSurf<T>::IntegrateU() const
{
	// create knot sets
	KnotSet kset = (*ksetu).CreateKnotSetIntegrate();
	
	return BspSurf<T>(IntegrateCPointsU(),kset.GetKnots(),*ktsv,ordu+1,ordv,kset.GetNum()-(ordu+1),numv);
}	


// compute the indefinite integral of the BezSurf and represent
// it as a BezSurf of one higher degree
template<class T>
BspSurf<T> BspSurf<T>::IntegrateV() const
{
	// create knot sets
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();

	return BspSurf<T>(IntegrateCPointsV(),*ktsu,kset.GetKnots(),ordu,ordv+1,numu,kset.GetNum()-(ordv+1));
}	


// integrate the BezSurf between the limits x1 and x2. Computes
// the indefinite integral as a BezSurf and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T BspSurf<T>::IntegrateUV1(double u1, double u2, double v1, double v2) const
{
	Vector<T> v(numv);

	// should be integrate CPoints
	for (int j=0; j<numv; j++)
		v[j] = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).Integrate(u1,u2);
	return BspCurv<T>(v,*ktsv,ordv,numv).Integrate(v1,v2);
}


// compute the indefinite integral of the BezSurf as a BezSurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> BspSurf<T>::IntegrateCPointsU() const
{
	KnotSet kset = (*ksetu).CreateKnotSetIntegrate();
	Matrix<T> ncpts(kset.GetNum()-(ordu+1),numv);

	for (int j=0; j<numv; j++) {
		Vector<T> v = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).IntegrateCPoints();			
		for (int i=0; i<kset.GetNum()-(ordu+1); i++) ncpts[i][j]=v[i];
	}
	
	return ncpts;
}	

// compute the indefinite integral of the BezSurf as a BezSurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> BspSurf<T>::IntegrateCPointsV() const
{
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();
	Matrix<T> ncpts(numu,kset.GetNum()-(ordv+1));

	// use create cpoints integrate 
	for (int i=0; i<numu; i++)  {
		Vector<T> v = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).IntegrateCPoints();			
		for (int j=0; j<kset.GetNum()-(ordv+1); j++) ncpts[i][j]=v[j];
	}

	return ncpts;
}

// inserts a vector of knots into the BspSurf and returns the new
// BspSurf
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotU(const Vector<double>& Kts, int n) const
{
	// find size of array
	if (n <= 0) return *this;
	KnotSet kset = (*ksetu).CreateKnotSetInsert(Kts,n);

	return BspSurf<T>(InsertKnotCPointsU(Kts,n),kset.GetKnots(), *ktsv, ordu, ordv, kset.GetNum()-ordu, numv);
}
    


// inserts a vector of knots into the BspSurf and returns the new
// BspSurf
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsU(const Vector<double>& Kts, int n) const
{
	if (n <= 0) return *cpts;
	Matrix<T> ncpts(numu+n,numv);
	// insert knots into u curves
	for (int j=0; j<numv; j++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).InsertKnotCPoints(Kts,n);
		for (int i=0; i<numu+n; i++) ncpts[i][j]=temp[i];
	}
	return ncpts;
}

         
// inserts a knot x into the BspSurf  and returns the new BspSurf
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotU(double x) const
{       
	// find size of array
	KnotSet kset = (*ksetu).CreateKnotSetInsert(x);
	
	return BspSurf<T>(InsertKnotCPointsU(x), kset.GetKnots(), *ktsv, ordu, ordv, kset.GetNum()-ordu, numv);
}                              
      

// inserts a knot x into the BspSurf  and returns the new BspSurf
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsU(double x) const
{       
	Matrix<T> ncpts(numu+1,numv);
	// insert knots into v curves
	for (int j=0; j<numv; j++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).InsertKnotCPoints(x);
		for (int i=0; i<numu+1; i++) ncpts[i][j]=temp[i];
	}
	return ncpts;
}                              
      



// inserts a knot x 'level' times into the BspSurf and returns
// the new BspSurf
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotU(double x, int level) const
{
	// create knot set 
	if (level <= 0) return *this;
	KnotSet kset = (*ksetu).CreateKnotSetInsert(x,level);
	
	return BspSurf<T>(InsertKnotCPointsU(x,level),kset.GetKnots(), *ktsv, ordu, ordv, kset.GetNum()-ordu, numv);
}



// inserts a knot x 'level' times into the BspSurf and returns
// the new BspSurf
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsU(double x, int level) const
{
	if (level <= 0) return *cpts;
	Matrix<T> ncpts(numu+level,numv);
	// insert into u curves
	for (int j=0; j<numv; j++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).InsertKnotCPoints(x,level);
		for (int i=0; i<numu+level; i++) ncpts[i][j]=temp[i];
	}
	return ncpts;
}

 
// inserts a vector of knots into the BspSurf and returns the new
// BspSurf
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotV(const Vector<double>& Kts, int N) const
{
	if (N <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(Kts,N);

	return BspSurf<T>(InsertKnotCPointsV(Kts,N),*ktsu, kset.GetKnots(), ordu, ordv,numu,kset.GetNum()-ordv);
}
    

// inserts a vector of knots into the BspSurf and returns the new
// BspSurf
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsV(const Vector<double>& Kts, int m) const
{
	if (m <= 0) return *cpts;
	Matrix<T> ncpts(numu,numv+m);
	// insert into v curves
	for (int i=0; i<numu; i++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).InsertKnotCPoints(Kts,m);
		for (int j=0; j<numv+m; j++) ncpts[i][j]=temp[j];
	}
	return ncpts;
}

         
// inserts a knot x into the BspSurf  and returns the new BspSurf
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsV(double x) const
{       
	Matrix<T> ncpts(numu,numv+1);

	// insert into v curves
	for (int i=0; i<numu; i++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).InsertKnotCPoints(x);
		for (int j=0; j<numv+1; j++) ncpts[i][j]=temp[j];
	}
	return ncpts;
}                              
      
         
// inserts a knot x into the BspSurf  and returns the new BspSurf
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotV(double x) const
{       
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(x);

	return BspSurf<T>(InsertKnotCPointsV(x),*ktsu, kset.GetKnots(), ordu, ordv,numu,kset.GetNum()-ordv);
}                              


// inserts a knot x 'level' times into the BspSurf and returns
// the new BspSurf
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsV(double x, int level) const
{
	if (level <= 0) return *cpts;
	Matrix<T> ncpts(numu,numv+level);
	// insert into v curves
	for (int i=0; i<numu; i++) {
		Vector<T> temp = BspCurv<T>((*cpts).Row(i),*ktsv,ordv,numv).InsertKnotCPoints(x,level);
		for (int j=0; j<numv+level; j++) ncpts[i][j]=temp[j];
	}
	return ncpts;
}


// inserts a knot x 'level' times into the BspSurf and returns
// the new BspSurf
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotV(double x, int level) const
{
	if (level <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(x,level);

	return BspSurf<T>(InsertKnotCPointsV(x,level),*ktsu, kset.GetKnots(), ordu, ordv,numu,kset.GetNum()-ordv);
}

template<class T>
BspSurf<T> BspSurf<T>::InsertKnotUV(const Vector<double>& Ktsu, int n, const Vector<double>& Ktsv, int m) const
{
	// insert into u then v
	return InsertKnotU(Ktsu,n).InsertKnotV(Ktsv,m);
}
	
             
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotUV(double u, double v) const
{
	// insert into u then v
	return InsertKnotU(u).InsertKnotV(v);
}
	
template<class T>
BspSurf<T> BspSurf<T>::InsertKnotUV(double u, double v, int levu, int levv) const
{
	// insert into u then v
	return InsertKnotU(u,levu).InsertKnotV(v,levv);
}

             
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsUV(double u, double v) const
{
	// insert into u then v
	Matrix<T> mat = InsertKnotCPointsV(v);

	KnotSet kset = (*ksetv).CreateKnotSetInsert(v);

	return BspSurf<T>(mat,*ktsu,kset.GetKnots(),ordu,ordv,numu,kset.GetNum()-ordv).InsertKnotCPointsU(u);
}
	
template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsUV(double u, double v, int levu, int levv) const
{
	// insert into u then v
	Matrix<T> mat = InsertKnotCPointsV(v,levv);

	KnotSet kset = (*ksetv).CreateKnotSetInsert(v,levv);

	return BspSurf<T>(mat,*ktsu,kset.GetKnots(),ordu,ordv,numu,kset.GetNum()-ordv).InsertKnotCPointsU(u,levu);
}

template<class T>
Matrix<T> BspSurf<T>::InsertKnotCPointsUV(const Vector<double>& Ktsu, int n, const Vector<int>& Ktsv, int m) const
{
	// insert into u then v
	Matrix<T> mat = InsertKnotCPointsV(Ktsv,m);

	KnotSet kset = (*ksetv).CreateKnotSetInsert(Ktsv,m);

	return BspSurf<T>(mat,*ktsu,kset.GetKnots(),ordu,ordv,numu,kset.GetNum()-ordv).InsertKnotCPointsU(Ktsu,n);
}

template<class T>
BspSurf<T> BspSurf<T>::SubdivideU(int level) const
{
	if (level <= 0) return *this;
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(level);

	return BspSurf<T>(SubdivideCPointsU(level),kset.GetKnots(),*ktsv,ordu,ordv,kset.GetNum()-ordu,numv);
	
}

template<class T>
Matrix<T> BspSurf<T>::SubdivideCPointsU(int level) const
{
	if (level <= 0) return *cpts;
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(level);
	int count = kset.GetNum()-ordu;
	Matrix<T> ncpts(count,numv);
	// insert into u curves
	for (int j=0; j<numv; j++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).SubdivideCPoints(level);
		for (int i=0; i<count; i++) ncpts[i][j]=temp[i];
	}
	return ncpts;

}


template<class T>
BspSurf<T> BspSurf<T>::SubdivideV(int level) const
{
	if (level <= 0) return *this;
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(level);

	return BspSurf<T>(SubdivideCPointsV(level),*ktsu,kset.GetKnots(),ordu,ordv,numu,kset.GetNum()-ordv);
}

template<class T>
Matrix<T> BspSurf<T>::SubdivideCPointsV(int level) const
{
	if (level <= 0) return cpts;
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(level);
	int count=kset.GetNum()-ordv;
	
	Matrix<T> ncpts(numu,count);

	// insert into v curves
	for (int i=0; i<numu; i++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).SubdivideCPoints(level);
		for (int j=0; j<count; j++) ncpts[i][j]=temp[j];
	}
	return ncpts;
}


// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
BspSurf<T> BspSurf<T>::SubdivideU(double u1, double u2) const
{
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(u1,u2);
	
	int count = kset.GetNum() - ordu;
	// return the new Bsurface
	return BspSurf<T>(SubdivideCPointsU(u1,u2), kset.GetKnots(), *ktsv, ordu, ordv, count, numv);
}

// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
Matrix<T> BspSurf<T>::SubdivideCPointsU(double u1, double u2) const
{
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(u1,u2);
	int count = kset.GetNum()-ordu;
	Matrix<T> ncpts(count,numv);

	for (int j=0; j<numv; j++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetCol(j),*ktsu,ordu,numu).SubdivideCPoints(u1,u2);
		// extract control points
		for (int i=0; i<count; i++) ncpts[i][j]=temp[i];
	}
	// return the new control points
	return ncpts;
}


// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
BspSurf<T> BspSurf<T>::SubdivideV(double v1, double v2) const
{
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);

	// return the new Bsurface
	return BspSurf<T>(SubdivideCPointsV(v1,v2), *ktsu, kset.GetKnots(), ordu, ordv,numu,kset.GetNum()-ordv);

}


// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
Matrix<T> BspSurf<T>::SubdivideCPointsV(double v1, double v2) const
{
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);
	int count = kset.GetNum()-ordv;
	Matrix<T> ncpts(numu,count);

	for (int i=0; i<numu; i++) {
		Vector<T> temp = BspCurv<T>((*cpts).GetRow(i),*ktsv,ordv,numv).SubdivideCPoints(v1,v2);
		// extract control points
		for (int j=0; j<count; j++) ncpts[i][j]=temp[j];
	}
	// return the new Bsurface
	return ncpts;

}



#endif
 
