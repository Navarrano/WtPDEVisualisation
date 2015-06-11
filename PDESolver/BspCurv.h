#ifndef BSPCURV
#define BSPCURV

#include "matrix.h"
#include "cpoints.h"
#include "compbezcurv.h"


class BspCurvBasisFunc : public Curve<double> {

	// data
	int ord;	 // order of basis function
	Ptr<Vector<double> > kts;		// knots 
	Ptr<BspCurv<double> > b;	// BspCurv representation of basis function

	// private functions
	BspCurv<double> CreateBspCurv() const;	// creates the BspCurv
	virtual ObjectID Identity() const { return std::string("class BspCurvBasisFunc"); }
public:
	// constructors
	BspCurvBasisFunc();
	BspCurvBasisFunc(const Vector<double>& Kts, int Ord);
	
	// access functions
	int ComputeDim() const;
	BspCurv<double> GetBspCurv() const;
	Vector<double> GetKnots() const;
	int GetOrd() const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;

	// evaluators
	virtual double operator()(double x) const;
	virtual double operator()(int, double x) const;
	double Eval(double x) const;
	virtual double Derive(int, double) const;

	// integration
	double CreateIntegral(double x1, double x2) const;
	
	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

class BspCurvBasisFuncSet : public TextObject, public FTextObject {

	int ord;
	int num;
	Ptr<Vector<double> > kts;
	Ptr<Vector<BspCurvBasisFunc> > b;
	virtual ObjectID Identity() const { return std::string("class BspCurvBasisFuncSet"); }
public:
	// constructors
	BspCurvBasisFuncSet();
	BspCurvBasisFuncSet(const Vector<double>& Kts, int Ord, int Num);
	BspCurvBasisFuncSet(const KnotSet& kset);
	
	// evaluators
	Vector<double> operator()(double x) const;
	Vector<double> operator()(int, double x) const;
	Vector<double> Derive(int, double x) const;
	Vector<double> Eval(double x) const;
	Vector<double> EvalNonZero(double x) const;

	// get functions
	Vector<double> GetKnots() const;
	int GetOrd() const;
	BspCurvBasisFunc GetBasisFunc(int i) const;


	// variational stuff
	Vector<double> CreateIntegralNew(const BspCurv<double>& c, double x1, double x2) const;
	Vector<double> CreateIntegral(double x1, double x2) const;
	Vector<double> CreateIntegralIntegral(double x1, double x2) const;
	Vector<double> CreateIntegralProduct(const Curve<double>& p, double x1, double x2) const;
	Vector<double> CreateIntegralIntegralIntegral(double x1, double x2) const;
	Vector<double> CreateIntegralIntegralProduct(const Curve<double>& p, double x1, double x2) const;
	Vector<double> ComputeVecBasis(double x1) const;
	Matrix<double> ComputeNMatrix(double x) const;
	Matrix<double> CreateMatrixIntegral(int lev, double x1, double x2) const;
	Matrix<double> CreateMatrixIntegral(int lev) const;
	Matrix<double> CreateMatrixMinimisation(int lev, double x1, double x2) const;
	Matrix<double> CreateMatrixMinimisation(int lev) const;
	Matrix<double> CreateMatrixIntegral(int lev1, int lev2, double x1, double x2) const;
	Matrix<double> CreateMatrixIntegral(int lev1, int lev2) const;
	Matrix<double> CreateMatrixMinimisation(int lev1, int lev2, double x1, double x2) const;
	Matrix<double> CreateMatrixMinimisation(int lev1, int lev2) const;
	void ComputeIntersection(Matrix<CompBezCurv<double> >& m1, Matrix<CompBezCurv<double> >& m2) const;
	Matrix<double> CreateMatrixKnotAverages() const;
	Matrix<double> CreateMatrixKnotAverages(Vector<double>& vec) const;
	Matrix<double> ComputeBasisMatrix() const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};


// class for a Bspline curve 
template<class T>
class BspCurv : public Curve<T> {
private:
	// data
	int ord;
	int num;
	Ptr<Vector<T> > cpts;
	Ptr<Vector<double> > kts;
	Ptr<KnotSet> kset;
	
	// private functions
	
	BspCurv<T> Elevate2(int level) const;
	BspCurv<T> OldProduct(const BspCurv<T>& c) const;
	Vector<T> ConvertCompBezCurvCPoints2() const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class BspCurv<double>");
		else {
			std::string s(typeid(T).name()), s1("class BspCurv<");
			return s1.append(s,0,s.size())+">";
		}
		/*if (std::string(typeid(T).name()) == std::string("class double"))
			return std::string("class BspCurvDouble"); 
		else if (std::string(typeid(T).name()) == std::string("class Point1D"))
			return std::string("class FBspCurv");
		else if (std::string(typeid(T).name()) == std::string("class Point2D"))
			return std::string("class BspCurv2D");
		else return std::string("class BspCurv3D"); */
	}
public:
	// constructors
	BspCurv();
	BspCurv(const Vector<T>& Cpts, const Vector<double>& Kts, int Ord, int Num);
	BspCurv(const Vector<T>& Cpts, int Ord, int Num);
	BspCurv(CPoints<T>& CPt, int Num);
	BspCurv(const PolyCurv<T>& p, const KnotSet& kset);
	BspCurv(const PolyCurv<T>& p, int Ord, const KnotSet& Kset);

	// access functions
	int GetOrd() const;
	int GetNum() const;
	Vector<T> GetCPoints() const;
	Vector<double> GetKnots() const;
	KnotSet GetKnotSet() const;
	double fabsnew(T t) const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
	Vector<double> GetLimits() const;
	BezCurv<T> GetSegment(int i) const;	
	BspCurvBasisFuncSet GetBasisFuncSet() const;
	BspCurvBasisFunc GetBasisFunc(int i) const;

	// evaluators
	virtual T operator()(double x) const;
	T Eval(double x) const;
	T Eval1(double u) const;
	Vector<T> EvalDeriv(int levu, double u) const;
	Vector<T> Eval2(double u) const;

	Vector<T> ComputePoints(int n) const;
	Matrix<T> ComputeTable(int n) const;
	
	// derivatives
	BspCurv<T> Derive(int level) const; 
	virtual T Derive(int level,double x) const;
	virtual T operator()(int, double x) const;
	Vector<T> DeriveCPoints(int deriv) const;
	T Derive1(int level, double x) const;
	Vector<double> Derive1(double x) const;
	Vector<double> Derive2(double x) const; 

	// make comptability with curve
	template<class T1>
	BspCurv<T> MakeCompatable(const BspCurv<T1>& b) const
	{
		// assume they are both of the same degree
		// create knot set for current object
	
		// create normalised knot set for current object relative to b
		KnotSet norm = (*kset).Normalise(b.GetKnotSet());
	
		// add knots to current object in knot set for b not present 
		// in knot set for Current object
		KnotSet add1 = norm.KnotUnion(b.GetKnotSet());
	
		// find differences in knots
		KnotSet kdiff1 = add1.KnotDifference(norm);

		if (kdiff1.GetNum() > 0) {
			// insert knots into current object to get compatibility with b
			return 	BspCurv<T>(*cpts,norm.GetKnots(),ord,num).InsertKnot(kdiff1.GetDistinctKnots(),kdiff1.GetMult(),kdiff1.GetNumDistinct());
		} else return BspCurv<T>(*cpts,norm.GetKnots(),ord,num);
	}


	BspCurv<T> MakeBreakCompatable(const KnotSet& Kset) const;
	
	template<class T1>
	BspCurv<T> MakeBreakCompatable(const BspCurv<T1>& b) const
	{
		// assume they are both of the same degree
		// create knot set for current object
	
		// create normalised knot set for current object relative to b
		KnotSet norm = (*kset).Normalise(b.GetKnotSet());
	
		// add knots to current object in knot set for b not present 
		// in knot set for Current object
		KnotSet add1 = norm.KnotUnion(b.GetKnotSet());
	
		// find differences in knots
		KnotSet kdiff1 = add1.KnotDifference(norm);
	
		if (kdiff1.GetNum() > 0) {
			// insert knots into current object to get compatibility with b
			return 	BspCurv<T>(*cpts,norm.GetKnots(),ord,num).InsertKnot(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
		} else return BspCurv<T>(*cpts,norm.GetKnots(),ord,num);
	
		/*if (kdiff1.GetNumDistinct() > 2) {
			Vector<double> nkts(kdiff1.GetNumDistinct()-2);
			Vector<int> nmult(kdiff1.GetNumDistinct()-2,1);
			for (int i=1; i<kdiff1.GetNumDistinct()-1; i++) nkts[i-1] = kdiff1.GetDistinctKnots()[i];
			// insert knots into current object to get compatibility with b
			return 	BspCurv<T>(*cpts,norm.GetKnots(),ord,num).InsertKnot(nkts,nmult,kdiff1.GetNumDistinct()-2);
		} else return BspCurv<T>(*cpts,norm.GetKnots(),ord,num);*/
	}

	BspCurv<T> MakeKnotSetCompatable(const KnotSet& kset) const;

	// knot removal
	BspCurv<T> KnotRemoval(double knot) const;
	bool IsKnotRemovable(double knot, double tol=cpntTol) const;
	BspCurv<T> RemovePossibleKnots(double tol=cpntTol) const;

	// degree elevation
	BspCurv<T> Elevate(int level) const;
	Vector<T> ElevateCPoints(int level) const;

	// integration
	BspCurv<T> IntegrateProduct(const BspCurv<T>& b) const;
	BspCurv<T> Integrate() const;
	BspCurv<T> Integrate(int level) const;
	T Integrate(double x1, double x2) const;
	Vector<T> IntegrateCPoints() const;
	Matrix<double> CreateMatrixIntegral(int deriv, double x1, double x2) const;

	// product
	template<class T1>
	BspCurv<T> Product(const BspCurv<T1>& b) const
	{
		BspCurv<T> d(*this);
		BspCurv<T1> e(b);

		// test to see whether knot sets are the same
		if (!(*kset).IsSameDistinctKnotSet(b.GetKnotSet())) {
			// needs to add distinct knots from b not in current object
			// into current object (independent of order)
			d = MakeBreakCompatable(b);
			e = b.MakeBreakCompatable(d);
		}
	
		// convert to composite Bezier 
		CompBezCurv<T> b3 = d.ConvertCompBezCurv();
		CompBezCurv<T1> b4 = e.ConvertCompBezCurv();
		CompBezCurv<T> b5 = b3.Product(b4);

		// form product knot set
		KnotSet kset3 = (*kset).CreateKnotSetProduct(b.GetKnotSet());

		// form knot set for product based on Bezier form
		KnotSet kset4(b5.GetKnots(),b5.GetOrd(),(b5.GetNum()+1)*b5.GetOrd());
		// form difference knot set
		KnotSet kdiff3 = kset4.KnotDifference(kset3);
		// remove knots according to difference knot set
		BspCurv<T> b6 = b5.ConvertBspCurv();
		for (int i=0; i<kdiff3.GetNumDistinct(); i++)
			for (int j=0; j<kdiff3.GetMult()[i]; j++) 
				b6 = b6.KnotRemoval(kdiff3.GetDistinctKnots()[i]);
	
		return b6;
	}
	template<class T1>
	BspCurv<T> ProductWithRemoval(const BspCurv<T1>& c) const
	{
		// needs to add distinct knots from b not in current object
		// into current object (independent of order)
		BspCurv<T> d = MakeBreakCompatable(b);
		BspCurv<T> e = b.MakeBreakCompatable(d);
		// convert to composite Bezier 
		CompBezCurv<T> b3 = d.ConvertCompBezCurv();
		CompBezCurv<T> b4 = e.ConvertCompBezCurv();
		CompBezCurv<T> b5 = b3.Product(b4);

		return b5.ConvertBspCurv().RemovePossibleKnots();
	}
	template<class T1>
	BspCurv<T> ProductWithoutRemoval(const BspCurv<T1>& b) const
	{
		// needs to add distinct knots from b not in current object
		// into current object (independent of order)
		BspCurv<T> d = MakeBreakCompatable(b);
		BspCurv<T> e = b.MakeBreakCompatable(d);
		// convert to composite Bezier 
		CompBezCurv<T> b3 = d.ConvertCompBezCurv();
		CompBezCurv<T> b4 = e.ConvertCompBezCurv();
		CompBezCurv<T> b5 = b3.Product(b4);

		return b5.ConvertBspCurv();
	}

	// subdivision
	BspCurv<T> Subdivide(double x1, double x2) const;
	BspCurv<T> Subdivide(int level) const;
	Vector<T> SubdivideCPoints(int level) const;
	Vector<T> SubdivideCPoints(double x1, double x2) const;

	// addition and subtraction
	BspCurv<T> Add(const BspCurv<T>& b) const;
	BspCurv<T> Subtract(const BspCurv<T>& b) const;

	// derivative matrix
	Matrix<T> Dmatrix(int deriv) const;

	// knot insertion
	Vector<T> InsertKnotCPoints(const Vector<double>& t,const Vector<int>& mult,int n) const;
	BspCurv<T> InsertKnot(const Vector<double>& Kts, const Vector<int>& mult, int N) const;
	BspCurv<T> InsertKnot(double x) const;
	Vector<T> InsertKnotCPoints(double x) const;
	Vector<T> InsertKnotCPoints(double x, int level) const;
	Vector<T> InsertKnotCPoints(const Vector<double>& Kts, int n) const;
	BspCurv<T> InsertKnot(double x, int level) const;
	BspCurv<T> InsertKnot(const Vector<double>& Kts, int n) const;

	Vector<T> ComputeKnotSetAverValues() const;
	BspCurv<T> CreateApproximation(const FnCurvObject<T>& curve) const;

	// conversion
	CompBezCurv<T> ConvertCompBezCurv() const;
	CompPolyCurv<T> ConvertCompPolyCurv() const;
	Vector<T> ConvertCompBezCurvCPoints() const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};


// CONSTRUCTORS

// constructor builds a BspCurv from a Vector of control points, knots
// an order and number of control points
template<class T>
BspCurv<T>::BspCurv(const Vector<T>& Cpts, const Vector<double>& Kts, int Ord, int Num) :
			ord(Ord), num(Num), cpts(new Vector<T>(Cpts)), kts(new Vector<double>(Kts)),
			kset(new KnotSet(*kts,ord,ord+num))
{
	// check to see that end knots have multiplicity ord
	// if not, curve is converted to this type using
	// de Boor triangular table evaluation algorithm
	if ((*kset).CheckKnotSet()) {

		// create CPoints object
		CPoints<T> cpt(*cpts, *kset,num);
		// update control points according to multipliicty ord at both ends
		*cpts = cpt.UpdateCPoints().GetCPoints();
		*kset = cpt.GetKnotSet();
		// assign new knot set with multipliicty ord at both ends
		*kts = (*kset).GetKnots();
	}
}


// constructor, builds a BspCurv from a Vector object 
//of control points an ord and num
// creates the knot set using the CreateKnotSetCPoints
// method in KnotSet class (uses method from Farin's book)
template<class T>
BspCurv<T>::BspCurv(const Vector<T>& Cpts, int Ord, int Num) : ord(Ord), 
	num(Num), cpts(new Vector<T>(Cpts))
{
	// initialise kset
	kts = new Vector<double>(Math::CreateKnots(num,ord));
	kset = new KnotSet(*kts,ord,ord+num);
}

// constructor, builds a BspCurv from a CPoints object 
// assume that the CPoint object has the KnotSet member 
// initialised (not necessarily true!!)
template<class T>
BspCurv<T>::BspCurv(CPoints<T>& Cpt, int Num) : ord(Cpt.GetKnotSet().GetOrd()), 
	num(Num), cpts(new Vector<T>(Num)), kts(new Vector<double>(Cpt.GetKnotSet().GetNum())), 
	kset(new KnotSet(Cpt.GetKnotSet()))
{
	for (int i=0; i<Num; i++) (*cpts)[i] = Cpt.GetCPoints()[i];
	for (int i=0; i<Cpt.GetKnotSet().GetNum(); i++) (*kts)[i] = Cpt.GetKnotSet().GetKnots()[i];
	// check knot set to see that it has multiplicity ord at both ends
	// update knot set and control points if not

	if (Cpt.GetKnotSet().CheckKnotSet()) {
		Cpt.UpdateCPoints();
		*cpts = Cpt.GetCPoints();
		*kts = Cpt.GetKnotSet().GetKnots();
	}
	ord = Cpt.GetKnotSet().GetOrd();
}


template<class T>
BspCurv<T>::BspCurv(const PolyCurv<T>& p, const KnotSet& Kset) : ord(Kset.GetOrd()), num(Kset.GetNum()-Kset.GetOrd()), 
					kts(new Vector<double>(Kset.GetKnots())), kset(new KnotSet(Kset))
{
	if (Kset.GetOrd() >= p.GetOrd())
		cpts = new Vector<T>(p.Elevate(Kset.GetOrd()-p.GetOrd()).ConvertBspCurv().MakeKnotSetCompatable(Kset).GetCPoints());
	//else cpts = new Vector<T>(ApproximationCurv<T>(Kset.GetOrd()-1,Kset.GetDistinctKnots(),2,100,FnCurvObject<T>(new PolyCurv<T>(p)),Kset.GetDistinctKnots().GetNum()-1,-1).ComputeApprox().GetCPoints());
	else {
		Matrix<double> mat = BspCurvBasisFuncSet(Kset.GetKnots(),Kset.GetOrd(),Kset.GetNum()).ComputeBasisMatrix();
		Vector<double> vec = Kset.ComputeKnotSetAver();
		Vector<T> vec1(Kset.GetNum()-Kset.GetOrd());
		for (int i=0; i<vec.GetNum(); i++) vec1[i] = p(vec[i]);
		cpts = new Vector<T>(Math::Solve(mat,vec1));
	}
}

template<class T>
BspCurv<T>::BspCurv(const PolyCurv<T>& p, int Ord, const KnotSet& Kset) 
{
	*this = p.ConvertBspCurv().MakeBreakCompatable(Kset);
}




// default constructor
template<class T>
BspCurv<T>::BspCurv() : ord(0), num(0), cpts(), kts(), kset() { }


// ACCESS FUNCTIONS

// get the order of the BspCurv
template<class T>
inline int BspCurv<T>::GetOrd() const { return ord; }


// get the number of control points
template<class T>
inline int BspCurv<T>::GetNum() const { return num; }


// get the control points
template<class T>
inline Vector<T> BspCurv<T>::GetCPoints() const { return *cpts; }


// get the knot vector
template<class T>
inline Vector<double> BspCurv<T>::GetKnots() const { return *kts; }

// get the knot vector
template<class T>
inline KnotSet BspCurv<T>::GetKnotSet() const { return *kset; }



template<class T>
inline double BspCurv<T>::GetLeftLimit() const
{
	return (*kts)[ord-1];
}

template<class T>
inline double BspCurv<T>::GetRightLimit() const
{
	return (*kts)[num];
}

template<class T>
inline Vector<double> BspCurv<T>::GetLimits() const
{
	return (*kset).GetDistinctKnots();
}

// get the ith sgement of the BspCurv and return it as a BezCurv
// assume that i starts at 1 (for first segment)
template<class T>	
BezCurv<T> BspCurv<T>::GetSegment(int i) const
{
	// subdivide the curve between i-1 and i
	BspCurv<T> b = Subdivide((*kset).GetDistinctKnots()[i-1],(*kset).GetDistinctKnots()[i]);
	// create and return the BezCurv
	return BezCurv<T>(b.GetCPoints(), b.GetKnots(), ord);
}


template<class T>
BspCurvBasisFuncSet BspCurv<T>::GetBasisFuncSet() const
{
	return BspCurvBasisFuncSet(*kts,ord,ord+num);
}

template<class T>
BspCurvBasisFunc BspCurv<T>::GetBasisFunc(int i) const
{
	Vector<double> v(ord+1);

	for (int j=0; j<=ord; j++) v[j] = (*kts)[i-1+j]; 
	
	return BspCurvBasisFunc(v,ord);
}

// EVALUATORS

// evaluate the BspCurv at the point x using de Boor algorithm
template<class T>
T BspCurv<T>::operator()(double x) const
{
	if (x < (*kts)[ord-1] || x > (*kts)[num]) return 0.0;

	Vector<T> temp(ord);
	double lbd;

	// find the interval for x 
	int ind = (*kset).Find_index(x); 

	//compute the point according to de Boor algorithm 
	// extract relevant control points
	for (int j=0; j<ord; j++) temp[j]=(*cpts)[ind-ord+j+1];

	// compute point on Curve according to triangular table
	for (int i=0; i<ord-1; i++) {
		for (int j=0; j<ord-i-1; j++) {
			lbd = (x - (*kts)[ind-ord+i+j+2])/((*kts)[ind+j+1]-(*kts)[ind-ord+i+j+2]);
			temp[j] = lbd*temp[j+1]+(1.0-lbd)*temp[j];
		
			
		}
	}
	return temp[0];
}

template<class T>
T BspCurv<T>::operator()(int lev, double x) const
{
	return Derive(lev)(x);
}

// evaluate the BspCurv using matrix method
// identical code to that of overloaded operator
template<class T>
T BspCurv<T>::Eval(double x) const
{
//	if (x < (*kts)[ord-1] || x > (*kts)[num]) return 0.0;
	Vector<double> v = (*kset).CreateVectorInterp(x);

	T sum=0.0;
	for (int i=0; i<num; i++) sum = sum+v[i]*(*cpts)[i];
	return sum;
	//return (*this)(x);

}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
T BspCurv<T>::Eval1(double u) const
{
	if (u < (*kts)[ord-1] || u > (*kts)[num]) return 0.0;
	Vector<double> uvec = kset->CreateVectorInterp(u);

	return Math::mult0(uvec,*cpts);
}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Vector<T> BspCurv<T>::Eval2(double u) const
{
	if (u < (*kts)[ord-1] || u > (*kts)[num]) return 0.0;
	Matrix<double> mu(kset->CreateVectorInterp(u));

	return Math::mult3(mu,*cpts);
}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Vector<T> BspCurv<T>::EvalDeriv(int levu, double u) const
{
	Matrix<double> mu(kset->CreateVectorInterpDeriv(levu,u));

	return Math::mult3(Math::mult1(mu,kset->CreateMatrixDeriv(levu)),*cpts);
}



template<class T>
Vector<T> BspCurv<T>::ComputePoints(int n) const
{
	Vector<T> v(n+1);

	double step = ((*kts)[num]-(*kts)[ord-1])/(double)n;
		// evaluate
	double val;
	for (int i=0; i<=n; i++) {
		val = (*kts)[ord-1]+i*step;
		v[i] = (*this)(val);
	}
	return v;
}

template<class T>
Matrix<T> BspCurv<T>::ComputeTable(int n) const
{
	Matrix<T> v(n+1,2);

	double step = ((*kts)[num]-(*kts)[ord-1])/(double)n;
		// evaluate
	double val;
	for (int i=0; i<=n; i++) {
		val = (*kts)[ord-1]+i*step;
		v[i][0] = i*step;
		v[i][1] = (*this)(val);
	}
	return v;
}


template<class T>
double BspCurv<T>::fabsnew(T t) const
{
	return sqrt(t*t);
}


template<class T>
BspCurv<T> BspCurv<T>::CreateApproximation(const FnCurvObject<T>& curve) const
{
	BspCurvBasisFuncSet b(*kts,ord,ord+num);
	Matrix<double> mat = b.ComputeBasisMatrix();
	Vector<double> vec = (*kset).ComputeKnotSetAver();
	Vector<T> vec1(num);
	for (int i=0; i<vec.GetNum(); i++) vec1[i] = curve(vec[i],0);
	return BspCurv<T>(Math::Solve(mat,vec1),*kts,ord,num);
}


// KNOT REMOVAL

// remove all possible knots from the curve 
// uses knot removal algorithm from Tiller paper
// requires a tolerance
template<class T>
BspCurv<T> BspCurv<T>::RemovePossibleKnots(double tol) const
{
	// make a local copy
	BspCurv<T> temp(*this);

	// try to remove each distinct knot up to its multiplicity
	// updates curve after each single removal
	for (int i=1; i<(*kset).GetNumDistinct()-1; i++)
		for (int j=0; j<(*kset).GetMult()[i]; j++)
			if (temp.IsKnotRemovable((*kset).GetDistinctKnots()[i],tol))
				temp = temp.KnotRemoval((*kset).GetDistinctKnots()[i]);
		
	return temp;
}


// removes a single knot from the curve
// need to check on the range for the knot
template<class T>
BspCurv<T> BspCurv<T>::KnotRemoval(double knot) const
{

	Vector<T> ncpts(*cpts);

	// find position of knot for removal
	int p = (*kset).Fndint1(knot);
	int s = (*kset).FindMultiplicity(knot);

	// compute i, j values 
	int i=p+s-ord+1;
	int j=p;

	// compute number of equations
	int num_eqns=ord-s;
	
	double alpha_i, alpha_j;

	// claculate new control point values
	while (j-i > 0)
	{
		alpha_i = (knot-(*kts)[i])/((*kts)[i+ord]-(*kts)[i]);
		alpha_j = (knot-(*kts)[j])/((*kts)[j+ord]-(*kts)[j]);

		ncpts[i] = ((*cpts)[i]-(1-alpha_i)*ncpts[i-1])/alpha_i;
		ncpts[j] = ((*cpts)[j]-alpha_j*ncpts[j+1])/(1-alpha_j);
		i++;
		j--;
	}
	// create CPoints object with new control points
	CPoints<T> cpts1(ncpts,*kset,num);

	// remove according to whether number of equations is even or odd
	div_t div_result;
	div_result = div(num_eqns,2);
	if (div_result.rem == 0) {
		i--;
		j++;
		cpts1= cpts1.Removal(j,p+s);
	} else {
		cpts1 = cpts1.Removal(i,p+s);
	}
	// return new BspCurv
	return BspCurv<T>(cpts1,num-1);
}


// checks to see if a knot is removable or not
// within a tolerance tol
template<class T>
bool BspCurv<T>::IsKnotRemovable(double knot, double tol) const
{

	Vector<T> ncpts(*cpts);


	int p = (*kset).Fndint1(knot);
	int s = (*kset).FindMultiplicity(knot);

	int i=p+s-ord+1;
	int j=p;

	// compute number of equations
	int num_eqns=ord-s;

	// compute new control points
	double alpha_i, alpha_j, alpha;
	while (j-i > 0)
	{
		alpha_i = (knot-(*kts)[i])/((*kts)[i+ord]-(*kts)[i]);
		alpha_j = (knot-(*kts)[j])/((*kts)[j+ord]-(*kts)[j]);

		ncpts[i] = ((*cpts)[i]-(1-alpha_i)*ncpts[i-1])/alpha_i;
		ncpts[j] = ((*cpts)[j]-alpha_j*ncpts[j+1])/(1-alpha_j);
		i++;
		j--;
	}

	// check depends on whether number of equations is even or odd
	div_t div_result;
	div_result = div(num_eqns,2);
	if (div_result.rem == 0) {
		i--;
		j++;
		if (fabs(ncpts[i]-ncpts[j]) < tol) return true;
		else {return false;}
	} else {
	   i--;
	   j++;
	   alpha = (knot-(*kts)[i+1])/((*kts)[i+1+ord]-(*kts)[i+1]);
	   if (fabs((*cpts)[i+1]-(alpha*ncpts[i]+(1-alpha)*ncpts[j])) < tol || 
		   fabs((*cpts)[i+1]-(alpha*ncpts[j]+(1-alpha)*ncpts[i])) < tol) return true;
	   else { return false;}
	}
}



template<class T>
BspCurv<T> BspCurv<T>::MakeBreakCompatable(const KnotSet& Kset) const
{	
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*kset).Normalise(Kset);
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(Kset);
	
	// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);
	
	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspCurv<T>(*cpts,norm.GetKnots(),ord,num).InsertKnot(kdiff1.GetDistinctKnots(),kdiff1.GetNumDistinct());
	} else return BspCurv<T>(*cpts,norm.GetKnots(),ord,num);
	
}


// ADDITION and SUBTRACTION


// add two BspCurv curves together by adding their control points
// curves are not necessarily compatable
template<class T>
BspCurv<T> BspCurv<T>::Add(const BspCurv<T>& b) const
{
	BspCurv<T> a(*this), c(b);
	// make the two curves degree compatable
	if (ord > b.GetOrd()) c = b.Elevate(ord-b.GetOrd());
	if (b.GetOrd() > ord) a = Elevate(b.GetOrd()-ord);

	BspCurv<T> d = a.MakeCompatable(c);
	BspCurv<T> e = c.MakeCompatable(d);

	// now add compatable curves together
	
	// add control points
	Vector<T> temp(d.GetNum());
	for (int i=0; i<d.GetNum(); i++)
		temp[i]=d.GetCPoints()[i]+e.GetCPoints()[i];

	return BspCurv<T>(temp,d.GetKnots(),d.GetOrd(),d.GetNum());
}

template<class T>
BspCurv<T> BspCurv<T>::MakeKnotSetCompatable(const KnotSet& nkset) const
{
	// assume they are both of the same degree
	// create knot set for current object
	
	// create normalised knot set for current object relative to b
	KnotSet norm = (*kset).Normalise(nkset);
	
	// add knots to current object in knot set for b not present 
	// in knot set for Current object
	KnotSet add1 = norm.KnotUnion(nkset);
	
		// find differences in knots
	KnotSet kdiff1 = add1.KnotDifference(norm);

	if (kdiff1.GetNum() > 0) {
		// insert knots into current object to get compatibility with b
		return 	BspCurv<T>(*cpts,norm.GetKnots(),ord,num).InsertKnot(kdiff1.GetDistinctKnots(),kdiff1.GetMult(),kdiff1.GetNumDistinct());
	} else return BspCurv<T>(*cpts,norm.GetKnots(),ord,num);
}


// subtract two BspCurv curves together by adding their control points
// curves should be compatable and have same order!
template<class T>
BspCurv<T> BspCurv<T>::Subtract(const BspCurv<T>& b) const
{
	BspCurv<T> a(*this), c(b);
	// make the two curves degree compatable
	if (ord > b.GetOrd()) c = b.Elevate(ord-b.GetOrd());
	if (b.GetOrd() > ord) a = Elevate(b.GetOrd()-ord);

	BspCurv<T> d = a.MakeCompatable(c);
	BspCurv<T> e = c.MakeCompatable(d);

	// now add compatable curves together
	
	// add control points
	Vector<T> temp(d.GetNum());
	for (int i=0; i<d.GetNum(); i++)
		temp[i]=d.GetCPoints()[i]-e.GetCPoints()[i];

	return BspCurv<T>(temp,d.GetKnots(),d.GetOrd(),d.GetNum());
}



// DEGREE ELEVATION

// elevate the degree of the BspCurv by level
template<class T>
BspCurv<T> BspCurv<T>::Elevate(int level) const
{
	if (level <= 0) return *this;
	// create knot set
	KnotSet nkset = (*kset).CreateKnotSetElevate(level);
	
	// create and return curve
	return BspCurv<T>(ElevateCPoints(level),nkset.GetKnots(),ord+level,nkset.GetNum()-(ord+level));
}


// elevate the degree of the BspCurv by level
template<class T>
Vector<T> BspCurv<T>::ElevateCPoints(int level) const
{
	if (level <= 0) return *cpts;
	// create BspCurv of order level+1 representing identity 
	Vector<double> knots(2*(level+1));

	// create knot set for identity
	for (int i=0; i<level+1; i++) { 
		knots[i]=(*kts)[ord-1]; 
		knots[level+1+i]=(*kts)[num];
	}
	// create control points for identity
	Vector<T> cpts1(level+1,1.0);

	// use product algorithm to elevate degree
	// multiply current object by identity
	return Product(BspCurv<T>(cpts1,knots,level+1,level+1)).GetCPoints();
}

// DERIVATIVES


// compute the derivative of the BspCurv of order deriv and
// represent the result as another BspCurv
template<class T>
BspCurv<T> BspCurv<T>::Derive(int level) const
{   
	// if deriv == 0 return current object	
	if (level <= 0) return *this;

	else if (level >= ord) {
		Vector<double> knots((*kset).GetDistinctKnots());
		T t = 0.0;
		Vector<T> cpts1((*kset).GetNumDistinct()-1,t);
		return BspCurv<T>(cpts1,knots,1,(*kset).GetNumDistinct()-1);
	}
	// create knot set object
	KnotSet nkset = (*kset).CreateKnotSetDeriv(level);
	
	// create the new BspCurv object
	return BspCurv<T>(DeriveCPoints(level),nkset.GetKnots(),ord-level,nkset.GetNum()-(ord-level));
}

// compute the derivative of the BspCurv of order deriv and
// return the Vector of control points for the derivative curve
template<class T>
Vector<T> BspCurv<T>::DeriveCPoints(int deriv) const
{   
	// return control points if deriv == 0
	if (deriv <= 0) return *cpts;
	else if (deriv >= ord) {
		T t = 0.0;
		return Vector<T>((*kset).GetNumDistinct()-1,t);
	}

	// create control points for derivative
	// note this computes the knot set as well
	return (CPoints<T>(*cpts,*kset,num).CreateCPointsDeriv(deriv)).GetCPoints();
}

// evaluate the derivative of the BspCurv of order deriv at a point
// x. Computes the derivative as a BspCurv and then evaluates this at x
template<class T>
T BspCurv<T>::Derive(int level, double x) const
{
	if (level <= 0) return (*this)(x);
	else if (level >= ord) return (T)0;
	return Derive(level)(x);
}

// computes level derivative at x
template<class T>
T BspCurv<T>::Derive1(int level, double x) const
{
	if (level <= 0) return (*this)(x);
	if (level >= ord) return (T)0;

	Matrix<double> a(ord,ord,0.0);
	Vector<double> dp(ord), dm(ord);

	int ind = (*kset).Find_index(x);

	for (int i=0; i<ord; i++) {
		a[i][0] = (*cpts)[ind+1-ord+i];
		dp[i] = (*kts)[ind+1+i]-x;
		dm[i] = x - (*kts)[ind+1-ord+i];	
	}

	for (int i=0; i<ord-1; i++)
		for (int j=i+1; j<ord; j++)
			a[j][i+1] = (a[j][i]-a[j-1][i])/(dm[j]+dp[j-i-1]);

	Matrix<double> mat = BspCurvBasisFuncSet(*kts,ord,ord+num).ComputeNMatrix(x);

	Vector<double> v(ord,0.0);

	for (int i=0; i<ord-level; i++) v[i+1]=v[i]+mat[i][ord-level-1]*a[i+level][level];
	
	for (int i=0; i<level; i++) v[ord-level] = (ord-i-1)*v[ord-level];

	return v[ord-level]; 
}

// computes all derivatives at x
template<class T>
Vector<double> BspCurv<T>::Derive1(double x) const
{ 
	if (ord == 1) return Vector<T>(ord,0.0);
	Matrix<double> a(ord,ord,0.0);
	Vector<double> dp(ord), dm(ord);

	int ind = (*kset).Find_index(x);

	for (int i=0; i<ord; i++) {
		a[i][0] = (*cpts)[ind+1-ord+i];
		dp[i] = (*kts)[ind+1+i]-x;
		dm[i] = x - (*kts)[ind+1-ord+i];
	}
	for (int i=0; i<ord-1; i++)
		for (int j=i+1; j<ord; j++)
			a[j][i+1] = (a[j][i]-a[j-1][i])/(dm[j]+dp[j-i-1]);

	Matrix<double> mat = BspCurvBasisFuncSet(*kts,ord,ord+num).ComputeNMatrix(x);
	Vector<double> v(ord,0.0);
	Vector<double> res(ord,0.0);

	for (int i=0; i<ord-1; i++) {
		for (int j=0; j<ord-i-1; j++) v[j+1]=v[j]+mat[j][ord-i-2]*a[j+i+1][i+1];
		for (int k=0; k<=i; k++) v[ord-i-1] = (ord-k-1)*v[ord-i-1];
		res[i] = v[ord-i-1];
	}
	return res;
}


// computes all derivatives at x using knot insertion
template<class T>
Vector<double> BspCurv<T>::Derive2(double x) const
{
	if (ord == 1) return Vector<T>(ord,0.0);
	std::multiset<double>::iterator r3 = (*kset).GetKnotSet().upper_bound(x);


	int segno = (*kset).Find_segment(x);
	// take two consectutive knots
	double x1, x2;

	if (r3 == (*kset).GetKnotSet().end()) {
		x1 = (*kset).GetDistinctKnots()[segno-1];
		x2 = (*kset).GetDistinctKnots()[segno];
	} else {
		x1 = x;
		x2 = (*kset).GetDistinctKnots()[segno];
	}
	int s1 = (*kset).FindMultiplicity(x1);
	int s2 = (*kset).FindMultiplicity(x2);
	
	BspCurv<T> p1 = InsertKnot(x1,ord-s1);
	BspCurv<T> p2 = p1.InsertKnot(x2,ord-s2);

	int nsegno = p2.GetKnotSet().Find_segment(x1);
	// can do better by using derivative of Bez at endpoints
	return p2.GetSegment(nsegno).Derive(x);
}



template<class T>
Vector<T> BspCurv<T>::ComputeKnotSetAverValues() const
{
	Vector<T> m(num);

	Vector<double> uval((*kset).ComputeKnotSetAver());

	for (int i=0; i<num; i++) m[i] = (*this)(uval[i]);

	return m;
}



// INTEGRATION
// integrate the BspCurv between the limits x1 and x2. Computes
// the indefinite integral as a BspCurv and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T BspCurv<T>::Integrate(double x1, double x2) const
{
	// create the indefinite integral;
	BspCurv<T> intCurve = Integrate(); 

	// evaluate and subtract
	KnotSet kset = intCurve.GetKnotSet();

	if (x2 < kset.GetKnots()[intCurve.GetOrd()-1] || x1 > kset.GetKnots()[intCurve.GetNum()]) return 0;

	if (x1 < kset.GetKnots()[intCurve.GetOrd()-1]) x1 = kset.GetKnots()[intCurve.GetOrd()-1];
	if (x2 > kset.GetKnots()[intCurve.GetNum()]) x2 = kset.GetKnots()[intCurve.GetNum()];

	return (intCurve(x2) - intCurve(x1));
}


// compute the indefinite integral of the BspCurv and represent
// it as a BspCurv of one higher degree
template<class T>
BspCurv<T> BspCurv<T>::Integrate(int level) const
{
	BspCurv<T> b = Integrate();

	if (level >= 0) {
		for (int i=1; i<level; i++) b = b.Integrate();
		return b;
	}
	else return *this;
}	


// compute the indefinite integral of the BspCurv and represent
// it as a BspCurv of one higher degree
template<class T>
BspCurv<T> BspCurv<T>::Integrate() const
{
	// create knot set for integral
	KnotSet nkset = (*kset).CreateKnotSetIntegrate();
	// create the Curve
	return BspCurv<T>(IntegrateCPoints(),nkset.GetKnots(),ord+1,nkset.GetNum()-(ord+1));
}	


// compute the indefinite integral of the BspCurv as a BspCurv
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Vector<T> BspCurv<T>::IntegrateCPoints() const
{
	// return the control points as a CPoints object
	return 	(CPoints<T>(*cpts,*kset,num).CreateCPointsIntegrate()).GetCPoints();
}	



// SUBDIVISION

// subdivide the BspCurv upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspCurv after the knot refinement
template<class T>
BspCurv<T> BspCurv<T>::Subdivide(int level) const
{
	if (level <= 0) return *this;
	// create knot set for subdivision
	KnotSet nkset = (*kset).CreateKnotSetSubdivide(level);

	// create and return the BspCurv representing the refinement
	return BspCurv<T>(SubdivideCPoints(level),nkset.GetKnots(),ord,nkset.GetNum()-ord);
}


// subdivide the BspCurv upto level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns only the control points after the knot refinement
// as a CPoints object
template<class T>
Vector<T> BspCurv<T>::SubdivideCPoints(int level) const
{
	if (level <= 0) return *cpts;
	// create and return the control points as a CPoints object
	return (CPoints<T>(*cpts,*kset,num).CreateCPointsSubdivide(level)).GetCPoints();
}
    

// subdivide the BspCurv between the limits x1 and x2
// and return the new BspCurv
template<class T>
BspCurv<T> BspCurv<T>::Subdivide(double x1, double x2) const
{	
	// create subdivision knot set
	//if ((x1 - (*kts)[ord-1]) >  && x2 <= (*kts)[num-1]) {
		KnotSet nkset = (*kset).CreateKnotSetSubdivide(x1,x2);
	
	// create and return the subdivided BspCurv 
		return BspCurv<T>(SubdivideCPoints(x1,x2),nkset.GetKnots(),ord,nkset.GetNum()-ord);
	//} else return *this;
}

// subdivide the BspCurv between the limits x1 and x2
// and return the new BspCurv
template<class T>
Vector<T> BspCurv<T>::SubdivideCPoints(double x1, double x2) const
{		
	//if (x1 >= (*kts)[ord-1] && x2 <= (*kts)[num-1]) {
	// create new control points
		return (CPoints<T>(*cpts,*kset,num).CreateCPointsSubdivide(x1,x2)).GetCPoints();
//	} else return *cpts;
}


// KNOT INSERTION

 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspCurv<T> BspCurv<T>::InsertKnot(const Vector<double>& t,const Vector<int>& mult,int n) const
{
		// create the knot set
	if ( n <= 0) return *this;
	KnotSet nkset = (*kset).CreateKnotSetInsert(t,mult,n);
	// create and return the new BspCurv
	return BspCurv<T>(InsertKnotCPoints(t,mult,n),nkset.GetKnots(),ord,nkset.GetNum()-ord);
}
             
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Vector<T> BspCurv<T>::InsertKnotCPoints(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// create a local copy of the BspCurv
	if (n <= 0) return *cpts;

	BspCurv<T> temp(*this);
	for (int i=0; i<n; i++) temp = temp.InsertKnot(t[i],mult[i]);	
    return temp.GetCPoints();
}


// MATRIX DERIVATIVE

// find the derivative of the BspCurv of level given by deriv and
// represent the resulting control points of the derivative BspCurv
// in terms of the original control points of the Curve. Returns
// the matrix where each row give each new control point as a linear
// combination of the old
template<class T>
Matrix<T> BspCurv<T>::Dmatrix(int deriv) const
{	
	KnotSet nkset(*kset);

	// create matrix for first derivative
	Matrix<double> d = (*kset).CreateMatrixDeriv();

	// create knot sets for each subsequent derivative up to deriv
	// find new matrix and multiply by previous matrix.
	for (int i=1; i<deriv; i++) {
		nkset = nkset.CreateKnotSetDeriv(1);
		d = nkset.CreateMatrixDeriv()*d;
	}	
	// return matrix of linear combinations
	return d;
}


// CONVERSIONS

// convert the curve to composite Bezier form
template<class T>
CompBezCurv<T> BspCurv<T>::ConvertCompBezCurv() const
{
	KnotSet kset1 = (*kset).CreateKnotSetCompBezCurv();
	KnotSet kset2 = kset1.KnotDifference(*kset);
	if (kset2.GetNum() > 0) {
		return CompBezCurv<T>(kset1.GetNumDistinct()-1,ord,InsertKnotCPoints(kset2.GetDistinctKnots(),kset2.GetMult(),kset2.GetNumDistinct()),kset1.GetKnots());
	} else {
		return CompBezCurv<T>(kset1.GetNumDistinct()-1,ord,*cpts,kset1.GetKnots());
	}
}



// convert the curve to composite polynomial form
template<class T>
CompPolyCurv<T> BspCurv<T>::ConvertCompPolyCurv() const
{
	// convert to composite Bezier and then to composite poly
	return ConvertCompBezCurv().ConvertCompPolyCurv();
}


// READ and WRITE
template <class T>
void BspCurv<T>::write(std::ostream& os) 
{
	os << "Bspline Curve\n";
	os << "order is " << ord << "\n";
	os << "number of control points is " << num << "\n";
	os << "knots are\n";
	os << *kts;
	os << "\ncontrol points are";
	os << *cpts;
}

template <class T>
void BspCurv<T>::read(std::istream& is)
{
	int Ord;
	std::cout << "order of Bspline curve\n";
	is >> Ord;
	int Num;
	std::cout << "number of control points\n";
	is >> Num;
	Vector<double> Kts(Ord+Num);
	std::cout << "input knots\n";
	is >> Kts;
	Vector<T> Cpts(Num);
	std::cout << "input control points\n";
	is >> Cpts;
	*this = BspCurv<T>(Cpts,Kts,Ord,Num);
} 


template <class T>
void BspCurv<T>::writefile(std::ofstream& ofs) 
{
	ofs << ord << " " << num << "\n";
	ofs << *kts;
	ofs << *cpts;
}

template <class T>
void BspCurv<T>::readfile(std::ifstream& ifs)
{
	int Ord, Num;
	ifs >> Ord >> Num;
	Vector<double> Kts(Ord+Num);
	Vector<T> Cpts(Num);
	ifs >> Kts;
	ifs >> Cpts;
	*this = BspCurv<T>(Cpts,Kts,Ord,Num);
} 




// elevate the degree of the BspCurv by level
template<class T>
BspCurv<T> BspCurv<T>::Elevate2(int level) const
{
	if (level <= 0) return *this;
	// create BspCurv of order level+1 representing identity 
	Vector<double> knots(2*(level+1));

	// create knot set for identity
	for (int i=0; i<level+1; i++) { 
		knots[i]=(*kts)[ord-1]; 
		knots[level+1+i]=(*kts)[num];
	}
	// create control points for identity
	Vector<T> cpts1(level+1,1.0);

	// use product algorithm to elevate degree
	// multiply current object by identity
	return Product(BspCurv<T>(cpts1,knots,level+1,level+1));
}



// inserts a vector of n knots into the BspCurv and returns the new
// BspCurv
template<class T>
BspCurv<T> BspCurv<T>::InsertKnot(const Vector<double>& Kts, int n) const
{
	if (n <= 0) return *this;
	// create knot set
	KnotSet nkset = (*kset).CreateKnotSetInsert(Kts,n);
	
	// create and return the new BspCurv
	return BspCurv<T>(InsertKnotCPoints(Kts,n),nkset.GetKnots(),ord,nkset.GetNum()-ord);
}
   

// inserts a vector of n knots into the BspCurv and returns the new
// BspCurv
template<class T>
Vector<T> BspCurv<T>::InsertKnotCPoints(const Vector<double>& Kts, int n) const
{
	// create a local copy of the BspCurv
	if (n<=0) return *cpts;
    BspCurv temp(*this);
    
	// perform knot insertion 
	for(int i=0;i<n;i++)
		temp=temp.InsertKnot(Kts[i]);

	// return the new Curve
	return temp.GetCPoints();
}


         
// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
BspCurv<T> BspCurv<T>::InsertKnot(double x) const
{       
	// create knot set
	KnotSet nkset = (*kset).CreateKnotSetInsert(x);
	
	// create and return the new BspCurv
	return BspCurv<T>(InsertKnotCPoints(x),nkset.GetKnots(),ord,nkset.GetNum()-ord);
}                              
      
// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
Vector<T> BspCurv<T>::InsertKnotCPoints(double x) const
{       
	// calculate the new control points
	// also updates the knot set
	return (CPoints<T>(*cpts,*kset,num).CreateCPointsInsert(x)).GetCPoints();
}              


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
BspCurv<T> BspCurv<T>::InsertKnot(double x, int level) const
{	
	// create the knot set
	if (level <= 0) return *this;
	KnotSet nkset = (*kset).CreateKnotSetInsert(x,level);
		// create and return the new BspCurv
	return BspCurv<T>(InsertKnotCPoints(x,level),nkset.GetKnots(),ord,nkset.GetNum()-ord);
}


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
Vector<T> BspCurv<T>::InsertKnotCPoints(double x, int level) const
{
	// calculate the new control points
	// also updates the knot set
	if (level <= 0) return *cpts;
	return (CPoints<T>(*cpts,*kset,num).CreateCPointsInsert(x,level)).GetCPoints();
}

// convert the curve to composite Bezier form
template<class T>
Vector<T> BspCurv<T>::ConvertCompBezCurvCPoints() const
{	
	// check to see it is not already Bezier or CompBezCurv

	Vector<double> ndts((*kset).GetNumDistinct());
	Vector<int> nmult((*kset).GetNumDistinct());

	int count=0;
	// insert knots up to multiplicity ord 
	for (int i=1; i<(*kset).GetNumDistinct()-1; i++) 
		if ((*kset).GetMult()[i] < ord) {
			ndts[count]=(*kset).GetDistinctKnots()[i];
			nmult[count]=ord-(*kset).GetMult()[i];
			count++;
		}

	if (count > 0) {
		// insert the knots, what if nmult[i]=0?
		return InsertKnotCPoints(ndts,nmult,count);
	} else return *cpts;
}



// convert the curve to composite Bezier form
template<class T>
Vector<T> BspCurv<T>::ConvertCompBezCurvCPoints2() const
{	
	// check to see it is not already Bezier or CompBezCurv
	Vector<double> ndts((*kset).GetNumDistinct());
	Vector<int> nmult((*kset).GetNumDistinct());

	int count=0;
	// insert knots up to multiplicity ord - 1 
	for (int i=1; i<(*kset).GetNumDistinct()-1; i++) 
		if ((*kset).GetMult()[i] < ord-1) {
			ndts[count]=(*kset).GetDistinctKnots()[i];
			nmult[count]=ord-(*kset).GetMult()[i];
			count++;
		}

	if (count > 0) {
		// insert the knots, what if nmult[i]=0?
		return InsertKnotCPoints(ndts,nmult,count);
	} else return *cpts;
}




#endif
 
