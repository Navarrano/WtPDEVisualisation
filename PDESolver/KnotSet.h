


#ifndef KNOT
#define KNOT

#include "mathfunctions.h"
#include "matrix.h"


struct SortCollection : public TextObject, public FTextObject {
	double kts;
	double dts;
	int mult;

	SortCollection() { }

	SortCollection(double Kts, double Dts, int Mult) :
	kts(Kts), dts(Dts), mult(Mult)	{ }
	
	bool operator<(const SortCollection& s) const
	{
		return dts < s.dts;
	}
	
	bool operator==(const SortCollection& s) const 
	{
		return ((kts == s.kts) && (dts == s.dts) && (mult == s.mult));
	}
	
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

class KnotSet : public TextObject, public FTextObject {

	// data 
	int ord; // order of associated curve
	int num; // number of knots
	int n;	// number of distinct knots
	
	Ptr<Vector<double> > kts;	// knots
	Ptr<std::multiset<double> > mset;
	Ptr<Vector<double> > dts;	// distinct knots
	Ptr<Vector<int> > mult; // multiplicities
	

	// internal functions	
	// alternate constructor
	KnotSet SortKnotSet() const;
	KnotSet Sort(int N) const;
	KnotSet(const Vector<double>& Kts, const Vector<double>& Dts, 
		const Vector<int>& Mult, int N); // constructor used in Sort and product
	virtual ObjectID Identity() const { return std::string("class KnotSet"); }
public:
	// constructors
	KnotSet();
	KnotSet(int Num);
	KnotSet(int Ord, int Num);
	KnotSet(const Vector<double>& Kts, int Ord, int Num);
	KnotSet(const Vector<double>& Dts, const Vector<int>& Mult, int Ord, int N);
	
	// check functions
	bool CheckKnotSet() const; // checks to see if end knots have multiplicity ord
	KnotSet UpdateKnotSet();
	KnotSet& UpdateKnotsLeft(); // converts knot set to have multiplicity ord at both ends
	KnotSet& UpdateKnotsRight();

	// test for equality of two knot sets

	bool IsSameKnotSet(const KnotSet& kset) const;
	bool IsSameDistinctKnotSet(const KnotSet& kset) const;

	// knot averages
	Vector<double> ComputeKnotSetAver() const;

	KnotSet AddDistinctKnots(const KnotSet& kset) const;
	
	// matrices for derivatives
	Matrix<double> CreateMatrixDeriv() const;
	Matrix<double> CreateMatrixDeriv(int level) const;

	// interpolation vectors
	Vector<double> CreateVectorInterp(double val) const;
	Vector<double> CreateVectorInterpDeriv(int level, double val) const;

	// access functions
	Vector<int> GetMult() const;
	Vector<double> GetDistinctKnots() const;
	Vector<double> GetKnots() const;
	int GetOrd() const;
	int GetNum() const;
	int GetNumDistinct() const;
	int Fndint1(double x) const;
	int Find_index(double x) const;
	int Find_segment(double x) const;
	std::multiset<double> GetKnotSet() const;
	Vector<double> ComputeDistinctKnots() const;
	Vector<double> ComputeKnots() const;
	Vector<int> ComputeMultiplicities() const;
	int FindNumKnots() const;
	int FindNumDistinctKnots() const;
	int FindMultiplicity(double x) const;
	int ComputeParameterAverUpper(double x) const;
	int ComputeParameterAverLower(double x) const;

	// find the difference between two knot sets, one contains the other
	KnotSet KnotDifference(const KnotSet& kset) const; // for Bcurve product
	KnotSet KnotIntersection(const KnotSet& kset) const;
	KnotSet KnotUnion(const KnotSet& kset) const;
	bool Contains(const KnotSet& kset) const;

	// knot set normalisation
	KnotSet Normalise() const;
	KnotSet Normalise(double a, double b) const;
	KnotSet Normalise(const KnotSet& kset) const;
	
	// particular knot set creation
	KnotSet CreateKnotSetProduct(const KnotSet& kset) const;
	KnotSet CreateKnotSetDeriv(int deriv) const;
	KnotSet CreateKnotSetDeriv() const;
	KnotSet CreateKnotSetRemoval(double knot) const;
	KnotSet CreateKnotSetCompBezCurv() const;
	KnotSet CreateKnotSetIntegrate() const;
	KnotSet CreateKnotSetCurveBasis(int index) const;
	KnotSet CreateKnotSetElevate(int level) const;
	KnotSet CreateKnotSetSubdivide(int level) const;
	KnotSet CreateKnotSetSubdivide(double x1, double x2) const;
	KnotSet CreateKnotSetInsert(double x) const;
	KnotSet CreateKnotSetInsert(double x, int level) const;
	KnotSet CreateKnotSetInsert(const Vector<double>& Kts, int N) const;
	KnotSet CreateKnotSetInsert(const Vector<double>& Kts, const Vector<int>& mult, int N) const;


	// matrix of integrals of products of basis functions
	Matrix<double> CreateMatrixIntegral(int deriv, double x1, double x2) const;

	Matrix<double> ComputeLeastSquaresMatrix(const Vector<double>& tau, int m) const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};


#endif