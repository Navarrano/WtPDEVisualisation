
#ifndef BEZCURV
#define BEZCURV


#include "knotset.h"
#include "polycurv.h"


// class for a basis function of a BezCurv (a Bezier basis function)
class BezCurvBasisFunc : public Curve<double> {

	// data
	int ord;	 // order of basis function
	Ptr<Vector<double> > kts;		// knots 
	Ptr<BezCurv<double> > b;	// BezCurv representation of basis function
	BezCurv<double> CreateBezCurv() const;	// creates the BezCurv
	virtual ObjectID Identity() const { return std::string("class BezCurvBasisFunc"); }
public:

	// constructors
	BezCurvBasisFunc();		// default constructor
	BezCurvBasisFunc(const Vector<double>& Kts, int Ord);

	// access functions
	int ComputeDim() const;		// computes the dimensio
	BezCurv<double> GetBezCurv() const;  // gets the BezCurv
	Vector<double> GetKnots() const;
	int GetOrd() const;

	// evaluators
	double Eval(double x) const;	// evaluates the basis function
	virtual double operator()(double x) const;
	virtual double operator()(int val, double x) const;
	
	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual double Derive(int, double) const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
};


class BezCurvBasisFuncSet : public TextObject, public FTextObject {

	int ord;
	Ptr<Vector<double> > kts;
	Ptr<Vector<BezCurvBasisFunc> > b;
	virtual ObjectID Identity() const { return std::string("class BezCurvBasisFuncSet"); }
public:
	// constructors
	BezCurvBasisFuncSet();
	BezCurvBasisFuncSet(const Vector<double>& Kts, int Ord);

	// get functions 
	Vector<double> GetKnots() const;
	int GetOrd() const;
	BezCurvBasisFunc GetBasisFunc(int i) const;

	
	// evaluators
	Vector<double> operator()(double x) const;
	Vector<double> Eval(double x) const;
	Vector<double> operator()(int val, double x) const;
	Vector<double> Derive(int, double) const;

	// integration
	Matrix<double> CreateMatrixIntegral(int deriv, double x1, double x2) const;

	// read and write
	virtual void read(std::istream& is);
	virtual void readfile(std::ifstream& ifs);
	virtual void write(std::ostream& os);
	virtual void writefile(std::ofstream& ofs) ;
};


template<class T>
class BezCurv : public Curve<T> {
private:
	// data
	int ord;
	Ptr<Vector<T> > cpts;
	Ptr<Vector<double> > kts;
	Ptr<KnotSet> kset;
	double leftLimit;
	double rightLimit;

	// private functions
	BspCurv<T> InsertKnot(const Vector<double>& Kts, int n) const;
	BspCurv<T> InsertKnot(double x) const;
	Vector<T> InsertKnotCPoints(double x) const;
	Vector<T> InsertKnotCPoints(const Vector<double>& t,const Vector<int>& mult,int n) const;
	BspCurv<T> InsertKnot(double x, int level) const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class BezCurv<double>");
		else {
			std::string s(typeid(T).name()), s1("class BezCurv<");
			return s1.append(s,0,s.size())+">";
		}
		/*if (std::string(typeid(T).name()) == std::string("class double"))
			return std::string("class BezCurvDouble"); 
		else if (std::string(typeid(T).name()) == std::string("class Point1D"))
			return std::string("class FBezCurv");
		else if (std::string(typeid(T).name()) == std::string("class Point2D"))
			return std::string("class BezCurv2D");
		else return std::string("class BezCurv3D"); */
	}
public:
	// constructors
	BezCurv();
	BezCurv(const Vector<T>& Cpts, const Vector<double>& Kts, int Ord);
	BezCurv(const Vector<T>& Cpts, int Ord);
	BezCurv(const Vector<T>& Cpts, int Ord, double LeftLimit, double RightLimit);
	BezCurv(int Ord);


	// access functions
	int GetOrd() const;
	int GetNum() const;
	Vector<T> GetCPoints() const;
	Vector<double> GetKnots() const;
	KnotSet GetKnotSet() const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
	BezCurvBasisFuncSet GetBasisFuncSet() const;
	BezCurvBasisFunc GetBasisFunc(int i) const;

	

	// evaluation
	virtual T operator()(double x) const;
	T Eval(double x) const;

	// addition and subtraction
	BezCurv<T> Add(const BezCurv<T>& b) const;
	BezCurv<T> Subtract(const BezCurv<T>& b) const;

	// derivatives
	Vector<T> Derive(double x) const;
	BezCurv<T> Derive(int level) const;
	Vector<T> DeriveCPoints(int level) const;
	virtual T Derive(int level, double x) const;
	virtual T operator()(int, double) const;

	// degree elevation
	BezCurv<T> Elevate(int level) const;
	Vector<T> ElevateCPoints(int level) const;

	// integration
	BezCurv<T> Integrate() const;
	T Integrate(double x1, double x2) const;
	Vector<T> IntegrateCPoints() const;

	// product 
	template<class T1>
	BezCurv<T> Product(const BezCurv<T1>& b) const 
	{
		// Create the product knot set
		KnotSet nkset = (*kset).CreateKnotSetProduct(KnotSet(b.GetKnots(),b.GetOrd(),2*b.GetOrd()));

		// create the new BezCurv Curve
		return BezCurv<T>(ProductCPoints(b),nkset.GetKnots(),nkset.GetOrd());
	}

	template<class T1>
	Vector<T> ProductCPoints(const BezCurv<T1>& b) const 
	{ 
		int ordb = b.GetOrd();

		// create the control points
		// algorithm for product from Piegl/Tiller Paper
		Vector<T> temp(ord+ordb-1);
		T sum = 0.0;
		for (int k=0; k<ord+ordb-1; k++) {
			for (int l=Math::max1(0,k-ordb+1); l<=Math::min1(ord-1,k); l++) {
				sum = sum + (Math::Combin(ord-1,l)*Math::Combin(ordb-1,k-l)/Math::Combin(ord+ordb-2,k))*(*cpts)[l]*b.GetCPoints()[k-l];
			}
			temp[k] = sum;
			sum=0.0;
		}
		// create the new BezCurv Curve
		return temp;
	}

	// subdivision
	BezCurv<T> Subdivide(double x1, double x2) const;
	Vector<T> SubdivideCPoints(double x1, double x2) const;
	Vector<T> SubdivideCPoints(int level) const;
	BspCurv<T> Subdivide(int level) const;

	// insert knots
	BspCurv<T> InsertKnot(const Vector<double>& t,const Vector<int>& mult,int n) const;
	Vector<T> InsertKnotCPoints(double x, int level) const;

	// derivative matrix
	Matrix<T> Dmatrix(int deriv) const;

	// conversion 
	PolyCurv<T> ConvertPolyCurv() const;
	BspCurv<T> ConvertBspCurv() const;
	CompPolyCurv<T> ConvertCompPolyCurv() const;
	CompBezCurv<T> ConvertCompBezCurv() const;



	// read and write functions
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};


// CONSTRUCTORS

// constructor builds a BezCurv from a Vector of control points, knots
// and an order
template<class T>
BezCurv<T>::BezCurv(const Vector<T>& Cpts, const Vector<double>& Kts, int Ord) :
			ord(Ord), cpts(new Vector<T>(Cpts)), kts(new Vector<double>(Kts)),
			kset(new KnotSet(*kts,ord,2*ord))
{
	leftLimit = (*kts)[ord-1];
	rightLimit = (*kts)[ord];
}

			// constructor, builds a BezCurv from a Vector of control points
// Num in size
template<class T>
BezCurv<T>::BezCurv(const Vector<T>& Cpts, int Ord, double LeftLimit, double RightLimit) : ord(Ord), 
	cpts(new Vector<T>(Cpts)), leftLimit(LeftLimit), rightLimit(RightLimit)
{
	Vector<double> Limits(2);
	Limits[0] = leftLimit;
	Limits[1] = rightLimit;
	kts = new Vector<double>(Math::CreateKnots(1,ord,Limits));	
	kset = new KnotSet(*kts,ord,2*ord);
}


// constructor, builds a BezCurv from a Vector of control points
// Num in size
template<class T>
BezCurv<T>::BezCurv(const Vector<T>& Cpts, int Ord) : ord(Ord), 
	cpts(new Vector<T>(Cpts))
{
	kts = new Vector<double>(Math::CreateKnots(ord));	
	kset = new KnotSet(*kts,ord,2*ord);
	leftLimit = (*kts)[ord-1];
	rightLimit = (*kts)[ord];
}


// default constructor
template<class T>
BezCurv<T>::BezCurv() : ord(0), cpts(), kts(), kset() { }

// builds a BezCurv from an order
template<class T>
BezCurv<T>::BezCurv(int Ord) : ord(Ord), cpts(new Vector<T>(Ord))
{ 
	// assign knots 0, 1
	kts= new Vector<double>(Math::CreateKnots(ord));
	kset = new KnotSet(*kts,ord,2*ord);
	leftLimit = (*kts)[ord-1];
	rightLimit = (*kts)[ord];
}



// ACCESS FUNCTIONS

// get the order of the BezCurv
template<class T>
inline int BezCurv<T>::GetOrd() const { return ord; }

// get the number of control points
template<class T>
inline int BezCurv<T>::GetNum() const { return ord; }

// get the control points
template<class T>
inline Vector<T> BezCurv<T>::GetCPoints() const { return *cpts; }

// get the knot vector
template<class T>
inline Vector<double> BezCurv<T>::GetKnots() const { return *kts; }

// get the knot vector
template<class T>
inline KnotSet BezCurv<T>::GetKnotSet() const { return *kset; }

template<class T>
inline double BezCurv<T>::GetLeftLimit() const { return leftLimit; }

template<class T>
inline double BezCurv<T>::GetRightLimit() const { return rightLimit; }


template<class T>
BezCurvBasisFuncSet BezCurv<T>::GetBasisFuncSet() const
{
	return BezCurvBasisFuncSet(*kts,ord);
}

template<class T>
BezCurvBasisFunc BezCurv<T>::GetBasisFunc(int i) const
{
	Vector<double> v(ord+1);

	for (int j=0; j<=ord; j++) v[j] = (*kts)[i-1+j]; 
	
	return BezCurvBasisFunc(v,ord);
}

// ADD and SUBTRACT

// add two BezCurv's together, should have the same knot range
// addtion has knot vector from current object
template<class T>
BezCurv<T> BezCurv<T>::Add(const BezCurv<T>& b) const
{
	BezCurv<T> c(*this), d(b);

	if (ord > b.GetOrd()) d = d.Elevate(ord-b.GetOrd());
	else if (b.GetOrd() > ord) c = c.Elevate(b.GetOrd()-ord); 

	Vector<T> temp(c.GetOrd());	
	// add the control points together
	for (int i=0; i<c.GetOrd(); i++)
		temp[i]=c.GetCPoints()[i]+d.GetCPoints()[i];

	return BezCurv<T>(temp,d.GetKnots(),c.GetOrd());
}


// subtract two BezCurv's together (*this-b), 
// should have the same knot range
// subtraction has knot vector from current object
template<class T>
BezCurv<T> BezCurv<T>::Subtract(const BezCurv<T>& b) const
{
	BezCurv<T> c(*this), d(b);

	if (ord > b.GetOrd()) d = d.Elevate(ord-b.GetOrd());
	else if (b.GetOrd() > ord) c = c.Elevate(b.GetOrd()-ord); 

	Vector<T> temp(c.GetOrd());	
	// add the control points together
	for (int i=0; i<c.GetOrd(); i++)
		temp[i]=c.GetCPoints()[i]-d.GetCPoints()[i];

	return BezCurv<T>(temp,d.GetKnots(),c.GetOrd());
}

// EVALUATION 

// evaluate the BezCurv at the point x using de Boor algorithm
template<class T>
T BezCurv<T>::operator()(double x) const
{
	Vector<T> temp(*cpts);
  
	double xnew = (x-leftLimit)/(rightLimit-leftLimit);
	//compute the point according to de Boor algorithm 
	// compute point on Curve according to triangular table
	for (int i=0; i<ord-1; i++) 
		for (int j=0; j<ord-i-1; j++) 
			temp[j] = xnew*temp[j+1]+(1.0-xnew)*temp[j];

	return temp[0];
}


// evaluate the BezCurv using member function rather 
// than overloaded operator
// identical code to that of overloaded operator
template<class T>
T BezCurv<T>::Eval(double x) const
{
  	// use matrix form of multiplication
	return ConvertPolyCurv()(x);
}

// DEGREE ELEVATION

// elevate the degree of the BezCurv by level
template<class T>
BezCurv<T> BezCurv<T>::Elevate(int level) const
{
	if (level <= 0) return *this;
	KnotSet nkset = (*kset).CreateKnotSetElevate(level);

	// return degree elevated curve
	return BezCurv<T>(ElevateCPoints(level),nkset.GetKnots(),ord+level);
}


// elevate the degree of the BezCurv by level
template<class T>
Vector<T> BezCurv<T>::ElevateCPoints(int level) const
{
	if (level <= 0) return *cpts;
	// create control points
	Vector<T> temp(ord+level,0.0);
	T sum = 0.0; // need a conversion from 0.0 to T 

	// BezCurv elevation algorithm
	for (int i=0; i<ord+level; i++) {
		for (int j=0; j<ord; j++)
			sum = sum + (*cpts)[j]*(Math::Combin(ord-1,j)/Math::Combin(ord-1+level,i))*Math::Combin(level,i-j);
		temp[i]=sum;
		sum =0.0;
	}
	return temp;
}


// DERIVATIVES

template<class T>
T BezCurv<T>::operator()(int lev, double x) const
{
	return Derive(lev)(x);
}


// compute the derivative of the BezCurv of order deriv and
// represent the result as another BezCurv
template<class T>
Vector<T> BezCurv<T>::Derive(double x) const
{   
	Vector<T> res(ord);
	double x1 = leftLimit;
	double x2 = rightLimit;
	Vector<double> dts(2);
	Vector<int> mult(2);
	dts[0] = x1;
	dts[1] = x2;
	res[0] = (*this)(x);

	for (int i=0; i<ord-1; i++) {
		Vector<T> deriv(ord-i-1,0.0);
		T sum = 0.0;
		mult[0] = ord-i-1;
		mult[1] = ord-i-1;
		for (int k=0; k<ord-1-i; k++) {
			for (int l=k; l<=k+i+1; l++) 
				sum = sum+ (*cpts)[l]*Math::Combin(i+1,i+1-l+k)*pow(-1.0,(double)(i+1-l+k));
			deriv[k] = sum*Math::Combin(ord-1,i+1)*Math::Fact(i+1)/pow(x2-x1,(double)(i+1));
		}
		res[i+1] = BezCurv<T>(deriv,KnotSet(dts,mult,ord-i-1,2).GetKnots(),ord-i-1)(x);
	}
	return res;
}



// compute the derivative of the BezCurv of order deriv and
// represent the result as another BezCurv
template<class T>
BezCurv<T> BezCurv<T>::Derive(int level) const
{   
	if (level <= 0) return *this;
	if (level >= ord) {
		// perhaps create a default degenerate curve in this case
		Vector<double> knots((*kset).GetDistinctKnots());
		Vector<T> cpts1((*kset).GetNumDistinct()-1,0.0);
		return BezCurv<T>(cpts1,knots,1);
	}
	// create knot set for curve
	KnotSet nkset = (*kset).CreateKnotSetDeriv(level);

	// create and return the derivative BezCurv
	return BezCurv<T>(DeriveCPoints(level),nkset.GetKnots(),ord-level);
}


// compute the derivative of the BezCurv of order deriv and
// represent the result as another BezCurv
template<class T>
Vector<T> BezCurv<T>::DeriveCPoints(int level) const
{   
	if (level <= 0) return *cpts;
	else if (level > ord-1) return Vector<T>(1,0.0);
	else {
	// create control points for derivative
	// according to BezCurv derivative formula
		Vector<T> temp(ord-level);
		T sum=0.0; // conversion of 0.0 to T required

		// derivative algorithm for new control points
		for (int i=0; i<ord-level; i++) {
			for (int j=0; j<=level; j++)
				sum = sum + Math::Combin(level,j)*pow(-1.0,(double)(level-j))*(*cpts)[i+j];
			temp[i] = sum * Math::Fact(ord-1)/Math::Fact(ord-1-level);
			temp[i] = temp[i]/pow(((*kts)[ord]-(*kts)[ord-1]),(double)level);
			sum = 0.0;
		}
		// create and return the derivative BezCurv
		return temp;
	}
}


// evaluate the derivative of the BezCurv of order deriv at a point
// x. Computes the derivative as a BezCurv and then evaluates this at x
template<class T>
T BezCurv<T>::Derive(int level, double x) const
{
	return Derive(level)(x);
}


// INTEGRATION

// integrate the BezCurv between the limits x1 and x2. Computes
// the indefinite integral as a BezCurv and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T BezCurv<T>::Integrate(double x1, double x2) const
{
	// create the indefinite integral;
	BezCurv<T> intCurve = Integrate(); 
	KnotSet kset = intCurve.GetKnotSet();

	if (x2 < kset.GetKnots()[ord] || x1 > kset.GetKnots()[ord+1]) return 0;

	if (x1 < kset.GetKnots()[ord]) x1 = kset.GetKnots()[ord];
	if (x2 > kset.GetKnots()[ord+1]) x2 = kset.GetKnots()[ord+1];

	// evaluate and subtract
	return (intCurve(x2) - intCurve(x1));
}


// compute the indefinite integral of the BezCurv and represent
// it as a BezCurv of one degree higher
template<class T>
BezCurv<T> BezCurv<T>::Integrate() const
{
	// create knot set for integral
	KnotSet nkset = (*kset).CreateKnotSetIntegrate();

	// create the Curve
	return BezCurv<T>(IntegrateCPoints(),nkset.GetKnots(),nkset.GetOrd());
}	


// compute the indefinite integral of the BezCurv as a BezCurv
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Vector<T> BezCurv<T>::IntegrateCPoints() const
{
	return (CPoints<T>(*cpts,*kset,ord).CreateCPointsIntegrate()).GetCPoints();
}	


// PRODUCT


// compute the product of the BezCurv with another BezCurv and 
// represent the result as a new BezCurv 
/*template <class T>
BezCurv<T> BezCurv<T>::Product(const BezCurv<T>& b) const
{
	// Create the product knot set
	KnotSet nkset = (*kset).CreateKnotSetProduct(KnotSet(b.GetKnots(),b.GetOrd(),2*b.GetOrd()));

	// create the new BezCurv Curve
	return BezCurv<T>(ProductCPoints(b),nkset.GetKnots(),nkset.GetOrd());
}   


// compute the product of the BezCurv with another BezCurv and 
// represent the result as a new BezCurv 
template<class T>
Vector<T> BezCurv<T>::ProductCPoints(const BezCurv<T>& b) const
{
	int ordb = b.GetOrd();

	// create the control points
	// algorithm for product from Piegl/Tiller Paper
	Vector<T> temp(ord+ordb-1);
	T sum = 0.0;
	for (int k=0; k<ord+ordb-1; k++) {
		for (int l=Math::max1(0,k-ordb+1); l<=Math::min1(ord-1,k); l++) {
			sum = sum + (Math::Combin(ord-1,l)*Math::Combin(ordb-1,k-l)/Math::Combin(ord+ordb-2,k))*(*cpts)[l]*b.GetCPoints()[k-l];
		}
		temp[k] = sum;
		sum=0.0;
	}
	// create the new BezCurv Curve
	return temp;
} */  


// KNOT INSERTION

 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspCurv<T> BezCurv<T>::InsertKnot(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	return BspCurv<T>(*cpts,*kts,ord,ord).InsertKnot(t,mult,n);
}

// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Vector<T> BezCurv<T>::InsertKnotCPoints(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	return BspCurv<T>(*cpts,*kts,ord,ord).InsertKnotCPoints(t,mult,n);
}

// SUBDIVISION      
     
// subdivide the BezCurv between the limits x1 and x2
// and return the new BezCurv
template<class T>
BezCurv<T> BezCurv<T>::Subdivide(double x1, double x2) const
{
	// create knot set for new Curve	
	// create subdivision knot set
//	if (x1 >= leftLimit && x2 <= rightLimit) {
		KnotSet nkset = (*kset).CreateKnotSetSubdivide(x1,x2);

	// create and return the subdivided BezCurv 
		return BezCurv<T>(SubdivideCPoints(x1,x2),nkset.GetKnots(),ord);
//	} else return *this;
}

// subdivide the BezCurv between the limits x1 and x2
// and return the new BezCurv
template<class T>
Vector<T> BezCurv<T>::SubdivideCPoints(double x1, double x2) const
{
	// create control points for new Curve
//	if (x1 >= leftLimit && x2 <= rightLimit) {
		return (CPoints<T>(*cpts,*kset,ord).CreateCPointsSubdivide(x1,x2)).GetCPoints();
//	} else return *cpts;
}

// subdivide the BezCurv between the limits x1 and x2
// and return the new BezCurv
template<class T>
Vector<T> BezCurv<T>::SubdivideCPoints(int level) const
{
	// create control points for new Curve
	return BspCurv<T>(*cpts,*kts,ord,ord).SubdivideCPoints(level);
}

// subdivide the BezCurv between the limits x1 and x2
// and return the new BezCurv
template<class T>
BspCurv<T> BezCurv<T>::Subdivide(int level) const
{
	// create control points for new Curve
	return BspCurv<T>(*cpts,*kts,ord,ord).Subdivide(level);
}


// DERIVATIVE MATRIX

// find the derivative of the BezCurv of level given by deriv and
// represent the resulting control points of the derivative BezCurv
// in terms of the original control points of the Curve. Returns
// the matrix where each row give each new control point as a linear
// combination of the old
template<class T>
Matrix<T> BezCurv<T>::Dmatrix(int deriv) const
{
	// create matrix for first derivative
	if (deriv <= 0) return Math::CreateIdentity(ord);
	else if (deriv >= ord) return Matrix<T>(4,4,0.0);
	
	Matrix<double> d = (*kset).CreateMatrixDeriv();

	KnotSet nkset(*kset);
	// create knot sets for each subsequent derivative up to deriv
	// find new matrix and multiply by previous matrix.
	for (int i=1; i<deriv; i++) {
		nkset = nkset.CreateKnotSetDeriv(1);

		d = nkset.CreateMatrixDeriv()*d;
	}	
	// return matrix of linear combinations
	return d;
}


// CONVERSION

// Convert the BezCurv to a PolyCurv
// use matrix conversion
template<class T>
PolyCurv<T> BezCurv<T>::ConvertPolyCurv() const
{
	return PolyCurv<T>(Math::mult3(Math::ComputePolyCurvMatrix(ord),*cpts),ord,0.0,1.0);
}


// Convert BezCurv to BspCurv
// nothing to do just include number of control point
template<class T>
BspCurv<T> BezCurv<T>::ConvertBspCurv() const
{
	return BspCurv<T>(*cpts, *kts, ord, ord);
}

// READ and WRITE

template <class T>
void BezCurv<T>::write(std::ostream& os) 
{
	os << "Bezier Curve\n";
	os << "order is " << ord << "\n";
	os << "knots are\n";
	os << *kts;
	os << "\ncontrol points are\n";
	os << *cpts;
}
	
	
template <class T>
void BezCurv<T>::read(std::istream& is)
{
	int Ord;
	std::cout << "order of Bezier curve\n";
	is >> Ord;
	Vector<T> Cpts(Ord);
	std::cout << "input control points\n";
	is >> Cpts;
	*this = BezCurv<T>(Cpts,Ord);
} 


template <class T>
void BezCurv<T>::writefile(std::ofstream& ofs)
{
	ofs << "Bezier Curve\n";
	ofs << "order is " << ord << "\n";
	ofs << "knots are\n";
	ofs << *kts;
	ofs << "\ncontrol points are\n";
	ofs << *cpts;
}

template <class T>
void BezCurv<T>::readfile(std::ifstream& ifs)
{
	int Ord;
	ifs >> Ord;
	Vector<T> Cpts(Ord);
	ifs >> Cpts;
	*this = BezCurv<T>(Cpts,Ord);
} 


// PRIVATE FUNCTIONS

         
// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
BspCurv<T> BezCurv<T>::InsertKnot(double x) const
{
	return BSpCurv<T>(*cpts,*kts,ord,ord).InsertKnot(x);
}                              
      
// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
Vector<T> BezCurv<T>::InsertKnotCPoints(double x) const
{       
	// calculate the new control points
	// also updates the knot set
	return BspCurv<T>(*cpts,*kts,ord,ord).InsertKnotCPoints(x);
}              


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
BspCurv<T> BezCurv<T>::InsertKnot(double x, int level) const
{	
	return BspCurv<T>(*cpts,*kts,ord,ord).InsertKnot(x,level);
}


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
Vector<T> BezCurv<T>::InsertKnotCPoints(double x, int level) const
{	
	return BspCurv<T>(*cpts,*kts,ord,ord).InsertKnotCPoints(x,level);
}

// inserts a vector of n knots into the BspCurv and returns the new
// BspCurv
template<class T>
BspCurv<T> BezCurv<T>::InsertKnot(const Vector<double>& Kts, int n) const
{
	return BspCurv<T>(*cpts,*kts,ord,ord).InsertKnot(Kts,n);
}
   

#endif

