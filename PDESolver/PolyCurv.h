
#ifndef POLYCURV
#define POLYCURV

#include "matrix.h"
#include "bspcurv.h"

class PolyCurvBasisFunc : public Curve<double> {

	// data
	int ord;	 // order of basis function
	double leftlimit;
	double rightlimit;
	int index;
	Ptr<PolyCurv<double> > b;	// PolyCurv representation of basis function

	// private functions
	PolyCurv<double> CreatePolyCurv() const;	// creates the PolyCurv
	virtual ObjectID Identity() const { return std::string("class PolyCurvBasisFunc"); }
public:
	// constructors
	PolyCurvBasisFunc();
	PolyCurvBasisFunc(int Ord, int Index, double LeftLimit=0.0, double RightLimit=1.0);

	// evaluators
	double Eval(double x) const;
	virtual double operator()(double x) const;
	virtual double operator()(int, double x) const;


	// access functions
	PolyCurv<double> GetPolyCurv() const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	
	virtual double Derive(int, double) const;
	double GetLeftLimit() const;
	double GetRightLimit() const;
};


class PolyCurvBasisFuncSet : public TextObject, public FTextObject {

	// data
	int ord;	 // order of basis function
	double leftlimit;
	double rightlimit;
	Ptr<Vector<PolyCurvBasisFunc> > b;	// PolyCurv representation of basis function
	virtual ObjectID Identity() const { return std::string("class PolyCurvBasisFuncSet"); }
public:
	// constructors
	PolyCurvBasisFuncSet();
	PolyCurvBasisFuncSet(int Ord, double LeftLimit, double RightLimit);

	// evaluators
	Vector<double> Eval(double x) const;
	Vector<double> operator()(double x) const;

	// read and write
	void read(std::istream& is);
	void readfile(std::ifstream& ifs);
	void write(std::ostream& os);
	void writefile(std::ofstream& ofs);
};



template<class T>
class PolyCurv : public Curve<T> {
private:
	// data
	int ord;  // order 
	Ptr<Vector<T> > coeffs;  // coeffs
	double t1, t2; // range of curve
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class PolyCurv<double>");
		else {
			std::string s(typeid(T).name()), s1("class PolyCurv<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:

	// constructors
	PolyCurv();
	PolyCurv(const Vector<T>& Coeffs, int Ord, double t1=0.0, double t2=1.0);
	PolyCurv(int Ord, double t1=0.0, double t2=1.0);

	// access functions
	int GetOrd() const;
	int GetNum() const;
	double GetLeftLimit() const;
	double GetRightLimit() const;
	Vector<T> GetCoeffs() const;

	// evaluators
	virtual T operator()(double x) const;
	virtual T operator()(int, double x) const;
	T Eval(double x) const;
	Vector<T> ComputePoints(int n) const;

	// derivatives
	PolyCurv<T> Derive() const;
	PolyCurv<T> Derive(int level) const; 
	Vector<T> DeriveCPoints(int level) const;
	virtual T Derive(int level, double x) const;
	Vector<T> Derive2(double x) const;

	// degree elevation
	PolyCurv<T> Elevate(int level) const;
	Vector<T> ElevateCPoints(int level) const;

	// integration
	PolyCurv<T> Integrate() const;
	Vector<T> IntegrateCPoints() const;
	double Integrate(double x1, double x2) const;

	// reparameterisation
	PolyCurv<T> Reparameterise1() const;
	Vector<T> Reparameterise1CPoints() const;
	PolyCurv<T> Reparameterise2(double T1, double T2) const;
	PolyCurv<T> Reparameterise3(double T1, double T2) const;
	PolyCurv<T> Reparameterise4(double T1, double T2) const;


	// product
	template<class T1>
	PolyCurv<T> Product(const PolyCurv<T1>& b) const
	{
		// create the new PolyCurv Curve
		// should be over the same range, reparameterise first!
		return PolyCurv<T>(ProductCPoints(b),ord+b.GetOrd()-1,b.GetLeftLimit(),b.GetRightLimit());
	}
	template<class T1>
	Vector<T> ProductCPoints(const PolyCurv<T1>& b) const
	{
		PolyCurv<T> b1 = Reparameterise3(b.GetLeftLimit(),b.GetRightLimit());
		int ord1 = b.GetOrd();

		// create the new coeffs
		Vector<T> temp(ord+ord1-1);
		T sum = 0.0;
		for (int k=0; k<ord+ord1-1; k++) {
			for (int l=Math::max1(0,k-ord1+1); l<=Math::min1(ord-1,k); l++)
				sum = sum + (*coeffs)[l]*b.GetCoeffs()[k-l];
			temp[k] = sum;
			sum=0.0;
		}
		// create the new PolyCurv Curve
		return temp;
	}

	// subdivision
	PolyCurv<T> Subdivide(double x1, double x2) const;
	Vector<T> SubdivideCPoints(double x1, double x2) const;

	// matrix derivative
	Matrix<double> CreateMatrixDeriv() const;
	Matrix<T> Dmatrix(int deriv) const;

	// conversion
	BezCurv<T> ConvertBezCurv() const;
	Vector<T> ConvertBezCurvCPoints() const;
	BspCurv<T> ConvertBspCurv() const;

	// add and subtract
	PolyCurv<T> Add(const PolyCurv<T>& p) const;
	PolyCurv<T> Subtract(const PolyCurv<T>& p) const;
	
	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

// CONSTRUCTORS

// constructor builds a PolyCurv from a Vector of control points, knots
// an order and number of control points
template<class T>
PolyCurv<T>::PolyCurv(const Vector<T>& Coeffs, int Ord, double T1, double T2) :
			ord(Ord), coeffs(new Vector<T>(Coeffs)), t1(T1), t2(T2)
{
}

// constructor taking just ord
template<class T>
PolyCurv<T>::PolyCurv(int Ord, double T1, double T2) : ord(ord), coeffs(new Vector<T>(Ord,0.0)), t1(T1), t2(T2) 
{ 
}

// default constructor
template<class T>
PolyCurv<T>::PolyCurv() : ord(0), coeffs(), t1(0.0), t2(0.0) { }

// ACCESS FUNCTIONS

// get the order of the PolyCurv
template<class T>
inline int PolyCurv<T>::GetOrd() const { return ord; }

// get the number of coeffs
template<class T>
inline int PolyCurv<T>::GetNum() const { return ord; }

// get the control points
template<class T>
inline Vector<T> PolyCurv<T>::GetCoeffs() const { return *coeffs; }


// get the left limit
template<class T>
inline double PolyCurv<T>::GetLeftLimit() const { return t1; }

// get the left limit
template<class T>
inline double PolyCurv<T>::GetRightLimit() const { return t2; }


// EVALUATORS

// evaluate the PolyCurv at the point x using Nested multiplication
template<class T>
T PolyCurv<T>::operator()(double x) const
{
	T b = (*coeffs)[ord-1];
	for (int i=ord-2; i>=0; i--) b = (*coeffs)[i] + x*b;
	return b;
}


// evaluate the PolyCurv using standard method rather than overloaded 
// operator identical code to that of overloaded operator
template<class T>
T PolyCurv<T>::Eval(double x) const
{
	T b = (*coeffs)[ord-1];
	for (int i=ord-2; i>=0; i--) b = (*coeffs)[i] + x*b;
	return b;
}


// REPARAMETERISATION

// reparameterise from (t1,t2) to (0,1)
// formula in Sherar/Goult book
template<class T>
PolyCurv<T> PolyCurv<T>::Reparameterise1() const
{
	Vector<T> ncoeffs(ord);
	T sum=0.0;

	for (int i=0; i<ord; i++) {
		for (int j=0; j<ord-i; j++) sum = sum+(*coeffs)[i+j]*Math::Combin(i+j,i)*pow(t1,(double)j);
		sum = sum*pow(t2-t1,(double)i);
		ncoeffs[i]=sum;
		sum=0.0;
	}
	return PolyCurv<T>(ncoeffs,ord,0.0,1.0);
}

// reparameterise from (t1,t2) to (0,1)
// formula in Sherar/Goult book
template<class T>
Vector<T> PolyCurv<T>::Reparameterise1CPoints() const
{
	Vector<T> ncoeffs(ord);
	T sum=0.0;

	for (int i=0; i<ord; i++) {
		for (int j=0; j<ord-i; j++) sum = sum+(*coeffs)[i+j]*Math::Combin(i+j,i)*pow(t1,(double)j)*pow(t2-t1,(double)i);
		ncoeffs[i]=sum;
		sum=0.0;
	}
	return ncoeffs;
}


// reparameterise from (0,1) to (T1, T2)
// formula in Sherar/Goult book
template<class T>
PolyCurv<T> PolyCurv<T>::Reparameterise2(double T1, double T2) const
{
	Vector<T> ncoeffs(ord);
	T sum=0.0;

	for (int i=0; i<ord; i++) {
		for (int j=0; j<ord-i; j++) sum = sum+(*coeffs)[i+j]*Math::Combin(i+j,i)*pow(-T1,(double)j)/pow((T2-T1),(double)(i+j));
		ncoeffs[i]=sum;
		sum=0.0;
	}
	return PolyCurv<T>(ncoeffs,ord,T1,T2);
}

// reparameterise from (t1,t2) to (T1, T2)
template<class T>
PolyCurv<T> PolyCurv<T>::Reparameterise3(double T1, double T2) const
{
	// go from t1,t2 to 0,1 then 0,1 to T1, T2
	return (Reparameterise1().Reparameterise2(T1,T2));
}

template<class T>
PolyCurv<T> PolyCurv<T>::Reparameterise4(double T1, double T2) const
{
	T sum = 0.0;
	Vector<T> ncoeffs(ord);
	for (int i=0; i<ord; i++) {
		for (int j=i; j<ord; j++) 
			sum = sum + (*coeffs)[j]*pow(T1,(double)(j-i))*pow((T2-T1),(double)(i))*Math::Combin(j,i);
		ncoeffs[i] = sum;
		sum = 0.0;
	}
	return PolyCurv<T>(ncoeffs,ord,T1,T2);
}


// DEGREE ELEVATION

// elevate the degree of the PolyCurv by level
template<class T>
PolyCurv<T> PolyCurv<T>::Elevate(int level) const
{
	if (level <= 0) return *this; 
	Vector<T> ncoeff(ord+level);
	for (int i=0; i<ord; i++) ncoeff[i]=(*coeffs)[i];
	// set max coeff to 0
	for (int i=ord; i<ord+level; i++) ncoeff[i] = 0.0;
	return PolyCurv<T>(ncoeff,ord+level,t1,t2);
}


// elevate the degree of the PolyCurv by level
template<class T>
Vector<T> PolyCurv<T>::ElevateCPoints(int level) const
{
	if (level <= 0) return *coeffs;
	Vector<T> ncoeff(ord+level);
	for (int i=0; i<ord; i++) ncoeff[i]=(*coeffs)[i];
	// set max coeff to 0
	for (int i=ord; i<ord+level; i++) ncoeff[i] = 0.0;
	return ncoeff;
}


// DERIVATIVES

// compute the derivative of the PolyCurv of order deriv and
// represent the result as another PolyCurv
template<class T>
PolyCurv<T> PolyCurv<T>::Derive() const
{   
	if (ord <= 0) return *this;
	return PolyCurv<T>(DeriveCPoints(1), ord-1, t1, t2);
}



// compute the derivative of the PolyCurv of order deriv and
// represent the result as another PolyCurv
template<class T>
PolyCurv<T> PolyCurv<T>::Derive(int level) const
{   
	if (level <= 0) return *this;
	return PolyCurv<T>(DeriveCPoints(level), ord-level, t1, t2);
}

// compute the derivative of the PolyCurv of order deriv and
// represent the result as another PolyCurv
template<class T>
Vector<T> PolyCurv<T>::DeriveCPoints(int level) const
{   
	if (level <= 0) return *coeffs;
	Vector<T> ncoeff(ord-level);
	for (int i=level; i<ord; i++)
		ncoeff[i-level] = (*coeffs)[i]*Math::Fact(i)/Math::Fact(i-level);
	
	return ncoeff;
}


// to finish
template<class T>
Vector<T> PolyCurv<T>::Derive2(double x) const
{
	Vector<T> ncoeff(ord);
	ncoeff[0] = (*this)(x);
	for (int j=1; j<=ord-1; j++)
		for (int i=level; i<ord; i++)
			ncoeff[i-level] = (*coeffs)[i]*Math::Fact(i)/Math::Fact(i-level);
	
	return ncoeff;
}	


template<class T>
T PolyCurv<T>::operator()(int lev, double x) const
{
	return Derive(lev)(x);
}

// evaluate the derivative of the PolyCurv of order deriv at a point
// x. Computes the derivative as a PolyCurv and then evaluates this at x
template<class T>
T PolyCurv<T>::Derive(int level, double x) const
{
	return Derive(level)(x);
}


// INTEGRATION

// integrate the PolyCurv between the limits x1 and x2. Computes
// the indefinite integral as a PolyCurv and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
double PolyCurv<T>::Integrate(double x1, double x2) const
{
	// create the indefinite integral;

	if (x2 < t1 || x1 > t2) return 0.0;

	if (x1 < t1) x1 = t1;
	if (x2 > t2) x2 = t2;

	PolyCurv<T> intCurve = Integrate(); 
	
		// evaluate and subtract
	return (intCurve(x2) - intCurve(x1));
}


// compute the indefinite integral of the PolyCurv and represent
// it as a PolyCurv of one higher degree
template<class T>
PolyCurv<T> PolyCurv<T>::Integrate() const
{
	return PolyCurv<T>(IntegrateCPoints(),ord+1,t1,t2);
}	


// compute the indefinite integral of the PolyCurv as a PolyCurv
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Vector<T> PolyCurv<T>::IntegrateCPoints() const
{
	Vector<T> temp(ord+1,0.0);

	for (int i=1; i<ord+1; i++) temp[i] = (*coeffs)[i-1]/(double)i;
	return temp;
}

/*
// PRODUCT

// compute the product of the PolyCurv with another PolyCurv and 
// represent the result as a new PolyCurv 
template<class T>  
PolyCurv<T> PolyCurv<T>::Product(const PolyCurv<T>& b) const
{
	// create the new PolyCurv Curve
	// should be over the same range, reparameterise first!
	return PolyCurv<T>(ProductCPoints(b),ord+b.GetOrd()-1,b.GetLeftLimit(),b.GetRightLimit());
}   

// compute the product of the PolyCurv with another PolyCurv and 
// represent the result as a new PolyCurv 
template<class T>  
Vector<T> PolyCurv<T>::ProductCPoints(const PolyCurv<T>& b) const 
{
	PolyCurv<T> b1 = Reparameterise3(b.GetLeftLimit(),b.GetRightLimit());
	int ord1 = b.GetOrd();

	// create the new coeffs
	Vector<T> temp(ord+ord1-1);
	T sum = 0.0;
	for (int k=0; k<ord+ord1-1; k++) {
		for (int l=Math::max1(0,k-ord1+1); l<=Math::min1(ord-1,k); l++)
			sum = sum + (*coeffs)[l]*b.GetCoeffs()[k-l];
		temp[k] = sum;
		sum=0.0;
	}
	// create the new PolyCurv Curve
	return temp;
}   
*/

// ADD and SUBTRACT

// add two BezCurv's together, should have the same knot range
// addtion has knot vector from current object
template<class T>
PolyCurv<T> PolyCurv<T>::Add(const PolyCurv<T>& b) const
{
	PolyCurv<T> b1 = *this;//Reparameterise3(b.GetLeftLimit(),b.GetRightLimit());
	PolyCurv<T> c(b1), d(b);

	if (ord > b.GetOrd()) d = d.Elevate(ord-b.GetOrd());
	else if (b.GetOrd() > ord) c = c.Elevate(b.GetOrd()-ord); 

	Vector<T> temp(c.GetOrd());	
	// add the control points together
	for (int i=0; i<c.GetOrd(); i++)
		temp[i]=c.GetCoeffs()[i]+d.GetCoeffs()[i];

	return PolyCurv<T>(temp,c.GetOrd(),t1,t2);
}


// subtract two BezCurv's together (*this-b), 
// should have the same knot range
// subtraction has knot vector from current object
template<class T>
PolyCurv<T> PolyCurv<T>::Subtract(const PolyCurv<T>& b) const
{
	PolyCurv<T> b1 = *this;//Reparameterise3(b.GetLeftLimit(),b.GetRightLimit());
	PolyCurv<T> c(b1), d(b);

	if (ord > b.GetOrd()) d = d.Elevate(ord-b.GetOrd());
	else if (b.GetOrd() > ord) c = c.Elevate(b.GetOrd()-ord); 

	Vector<T> temp(c.GetOrd());	
	// add the control points together
	for (int i=0; i<c.GetOrd(); i++)
		temp[i]=c.GetCoeffs()[i]-d.GetCoeffs()[i];

	return PolyCurv<T>(temp,c.GetOrd(),t1,t2);
}


// SUBDIVISION

// subdivide the PolyCurv between the limits x1 and x2
// and return the new PolyCurv
template<class T>
PolyCurv<T> PolyCurv<T>::Subdivide(double x1, double x2) const
{
	// subdivide and represent segment over 0,1
	if (x1 >= t1 && x2 <= t2) return (PolyCurv<T>(*coeffs, ord, x1, x2).Reparameterise1());
	else return *this;
}

// subdivide the PolyCurv between the limits x1 and x2
// and return the new PolyCurv
template<class T>
Vector<T> PolyCurv<T>::SubdivideCPoints(double x1, double x2) const
{
	// subdivide and represent segment over 0,1
	if (x1 >= t1 && x2 <= t2) return (PolyCurv<T>(*coeffs, ord, x1, x2).Reparameterise1CPoints());
	else return *coeffs;
}


// MATRIX DERIVATIVE

// find the derivative of the PolyCurv of level given by deriv and
// represent the resulting control points of the derivative PolyCurv
// in terms of the original control points of the Curve. Returns
// the matrix where each row give each new control point as a linear
// combination of the old
template<class T>
Matrix<T> PolyCurv<T>::Dmatrix(int deriv) const
{
	// create matrix for first derivative
	Matrix<double> d = CreateMatrixDeriv();

	// find new matrix and multiply by previous matrix.
	PolyCurv<T> p(*this);
	for (int i=1; i<deriv; i++) {
		p = p.Derive(1);
		d = p.CreateMatrixDeriv()*d;	
	}
	// return matrix of linear combinations
	return d;
}

// create matrix representing derivative of PolyCurv
template<class T>
Matrix<double> PolyCurv<T>::CreateMatrixDeriv() const
{
	Matrix<double> res(ord-1,ord,0.0);

	for (int i=0; i<ord-1; i++) 
		res[i][i+1] = (double)(i+1); 
	return res;
}


// CONVERSION

// computes the BezCurv form of the PolyCurv
template<class T>
BezCurv<T> PolyCurv<T>::ConvertBezCurv() const
{
	// convert to 0,1 form
	PolyCurv<T> p = Reparameterise1();
	return BezCurv<T>(p.ConvertBezCurvCPoints(), ord, t1, t2);
}



// computes the BezCurv form of the PolyCurv
template<class T>
Vector<T> PolyCurv<T>::ConvertBezCurvCPoints() const
{
	// convert to 0,1 form
//	Reparameterise1();
	return Math::mult3(Math::ComputePolyCurvMatrixInverse(ord),*coeffs);
}


// computes the BezCurv form of the PolyCurv
template<class T>
BspCurv<T> PolyCurv<T>::ConvertBspCurv() const
{
	return ConvertBezCurv().ConvertBspCurv();
}

// READ and WRITE

template <class T>
void PolyCurv<T>::write(std::ostream& os)
{
	os << "Poly Curve\n";
	os << "order is " << ord << "\n";
	os << "limits are\n";
	os << GetLeftLimit() << " " << GetRightLimit();
	os << "\ncoefficients are\n";
	(*coeffs).write(os);
}
	
	
template <class T>
void PolyCurv<T>::read(std::istream& is)
{
	int Ord;
	std::cout << "order of Poly curve\n";
	is >> Ord;
	double left, right;
	std::cout << "input limits\n";
	is >> left >> right;
	Vector<T> Coeffs(Ord);
	std::cout << "input coeffs\n";
	Coeffs.read(is);
	*this = PolyCurv<T>(Coeffs,Ord,left,right);
} 


template <class T>
void PolyCurv<T>::writefile(std::ofstream& ofs) 
{
	ofs << "Poly Curve\n";
	ofs << "order is " << ord << "\n";
	ofs << "limits are\n";
	ofs << GetLeftLimit() << " " << GetRightLimit();
	ofs << "\ncoefficients are\n";
	(*coeffs).writefile(ofs);
}

template <class T>
void PolyCurv<T>::readfile(std::ifstream& ifs)
{
	int Ord;
	ifs >> Ord;
	double left, right;
	ifs >> left >> right;
	Vector<T> Coeffs(Ord);
	Coeffs.readfile(ifs);
	*this = PolyCurv<T>(Coeffs,Ord,left,right);
} 


#endif  
 