#ifndef POLYSURF
#define POLYSURF


#include "MathFunctions.h"
#include "Matrix.h"
#include "BezSurf.h"
#include "PolyCurv.h"


class PolySurfBasisFunc : public Surf<double> {

	// data
	int ordu;	 // order of basis function
	int ordv;	// order in v
	int indexu;
	int indexv;
	double leftlimitu;  // ranges in u
	double rightlimitu;
	double leftlimitv; // rnages in v
	double rightlimitv;
	Ptr<PolySurf<double> > b;	// Polysurf representation of basis function

	// private functions
	PolySurf<double> CreatePolySurf() const;	// creates the Polysurf
	virtual ObjectID Identity() const { return std::string("class PolySurfBasisFunc"); }
public:
	// constructors
	PolySurfBasisFunc();
	PolySurfBasisFunc(int Ordu, int Ordv, int Indexu, int Indexv, double LeftLimitU=0.0, double RightLimitU=1.0, double LeftLimitV=0.0, double RightLimitV=1.0);

	// evaluators
	double Eval(double u, double v) const;
	virtual double operator()(double u, double v) const;
	virtual double operator()(int, int, double u, double v) const;
	
	// access fuynctions
	PolySurf<double> GetPolySurf() const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);

	virtual double Derive(int, int , double, double) const;
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
};

class PolySurfBasisFuncSet : public TextObject, public FTextObject {

	// data
	int ordu;	 // order of basis function
	int ordv;
	double leftlimitu;  // ranges in u
	double rightlimitu;
	double leftlimitv; // ranges in v
	double rightlimitv;

	Ptr<Matrix<PolySurfBasisFunc> > b;	// PolySurf representation of basis function
	virtual ObjectID Identity() const { return std::string("class PolySurfBasisFuncSet"); }
public:
	// constructors
	PolySurfBasisFuncSet();
	PolySurfBasisFuncSet(int Ordu, int Ordv, double LU, double LV, double RU, double RV);
	
	// evaluators
	Matrix<double> Eval(double u, double v) const;
	Matrix<double> operator()(double u, double v) const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};



template<class T>
class PolySurf : public Surf<T> {
private:
	// data
	int ordu;  // order in u
	int ordv;	// order in v
	Ptr<Matrix<T> > coeffs;  // coeffs
	double u1, u2;  // limits in u
	double v1, v2;  // limits in v

	// private functions

	T DeriveU(int levu, double u, double v) const;
	T DeriveV(int levv, double u, double v) const;
	
	Matrix<T> ElevateCPointsU(int level) const;
	Matrix<T> ElevateCPointsV(int level) const;
	T IntegrateUV2(double u1, double u2, double v1, double v2) const;
	PolySurf<T> IntegrateU() const;
	PolySurf<T> IntegrateV() const;
	Matrix<T> IntegrateCPointsU() const;
	Matrix<T> IntegrateCPointsV() const;
	PolySurf<T> Product2(const PolySurf<T>& c) const;
	PolySurf<T> Product1(const PolySurf<T>& b) const;
	PolySurf<T> SubdivideU(double u1, double u2) const;
	PolySurf<T> SubdivideV(double v1, double v2) const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class PolySurf<double>");
		else {
			std::string s(typeid(T).name()), s1("class PolySurf<");
			return s1.append(s,0,s.size())+">";
		}
		/*if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class PolySurfDouble"); 
		else if (std::string(typeid(T).name()) == std::string("Point1D"))
			return std::string("class FPolySurf");
		else return std::string("class PolySurf3D"); */
	}
public:

	// constructors
	PolySurf();
	PolySurf(const Matrix<T>& Coeffs, int Ordu, int Ordv, double U1=0.0, double U2=1.0, double V1=0.0, double V2=1.0);
	PolySurf(int Ordu, int Ordv);

	// access functions
	int GetOrdU() const;
	int GetOrdV() const;
	Matrix<T> GetCoeffs() const;
	virtual double GetRightLimitU() const;
	virtual double GetRightLimitV() const;
	virtual double GetLeftLimitU() const;
	virtual double GetLeftLimitV() const;

	// isoparametric curves
	PolyCurv<T> GetIsoparametricU(double u) const;
	PolyCurv<T> GetIsoparametricV(double v) const;

	// evaluators
	virtual T operator()(double u, double v) const;
	virtual T operator()(int, int, double u, double v) const;
	T Eval(double u, double v) const;

	// reparameterisation
	PolySurf<T> Reparameterise1() const;
	PolySurf<T> Reparameterise2(double U1, double U2, double V1, double V2) const;
	PolySurf<T> Reparameterise3(double U1, double U2, double V1, double V2) const;

	// addition and subtraction
	PolySurf<T> Add(const PolySurf<T>& b) const;
	PolySurf<T> Subtract(const PolySurf<T>& b) const;

	// derivatives
	PolySurf<T> Derive(int levu, int levv) const; 
	virtual T Derive(int levu, int levv, double u, double v) const;
	PolySurf<T> DeriveU(int levu) const; 
	PolySurf<T> DeriveV(int levv) const; 
	Matrix<T> DeriveCPointsU(int level) const;
	Matrix<T> DeriveCPointsV(int level) const;

	// degree elevation
	PolySurf<T> ElevateU(int levu) const;
	PolySurf<T> ElevateV(int levu) const;
	PolySurf<T> Elevate(int levu, int levv) const;

	// integration
	PolySurf<T> Integrate() const;
	T Integrate(double u1, double u2, double v1, double v2) const;
	Matrix<T> IntegrateCPoints() const;

	// product
	template<class T1>
	PolySurf<T> Product(const PolySurf<T1>& b) const
	{
		return PolySurf<T>(ProductCPoints(b),ordu+b.GetOrdU()-1,ordv+b.GetOrdV()-1,b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV());
	}

	template<class T1>
	Matrix<T> ProductCPoints(const PolySurf<T1>& b) const
	{
		PolySurf<T> b1 = Reparameterise3(b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV());
		int ordbu = b.GetOrdU();
		int ordbv = b.GetOrdV();

		Matrix<T> temp(ordbu+ordu-1,ordbv+ordv-1);
	
		double sum =0.0;
		for (int k=0; k<ordbu+ordu-1; k++) {
			for (int l=0; l<ordbv+ordv-1; l++) {
				for (int i=Math::max1(0,k-ordbu+1); i<=Math::min1(ordu-1,k); i++) {
					for (int j=Math::max1(0,l-ordbv+1); j<=Math::min1(ordv-1,l); j++) {
						sum = sum + (Math::Combin(ordu-1,i)*Math::Combin(ordbu-1,k-i)*Math::Combin(ordv-1,j)*Math::Combin(ordbv-1,l-j)*b1.GetCoeffs()[i][j]*b.GetCoeffs()[k-i][l-j])/(Math::Combin(ordu+ordbu-2,k)*Math::Combin(ordv+ordbv-2,l));
					}
				}
				temp[k][l]=sum;
				sum=0.0;
			}
		}
		return temp;
	}
	
	// subdivision
	PolySurf<T> Subdivide(double u1, double u2, double v1, double v2) const;

	// conversion
	BezSurf<T> ConvertBezSurf() const;
	Matrix<T> ConvertBezSurfCPoints() const;
	BspSurf<T> ConvertBspSurf() const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

// CONSTRUCTORS

// constructor builds a PolySurf from a Vector of control points, knots
// an order and number of control points
template<class T>
PolySurf<T>::PolySurf(const Matrix<T>& Coeffs, int Ordu, int Ordv, double U1, double U2, double V1, double V2) :
			ordu(Ordu), ordv(Ordv), coeffs(new Matrix<T>(Coeffs)), u1(U1), u2(U2), v1(V1), v2(V2)
{
}

// default constructor
template<class T>
PolySurf<T>::PolySurf() : ordu(0), ordv(0), coeffs() { }

// default constructor
template<class T>
PolySurf<T>::PolySurf(int Ordu, int Ordv) : ordu(Ordu), ordv(Ordv), coeffs(new Matrix<T>(ordu,ordv)), u1(0), u2(1), v1(0), v2(1) 
{
}

// ACCESS FUNCTIONS

// get the order of the PolySurf
template<class T>
inline int PolySurf<T>::GetOrdU() const { return ordu; }

// get the order of the PolySurf
template<class T>
inline int PolySurf<T>::GetOrdV() const { return ordv; }


// get the control points
template<class T>
inline Matrix<T> PolySurf<T>::GetCoeffs() const { return *coeffs; }

// get the knot vector
template<class T>
inline double PolySurf<T>::GetLeftLimitU() const { return u1; }

// get the knot vector
template<class T>
inline double PolySurf<T>::GetRightLimitU() const { return u2; }

// get the knot vector
template<class T>
inline double PolySurf<T>::GetLeftLimitV() const { return v1; }

// get the knot vector
template<class T>
inline double PolySurf<T>::GetRightLimitV() const { return v2; }


// ADDITION & SUBTRACTION

// subtract two PolySurf's together
// assuming same range in u and v different orders
template<class T>
PolySurf<T> PolySurf<T>::Subtract(const PolySurf<T>& b) const
{
	PolySurf<T> b1 = Reparameterise3(b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV());
	PolySurf<T> c(b1), d(b);

	if (b.GetOrdU() > ordu) c = c.ElevateU(b.GetOrdU()-ordu);
	else if (ordu > b.GetOrdU()) d = d.ElevateU(ordu-b.GetOrdU());
	if (b.GetOrdV() > ordv) c = c.ElevateV(b.GetOrdV()-ordv);
	else if (ordv > b.GetOrdV()) d = d.ElevateV(ordv-b.GetOrdV());

	Matrix<T> temp(d.GetOrdU(),d.GetOrdV());
	for (int i=0; i<d.GetOrdU(); i++)
		for (int j=0; j<d.GetOrdV(); j++)
			temp[i][j]=c.GetCoeffs()[i][j]-d.GetCoeffs()[i][j];

	return PolySurf<T>(temp,d.GetOrdU(),d.GetOrdV(),u1,u2,v1,v2);
}

// add two PolySurf's together
// assuming same range in u and v different orders
template<class T>
PolySurf<T> PolySurf<T>::Add(const PolySurf<T>& b) const
{
	PolySurf<T> b1 = Reparameterise3(b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV());
	PolySurf<T> c(b1), d(b);

	if (b.GetOrdU() > ordu) c = c.ElevateU(b.GetOrdU()-ordu);
	else if (ordu > b.GetOrdU()) d = d.ElevateU(ordu-b.GetOrdU());
	if (b.GetOrdV() > ordv) c = c.ElevateV(b.GetOrdV()-ordv);
	else if (ordv > b.GetOrdV()) d = d.ElevateV(ordv-b.GetOrdV());

	Matrix<T> temp(d.GetOrdU(),d.GetOrdV());
	for (int i=0; i<d.GetOrdU(); i++)
		for (int j=0; j<d.GetOrdV(); j++)
			temp[i][j]=d.GetCoeffs()[i][j]+c.GetCoeffs()[i][j];

	return PolySurf<T>(temp,d.GetOrdU(),d.GetOrdV(),u1,u2,v1,v2);
}

// REPARAMETERISATION

// reparameterise from (t1,t2) to (0,1)
// formula in Sherar/Goult book
template<class T>
PolySurf<T> PolySurf<T>::Reparameterise1() const
{
	Matrix<T> mat1(ordu,ordv);
	Matrix<T> mat2(ordu,ordv);
	// evaluate the rows

	for (int i=0; i<ordu; i++) {
		PolyCurv<T> b = PolyCurv<T>((*coeffs).GetRow(i),ordv,v1,v2).Reparameterise1();
		for (int j=0; j<ordv; j++) mat1[i][j]=b.GetCoeffs()[j];
	}
	for (int j=0; j<ordv; j++) {
		PolyCurv<T> b = PolyCurv<T>(mat1.GetCol(j),ordu,u1,u2).Reparameterise1();
		for (int i=0; i<ordu; i++) mat2[i][j]=b.GetCoeffs()[i];
	}
	return PolySurf<T>(mat2,ordu,ordv,0.0,1.0,0.0,1.0);
}

// reparameterise from (0,1) to (T1, T2)
// formula in Sherar/Goult book
template<class T>
PolySurf<T> PolySurf<T>::Reparameterise2(double U1, double U2, double V1, double V2) const
{
	Matrix<T> mat1(ordu,ordv);
	Matrix<T> mat2(ordu,ordv);
	// evaluate the rows
	for (int i=0; i<ordu; i++) {
		PolyCurv<T> b = PolyCurv<T>((*coeffs).GetRow(i),ordv,v1,v2).Reparameterise2(V1,V2);
		for (int j=0; j<ordv; j++) mat1[i][j]=b.GetCoeffs()[j];
	}
	for (int j=0; j<ordv; j++) {
		PolyCurv<T> b = PolyCurv<T>(mat1.GetCol(j),ordu,u1,u2).Reparameterise2(U1,U2);
		for (int i=0; i<ordu; i++) mat2[i][j]=b.GetCoeffs()[i];
	}
	return PolySurf<T>(mat2,ordu,ordv,U1,U2,V1,V2);
}

// reparameterise from (t1,t2) to (T1, T2)
template<class T>
PolySurf<T> PolySurf<T>::Reparameterise3(double U1, double U2, double V1, double V2) const
{
	// go from t1,t2 to 0,1 then 0,1 to T1, T2
	return Reparameterise1().Reparameterise2(U1,U2,V1,V2);
}


// EVALUATORS

// evaluate the PolySurf at the point x using de Boor algorithm
template<class T>
T PolySurf<T>::operator()(double u, double v) const
{
  	Vector<T> temp(ordu);
	// evaluate the rows
	for (int i=0; i<ordu; i++) temp[i]= PolyCurv<T>((*coeffs).GetRow(i),ordv,v1,v2)(v);
	// evaluate resulting column
	return PolyCurv<T>(temp,ordu,u1,u2)(u);
}


// evaluate the PolySurf using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T PolySurf<T>::Eval(double u, double v) const
{
  	Vector<T> temp(ordv);

	// evaluate rows
	for (int j=0; j<ordv; j++) temp[j]= PolyCurv<T>((*coeffs).GetCol(j),ordu,u1,u2).Eval(u);
	// evaluate resulting columns
	return PolyCurv<T>(temp,ordv,v1,v2).Eval(v);
}

// ISOPARAMETRIC CURVES

template<class T>
PolyCurv<T> PolySurf<T>::GetIsoparametricV(double v) const
{
	Vector<T> ncpts(ordu);
	// evaluate the v curves at v
	for (int i=0; i<ordu; i++)
		ncpts[i] = PolyCurv<T>((*coeffs).GetRow(i),ordv).Eval(v);
	// return the u curve
	return PolyCurv<T>(ncpts,ordu);
}


template<class T>
PolyCurv<T> PolySurf<T>::GetIsoparametricU(double u) const
{
	Vector<T> ncpts(ordv);
	// evaluate the v curves at v
	for (int i=0; i<ordv; i++)
		ncpts[i] = PolyCurv<T>((*coeffs).GetCol(i),ordu).Eval(u);
	// return the u curve
	return PolyCurv<T>(ncpts,ordv);
}


// DEGREE ELEVATION


// elevate the degree of the PolySurf by level
template<class T>
PolySurf<T> PolySurf<T>::Elevate(int levu, int levv) const
{
	return ElevateU(levu).ElevateV(levv);
}


// DERIVATIVES

// compute the derivative of the PolySurf of order deriv and
// represent the result as another PolySurf
template<class T>
PolySurf<T> PolySurf<T>::Derive(int levu, int levv) const
{
	return DeriveU(levu).DeriveV(levv);
}

template<class T>
T PolySurf<T>::operator() (int valu, int valv, double u, double v) const
{
	return Derive(valu,valv,u,v);
}


// evaluate the derivative of the PolySurf of order deriv at a point
// x. Computes the derivative as a PolySurf and then evaluates this at x
template<class T>
T PolySurf<T>::Derive(int levu, int levv, double u, double v) const
{
	return (DeriveU(levu).DeriveV(levv))(u,v);
}


// INTEGRATION

// integrate the PolySurf between the limits x1 and x2. Computes
// the indefinite integral as a PolySurf and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T PolySurf<T>::Integrate(double U1, double U2, double V1, double V2) const
{
	PolyCurv<T> surf = IntegrateU().IntegrateV();

	return surf(U2,V2)-surf(U1,V2)-surf(U2,V1)+surf(U1,V1);
}

// compute the indefinite integral of the PolySurf and represent
// it as a PolySurf of one higher degree
template<class T>
PolySurf<T> PolySurf<T>::Integrate() const
{
	return IntegrateU().IntegrateV();
}	

// compute the indefinite integral of the PolySurf as a PolySurf
// and return just the control points of this Curvve as a CPoints
// object
template<class T>
Matrix<T> PolySurf<T>::IntegrateCPoints() const
{
	Matrix<T> nceoffs1(ordu,ordv+1);
	Matrix<T> ncoeffs2(ordu+1,ordv+1);
	Vector<T> v;

	for (int i=0; i<ordu; i++) {
		Vector<T> v = PolyCurv<T>(cpts.GetRow(i),ordv,v1,v2).IntegrateCPoints();	
		for (int j=0; j<ordv+1; j++) ncoeffs1[i][j]=v1[j];
	}

	for (int j=0; j<ordv+1; j++) {
		Vector<T> v = PolyCurv<T>(ncoeffs1.GetCol(j),ordu,u1,u2).IntegrateCPoints();		
		for (int i=0; i<ordu+1; i++) ncoeffs2[i][j]=v2[i];
	}
	return ncoeffs2;
}	

/*
// PRODUCT

// compute the product of the PolySurf with another PolySurf and 
// represent the result as a new PolySurf 
template<class T>  
Matrix<T> PolySurf<T>::ProductCPoints(const PolySurf<T>& b) const
{
	PolySurf<T> b1 = Reparameterise3(b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV());
	int ordbu = b.GetOrdU();
	int ordbv = b.GetOrdV();

	Matrix<T> temp(ordbu+ordu-1,ordbv+ordv-1);
	
	double sum =0.0;
	for (int k=0; k<ordbu+ordu-1; k++) {
		for (int l=0; l<ordbv+ordv-1; l++) {
			for (int i=Math::max1(0,k-ordbu+1); i<=Math::min1(ordu-1,k); i++) {
				for (int j=Math::max1(0,l-ordbv+1); j<=Math::min1(ordv-1,l); j++) {
					sum = sum + (Math::Combin(ordu-1,i)*Math::Combin(ordbu-1,k-i)*Math::Combin(ordv-1,j)*Math::Combin(ordbv-1,l-j)*b1.GetCoeffs()[i][j]*b.GetCoeffs()[k-i][l-j])/(Math::Combin(ordu+ordbu-2,k)*Math::Combin(ordv+ordbv-2,l));
				}
			}
			temp[k][l]=sum;
			sum=0.0;
		}
	}
	return temp;
}


// compute the product of the PolySurf with another PolySurf and 
// represent the result as a new PolySurf 
template<class T>  
PolySurf<T> PolySurf<T>::Product(const PolySurf<T>& b) const
{
	return PolySurf<T>(ProductCPoints(b),ordu+b.GetOrdU()-1,ordv+b.GetOrdV()-1,b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV());
}
*/

// SUBDIVISION

// subdivide the BspSurf upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspSurf after the knot refinement
template<class T>
PolySurf<T> PolySurf<T>::Subdivide(double U1, double U2, double V1, double V2) const
{
	return SubdivideU(U1,U2).SubdivideV(V1,V2);
}


// CONVERSION

// convert to BezSurf form
template<class T>
BezSurf<T> PolySurf<T>::ConvertBezSurf() const
{
	PolySurf<T> p = Reparameterise1();
	return BezSurf<T>(p.ConvertBezSurfCPoints(),ordu,ordv,u1,u2,v1,v2);
}

// convert to BezSurf form
template<class T>
Matrix<T> PolySurf<T>::ConvertBezSurfCPoints() const
{
	Matrix<T> m = Math::mult2(Math::ComputePolyCurvMatrixInverse(ordu),*coeffs);
	
	return Math::mult1(m,Math::ComputePolyCurvMatrixInverseTranspose(ordv));
}

// convert to BspSurf form
template<class T>
BspSurf<T> PolySurf<T>::ConvertBspSurf() const
{
	return ConvertBezSurf().ConvertBspSurf();
}

// READ and WRITE

template <class T>
void PolySurf<T>::write(std::ostream& os) 
{
	os << "Poly Surface\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "limits in u and v are\n";
	os << GetLeftLimitU() << " " << GetRightLimitU();
	os << "\nlimits in v are\n";
	os << GetLeftLimitV() << " " << GetRightLimitV();
	os << "\ncoeffs are\n";
	os << *coeffs;
}

template <class T>
void PolySurf<T>::read(std::istream& is)
{
	int Ordu, Ordv;
	std::cout << "order of Poly Surface in u and v";
	is >> Ordu >> Ordv;
	double leftu, rightu, leftv, rightv;
	std::cout << "limits in u and v\n";
	is >> leftu >> rightu >> leftv >> rightv;
	Matrix<T> Coeffs(Ordu,Ordv);
	std::cout << "input control points";
	is >> Coeffs;
	*this = PolySurf<T>(Coeffs,Ordu,Ordv,leftu,rightu,leftv,rightv);
} 


template <class T>
void PolySurf<T>::writefile(std::ofstream& ofs) 
{
	ofs << "Poly Surface\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "limits in u and v are\n";
	ofs << GetLeftLimitU() << " " << GetRightLimitU();
	ofs << "\nlimits in v are\n";
	ofs << GetLeftLimitV() << " " << GetRightLimitV();
	ofs << "\ncoeffs are\n";
	ofs << *coeffs;
}

template <class T>
void PolySurf<T>::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv;
	ifs >> Ordu >> Ordv;
	double leftu, rightu, leftv, rightv;
	ifs >> leftu >> rightu >> leftv >> rightv;
	Matrix<T> Coeffs(Ordu,Ordv);
	ifs >> Coeffs;
	*this = PolySurf<T>(Coeffs,Ordu,Ordv,leftu,rightu,leftv,rightv);
} 

// PRIVATE FUNCTIONS

// elevate the degree of the PolySurf by level
template<class T>
PolySurf<T> PolySurf<T>::ElevateU(int level) const
{
	// return the new BspSurf
	if (level <= 0) return *this;
	else return PolySurf<T>(ElevateCPointsU(level), ordu+level, ordv, u1, u2, v1, v2);	
}

// elevate the degree of the PolySurf by level
template<class T>
PolySurf<T> PolySurf<T>::ElevateV(int level) const
{
	// return the new BspSurf
	if (level <= 0) return *this;
	else return PolySurf<T>(ElevateCPointsV(level), ordu, ordv+level, u1, u2, v1, v2);
}

// elevate the degree of the PolySurf by level
template<class T>
Matrix<T> PolySurf<T>::ElevateCPointsU(int level) const
{
	if (level <= 0) return *coeffs;
	Matrix<T> ncoeffs(ordu+level,ordv);

	// elevate the columns
	for (int j=0; j<ordv; j++) {
			Vector<T> temp = PolyCurv<T>((*coeffs).GetCol(j),ordu,u1,u2).ElevateCPoints(level);
			// extract control points
			for (int i=0; i<ordu+level; i++) ncoeffs[i][j]=temp[i];
	}
	// return the new BspSurf
	return ncoeffs;
}

// elevate the degree of the PolySurf by level
template<class T>
Matrix<T> PolySurf<T>::ElevateCPointsV(int level) const
{
	if (level <= 0) return *coeffs;
	Matrix<T> ncoeffs(ordu,ordv+level);

	// elevate the rows
	for (int i=0; i<ordu; i++) {
			Vector<T> temp = PolyCurv<T>((*coeffs).GetRow(i),ordv,v1,v2).ElevateCPoints(level);
			// extract control points
			for (int j=0; j<ordv+level; j++) ncoeffs[i][j]=temp[j];
	}
	// return the new BspSurf
	return ncoeffs;
}

// compute the derivative of the PolySurf of order deriv and
// represent the result as another PolySurf
template<class T>
PolySurf<T> PolySurf<T>::DeriveU(int level) const
{   
	if (level <= 0) return *this;
	// return the new BspSurf
	else return PolySurf<T>(DeriveCPointsU(level), ordu-level, ordv, u1, u2, v1, v2);
}

// compute the derivative of the PolySurf of order deriv and
// represent the result as another PolySurf
template<class T>
PolySurf<T> PolySurf<T>::DeriveV(int level) const
{
	if (level <= 0) return *this;
	// return the new BspSurf
	else return PolySurf<T>(DeriveCPointsV(level), ordu, ordv-level,u1,u2,v1,v2);
}


// compute the derivative of the PolySurf of order deriv and
// represent the result as another PolySurf
template<class T>
Matrix<T> PolySurf<T>::DeriveCPointsU(int level) const
{   
	if (level <= 0) return *coeffs;
	Matrix<T> ncoeffs(ordu-level,ordv);

	// derive the columns
	for (int j=0; j<ordv; j++) {
		Vector<T> temp = PolyCurv<T>((*coeffs).GetCol(j),ordu,u1,u2).DeriveCPoints(level);
			// extract control points
		for (int i=0; i<ordu-level; i++) ncoeffs[i][j]=temp[i];
	}
	// return the new BspSurf
	return ncoeffs;
}

// compute the derivative of the PolySurf of order deriv and
// represent the result as another PolySurf
template<class T>
Matrix<T> PolySurf<T>::DeriveCPointsV(int level) const
{
	if (level <= 0) return *coeffs;
	Matrix<T> ncoeffs(ordu,ordv-level);
	
	// derive the rows
	for (int i=0; i<ordu; i++) {
		Vector<T> temp = PolyCurv<T>((*coeffs).GetRow(i),ordv,v1,v2).DeriveCPoints(level);
			// extract control points
		for (int j=0; j<ordv; j++) ncoeffs[i][j]=temp[j];
	}
	// return the new BspSurf
	return ncoeffs;
}

// evaluate the derivative of the PolySurf of order deriv at a point
// x. Computes the derivative as a PolySurf and then evaluates this at x
template<class T>
T PolySurf<T>::DeriveU(int level, double u, double v) const
{
	return DeriveU(level)(u,v);
}

// evaluate the derivative of the PolySurf of order deriv at a point
// x. Computes the derivative as a PolySurf and then evaluates this at x
template<class T>
T PolySurf<T>::DeriveV(int level, double u, double v) const
{
	return DeriveV(level)(u,v);
}

// integrate the PolySurf between the limits x1 and x2. Computes
// the indefinite integral as a PolySurf and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T PolySurf<T>::IntegrateUV2(double U1, double U2, double V1, double V2) const
{
	Vector<T> v(ordv);

	for (int j=0; j<ordv; j++)
		v[j] = PolyCurv<T>(GetCol(j),ordu,u1,u2).Integrate(U1,U2);

	return PolyCurv<T>(v,ordv,v1,v2).Integrate(V1,V2);
}


// compute the indefinite integral of the PolySurf and represent
// it as a PolySurf of one higher degree
template<class T>
PolySurf<T> PolySurf<T>::IntegrateU() const
{
	return PolySurf<T>(IntegrateCPointsU(),ordu+1,ordv,u1,u2,v1,v2);
}	


// compute the indefinite integral of the PolySurf and represent
// it as a PolySurf of one higher degree
template<class T>
PolySurf<T> PolySurf<T>::IntegrateV() const
{
	return PolySurf<T>(IntegrateCPointsV(),ordu,ordv+1,u1,u2,v1,v2);
}	


// compute the indefinite integral of the PolySurf as a PolySurf
// and return just the control points of this Curvve as a CPoints
// object
template<class T>
Matrix<T> PolySurf<T>::IntegrateCPointsU() const
{
	Matrix<T> ncoeffs(ordu+1,ordv);

	for (int j=0; j<ordv; j++) {
		Vector<T> v = PolyCurv<T>((*coeffs).GetCol(j),ordu,u1,u2).IntegrateCPoints();		
		for (int i=0; i<ordu+1; i++) ncoeffs[i][j]=v[i];
	}
	return ncoeffs;
}	

// compute the indefinite integral of the PolySurf as a PolySurf
// and return just the control points of this Curvve as a CPoints
// object
template<class T>
Matrix<T> PolySurf<T>::IntegrateCPointsV() const
{
	Matrix<T> ncoeffs(ordu,ordv+1);
	Vector<T> v;
	// use create cpoints integrate to save time
	for (int i=0; i<ordu; i++) {
		Vector<T> v = PolyCurv<T>((*coeffs).GetRow(i),ordv,v1,v2).IntegrateCPoints();	
		for (int j=0; j<ordv+1; j++) ncoeffs[i][j]=v[j];
	}
	// return the matrix
	return ncoeffs;
}	

// compute the product of the PolySurf with another PolySurf and 
// represent the result as a new PolySurf 
template<class T>  
PolySurf<T> PolySurf<T>::Product2(const PolySurf<T>& b) const
{
	// Surfs should be compatible in terms of parameter limits
	Matrix<Vector<T> > mat1(ordu,ordu);	

	for (int i=0; i<ordu; i++) {
		temp1 = PolyCurv<T>(cpts.GetRow(i),ordv,v1,v2);
		for (int j=0; j<b.GetOrdU(); j++) {
			temp2 = PolyCurv<T>(b.GetCoeffs().GetRow(j),b.GetOrdV(),v1,v2);
			mat1[i][j] = temp1.ProductCPoints(temp2);	
		}
	}


	Matrix<Vector<T> > mat2(ordu,b.GetOrdU());
	Vector<double> v2(ordu,0.0);
	Vector<double> v3(b.GetOrdU(),0.0);

	for (int i=0; i<ordu; i++) {
		v2[i]=1.0;
		temp3 = PolyCurv<T>(v2, ordu, u1, u2);
		for (int j=0; j<b.GetOrdU(); j++) {
			v3[j] = 1.0;
			temp4 = PolyCurv<T>(v3, b.GetOrdU(),u1,u2);
			mat2[i][j] = temp3.ProductCPoints(temp4)
		}
	}

	int numcu = ordu+b.GetOrdU()-1;
	int numcv = ordv+b.GetOrdV()-1
	Matrix<T> ncoeffs(numcu,numcv);
	
	T sum=0.0;
	for (int i=0; i<numcu; i++)
		for (int j=0; j<numcv; j++) { 
			for (i1=0; i1<ordu; i1++) 
				for (j1=0; j1<b.GetOrdU(); j1++) sum = sum+(mat1[i1][j1])[j]*(mat2[i1][j1])[i];
			ncoeffs[i][j] = sum;
			sum = 0.0;
		}	
	return PolySurf<T>(ncpts,numcu,numcv,u1,u2,v1,v2);
}   


// subdivide the BspSurf upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspSurf after the knot refinement
template<class T>
PolySurf<T> PolySurf<T>::SubdivideU(double U1, double U2) const
{
	Matrix<T> ncoeffs(ordu,ordv);

	for (int j=0; j<ordv; j++) {
		Vector<T> temp = PolyCurv<T>((*coeffs).GetCol(j),ordu,u1,u2).SubdivideCPoints(U1,U2);
			// extract control points
		for (int i=0; i<ordu; i++) ncoeffs[i][j]=temp[i];
	}
	// return the new BspSurf
	return PolySurf<T>(ncoeffs, ordu, ordv, 0.0, 1.0, v1, v2);
}

// subdivide the BspSurf upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BspSurf after the knot refinement
template<class T>
PolySurf<T> PolySurf<T>::SubdivideV(double V1, double V2) const
{
	Matrix<T> ncoeffs(ordu,ordv);

	for (int i=0; i<ordu; i++) {
		Vector<T> temp = PolyCurv<T>((*coeffs).GetRow(i),ordv,v1,v2).SubdivideCPoints(V1,V2);
		// extract control points
		for (int j=0; j<ordv; j++) ncoeffs[i][j]=temp[j];
	}
	// return the new BspSurf
	return PolySurf<T>(ncoeffs, ordu, ordv, u1, u2, 0.0, 1.0);

}


#endif