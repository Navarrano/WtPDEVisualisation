#ifndef POLYVOL
#define POLYVOL

#include "BspVol.h"
#include "polysurf.h"

class PolyVolBasisFunc;

template<class T>
class BezVol;

template<class T>
class BspVol;

template<class T>
class PolyVol;

class PolyVolBasisFunc : public Vol<double> {

	// data
	int ordu;	 // order of basis function
	int ordv;	// order in v
	int ordw;
	int indexu;
	int indexv;
	int indexw;
	double leftlimitu;  // ranges in u
	double rightlimitu;
	double leftlimitv; // rnages in v
	double rightlimitv;
	double leftlimitw;
	double rightlimitw;
	Ptr<PolyVol<double> > b;	// PolyVol representation of basis function

	// private functions
	PolyVol<double> CreatePolyVol() const;	// creates the PolyVol
	virtual ObjectID Identity() const { return std::string("class PolyVolBasisFunc"); }
public:

	// constructors
	PolyVolBasisFunc();
	PolyVolBasisFunc(int Ordu, int Ordv, int Ordw, int Indexu, int Indexv, int Indexw, double LeftLimitU=0.0, double RightLimitU=1.0, 
	double LeftLimitV=0.0, double RightLimitV=1.0, double LeftLimitW=0.0, double RightLimitW=1.0);
	

	// evaluators
	double Eval(double u, double v, double w) const;
	virtual double operator()(double u, double v, double w) const;
	virtual double operator()(int, int, int, double u, double v, double w) const;
	// access functions
	PolyVol<double> GetPolyVol() const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	
	virtual double Derive(int, int , int, double, double, double) const;
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
	virtual double GetLeftLimitW() const;
	virtual double GetRightLimitW() const;
};


class PolyVolBasisFuncSet : public TextObject, public FTextObject {
	// data
	int ordu;	 // order of basis function
	int ordv;
	int ordw;
	double leftlimitu;
	double rightlimitu;
	double leftlimitv;
	double rightlimitv;
	double leftlimitw;
	double rightlimitw;

	Ptr<Matrix3D<PolyVolBasisFunc> > b;	// PolyVol representation of basis function
	virtual ObjectID Identity() const { return std::string("class PolyVolBasisFuncSet"); }
public:
	// constructors
	PolyVolBasisFuncSet();
	PolyVolBasisFuncSet(int Ordu, int Ordv, int Ordw, double LU, double RU, double LV, double RV, double LW, double RW);


	// evaluators
	Matrix3D<double> Eval(double u, double v, double w) const;
	Matrix3D<double> operator()(double u, double v, double w) const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

template<class T>
class PolyVol : public Vol<T> {
private:
	// data
	int ordu;  // order in u
	int ordv;	// order in v
	int ordw;
	Ptr<Matrix3D<T> > coeffs;  // coeffs
	double u1, u2;  // limits in u
	double v1, v2;  // limits in v
	double w1, w2;

	// private functions
	PolyVol<T> DeriveUV(int levu, int levv) const; 
	PolyVol<T> DeriveUW(int levu, int levw) const;
	PolyVol<T> DeriveVW(int levv, int levw) const;
	
	T DeriveU(int levu, double u, double v, double w) const;
	T DeriveV(int levv, double u, double v, double w) const;
	T DeriveW(int levw, double u, double v, double w) const;
	T DeriveUV(int levu, int levv, double u, double v, double w) const;
	T DeriveUW(int levu, int levw, double u, double v, double w) const;
	T DeriveVW(int levv, int levw, double u, double v, double w) const;
	Matrix3D<T> ElevateCPointsU(int levu) const;
	Matrix3D<T> ElevateCPointsV(int levu) const;
	Matrix3D<T> ElevateCPointsW(int levw) const;
	PolyVol<T> ElevateUV(int levu, int levv) const;
	PolyVol<T> ElevateUW(int levu, int levw) const;
	PolyVol<T> ElevateVW(int levv, int levw) const;
	PolyVol<T> IntegrateU() const;
	PolyVol<T> IntegrateV() const;
	PolyVol<T> IntegrateW() const;
	PolyVol<T> IntegrateUV() const;
	PolyVol<T> IntegrateUW() const;
	PolyVol<T> IntegrateVW() const;
	T IntegrateUVW2(double u1, double u2, double v1, double v2, double w1, double w2) const;
	PolyVol<T> SubdivideU(double u1, double u2) const;
	PolyVol<T> SubdivideV(double v1, double v2) const;
	PolyVol<T> SubdivideW(double w1, double w2) const;
	PolyVol<T> SubdivideUV(double u1, double u2, double v1, double v2) const;
	PolyVol<T> SubdivideUW(double u1, double u2, double w1, double w2) const;
	PolyVol<T> SubdivideVW(double v1, double v2, double w1, double w3) const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class PolyVol<double>");
		else {
			std::string s(typeid(T).name()), s1("class PolyVol<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	PolyVol();
	PolyVol(const Matrix3D<T>& Coeffs, int Ordu, int Ordv, int Ordw, 
		double U1=0.0, double U2=1.0, double V1=0.0, double V2=1.0, 
		double W1=0.0, double W2=1.0);
	PolyVol(int Ordu, int Ordv, int Ordw);

	// access functions
	int GetOrdU() const;
	int GetOrdV() const;
	int GetOrdW() const;
	Matrix3D<T> GetCoeffs() const;
	virtual double GetRightLimitU() const;
	virtual double GetRightLimitV() const;
	virtual double GetLeftLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitW() const;
	virtual double GetLeftLimitW() const;

	// reparameterise
	PolyVol<T> Reparameterise1() const;
	PolyVol<T> Reparameterise2(double U1, double U2, double V1, double V2, double W1, double W2) const;
	PolyVol<T> Reparameterise3(double U1, double U2, double V1, double V2, double W1, double W2) const;

	// evaluators
	virtual T operator()(double u, double v, double w) const;
	virtual T operator()(int, int, int, double u, double v, double w) const;
	T Eval(double u, double v, double w) const;

	// add and subtract
	PolyVol<T> Add(const PolyVol<T>& b) const;
	PolyVol<T> Subtract(const PolyVol<T>& b) const;

	// derivatives
	PolyVol<T> DeriveU(int levu) const; 
	PolyVol<T> DeriveV(int levv) const;
	PolyVol<T> DeriveW(int levw) const;
	PolyVol<T> Derive(int levu, int levv, int levw) const;
	virtual T Derive(int levu, int levv, int levw, double u, double v, double w) const;
	Matrix3D<T> DeriveCPointsU(int levu) const; 
	Matrix3D<T> DeriveCPointsV(int levv) const;
	Matrix3D<T> DeriveCPointsW(int levw) const;

	// degree elevation
	PolyVol<T> ElevateU(int levu) const;
	PolyVol<T> ElevateV(int levu) const;
	PolyVol<T> ElevateW(int levw) const;
	PolyVol<T> ElevateUVW(int levu, int levv, int levw) const;

	// integration
	PolyVol<T> IntegrateUVW() const;
	T Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<T> IntegrateCPoints() const;

	// product
	template<class T1>
	Matrix3D<T> ProductCPoints(const PolyVol<T1>& b) const
	{
		PolyVol<T> b1 = Reparameterise3(b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV(),b.GetLeftLimitW(),b.GetRightLimitW());
		int ordbu = b.GetOrdU();
		int ordbv = b.GetOrdV();
		int ordbw = b.GetOrdW();

		Matrix3D<T> temp(ordbu+ordu-1,ordbv+ordv-1,ordbw+ordw-1);
	
		double sum =0.0;

		for (int k=0; k<ordbu+ordu-1; k++) {
			for (int l=0; l<ordbv+ordv-1; l++) {
				for (int m=0; m<ordbw+ordw-1; m++) {
					for (int i=Math::max1(0,k-ordbu+1); i<=Math::min1(ordu-1,k); i++) {
						for (int j=Math::max1(0,l-ordbv+1); j<=Math::min1(ordv-1,l); j++) {
							for (int s=Math::max1(0,m-ordbw+1); s<=Math::min1(ordw-1,m); s++) {
								sum = sum + (Math::Combin(ordu-1,i)*Math::Combin(ordbu-1,k-i)*Math::Combin(ordv-1,j)*Math::Combin(ordbv-1,l-j)*Math::Combin(ordw-1,s)*Math::Combin(ordbw-1,m-s)*b1.GetCoeffs()[s][i][j]*b.GetCoeffs()[m-s][k-i][l-j])/(Math::Combin(ordu+ordbu-2,k)*Math::Combin(ordv+ordbv-2,l)*Math::Combin(ordw+ordbw-2,m));
							}
						}
					}
					temp[m][k][l]=sum;
					sum=0.0;	
				}
			}
		}
		return temp;
	}

	template<class T1>
	PolyVol<T> Product(const PolyVol<T1>& b) const
	{
		PolyVol<T> b1 = Reparameterise3(b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV(),b.GetLeftLimitW(),b.GetRightLimitW());
		int ordbu = b.GetOrdU();
		int ordbv = b.GetOrdV();
		int ordbw = b.GetOrdW();

		Matrix3D<T> temp(ordbu+ordu-1,ordbv+ordv-1,ordbw+ordw-1);
	
		double sum =0.0;

		for (int k=0; k<ordbu+ordu-1; k++) {
			for (int l=0; l<ordbv+ordv-1; l++) {
				for (int m=0; m<ordbw+ordw-1; m++) {
					for (int i=Math::max1(0,k-ordbu+1); i<=Math::min1(ordu-1,k); i++) {
						for (int j=Math::max1(0,l-ordbv+1); j<=Math::min1(ordv-1,l); j++) {
							for (int s=Math::max1(0,m-ordbw+1); s<=Math::min1(ordw-1,m); s++) {
								sum = sum + (Math::Combin(ordu-1,i)*Math::Combin(ordbu-1,k-i)*Math::Combin(ordv-1,j)*Math::Combin(ordbv-1,l-j)*Math::Combin(ordw-1,s)*Math::Combin(ordbw-1,m-s)*b1.GetCoeffs()[s][i][j]*b.GetCoeffs()[m-s][k-i][l-j])/(Math::Combin(ordu+ordbu-2,k)*Math::Combin(ordv+ordbv-2,l)*Math::Combin(ordw+ordbw-2,m));
							}
						}
					}
					temp[m][k][l]=sum;
					sum=0.0;	
				}
			}
		}
		// build the knots sets
		// build the knotset objects
	
		// return the new Surface
		return PolyVol<T>(temp, ordu+ordbu-1, ordv+ordbv-1, ordw+ordbw-1,b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV(),b.GetLeftLimitW(),b.GetRightLimitW());
	}

	PolyVol<T> Product2(const PolyVol<T>& b) const;

	// subdivision
	PolyVol<T> Subdivide(double u1, double u2, double v1, double v2, double w1, double w2) const;

	// conversion
	BezVol<T> ConvertBezVol() const;
	Matrix3D<T> ConvertBezVolCPoints() const;
	BspVol<T> ConvertBspVol() const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

// CONSTRUCTORS

// constructor builds a PolyVol from a Vector of control points, knots
// an order and number of control points
template<class T>
PolyVol<T>::PolyVol(const Matrix3D<T>& Coeffs, int Ordu, int Ordv, int Ordw, double U1, double U2, double V1, double V2, double W1, double W2) :
			ordu(Ordu), ordv(Ordv), ordw(Ordw), coeffs(new Matrix3D<T>(Coeffs)), u1(U1), u2(U2), v1(V1), v2(V2), w1(W1), w2(W2)
{
}

// default constructor
template<class T>
PolyVol<T>::PolyVol() : ordu(0), ordv(0), ordw(0), coeffs(), u1(0),u2(1),v1(0),v2(1),w1(0),w2(1) { }

// default constructor
template<class T>
PolyVol<T>::PolyVol(int Ordu, int Ordv, int Ordw) : ordu(Ordu), ordv(Ordv), ordw(Ordw), 
coeffs(new Matrix3D<T>(Ordu,Ordv,Ordw)), u1(0), u2(1), v1(0), v2(1), w1(0), w2(1) { }


// ACCESS FUNCTIONS

// get the order of the PolyVol
template<class T>
inline int PolyVol<T>::GetOrdU() const { return ordu; }

// get the order of the PolyVol
template<class T>
inline int PolyVol<T>::GetOrdV() const { return ordv; }

// get the order of the PolyVol
template<class T>
inline int PolyVol<T>::GetOrdW() const { return ordw; }
// get the control points

template<class T>
inline Matrix3D<T> PolyVol<T>::GetCoeffs() const { return *coeffs; }

// get the knot vector
template<class T>
inline double PolyVol<T>::GetLeftLimitU() const { return u1; }

// get the knot vector
template<class T>
inline double PolyVol<T>::GetRightLimitU() const { return u2; }

// get the knot vector
template<class T>
inline double PolyVol<T>::GetLeftLimitV() const { return v1; }

// get the knot vector
template<class T>
inline double PolyVol<T>::GetRightLimitV() const { return v2; }


// get the knot vector
template<class T>
inline double PolyVol<T>::GetLeftLimitW() const { return w1; }

// get the knot vector
template<class T>
inline double PolyVol<T>::GetRightLimitW() const { return w2; }



// REPARAMETERISE

// reparameterise from (t1,t2) to (0,1)
// formula in Sherar/Goult book
template<class T>
PolyVol<T> PolyVol<T>::Reparameterise1() const
{
	Matrix3D<T> ncoeffs(ordu,ordv,ordw);

	for (int k=0; k<ordw; k++) {
		PolySurf<T> b = PolySurf<T>((*coeffs).GetUV(k),ordu,ordv,u1,u2,v1,v2).Reparameterise1();
		for (int i=0; i<ordu; i++)
			for (int j=0; j<ordv; j++) ncoeffs[k][i][j] = b.GetCoeffs()[i][j];
	}

	Matrix3D<T> ncoeffs1(ordu,ordv,ordw);
	for (int i=0; i<ordu; i++)
		for (int j=0; j<ordv; j++) {
			PolyCurv<T> b = PolyCurv<T>(ncoeffs.GetW(i,j),ordw,w1,w2).Reparameterise1();
			for (int k=0; k<ordw; k++)  
				ncoeffs1[k][i][j] = b.GetCoeffs()[k];
		}		
    return PolyVol<T>(ncoeffs1,ordu,ordv,ordw,0.0,1.0,0.0,1.0,0.0,1.0);
}
	

// reparameterise from (0,1) to (T1, T2)
// formula in Sherar/Goult book
template<class T>
PolyVol<T> PolyVol<T>::Reparameterise2(double U1, double U2, double V1, double V2, double W1, double W2) const
{
	Matrix3D<T> ncoeffs(ordu,ordv,ordw);

	for (int k=0; k<ordw; k++) {
		PolySurf<T> b = PolySurf<T>((*coeffs).GetUV(k),ordu,ordv,u1,u2,v1,v2).Reparameterise2(U1,U2,V1,V2);
		for (int i=0; i<ordu; i++)
			for (int j=0; j<ordv; j++) ncoeffs[k][i][j] = b.GetCoeffs()[i][j];
	}

	Matrix3D<T> ncoeffs1(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++)
		for (int j=0; j<ordv; j++) {
			PolyCurv<T> b = PolyCurv<T>(ncoeffs.GetW(i,j),ordw,w1,w2).Reparameterise2(W1,W2);
			for (int k=0; k<ordw; k++)  
				ncoeffs1[k][i][j] = b.GetCoeffs()[k];
		}		
    return PolyVol<T>(ncoeffs1,ordu,ordv,ordw,U1,U2,V1,V2,W1,W2);
}

// reparameterise from (t1,t2) to (T1, T2)
template<class T>
PolyVol<T> PolyVol<T>::Reparameterise3(double U1, double U2, double V1, double V2, double W1, double W2) const
{
	// go from t1,t2 to 0,1 then 0,1 to T1, T2
	return Reparameterise1().Reparameterise2(U1,U2,V1,V2,W1,W2);
}


// ADD and SUBTRACT

// subtract two PolyVol's together
// assuming same range in u and v different orders
template<class T>
PolyVol<T> PolyVol<T>::Subtract(const PolyVol<T>& b) const
{
	PolyVol<T> b1 = Reparameterise3(b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV(),b.GetLeftLimitW(),b.GetRightLimitW());
	PolyVol<T> c(b1), d(b);

	if (b.GetOrdU() > ordu) c = c.ElevateU(b.GetOrdU()-ordu);
	else if (ordu > b.GetOrdU()) d = d.ElevateU(ordu-b.GetOrdU());
	if (b.GetOrdV() > ordv) c = c.ElevateV(b.GetOrdV()-ordv);
	else if (ordv > b.GetOrdV()) d = d.ElevateV(ordv-b.GetOrdV());
	if (b.GetOrdW() > ordw) c = c.ElevateW(b.GetOrdW()-ordw);
	else if (ordw > b.GetOrdW()) d = d.ElevateW(ordw-b.GetOrdW());

	Matrix3D<T> temp(d.GetOrdU(),d.GetOrdV(),d.GetOrdW());
	for (int i=0; i<d.GetOrdU(); i++)
		for (int j=0; j<d.GetOrdV(); j++)
			for (int k=0; k<d.GetOrdW(); k++)
				temp[k][i][j]=c.GetCoeffs()[k][i][j]-d.GetCoeffs()[k][i][j];

	return PolyVol<T>(temp,d.GetOrdU(),d.GetOrdV(),d.GetOrdW(),u1,u2,v1,v2,w1,w2);
}

// add two PolyVol's together
// assuming same range in u and v different orders
template<class T>
PolyVol<T> PolyVol<T>::Add(const PolyVol<T>& b) const
{
	PolyVol<T> b1 = Reparameterise3(b.GetLeftLimitU(),b.GetRightLimitU(),b.GetLeftLimitV(),b.GetRightLimitV(),b.GetLeftLimitW(),b.GetRightLimitW());
	PolyVol<T> c(b1), d(b);

	if (b.GetOrdU() > ordu) c = c.ElevateU(b.GetOrdU()-ordu);
	else if (ordu > b.GetOrdU()) d = d.ElevateU(ordu-b.GetOrdU());
	if (b.GetOrdV() > ordv) c = c.ElevateV(b.GetOrdV()-ordv);
	else if (ordv > b.GetOrdV()) d = d.ElevateV(ordv-b.GetOrdV());
	if (b.GetOrdW() > ordw) c = c.ElevateW(b.GetOrdW()-ordw);
	else if (ordw > b.GetOrdW()) d = d.ElevateW(ordw-b.GetOrdW());

	Matrix3D<T> temp(d.GetOrdU(),d.GetOrdV(),d.GetOrdW());
	for (int i=0; i<d.GetOrdU(); i++)
		for (int j=0; j<d.GetOrdV(); j++)
			for (int k=0; k<d.GetOrdW(); k++)
				temp[k][i][j]=c.GetCoeffs()[k][i][j]+d.GetCoeffs()[k][i][j];

	return PolyVol<T>(temp,d.GetOrdU(),d.GetOrdV(),d.GetOrdW(),u1,u2,v1,v2,w1,w2);			
}


// EVALUATORS

// evaluate the PolyVol at the point x using de Boor algorithm
template<class T>
T PolyVol<T>::operator()(double u, double v, double w) const
{
  	Vector<T> temp(ordw);

	for (int k=0; k<ordw; k++) temp[k]= PolySurf<T>((*coeffs).GetUV(k),ordu,ordv,u1,u2,v1,v2)(u,v);

    return PolyCurv<T>(temp,ordw,w1,w2)(w);
}


// evaluate the PolyVol using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T PolyVol<T>::Eval(double u, double v, double w) const
{
  	Vector<T> temp(ordw);

	for (int k=0; k<ordw; k++) temp[k]= PolySurf<T>((*coeffs).GetUV(k),ordu,ordv,u1,u2,v1,v2)(u,v);

    return PolyCurv<T>(temp,ordw,w1,w2)(w);	
}

// DEGREE ELEVATION


// elevate the degree of the PolyVol by level
template<class T>
PolyVol<T> PolyVol<T>::ElevateUVW(int levu, int levv, int levw) const
{
	return (ElevateU(levu).ElevateV(levv)).ElevateW(levw);
}


// DERIVATIVES

// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
PolyVol<T> PolyVol<T>::Derive(int levu, int levv, int levw) const
{
	return DeriveU(levu).DeriveV(levv).DeriveW(levw);
}

template<class T>
T PolyVol<T>::operator() (int valu, int valv, int valw, double u, double v, double w) const
{
	return Derive(valu,valv,valw,u,v,w);
}



// evaluate the derivative of the PolyVol of order deriv at a point
// x. Computes the derivative as a PolyVol and then evaluates this at x
template<class T>
T PolyVol<T>::Derive(int levu, int levv, int levw, double u, double v, double w) const
{
	return (DeriveU(levu).DeriveV(levv)).DeriveW(levw)(u,v,w);
}

// INTEGRATION

// integrate the PolyVol between the limits x1 and x2. Computes
// the indefinite integral as a PolyVol and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T PolyVol<T>::Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	PolyVol<T> temp = IntegrateUVW();
	
	return temp(u2,v2,w2)-temp(u2,v2,w1)-temp(u2,v1,w2)-temp(u1,v2,w2)+
		temp(u2,v1,w1)+temp(u1,v2,w1)+temp(u1,v1,w2)-temp(u1,v1,w1);
}




// SUBDIVISION


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
PolyVol<T> PolyVol<T>::Subdivide(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	return (SubdivideU(u1,u2).SubdivideV(v1,v2)).SubdivideW(w1,w2);
}


// CONVERSION

template<class T>
BezVol<T> PolyVol<T>::ConvertBezVol() const
{
	PolyVol<T> p = Reparameterise1();
	return BezVol<T>(p.ConvertBezVolCPoints(),ordu,ordv,ordw,u1,u2,v1,v2,w1,w2);
}

template<class T>
Matrix3D<T> PolyVol<T>::ConvertBezVolCPoints() const
{
	Matrix3D<T> m = Math::MM3DI_1(Math::ComputePolyCurvMatrixInverseTranspose(ordu), *coeffs);

	m = Math::MM3DJ_1(Math::ComputePolyCurvMatrixInverseTranspose(ordv), m);

	return Math::MM3DK_1(Math::ComputePolyCurvMatrixInverseTranspose(ordw), m);
}


template<class T>
BspVol<T> PolyVol<T>::ConvertBspVol() const
{
	return (ConvertBezVol().ConvertBspVol());
}

// READ and WRITE

template <class T>
void PolyVol<T>::write(std::ostream& os) 
{
	os << "Poly Volume\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "order in w is " << ordw << "\n";
	os << "limits in u are\n";
	os << GetLeftLimitU() << " " << GetRightLimitU();
	os << "\nlimits in v are \n";
	os << GetLeftLimitV() << " " << GetRightLimitV();
	os << "\nlimits in w are \n";
	os << GetLeftLimitW() << " " << GetRightLimitW();
	os << "\ncoeffs are\n";
	os << *coeffs;
}

template <class T>
void PolyVol<T>::read(std::istream& is)
{
	int Ordu, Ordv, Ordw;
	std::cout << "order of Poly Volume in u and v and w";
	is >> Ordu >> Ordv >> Ordw;
	Matrix3D<T> Coeffs(Ordu,Ordv,Ordw);
	double leftu, leftv, rightu, rightv, leftw, rightw;
	std::cout << "input limits in u\n";
	is >> leftu >> rightu;
	std::cout << "input limits in v\n";
	is >> leftv >> rightv;
	std::cout << "input limits in w\n";
	is >> leftw >> rightw;
	std::cout << "\ninput coeffs\n";
	is >> Coeffs;
	*this = PolyVol<T>(Coeffs,Ordu,Ordv,Ordw,leftu,rightu,leftv,rightv,leftw,rightw);
} 


template <class T>
void PolyVol<T>::writefile(std::ofstream& ofs) 
{
	ofs << "Poly Volume\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "order in w is " << ordw << "\n";
	ofs << "limits in u are\n";
	ofs << GetLeftLimitU() << " " << GetRightLimitU();
	ofs << "\nlimits in v are \n";
	ofs << GetLeftLimitV() << " " << GetRightLimitV();
	ofs << "\nlimits in w are \n";
	ofs << GetLeftLimitW() << " " << GetRightLimitW();
	ofs << "\ncoeffs are\n";
	ofs << *coeffs;
}


template <class T>
void PolyVol<T>::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw;
	ifs >> Ordu >> Ordv >> Ordw;
	Matrix3D<T> Coeffs(Ordu,Ordv,Ordw);
	double leftu, leftv, rightu, rightv, leftw, rightw;
	ifs >> leftu >> rightu;
	ifs >> leftv >> rightv;
	ifs >> leftw >> rightw;
	ifs >> Coeffs;
	*this = PolyVol<T>(Coeffs,Ordu,Ordv,Ordw,leftu,rightu,leftv,rightv,leftw,rightw);
} 

// PRIVATE FUNCTIONS

// elevate the degree of the PolyVol by level
template<class T>
PolyVol<T> PolyVol<T>::ElevateU(int level) const
{
	if (level <= 0) return *this;
	Matrix3D<T> ncpts(ordu+level,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = PolyCurv<T>((*coeffs).GetU(j,k),ordu,u1,u2).ElevateCPoints(level);
			// extract control points
			for (int i=0; i<ordu+level; i++) ncpts[k][i][j]=temp[i];
		}

	// return the new BVol
	return PolyVol<T>(ncpts, ordu+level, ordv, ordw, u1,u2,v1,v2,w1,w2);		
}


// elevate the degree of the PolyVol by level
template<class T>
Matrix3D<T> PolyVol<T>::ElevateCPointsU(int level) const
{
	if (level <= 0) return *coeffs;
	Matrix3D<T> ncpts(ordu+level,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = PolyCurv<T>((*coeffs).GetU(j,k),ordu,u1,u2).ElevateCPoints(level);
			// extract control points
			for (int i=0; i<ordu+level; i++) ncpts[k][i][j]=temp[i];
		}

	return ncpts;
}



// elevate the degree of the PolyVol by level
template<class T>
PolyVol<T> PolyVol<T>::ElevateV(int level) const
{
	if (level <= 0) return *this;
	Matrix3D<T> ncpts(ordu,ordv+level,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = PolyCurv<T>((*coeffs).GetV(i,k),ordv,v1,v2).ElevateCPoints(level);
			// extract control points
			for (int j=0; j<ordv+level; j++) ncpts[k][i][j]=temp[j];
	}
	// return the new BVol
	return PolyVol<T>(ncpts,ordu, ordv+level, ordw, u1,u2,v1,v2,w1,w2);
}


// elevate the degree of the PolyVol by level
template<class T>
Matrix3D<T> PolyVol<T>::ElevateCPointsV(int level) const
{
	if (level <= 0) return *coeffs;
	Matrix3D<T> ncpts(ordu,ordv+level,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = PolyCurv<T>((*coeffs).GetV(i,k),ordv,v1,v2).ElevateCPoints(level);
			// extract control points
			for (int j=0; j<ordv+level; j++) ncpts[k][i][j]=temp[j];
	}
	return ncpts;
}


// elevate the degree of the PolyVol by level
template<class T>
PolyVol<T> PolyVol<T>::ElevateW(int level) const
{
	if (level <= 0) return *this;
	Matrix3D<T> ncpts(ordu,ordv,ordw+level);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> temp = PolyCurv<T>((*coeffs).GetW(i,j),ordw,w1,w2).ElevateCPoints(level);
			// extract control points
			for (int k=0; k<ordw+level; k++) ncpts[k][i][j]=temp[k];
	}
	// return the new BVol
	return PolyVol<T>(ncpts, ordu, ordv, ordw+level,u1,u2,v1,v2,w1,w2);
}


// elevate the degree of the PolyVol by level
template<class T>
Matrix3D<T> PolyVol<T>::ElevateCPointsW(int level) const
{
	if (level <= 0) return *coeffs;
	Matrix3D<T> ncpts(ordu,ordv,ordw+level);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> temp = PolyCurv<T>((*coeffs).GetW(i,j),ordw,w1,w2).ElevateCPoints(level);
			// extract control points
			for (int k=0; k<ordw+level; k++) ncpts[k][i][j]=temp[k];
	}
	return ncpts;
}


// elevate the degree of the PolyVol by level
template<class T>
PolyVol<T> PolyVol<T>::ElevateUV(int levu, int levv) const
{
	return ElevateU(levu).ElevateV(levv);
}

// elevate the degree of the PolyVol by level
template<class T>
PolyVol<T> PolyVol<T>::ElevateUW(int levu, int levw) const
{
	return ElevateU(levu).ElevateW(levw);
}

// elevate the degree of the PolyVol by level
template<class T>
PolyVol<T> PolyVol<T>::ElevateVW(int levv, int levw) const
{
	return ElevateV(levv).ElevateW(levw);
}

// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
PolyVol<T> PolyVol<T>::DeriveU(int level) const
{   
	if (level <= 0) return *this;
	Matrix3D<T> ncpts(ordu-level,ordv,ordw);

	// only need to compute knot vector once, better to use version that
	// returns control points only
	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetU(j,k),ordu,u1,u2).DeriveCPoints(level);
			// extract control points
			for (int i=0; i<ordu-level; i++) ncpts[k][i][j]=v[i];
	}
	// return the new BVol
	return PolyVol<T>(ncpts,ordu-level, ordv, ordw,u1,u2,v1,v2,w1,w2);
}

// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
Matrix3D<T> PolyVol<T>::DeriveCPointsU(int level) const
{   
	if (level <= 0) return *coeffs;
	Matrix3D<T> ncpts(ordu-level,ordv,ordw);

	// only need to compute knot vector once, better to use version that
	// returns control points only
	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetU(j,k),ordu,u1,u2).DeriveCPoints(level);
			// extract control points
			for (int i=0; i<ordu-level; i++) ncpts[k][i][j]=v[i];
	}
	// return the new BVol
	return ncpts;
}


// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
PolyVol<T> PolyVol<T>::DeriveV(int level) const
{
	if (level <= 0) return *this;
	Matrix3D<T> ncpts(ordu,ordv-level,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetV(i,k),ordv,v1,v2).DeriveCPoints(level);
			// extract control points
			for (int j=0; j<ordv-level; j++) ncpts[k][i][j]=v[j];
	}
	// return the new BVol
	return PolyVol<T>(ncpts, ordu, ordv-level, ordw, u1,u2,v1,v2,w1,w2);
}


// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
Matrix3D<T> PolyVol<T>::DeriveCPointsV(int level) const
{
	if (level <= 0) return *coeffs;
	Matrix3D<T> ncpts(ordu,ordv-level,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetV(i,k),ordv,v1,v2).DeriveCPoints(level);
			// extract control points
			for (int j=0; j<ordv-level; j++) ncpts[k][i][j]=v[j];
	}
	// return the new BVol
	return ncpts;
}



// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
PolyVol<T> PolyVol<T>::DeriveW(int level) const
{
	if (level <= 0) return *this;
	Matrix3D<T> ncpts(ordu,ordv,ordw-level);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetW(i,j),ordw,w1,w2).DeriveCPoints(level);
			// extract control points
			for (int k=0; k<ordw-level; k++) ncpts[k][i][j]=v[k];
	}
	// return the new BVol
	return PolyVol<T>(ncpts, ordu, ordv, ordw-level, u1,u2,v1,v2,w1,w2);
}


// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
Matrix3D<T> PolyVol<T>::DeriveCPointsW(int level) const
{
	if (level <= 0) return *coeffs;
	Matrix3D<T> ncpts(ordu,ordv,ordw-level);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetW(i,j),ordw,w1,w2).DeriveCPoints(level);
			// extract control points
			for (int k=0; k<ordw-level; k++) ncpts[k][i][j]=v[k];
	}
	// return the new BVol
	return ncpts;
}


// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
PolyVol<T> PolyVol<T>::DeriveUV(int levu, int levv) const
{
	return DeriveU(levu).DeriveV(levv);
}

// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
PolyVol<T> PolyVol<T>::DeriveUW(int levu, int levw) const
{
	return DeriveU(levu).DeriveW(levw);
}


// compute the derivative of the PolyVol of order deriv and
// represent the result as another PolyVol
template<class T>
PolyVol<T> PolyVol<T>::DeriveVW(int levv, int levw) const
{
	return DeriveV(levv).DeriveW(levw);
}




// evaluate the derivative of the PolyVol of order deriv at a point
// x. Computes the derivative as a PolyVol and then evaluates this at x
template<class T>
T PolyVol<T>::DeriveU(int level, double u, double v, double w) const
{
	return DeriveU(level)(u,v,w);
}


// evaluate the derivative of the PolyVol of order deriv at a point
// x. Computes the derivative as a PolyVol and then evaluates this at x
template<class T>
T PolyVol<T>::DeriveV(int level, double u, double v, double w) const
{
	return DeriveV(level)(u,v,w);
}


// evaluate the derivative of the PolyVol of order deriv at a point
// x. Computes the derivative as a PolyVol and then evaluates this at x
template<class T>
T PolyVol<T>::DeriveW(int level, double u, double v, double w) const
{
	return DeriveW(level)(u,v,w);
}


// evaluate the derivative of the PolyVol of order deriv at a point
// x. Computes the derivative as a PolyVol and then evaluates this at x
template<class T>
T PolyVol<T>::DeriveUV(int levu, int levv, double u, double v, double w) const
{
	return (DeriveU(levu).DeriveV(levv))(u,v,w);
}


// evaluate the derivative of the PolyVol of order deriv at a point
// x. Computes the derivative as a PolyVol and then evaluates this at x
template<class T>
T PolyVol<T>::DeriveUW(int levu, int levw, double u, double v, double w) const
{
	return (DeriveU(levu).DeriveW(levw))(u,v,w);
}


// evaluate the derivative of the PolyVol of order deriv at a point
// x. Computes the derivative as a PolyVol and then evaluates this at x
template<class T>
T PolyVol<T>::DeriveVW(int levv, int levw, double u, double v, double w) const
{
	return (DeriveV(levv).DeriveW(levw))(u,v,w);
}
// integrate the PolyVol between the limits x1 and x2. Computes
// the indefinite integral as a PolyVol and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T PolyVol<T>::IntegrateUVW2(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Vector<T> w(ordw);

	// should be integrate CPoints
	for (int k=0; k<ordw; k++)
		w[k] = PolySurf<T>((*coeffs).GetUV(k),ordu,ordv,u1,u2,v1,v2).Integrate(u1,u2,v1,v2);
	return PolyCurv<T>(w,ordw,w1,w2).Integrate(w1,w2);
}


// compute the indefinite integral of the PolyVol and represent
// it as a PolyVol of one higher degree
template<class T>
PolyVol<T> PolyVol<T>::IntegrateU() const
{
	Matrix3D<T> ncpts(ordu+1,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetU(j,k),ordu,u1,u2).IntegrateCPoints();
			// extract control points
			for (int i=0; i<ordu+1; i++) ncpts[k][i][j]=v[i];
	}	

	return PolyVol<T>(ncpts,ordu+1,ordv,ordw,u1,u2,v1,v2,w1,w2);
}	


// compute the indefinite integral of the PolyVol and represent
// it as a PolyVol of one higher degree
template<class T>
PolyVol<T> PolyVol<T>::IntegrateV() const
{
	Matrix3D<T> ncpts(ordu,ordv+1,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetV(i,k),ordv,v1,v2).IntegrateCPoints();
			// extract control points
			for (int j=0; j<ordv+1; j++) ncpts[k][i][j]=v[j];
	}	

	return PolyVol<T>(ncpts,ordu,ordv+1,ordw,u1,u2,v1,v2,w1,w2);	
}	

// compute the indefinite integral of the PolyVol and represent
// it as a PolyVol of one higher degree
template<class T>
PolyVol<T> PolyVol<T>::IntegrateW() const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw+1);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetW(i,j),ordw,w1,w2).IntegrateCPoints();
			// extract control points
			for (int k=0; k<ordw+1; k++) ncpts[k][i][j]=v[k];
	}	
	
	return PolyVol<T>(ncpts,ordu,ordv,ordw+1,u1,u2,v1,v2,w1,w2);	
}	


// compute the indefinite integral of the PolyVol and represent
// it as a PolyVol of one higher degree
template<class T>
PolyVol<T> PolyVol<T>::IntegrateUV() const
{
	return IntegrateU().IntegrateV();
}	

// compute the indefinite integral of the PolyVol and represent
// it as a PolyVol of one higher degree
template<class T>
PolyVol<T> PolyVol<T>::IntegrateUW() const
{
	return IntegrateU().IntegrateW();
}	


// compute the indefinite integral of the PolyVol and represent
// it as a PolyVol of one higher degree
template<class T>
PolyVol<T> PolyVol<T>::IntegrateVW() const
{
	return IntegrateV().IntegrateW();
}	


// compute the product of the PolyVol with another PolyVol and 
// represent the result as a new PolyVol 
template<class T>  
PolyVol<T> PolyVol<T>::Product2(const PolyVol<T>& b) const
{
	return (ConvertBezVol().Product(b.ConvertBezVol())).ConvertPolyVol();
}   

// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
PolyVol<T> PolyVol<T>::SubdivideU(double u1, double u2) const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetU(j,k),ordu,u1,u2).SubdivideCPoints(u1,u2);
			// extract control points
		for (int i=0; i<ordu; i++) ncpts[k][i][j]=v[i];
	}
	// return the new BVolace
	return PolyVol<T>(ncpts, ordu, ordv, ordw,u1,u2,v1,v2,w1,w2);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVolace after the knot refinement
template<class T>
PolyVol<T> PolyVol<T>::SubdivideV(double v1, double v2) const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetV(i,k),ordv,v1,v2).SubdivideCPoints(v1,v2);
			// extract control points
		for (int j=0; j<ordv; j++) ncpts[k][i][j]=v[j];
	}
	// return the new BVol
	return PolyVol<T>(ncpts, ordu, ordv, ordw,u1,u2,v1,v2,w1,w2);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVolace after the knot refinement
template<class T>
PolyVol<T> PolyVol<T>::SubdivideW(double w1, double w2) const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> v = PolyCurv<T>((*coeffs).GetW(i,j),ordw,w1,w2).SubdivideCPoints(w1,w2);
			// extract control points
		for (int k=0; k<ordw; k++) ncpts[k][i][j]=v[k];
	}
	// return the new BVol
	return PolyVol<T>(ncpts, ordu, ordv, ordw,u1,u2,v1,v2,w1,w2);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
PolyVol<T> PolyVol<T>::SubdivideUV(double u1, double u2, double v1, double v2) const
{
	return SubdivideU(u1,u2).SubdivideV(v1,v2);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
PolyVol<T> PolyVol<T>::SubdivideUW(double u1, double u2, double w1, double w2) const
{
	return SubdivideU(u1,u2).SubdivideW(w1,w2);
}

// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
PolyVol<T> PolyVol<T>::SubdivideVW(double v1, double v2, double w1, double w2) const
{
	return SubdivideV(v1,v2).SubdivideW(w1,w2);
}


#endif

