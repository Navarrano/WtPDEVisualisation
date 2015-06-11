
#ifndef BEZSURF
#define BEZSURF


#include "Matrix.h"
#include "BezCurv.h"


class BezSurfBasisFunc : public Surf<double> {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<BezSurf<double> > b;	// BezSurf representation of basis function
	BezSurf<double> CreateBezSurf() const;	// creates the BezSurf
	virtual ObjectID Identity() const { return std::string("class BezSurfBasisFunc"); }
	int ComputeDimU() const;
	int ComputeDimV() const;
public:

	// constructors
	BezSurfBasisFunc();
	BezSurfBasisFunc(int Ordu, int Ordv, const Vector<double>& Ktsu, const Vector<double>& Ktsv);
	
	// evaluators
	double Eval(double u, double v) const;	// evaluate basis function
	virtual double operator()(double u, double v) const;
	virtual double operator() (int, int, double, double) const;
	virtual double Derive(int, int, double, double) const;

	// access functions
	BezSurf<double> GetBezSurf() const;		// get BezSurf representation
	int GetOrdU() const;
	int GetOrdV() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;


	// read and write
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;

	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};


class BezSurfBasisFuncSet : public TextObject, public FTextObject {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v

	// private functions
	Ptr<Matrix<BezSurfBasisFunc> > b;
	virtual ObjectID Identity() const { return std::string("class BezSurfBasisFuncSet"); }
public:

	// constructors
	BezSurfBasisFuncSet();
	BezSurfBasisFuncSet(int Ordu, int Ordv, const Vector<double>& Ktsu, const Vector<double>& Ktsv);
	
	// get functions
	int GetOrdU() const;
	int GetOrdV() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	BezSurfBasisFunc GetBezSurfBasisFunc(int i, int j) const;

	// evaluators
	Matrix<double> Eval(double u, double v) const;	// evaluate basis function
	Matrix<double> operator()(double u, double v) const;
	Matrix<double> operator() (int, int, double, double) const;
	Matrix<double> Derive(int, int, double, double) const;

	
	// read and write
	virtual void read(std::istream& is);
	virtual void readfile(std::ifstream& ifs);
	virtual void write(std::ostream& os);
	virtual void writefile(std::ofstream& ofs);
};



template<class T>
class BezSurf : public Surf<T> {
private:
	// data
	int ordu;
	int ordv;
	Ptr<Matrix<T> > cpts;
	Ptr<Vector<double> > ktsu;
	Ptr<Vector<double> > ktsv;
	Ptr<KnotSet> ksetu;
	Ptr<KnotSet> ksetv;
	double leftLimitU;
	double rightLimitU;
	double leftLimitV;
	double rightLimitV;

	// private functions
	// derivatives

	// degree elevation

	// integration
	Matrix<T> IntegrateCPointsUV() const;
	Matrix<T> IntegrateCPointsU() const;
	Matrix<T> IntegrateCPointsV() const;

	// subdivision
	BspSurf<T> SubdivideU(int level) const;
	BspSurf<T> SubdivideV(int level) const;
	Matrix<T> SubdivideCPointsU(int level) const;
	Matrix<T> SubdivideCPointsV(int level) const;
	BezSurf<T> SubdivideU(double u1, double u2) const;
	BezSurf<T> SubdivideV(double v1, double v2) const;
	Matrix<T> SubdivideCPointsU(double u1, double u2) const;
	Matrix<T> SubdivideCPointsV(double v1, double v2) const;

	// insert knots u
	Matrix<T> InsertKnotCPointsU(double x) const;
	BspSurf<T> InsertKnotU(double x) const;
	Matrix<T> InsertKnotCPointsU(double x,int level) const;
	BspSurf<T> InsertKnotU(double x,int level) const;
	Matrix<T> InsertKnotCPointsU(const Vector<double>& Kts,int n) const;
	BspSurf<T> InsertKnotU(const Vector<double>& Kts, int n) const;
	
	// insert knots v
	BspSurf<T> InsertKnotV(double x) const;
	Matrix<T> InsertKnotCPointsV(double x,int level) const;
	BspSurf<T> InsertKnotV(double x,int level) const;
	Matrix<T> InsertKnotCPointsV(const Vector<double>& Kts,int n) const;
	BspSurf<T> InsertKnotV(const Vector<double>& Kts, int n) const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class BezSurf<double>");
		else {
			std::string s(typeid(T).name()), s1("class BezSurf<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	BezSurf();
	BezSurf(const Matrix<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv);
	BezSurf(const Matrix<T>& Cpts, int Ordu, int Ordv);
	BezSurf(const Matrix<T>& Cpts, int Ordu, int Ordv, double Lu, double Ru, double Lv, double Rv);
	BezSurf(int Ordu, int Ordv);
	
	// access functions
	int GetOrdU() const;
	int GetNumU() const;
	int GetOrdV() const;
	int GetNumV() const;
	Matrix<T> GetCPoints() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	KnotSet GetKnotSetU() const;
	KnotSet GetKnotSetV() const;
	virtual double GetLeftLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitU() const;
	virtual double GetRightLimitV() const;


	// evaluation functions
	virtual T operator()(double u, double v) const;
	T Eval(double u, double v) const;
	
	// addition and subtraction
	BezSurf<T> Add(const BezSurf<T>& b) const;
	BezSurf<T> Subtract(const BezSurf<T>& b) const;

	// derivatives
	BezSurf<T> Derive(int levu, int levv) const; 
	Matrix<T> DeriveCPoints(int levu, int levv) const;
	virtual T Derive(int levu, int levv, double u, double v) const;
	BezSurf<T> DeriveU(int levu) const; 
	BezSurf<T> DeriveV(int levv) const; 
	Matrix<T> DeriveCPointsU(int levu) const;
	Matrix<T> DeriveCPointsV(int levv) const;
	T DeriveU(int levu, double u, double v) const;
	T DeriveV(int levv, double u, double v) const;
	virtual T operator() (int, int, double, double) const;
	
	// degree elevation
	BezSurf<T> ElevateUV(int levu, int levv) const;
	Matrix<T> ElevateCPointsUV(int levu, int levv) const;
	BezSurf<T> ElevateU(int levu) const;
	BezSurf<T> ElevateV(int levu) const;
	Matrix<T> ElevateCPointsU(int levu) const;
	Matrix<T> ElevateCPointsV(int levu) const;
	
	// integration
	BezSurf<T> IntegrateU() const;
	BezSurf<T> IntegrateV() const;
	BezSurf<T> IntegrateUV() const;
	T IntegrateUV1(double u1, double u2, double v1, double v2) const;
	T IntegrateUV(double u1, double u2, double v1, double v2) const;
	
	// product
	template<class T1>
	BezSurf<T> Product(const BezSurf<T1>& b) const
	{
		// form the product knot sets
		KnotSet produ=(*ksetu).CreateKnotSetProduct(b.GetKnotSetU());
		KnotSet prodv=(*ksetv).CreateKnotSetProduct(b.GetKnotSetV());

		// return the new Surface
		return BezSurf<T>(ProductCPoints(b), produ.GetKnots(), prodv.GetKnots(), produ.GetOrd(), prodv.GetOrd());
	}

	template<class T1>
	BezSurf<T> Product1(const BezSurf<T1>& b) const
	{
		// build the knotset objects
	
		// form the product knot sets
		KnotSet produ=(*ksetu).CreateKnotSetProduct(b.GteKnotSetU());
		KnotSet prodv=(*ksetv).CreateKnotSetProduct(b.GetKnotSetV());

		int numcu = produ.GetNum() - produ.GetOrd();
		int numcv = prodv.GetNum() - prodv.GetOrd();

		Matrix<Vector<T> > mat1(ordu,b.GetNumU());	

		int numu = b.GetNumU();
		for (int i=0; i<ordu; i++) {
			BezCurv<T> temp1 = BezCurv<T>((*cpts).GetRow(i), *ktsv, ordv);
			for (int j=0; j<numu; j++) {
				BezCurv<T> temp2 = BezCurv<T>((b.GetCPoints()).GetRow(j),b.GetKnotsV(),b.GetOrdV());
				mat1[i][j] = temp1.ProductCPoints(temp2);
			}
		}


		Matrix<Vector<T> > mat2(ordu,numu);
		Vector<double> v2(ordu,0.0);
		Vector<double> v3(numu,0.0);

		for (int i=0; i<ordu; i++) {
			v2[i]=1.0;
			BezCurv<T> temp3 = BezCurv<T>(v2, *ktsu, ordu);
			for (int j=0; j<numu; j++) {
				v3[j] = 1.0;
				BezCurv<T> temp4 = BezCurv<T>(v3, b.GetKnotsU(), b.GetOrdU());
				mat2[i][j] = temp3.ProductCPoints(temp4);
			}
		}

		Matrix<T> ncpts(numcu,numcv);
	
		T sum=0.0;	// conversion required
		for (int i=0; i<numcu; i++)
			for (int j=0; j<numcv; j++) { 
				for (int i1=0; i1<ordu; i1++) 
					for (int j1=0; j1<b.GetOrdU(); j1++) sum = sum+(mat1[i1][j1])[j]*(mat2[i1][j1])[i];
				ncpts[i][j] = sum;
				sum = 0.0;
			}	
		return BezSurf<T>(ncpts, produ.GetKnots(), prodv.GetKnots(), produ.GetOrd(), prodv.GetOrd());
	}

	template<class T1>
	Matrix<T> ProductCPoints(const BezSurf<T1>& b) const
	{
		int ordbu = b.GetOrdU();
		int ordbv = b.GetOrdV();

		Matrix<T> temp(ordbu+ordu-1,ordbv+ordv-1);
	
		double sum =0.0;
		for (int k=0; k<ordbu+ordu-1; k++) {
			for (int l=0; l<ordbv+ordv-1; l++) {
				for (int i=Math::max1(0,k-ordbu+1); i<=Math::min1(ordu-1,k); i++) {
					for (int j=Math::max1(0,l-ordbv+1); j<=Math::min1(ordv-1,l); j++) {
						sum = sum + (Math::Combin(ordu-1,i)*Math::Combin(ordbu-1,k-i)*Math::Combin(ordv-1,j)*Math::Combin(ordbv-1,l-j)*(*cpts)[i][j]*b.GetCPoints()[k-i][l-j]/(Math::Combin(ordu+ordbu-2,k)*Math::Combin(ordv+ordbv-2,l)));
					}
				}
				temp[k][l]=sum;
				sum=0.0;
			}
		}
		return temp;
	}

	// subdivision
	BspSurf<T> Subdivide(int levu, int levv) const;
	BspSurf<T> SubdivideCPoints(int levu, int levv) const;
	BezSurf<T> Subdivide(double u1, double u2, double v1, double v2) const;
	Matrix<T> SubdivideCPoints(double u1, double u2, double v1, double v2) const;

	// insert knots
	Matrix<T> InsertKnotCPointsU(const Vector<double>& Kts,const Vector<int>& mult, int n) const;
	BspSurf<T> InsertKnotU(const Vector<double>& Kts, const Vector<int>& mult, int n) const;
	Matrix<T> InsertKnotCPointsV(const Vector<double>& Kts,const Vector<int>& mult, int n) const;
	BspSurf<T> InsertKnotV(const Vector<double>& Kts, const Vector<int>& mult, int n) const;
	
	// conversion functions
	PolySurf<T> ConvertPolySurf() const;
	BspSurf<T> ConvertBspSurf() const;

	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};

// CONSTRUCTORS

// constructor builds a BezSurf from a Vector of control points, knots
// an order and number of control points
template<class T>
BezSurf<T>::BezSurf(const Matrix<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv) :
			ordu(Ordu), ordv(Ordv), cpts(new Matrix<T>(Cpts)), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)),
			ksetu(new KnotSet(Ktsu,ordu,2*ordu)), ksetv(new KnotSet(Ktsv,ordv,2*ordv))
{
	leftLimitU = (*ktsu)[ordu-1];
	rightLimitU = (*ktsu)[ordu];
	leftLimitV = (*ktsv)[ordv-1];
	rightLimitV = (*ktsv)[ordv-1];
}

// default constructor
template<class T>
BezSurf<T>::BezSurf() : ordu(0), ordv(0), cpts(), ktsu(), ktsv(), ksetu(), ksetv() { }


// constructor builds a BezSurf from a Vector of control points, knots
// an order and number of control points
template<class T>
BezSurf<T>::BezSurf(const Matrix<T>& Cpts, int Ordu, int Ordv) :
			ordu(Ordu), ordv(Ordv), cpts(new Matrix<T>(Cpts))
{
	ktsu = new Vector<double>(Math::CreateKnots(ordu));
	ktsv = new Vector<double>(Math::CreateKnots(ordv));
	ksetu = new KnotSet(*ktsu,ordu,2*ordu);
	ksetv = new KnotSet(*ktsv,ordv,2*ordv);
	leftLimitU = (*ktsu)[ordu-1];
	rightLimitU = (*ktsu)[ordu];
	leftLimitV = (*ktsv)[ordv-1];
	rightLimitV = (*ktsv)[ordv-1];
}


// constructor builds a BezSurf from a Vector of control points, knots
// an order and number of control points
template<class T>
BezSurf<T>::BezSurf(const Matrix<T>& Cpts, int Ordu, int Ordv, double Lu, double Ru, double Lv, double Rv) :
			ordu(Ordu), ordv(Ordv), cpts(new Matrix<T>(Cpts))
{
	Vector<double> Limitsu(2), Limitsv(2);
	Limitsu[0] = Lu;
	Limitsu[1] = Ru;
	Limitsv[0] = Lv;
	Limitsv[1] = Rv;
	ktsu = new Vector<double>(Math::CreateKnots(1,ordu,Limitsu));
	ktsv = new Vector<double>(Math::CreateKnots(1,ordv,Limitsv));
	ksetu = new KnotSet(*ktsu,ordu,2*ordu);
	ksetv = new KnotSet(*ktsv,ordv,2*ordv);
	leftLimitU = Lu;
	rightLimitU = Ru;
	leftLimitV = Lv;
	rightLimitV = Rv;
}

// builds a BezCurv from an order
template<class T>
BezSurf<T>::BezSurf(int Ordu, int Ordv) : ordu(Ordu), ordv(Ordv), cpts(new Matrix<T>(Ordu,Ordv))
{ 
	// assign knots 0, 1
	ktsu = new Vector<double>(Math::CreateKnots(ordu));
	ktsv = new Vector<double>(Math::CreateKnots(ordv));
	ksetu = new KnotSet(*ktsu,ordu,2*ordu);
	ksetv = new KnotSet(*ktsv,ordv,2*ordv);
	leftLimitU = (*ktsu)[ordu-1];
	rightLimitU = (*ktsu)[ordu];
	leftLimitV = (*ktsv)[ordv-1];
	rightLimitV = (*ktsv)[ordv];
}



// ACCESS FUNCTIONS

// get the order of the BezSurf
template<class T>
inline int BezSurf<T>::GetOrdU() const { return ordu; }

// get the order of the BezSurf
template<class T>
inline int BezSurf<T>::GetOrdV() const { return ordv; }

// get the number of control points
template<class T>
inline int BezSurf<T>::GetNumU() const { return ordu; }

// get the number of control points
template<class T>
inline int BezSurf<T>::GetNumV() const { return ordv; }


// get the control points
template<class T>
inline Matrix<T> BezSurf<T>::GetCPoints() const { return *cpts; }

// get the knot vector
template<class T>
inline Vector<double> BezSurf<T>::GetKnotsU() const { return *ktsu; }

// get the knot vector
template<class T>
inline Vector<double> BezSurf<T>::GetKnotsV() const { return *ktsv; }

// get the knot vector
template<class T>
inline KnotSet BezSurf<T>::GetKnotSetU() const { return *ksetu; }

// get the knot vector
template<class T>
inline KnotSet BezSurf<T>::GetKnotSetV() const { return *ksetv; }

template<class T>
inline double BezSurf<T>::GetLeftLimitU() const { return leftLimitU; }

template<class T>
inline double BezSurf<T>::GetLeftLimitV() const { return leftLimitV; }


template<class T>
inline double BezSurf<T>::GetRightLimitU() const { return rightLimitU; }

template<class T>
inline double BezSurf<T>::GetRightLimitV() const { return rightLimitV; }



// ADDITION AND SUBTRACTION


// add two BezSurf's
template<class T>
BezSurf<T> BezSurf<T>::Add(const BezSurf<T>& b) const
{
	BezSurf<T> c(*this), d(b);
	
	if (ordu > b.GetOrdU()) d = d.ElevateU(ordu-b.GetOrdU());
	else if (b.GetOrdU() > ordu) c = c.ElevateU(b.GetOrdU()-ordu);
	if (ordv > b.GetOrdV()) d = d.ElevateV(ordv-b.GetOrdV());
	else if (b.GetOrdV() > ordv) c = c.ElevateV(b.GetOrdV()-ordv);

	Matrix<T> temp(b.GetOrdU(),b.GetOrdV());

	for (int i=0; i<b.GetOrdU(); i++)
		for (int j=0; j<b.GetOrdV(); j++)
			temp[i][j]=c.GetCPoints()[i][j]+d.GetCPoints()[i][j];

	return BezSurf<T>(temp,d.GetKnotsU(),d.GetKnotsV(),d.GetOrdU(),d.GetOrdV());
}


// subtract two BezSurf's
template<class T>
BezSurf<T> BezSurf<T>::Subtract(const BezSurf<T>& b) const
{
	BezSurf<T> c(*this), d(b);
	
	if (ordu > b.GetOrdU()) d = d.ElevateU(ordu-b.GetOrdU());
	else if (b.GetOrdU() > ordu) c = c.ElevateU(b.GetOrdU()-ordu);
	if (ordv > b.GetOrdV()) d = d.ElevateV(ordv-b.GetOrdV());
	else if (b.GetOrdV() > ordv) c = c.ElevateV(b.GetOrdV()-ordv);


	Matrix<T> temp(b.GetOrdU(),b.GetOrdV());

	for (int i=0; i<b.GetOrdU(); i++)
		for (int j=0; j<b.GetOrdV(); j++)
			temp[i][j]=c.GetCPoints()[i][j]-d.GetCPoints()[i][j];

	return BezSurf<T>(temp,d.GetKnotsU(),d.GetKnotsV(),d.GetOrdU(),d.GetOrdV());
}



// EVALUATION

// evaluate the BezSurf at the point x using de Boor algorithm
template<class T>
T BezSurf<T>::operator()(double u, double v) const
{
  	Vector<T> v2(ordu);
	
	for (int i=0; i<ordu; i++) 
		v2[i]= BezCurv<T>((*cpts).GetRow(i),*ktsv,ordv)(v);

	return BezCurv<T>(v2,*ktsu,ordu)(u);
}


// evaluate the BezSurf using matrix method
template<class T>
T BezSurf<T>::Eval(double u, double v) const
{
	return operator()(u,v);
}


// DEGREE ELEVATION

// elevate the degree of the BezSurf by level
template<class T>
BezSurf<T> BezSurf<T>::ElevateUV(int levu, int levv) const
{
	return ElevateU(levu).ElevateV(levv);
}

// elevate the degree of the BezSurf by level
template<class T>
Matrix<T> BezSurf<T>::ElevateCPointsUV(int levu, int levv) const
{
	// elevate in v
	Matrix<T> ncpts = ElevateCPointsV(levv);
	
	KnotSet kset = (*ksetv).CreateKnotSetElevate(levv);

	return BezSurf<T>(ncpts,*ktsu,kset.GetKnots(),ordu,ordv+levv).ElevateCPointsU(levu);
}


// DERIVATIVES

// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
BezSurf<T> BezSurf<T>::Derive(int levu, int levv) const
{
	return DeriveU(levu).DeriveV(levv);
}

// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
Matrix<T> BezSurf<T>::DeriveCPoints(int levu, int levv) const
{	
	// derive in v
	Matrix<T> ncpts = DeriveCPointsV(levv);
	KnotSet kset = (*ksetv).CreateKnotSetDerive(levv);

	return BezSurf<T>(ncpts,*ktsu,kset.GetKnots(),ordu,ordv-levv).DeriveCPointsU(levu);
}

template<class T>
T BezSurf<T>::operator() (int valu, int valv, double u, double v) const
{
	return Derive(valu,valv,u,v);
}

// evaluate the derivative of the BezSurf of order deriv at a point
// x. Computes the derivative as a BezSurf and then evaluates this at x
template<class T>
T BezSurf<T>::Derive(int levu, int levv, double u, double v) const
{
	return (DeriveU(levu).DeriveV(levv))(u,v);
}


// INTEGRATION


// integrate the BezSurf between the limits x1 and x2. Computes
// the indefinite integral as a BezSurf and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T BezSurf<T>::IntegrateUV(double u1, double u2, double v1, double v2) const
{
	BezSurf<T> temp = IntegrateUV();
	
	return temp(u2,v2)-temp(u1,v2)-temp(u2,v1)+temp(u1,v1);
}


// compute the indefinite integral of the BezSurf and represent
// it as a BezSurf of one higher degree
template<class T>
BezSurf<T> BezSurf<T>::IntegrateUV() const
{
	return IntegrateU().IntegrateV();
}	


// compute the indefinite integral of the BezSurf as a BezSurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> BezSurf<T>::IntegrateCPointsUV() const
{
	Matrix<T> ncpts = IntegrateCPointsV();
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();
	return BezSurf<T>(ncpts,*ktsu,kset.GetKnots(),ordu,ordv+1).IntegrateCPointsU();
}	


// PRODUCT
/*
// compute the product of the BezSurf with another BezSurf and 
// represent the result as a new BezSurf 
template<class T>  
BezSurf<T> BezSurf<T>::Product(const BezSurf<T>& b) const
{
	// form the product knot sets
	KnotSet produ=(*ksetu).CreateKnotSetProduct(b.GetKnotSetU());
	KnotSet prodv=(*ksetv).CreateKnotSetProduct(b.GetKnotSetV());

	// return the new Surface
	return BezSurf<T>(ProductCPoints(b), produ.GetKnots(), prodv.GetKnots(), produ.GetOrd(), prodv.GetOrd());
}   


// compute the product of the BezSurf with another BezSurf and 
// represent the result as a new BezSurf 
template<class T>  
Matrix<T> BezSurf<T>::ProductCPoints(const BezSurf<T>& b) const
{
	int ordbu = b.GetOrdU();
	int ordbv = b.GetOrdV();

	Matrix<T> temp(ordbu+ordu-1,ordbv+ordv-1);
	
	double sum =0.0;
	for (int k=0; k<ordbu+ordu-1; k++) {
		for (int l=0; l<ordbv+ordv-1; l++) {
			for (int i=Math::max1(0,k-ordbu+1); i<=Math::min1(ordu-1,k); i++) {
				for (int j=Math::max1(0,l-ordbv+1); j<=Math::min1(ordv-1,l); j++) {
					sum = sum + (Math::Combin(ordu-1,i)*Math::Combin(ordbu-1,k-i)*Math::Combin(ordv-1,j)*Math::Combin(ordbv-1,l-j)*(*cpts)[i][j]*b.GetCPoints()[k-i][l-j]/(Math::Combin(ordu+ordbu-2,k)*Math::Combin(ordv+ordbv-2,l)));
				}
			}
			temp[k][l]=sum;
			sum=0.0;
		}
	}
	return temp;
} 
*/  


//SUBDIVISION

template<class T>
BspSurf<T> BezSurf<T>::Subdivide(int levu, int levv) const
{
	return BspSurf<T>(*cpts,*ktsu,*ktsv,ordu,ordv,ordu,ordv).SubdivideUV(levu,levv);
}



// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
BezSurf<T> BezSurf<T>::Subdivide(double u1, double u2, double v1, double v2) const
{
	return SubdivideU(u1,u2).SubdivideV(v1,v2);
}

// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
Matrix<T> BezSurf<T>::SubdivideCPoints(double u1, double u2, double v1, double v2) const
{
	Matrix<T> ncpts = SubdivideV(v1,v2);
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);
	
	return BezSurf<T>(ncpts,*ktsu,kset.GetKnots(),ordu,ordv).SubdivideCPointsU(u1,u2);
}


// CONVERSION

template<class T>
PolySurf<T> BezSurf<T>::ConvertPolySurf() const
{
	Matrix<T> m = Math::mult2(Math::ComputePolyCurvMatrix(ordu),*cpts);
	Matrix<T> ncpts = Math::mult1(m,Math::ComputePolyCurvMatrixTranspose(ordv));
	//return PolySurf<T>(ncpts,ordu,ordv,0.0,1.0,0.0,1.0);
	return PolySurf<T>(ncpts,ordu,ordv,0.0,1.0,0.0,1.0).Reparameterise2(leftLimitU,rightLimitU,leftLimitV,rightLimitV);
}


template<class T>
BspSurf<T> BezSurf<T>::ConvertBspSurf() const
{
	return BspSurf<T>(*cpts, *ktsu, *ktsv, ordu, ordv, ordu, ordv);
}


// KNOT INSERTION

 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspSurf<T> BezSurf<T>::InsertKnotU(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// calculate the new knot set
	KnotSet kset = (*ksetu).CreateKnotSetInsert(t,mult,n);

	// create and return the new BspSurf
	return BspSurf<T>(InsertKnotCPointsU(t,mult,n),kset.GetKnots(),*ktsv,ordu,ordv,kset.GetNum()-ordu,ordv);
}
          
 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Matrix<T> BezSurf<T>::InsertKnotCPointsU(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// create a local copy of the BspSurf
	BspSurf temp = BezSurf(*cpts,*ktsu,*ktsv,ordu,ordv).InsertKnotU(t[0],mult[0]);
    
	for (int i=1; i<n; i++) temp = temp.InsertKnotV(t[i],mult[i]);	
    return temp.GetCPoints();
}

 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspSurf<T> BezSurf<T>::InsertKnotV(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// calculate the new knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(t,mult,n);

	// create and return the new BspSurf
	return BspSurf<T>(InsertKnotCPointsV(t,mult,n),ktsu,kset.GetKnots(),ordu,ordv,ordu,kset.GetNum()-ordv);
}
          
 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Matrix<T> BezSurf<T>::InsertKnotCPointsV(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// create a local copy of the BspSurf
	BspSurf temp = BezSurf(*cpts,*ktsu,*ktsv,ordu,ordv).InsertKnotV(t[0],mult[0]);
    
	for (int i=1; i<n; i++) temp = temp.InsertKnotV(t[i],mult[i]);	
    return temp.GetCPoints();
}


// READ and WRITE FUNCTIONS

template <class T>
void BezSurf<T>::write(std::ostream& os) 
{
	os << "Bezier Surface\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are\n";
	os << *ktsv;
	os << "\ncontrol points are\n";
	os << *cpts;
}

template <class T>
void BezSurf<T>::read(std::istream& is)
{
	int Ordu, Ordv;
	std::cout << "order of Bezier Surface in u and v";
	is >> Ordu >> Ordv;
	Matrix<T> Cpts(Ordu,Ordv);
	std::cout << "input control points";
	is >> Cpts;
	*this = BezSurf<T>(Cpts,Ordu,Ordv);
} 


template <class T>
void BezSurf<T>::writefile(std::ofstream& ofs) 
{
	ofs << "Bezier Surface\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are\n";
	ofs << *ktsv;
	
	ofs << "\ncontrol points are\n";
	ofs << *cpts;
}

template <class T>
void BezSurf<T>::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv;
	ifs >> Ordu >> Ordv;
	Matrix<T> Cpts(Ordu,Ordv);
	ifs >> Cpts;
	*this = BezSurf<T>(Cpts,Ordu,Ordv);
} 


// PRIVATE FUNCTIONS

// DEGREE ELEVATION

// elevate the degree of the BezSurf by level
template<class T>
BezSurf<T> BezSurf<T>::ElevateU(int level) const
{
	if (level <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetu).CreateKnotSetElevate(level);

	// return the new Bsurf
	return BezSurf<T>(ElevateCPointsU(level), kset.GetKnots(), *ktsv, ordu+level, ordv);	
}

// elevate the degree of the BezSurf by level
template<class T>
Matrix<T> BezSurf<T>::ElevateCPointsU(int level) const
{
	if (level <= 0) return *cpts;
	Matrix<T> ncpts(ordu+level,ordv);
	Vector<T> temp;

	for (int j=0; j<ordv; j++) {
		temp = BezCurv<T>((*cpts).GetCol(j),*ktsu,ordu).ElevateCPoints(level);
		// extract control points
		for (int i=0; i<ordu+level; i++) ncpts[i][j]=temp[i];
	}
	// return the new control points
	return ncpts;
}


// elevate the degree of the BezSurf by level
template<class T>
BezSurf<T> BezSurf<T>::ElevateV(int level) const
{
	if (level <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetElevate(level);

	// return the new Bsurface
	return BezSurf<T>(ElevateCPointsV(level), *ktsu, kset.GetKnots(), ordu, ordv+level);
}

// elevate the degree of the BezSurf by level
template<class T>
Matrix<T> BezSurf<T>::ElevateCPointsV(int level) const
{
	if (level <= 0) return *cpts;
	Matrix<T> ncpts(ordu,ordv+level);

	for (int i=0; i<ordu; i++) {
		Vector<T> temp = BezCurv<T>((*cpts).GetRow(i),*ktsv,ordv).ElevateCPoints(level);
		// extract control points
		for (int j=0; j<ordv+level; j++) ncpts[i][j]=temp[j];
	}
	// return the new control points
	return ncpts;
}


// DERIVATIVES

// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
BezSurf<T> BezSurf<T>::DeriveU(int level) const
{   
	if (level <= 0) return *this;
	
	// create knot set
	KnotSet kset = (*ksetu).CreateKnotSetDeriv(level);

	// return the new Bsurf
	return BezSurf<T>(DeriveCPointsU(level), kset.GetKnots(), *ktsv, ordu-level, ordv);
}


// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
Matrix<T> BezSurf<T>::DeriveCPointsU(int level) const
{   
	if (level <= 0) return *cpts;
	Matrix<T> ncpts(ordu-level,ordv);

	// returns control points only
	for (int j=0; j<ordv; j++) {
		Vector<T> temp = BezCurv<T>((*cpts).GetCol(j),*ktsu,ordu).DeriveCPoints(level);
		// extract control points
		for (int i=0; i<ordu-level; i++) ncpts[i][j]=temp[i];
	}
	// return the new control points
	return ncpts;
}



// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
BezSurf<T> BezSurf<T>::DeriveV(int level) const
{
	if (level <= 0) return *this;
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetDeriv(level);

	// return the new Bsurf
	return BezSurf<T>(DeriveCPointsV(level), *ktsu, kset.GetKnots(), ordu, ordv-level);
}


// compute the derivative of the BezSurf of order deriv and
// represent the result as another BezSurf
template<class T>
Matrix<T> BezSurf<T>::DeriveCPointsV(int level) const
{
	if (level <= 0) return *cpts;
	Matrix<T> ncpts(ordu,ordv-level);

	// derive in v
	for (int i=0; i<ordu; i++) {
		Vector<T> temp = BezCurv<T>((*cpts).GetRow(i),*ktsv,ordv).DeriveCPoints(level);
		// extract control points
		for (int j=0; j<ordv; j++) ncpts[i][j]=temp[j];
	}

	// return the new control points
	return ncpts;
}


// evaluate the derivative of the BezSurf of order deriv at a point
// x. Computes the derivative as a BezSurf and then evaluates this at x
template<class T>
T BezSurf<T>::DeriveU(int level, double u, double v) const
{
	return DeriveU(level)(u,v);
}

// evaluate the derivative of the BezSurf of order deriv at a point
// x. Computes the derivative as a BezSurf and then evaluates this at x
template<class T>
T BezSurf<T>::DeriveV(int level, double u, double v) const
{
	return DeriveV(level)(u,v);
}


// INTEGRATION

// integrate the BezSurf between the limits x1 and x2. Computes
// the indefinite integral as a BezSurf and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T BezSurf<T>::IntegrateUV1(double u1, double u2, double v1, double v2) const
{
	Vector<T> v(ordv);

	// should be integrate CPoints
	for (int j=0; j<ordv; j++)
		v[j] = BezCurv<T>(GetCol(j),*ktsu,ordu).Integrate(u1,u2);
	return BezCurv<T>(v,*ktsv,ordv).Integrate(v1,v2);
}

// compute the indefinite integral of the BezSurf and represent
// it as a BezSurf of one higher degree
template<class T>
BezSurf<T> BezSurf<T>::IntegrateU() const
{
	// create knot sets
	KnotSet kset = (*ksetu).CreateKnotSetIntegrate();
	
	return BezSurf<T>(IntegrateCPointsU(),kset.GetKnots(),*ktsv,ordu+1,ordv);
}	


// compute the indefinite integral of the BezSurf and represent
// it as a BezSurf of one higher degree
template<class T>
BezSurf<T> BezSurf<T>::IntegrateV() const
{
	// create knot sets
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();

	return BezSurf<T>(IntegrateCPointsV(),*ktsu,kset.GetKnots(),ordu,ordv+1);
}	

// compute the indefinite integral of the BezSurf as a BezSurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> BezSurf<T>::IntegrateCPointsU() const
{
	Matrix<T> ncpts(ordu+1,ordv);

	for (int j=0; j<ordv; j++) {
		Vector<T> v = BezCurv<T>((*cpts).GetCol(j),*ktsu,ordu).IntegrateCPoints();			
		for (int i=0; i<ordu+1; i++) ncpts[i][j]=v[i];
	}
	
	return ncpts;
}	

// compute the indefinite integral of the BezSurf as a BezSurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> BezSurf<T>::IntegrateCPointsV() const
{
	Matrix<T> ncpts(ordu,ordv+1);

	// use create cpoints integrate 
	for (int i=0; i<ordu; i++)  {
		Vector<T> v = BezCurv<T>((*cpts).GetRow(i),*ktsv,ordv).IntegrateCPoints();			
		for (int j=0; j<ordv+1; j++) ncpts[i][j]=v[j];
	}

	return ncpts;
}
/*
// PRODUCT FUNCTIONS

// compute the product of the BezSurf with another BezSurf and 
// represent the result as a new BezSurf 
template<class T>  
BezSurf<T> BezSurf<T>::Product1(const BezSurf<T>& b) const
{
	// build the knotset objects
	
	// form the product knot sets
	KnotSet produ=(*ksetu).CreateKnotSetProduct(b.GteKnotSetU());
	KnotSet prodv=(*ksetv).CreateKnotSetProduct(b.GetKnotSetV());

	int numcu = produ.GetNum() - produ.GetOrd();
	int numcv = prodv.GetNum() - prodv.GetOrd();

	Matrix<Vector<T> > mat1(ordu,b.GetNumU());	

	int numu = b.GetNumU();
	for (int i=0; i<ordu; i++) {
		BezCurv<T> temp1 = BezCurv<T>((*cpts).GetRow(i), *ktsv, ordv);
		for (int j=0; j<numu; j++) {
			BezCurv<T> temp2 = BezCurv<T>((b.GetCPoints()).GetRow(j),b.GetKnotsV(),b.GetOrdV());
			mat1[i][j] = temp1.ProductCPoints(temp2);
		}
	}


	Matrix<Vector<T> > mat2(ordu,numu);
	Vector<double> v2(ordu,0.0);
	Vector<double> v3(numu,0.0);

	for (int i=0; i<ordu; i++) {
		v2[i]=1.0;
		BezCurv<T> temp3 = BezCurv<T>(v2, *ktsu, ordu);
		for (int j=0; j<numu; j++) {
			v3[j] = 1.0;
			BezCurv<T> temp4 = BezCurv<T>(v3, b.GetKnotsU(), b.GetOrdU());
			mat2[i][j] = temp3.ProductCPoints(temp4);
		}
	}

	Matrix<T> ncpts(numcu,numcv);
	
	T sum=0.0;	// conversion required
	for (int i=0; i<numcu; i++)
		for (int j=0; j<numcv; j++) { 
			for (int i1=0; i1<ordu; i1++) 
				for (int j1=0; j1<b.GetOrdU(); j1++) sum = sum+(mat1[i1][j1])[j]*(mat2[i1][j1])[i];
			ncpts[i][j] = sum;
			sum = 0.0;
		}	
	return BezSurf<T>(ncpts, produ.GetKnots(), prodv.GetKnots(), produ.GetOrd(), prodv.GetOrd());
}   
*/

// SUBDIVISION

template<class T>
BspSurf<T> BezSurf<T>::SubdivideU(int level) const
{
	return BspSurf<T>(*cpts,*ktsu,*ktsv,ordu,ordv,ordu,ordv).SubdivideU(level);
}

template<class T>
Matrix<T> BezSurf<T>::SubdivideCPointsU(int level) const
{
	return BspSurf<T>(*cpts,*ktsu,*ktsv,ordu,ordv,ordu,ordv).SubdivideCPointsU(level);
}


template<class T>
BspSurf<T> BezSurf<T>::SubdivideV(int level) const
{
	return BspSurf<T>(*cpts,*ktsu,*ktsv,ordu,ordv,ordu,ordv).SubdivideV(level);
}

template<class T>
Matrix<T> BezSurf<T>::SubdivideCPointsV(int level) const
{
	return BspSurf<T>(*cpts,*ktsu,*ktsv,ordu,ordv,ordu,ordv).SubdivideCPointsV(level);
}

// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
BezSurf<T> BezSurf<T>::SubdivideU(double u1, double u2) const
{
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(u1,u2);
	Matrix<T> ncpts(ordu,ordv);
	
	// return the new Bsurface
	return BezSurf<T>(SubdivideCPointsU(u1,u2), kset.GetKnots(), *ktsv, ordu, ordv);
}

// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
Matrix<T> BezSurf<T>::SubdivideCPointsU(double u1, double u2) const
{
	Matrix<T> ncpts(ordu,ordv);

	for (int j=0; j<ordv; j++) {
		Vector<T> temp = BezCurv<T>((*cpts).GetCol(j),*ktsu,ordu).SubdivideCPoints(u1,u2);
		// extract control points
		for (int i=0; i<ordu; i++) ncpts[i][j]=temp[i];
	}
	// return the new control points
	return ncpts;
}


// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
BezSurf<T> BezSurf<T>::SubdivideV(double v1, double v2) const
{
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);

	// return the new Bsurface
	return BezSurf<T>(SubdivideCPointsV(v1,v2), *ktsu, kset.GetKnots(), ordu, ordv);

}


// subdivide the Bsurface upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the Bsurface after the knot refinement
template<class T>
Matrix<T> BezSurf<T>::SubdivideCPointsV(double v1, double v2) const
{
	Matrix<T> ncpts(ordu,ordv);

	for (int i=0; i<ordu; i++) {
		Vector<T> temp = BezCurv<T>((*cpts).GetRow(i),*ktsv,ordv).SubdivideCPoints(v1,v2);
		// extract control points
		for (int j=0; j<ordv; j++) ncpts[i][j]=temp[j];
	}
	// return the new Bsurface
	return ncpts;

}


// INSERT KNOTS

// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
Matrix<T> BezSurf<T>::InsertKnotCPointsU(double x) const
{       
	Matrix<T> ncpts(ordu+1,ordv);

	for (int j=0; j<ordv; j++) {
		Vector<T> temp = BezCurv<T>((*cpts).GetCol(j),*ktsu,ordu).InsertKnotCPoints(x);
		// extract control points
		for (int i=0; i<ordu+1; i++) ncpts[i][j]=temp[i];
	}
	// return new points
	return ncpts;
}                              
      
// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
BspSurf<T> BezSurf<T>::InsertKnotU(double x) const
{       
	// calculate the new knot set
	KnotSet kset = KnotSet(*ktsu,ordu,2*ordu).CreateKnotSetInsert(x);

	// create and return the new BspSurf
	return BspSurf<T>(InsertKnotCPointsU(x),kset.GetKnots(),*ktsv,ordu,ordv,ordu+1,ordv);
}              


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
Matrix<T> BezSurf<T>::InsertKnotCPointsU(double x, int level) const
{	
	if (level <= 0) return cpts;
	Matrix<T> ncpts(ordu+level,ordv);

	for (int j=0; j<ordv; j++) {
		Vector<T> temp = BezCurv<T>((*cpts).GetCol(j),*ktsu,ordu).InsertKnotCPoints(x,level);
		// extract control points
		for (int i=0; i<ordu+level; i++) ncpts[i][j]=temp[i];
	}
	// return new points
	return ncpts;
}


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
BspSurf<T> BezSurf<T>::InsertKnotU(double x, int level) const
{	
	if (level <= 0) return *this;
	// calculate the new knot set
	KnotSet kset = (*ksetu).CreateKnotSetInsert(x,level);

	// create and return the new BspSurf
	return BspSurf<T>(InsertKnotCPointsU(x,level),kset.GetKnots(),*ktsv,ordu,ordv,ordu+level,ordv);
}


// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
BspSurf<T> BezSurf<T>::InsertKnotV(double x) const
{       
	// calculate the new knot set
	KnotSet kset = KnotSet(*ktsv,ordv,2*ordv).CreateKnotSetInsert(x);

	// create and return the new BspSurf
	return BspSurf<T>(InsertKnotCPointsV(x),*ktsu,kset.GetKnots(),ordu,ordv,ordu,ordv+1);
}              


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
Matrix<T> BezSurf<T>::InsertKnotCPointsV(double x, int level) const
{	
	if (level <= 0) return *this;
	Matrix<T> ncpts(ordu,ordv+level);

	for (int i=0; i<ordu; i++) {
		Vector<T> temp = BezCurv<T>((*cpts).GetRow(i),*ktsv,ordv).InsertKnotCPoints(x,level);
		// extract control points
		for (int j=0; j<ordv+level; j++) ncpts[i][j]=temp[j];
	}
	// return new points
	return ncpts;
}


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
BspSurf<T> BezSurf<T>::InsertKnotV(double x, int level) const
{	
	if (level <= 0) return *this;
	// calculate the new knot set
	KnotSet kset = KnotSet(*ktsv,ordv,2*ordv).CreateKnotSetInsert(x,level);

	// create and return the new BspSurf
	return BspSurf<T>(InsertKnotCPointsV(x,level),*ktsu,kset.GetKnots(),ordu,ordv,ordu,ordv+level);
}




#endif

