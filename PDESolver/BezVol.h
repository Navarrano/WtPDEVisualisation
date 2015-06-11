
#ifndef BEZVOL
#define BEZVOL

#include "BezSurf.h"
#include "bspvol.h"
template<class T>
class BezVol;

template<class T>
class BspVol;

template<class T>
class PolyVol;

class BezVolBasisFunc : public Vol<double> {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int ordw;	// order in w
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<BezVol<double> > b;	// BezVol representation of basis function
	
	// private functions
	BezVol<double> CreateBezVol() const;	// creates the BezVol
	virtual ObjectID Identity() const { return std::string("class BezVolBasisFunc"); }
public:

	// constructors
	BezVolBasisFunc();
	BezVolBasisFunc(int Ordu, int Ordv, int Ordw, const Vector<double>& Ktsu, const Vector<double>& Ktsv,
			const Vector<double>& Ktsw);

	// evaluators
	double Eval(double u, double v, double w) const;	// evaluate basis function
	virtual double operator()(double u, double v, double w) const;
	virtual double operator()(int, int, int, double u, double v, double w) const;

	// access functions
	BezVol<double> GetBezVol() const;		// get BezVol representation
	int ComputeDimU() const;
	int ComputeDimV() const;
	int ComputeDimW() const;

	
	virtual double Derive(int, int , int, double, double, double) const;
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
	virtual double GetLeftLimitW() const;
	virtual double GetRightLimitW() const;
	
	// read and write

	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};


class BezVolBasisFuncSet : public TextObject, public FTextObject {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int ordw;
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<Matrix3D<BezVolBasisFunc> > b;
	virtual ObjectID Identity() const { return std::string("class BezVolBasisFuncSet"); }
public:

	// constructors
	BezVolBasisFuncSet();
	BezVolBasisFuncSet(int Ordu, int Ordv, int Ordw, 
		const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw);

	// evaluators
	Matrix3D<double> Eval(double u, double v, double w) const;	// evaluate basis function
	Matrix3D<double> operator()(double u, double v, double w) const;
	Matrix3D<double> operator()(int, int, int, double u, double v, double w) const;
	Matrix3D<double> Derive(int, int, int, double u, double v, double w) const;
	
	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};


//#include "headers.h"

template<class T>
class BezVol : public Vol<T> {
private:
	// data
	int ordu;
	int ordv;
	int ordw;
	Ptr<Matrix3D<T> > cpts;
	Ptr<Vector<double> > ktsu;
	Ptr<Vector<double> > ktsv;
	Ptr<Vector<double> > ktsw;
	Ptr<KnotSet> ksetu;
	Ptr<KnotSet> ksetv;
	Ptr<KnotSet> ksetw;
	double leftLimitU;
	double rightLimitU;
	double leftLimitV;
	double rightLimitV;
	double leftLimitW;
	double rightLimitW;

	// private functions
	// derivatives
	BezVol<T> DeriveUV(int levu, int levv) const; 
	BezVol<T> DeriveUW(int levu, int levw) const; 
	BezVol<T> DeriveVW(int levv, int levw) const;


	
	T DeriveU(int level, double u, double v, double w) const;
	T DeriveV(int level, double u, double v, double w) const;
	T DeriveW(int level, double u, double v, double w) const;
	T DeriveUV(int levu, int levv, double u, double v, double w) const;
	T DeriveUW(int levu, int levw, double u, double v, double w) const;
	T DeriveVW(int levv, int levw, double u, double v, double w) const;
//	T Derive(int levu, int levv, int levw, double u, double v, double w) const;
	
	// degree elevation

	BezVol<T> ElevateUV(int levu, int levv) const;
	BezVol<T> ElevateUW(int levu, int levw) const;
	BezVol<T> ElevateVW(int levv, int levw) const;

	// integration
	BezVol<T> IntegrateUV() const;
	BezVol<T> IntegrateUW() const;
	BezVol<T> IntegrateVW() const;
	BezVol<T> IntegrateU() const;
	BezVol<T> IntegrateV() const;
	BezVol<T> IntegrateW() const;
	Matrix3D<T> IntegrateCPointsU() const;
	Matrix3D<T> IntegrateCPointsV() const;
	Matrix3D<T> IntegrateCPointsW() const;
	Matrix3D<T> IntegrateCPointsUV() const;
	Matrix3D<T> IntegrateCPointsUW() const;
	Matrix3D<T> IntegrateCPointsVW() const;
	T IntegrateUVW2(double u1, double u2, double v1, double v2, double w1, double w2) const;
	T IntegrateUVW1(double u1, double u2, double v1, double v2, double w1, double w2) const;

	// product
	BezVol<T> Product2(const BezVol<T>& b) const;

	// insert knots
	Matrix3D<T> InsertKnotCPointsU(double x) const;
	BspVol<T> InsertKnotU(double x) const;
	BspVol<T> InsertKnotU(double x, int level) const;
	Matrix3D<T> InsertKnotCPointsU(double x, int level) const;
	BspVol<T> InsertKnotU(const Vector<double>& Kts, int n) const;
	Matrix3D<T> InsertKnotCPointsU(const Vector<double>& Kts, int n) const;
	BspVol<T> InsertKnotV(double x) const;
	BspVol<T> InsertKnotV(double x, int level) const;
	Matrix3D<T> InsertKnotCPointsV(double x, int level) const;
	BspVol<T> InsertKnotV(const Vector<double>& Kts, int n) const;
	Matrix3D<T> InsertKnotCPointsV(const Vector<double>& Kts, int n) const;
	BspVol<T> InsertKnotW(double x) const;
	BspVol<T> InsertKnotW(double x, int level) const;
	Matrix3D<T> InsertKnotCPointsW(double x, int level) const;
	BspVol<T> InsertKnotW(const Vector<double>& Kts, int n) const;
	Matrix3D<T> InsertKnotCPointsW(const Vector<double>& Kts, int n) const;
	
	// subdivision
	BspVol<T> SubdivideU(int level) const;
	BspVol<T> SubdivideV(int level) const;
	BspVol<T> SubdivideW(int level) const;
	Matrix3D<T> SubdivideCPointsU(int level) const;
	Matrix3D<T> SubdivideCPointsV(int level) const;
	Matrix3D<T> SubdivideCPointsW(int level) const;
	Matrix3D<T> SubdivideCPointsUV(int levu, int levv) const;
	Matrix3D<T> SubdivideCPointsUW(int levu, int levw) const;
	Matrix3D<T> SubdivideCPointsVW(int levv, int levw) const;
	BspVol<T> SubdivideUV(int levu, int levv) const;
	BspVol<T> SubdivideUW(int levu, int levw) const;
	BspVol<T> SubdivideVW(int levv, int levw) const;
	BezVol<T> SubdivideU(double u1, double u2) const;
	BezVol<T> SubdivideV(double v1, double v2) const;
	BezVol<T> SubdivideW(double w1, double w2) const;
	Matrix3D<T> SubdivideCPointsU(double u1, double u2) const;
	Matrix3D<T> SubdivideCPointsV(double v1, double v2) const;
	Matrix3D<T> SubdivideCPointsW(double w1, double w2) const;
	BezVol<T> SubdivideUV(double u1, double u2, double v1, double v2) const;
	BezVol<T> SubdivideUW(double u1, double u2, double w1, double w2) const;
	BezVol<T> SubdivideVW(double v1, double v2, double w1, double w2) const;
	Matrix3D<T> SubdivideCPointsUV(double u1, double u2, double v1, double v2) const;
	Matrix3D<T> SubdivideCPointsUW(double u1, double u2, double w1, double w2) const;
	Matrix3D<T> SubdivideCPointsVW(double v1, double v2, double w1, double w2) const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class BezVol<double>");
		else {
			std::string s(typeid(T).name()), s1("class BezVol<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	BezVol();
	BezVol(const Matrix3D<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw);
	BezVol(const Matrix3D<T>& Cpts, int Ordu, int Ordv, int Ordw);
	BezVol(const Matrix3D<T>& Cpts, int Ordu, int Ordv, int Ordw, double Lu, double Ru, double Lv, double Rv, double Lw, double Rw);
	BezVol(int Ordu, int Ordv, int Ordw);

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
	
	
	// evaluation
	virtual T operator()(double u, double v, double w) const;
	virtual T operator()(int levu, int levv, int levw, double u, double v, double w) const;
	virtual T Derive(int levu, int levv, int levw, double u, double v, double w) const;
	T Eval(double u, double v, double w) const;
	T Eval1(double u, double v, double w) const;
	Matrix3D<T> Eval2(double u, double v, double w) const;
	Matrix3D<T> EvalDeriv(int levu, int levv, int levw, double u, double v, double w) const;



	// addition and subtraction
	BezVol<T> Add(const BezVol<T>& b) const;
	BezVol<T> Subtract(const BezVol<T>& b) const;

	// derivatives
	BezVol<T> DeriveU(int levu) const; 
	BezVol<T> DeriveV(int levv) const; 
	BezVol<T> DeriveW(int levw) const;
	BezVol<T> Derive(int levu, int levv, int levw) const;
	Matrix3D<T> DeriveCPointsU(int level) const;
	Matrix3D<T> DeriveCPointsV(int level) const;
	Matrix3D<T> DeriveCPointsW(int level) const;
	
	// elevation
	BezVol<T> ElevateUVW(int levu, int levv, int levw) const;
	Matrix3D<T> ElevateCPoints(int levu, int levv, int levw) const;	
	BezVol<T> ElevateU(int level) const;
	BezVol<T> ElevateV(int level) const;
	BezVol<T> ElevateW(int level) const;
	Matrix3D<T> ElevateCPointsU(int level) const;
	Matrix3D<T> ElevateCPointsV(int level) const;
	Matrix3D<T> ElevateCPointsW(int level) const;

	// integration
	BezVol<T> Integrate() const;
	Matrix3D<T> IntegrateCPoints() const;
	T Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const;
	
	// product
	template<class T1>
	BezVol<T> Product(const BezVol<T1>& b) const
	{
		// build the knots sets
		// build the knotset objects
	
		// form the product knot sets
		KnotSet produ=(*ksetu).CreateKnotSetProduct(b.GetKnotSetU());
		KnotSet prodv=(*ksetv).CreateKnotSetProduct(b.GetKnotSetV());
		KnotSet prodw=(*ksetw).CreateKnotSetProduct(b.GetKnotSetW());

		// return the new Surface
		return BezVol<T>(ProductCPoints(b), produ.GetKnots(), prodv.GetKnots(), prodw.GetKnots(), produ.GetOrd(), prodv.GetOrd(), prodw.GetOrd());
	}
	template<class T1>
	Matrix3D<T> ProductCPoints(const BezVol<T1>& b) const
	{
		int ordbu = b.GetOrdU();
		int ordbv = b.GetOrdV();
		int ordbw = b.GetOrdW();

		Matrix3D<T> temp(ordbu+ordu-1,ordbv+ordv-1,ordbw+ordw-1);
	
		T sum =0.0;
		for (int k=0; k<ordbu+ordu-1; k++) {
			for (int l=0; l<ordbv+ordv-1; l++) {
				for (int m=0; m<ordbw+ordw-1; m++) {
					for (int s=Math::max1(0,m-ordbw+1); s<=Math::min1(ordw-1,m); s++) {
						for (int i=Math::max1(0,k-ordbu+1); i<=Math::min1(ordu-1,k); i++) {
							for (int j=Math::max1(0,l-ordbv+1); j<=Math::min1(ordv-1,l); j++) {
								sum = sum + (*cpts)[s][i][j]*b.GetCPoints()[m-s][k-i][l-j]*(Math::Combin(ordu-1,i)*Math::Combin(ordbu-1,k-i)*Math::Combin(ordv-1,j)*Math::Combin(ordbv-1,l-j)*Math::Combin(ordw-1,s)*Math::Combin(ordbw-1,m-s)/(Math::Combin(ordu+ordbu-2,k)*Math::Combin(ordv+ordbv-2,l)*Math::Combin(ordw+ordbw-2,m)));
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

	virtual double GetLeftLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetLeftLimitW() const;

	virtual double GetRightLimitU() const;
	virtual double GetRightLimitV() const;
	virtual double GetRightLimitW() const;
	
	// insert knots
	BspVol<T> InsertKnotU(const Vector<double>& Kts, const Vector<int>& mult, int n) const;
	Matrix3D<T> InsertKnotCPointsU(const Vector<double>& Kts, const Vector<int>& mult, int n) const;	
	BspVol<T> InsertKnotV(const Vector<double>& Kts, const Vector<int>& mult, int n) const;
	Matrix3D<T> InsertKnotCPointsV(const Vector<double>& Kts, const Vector<int>& mult, int n) const;
	BspVol<T> InsertKnotW(const Vector<double>& Kts, const Vector<int>& mult, int n) const;
	Matrix3D<T> InsertKnotCPointsW(const Vector<double>& Kts, const Vector<int>& mult, int n) const;
	
	// subdivision
	Matrix3D<T> SubdivideCPoints(int levu, int levv, int levw) const;
	BspVol<T> Subdivide(int levu, int levv, int levw) const;	
	BezVol<T> Subdivide(double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<T> SubdivideCPoints(double u1, double u2, double v1, double v2, double w1, double w2) const;

	// conversion
	PolyVol<T> ConvertPolyVol() const;
	BspVol<T> ConvertBspVol() const;

	// read and write
	
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);

};

// CONSTRUCTORS

// constructor builds a BezVol from a Vector of control points, knots
// an order and number of control points
template<class T>
BezVol<T>::BezVol(const Matrix3D<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw) :
			ordu(Ordu), ordv(Ordv), ordw(Ordw), cpts(new Matrix3D<T>(Cpts)), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)), ktsw(new Vector<double>(Ktsw)),
			ksetu(new KnotSet(*ktsu,ordu,2*ordu)),ksetv(new KnotSet(*ktsv,ordv,2*ordv)),
			ksetw(new KnotSet(*ktsw,ordw,2*ordw))
{
}

// default constructor
template<class T>
BezVol<T>::BezVol() : ordu(0), ordv(0), ordw(0), cpts(), ktsu(), ktsv(), ktsw() { }


// constructor builds a BezVol from a Vector of control points, knots
// an order and number of control points
template<class T>
BezVol<T>::BezVol(const Matrix3D<T>& Cpts, int Ordu, int Ordv, int Ordw) :
			ordu(Ordu), ordv(Ordv), ordw(Ordw), cpts(new Matrix3D<T>(Cpts))
{
	ktsu = new Vector<double>(Math::CreateKnots(Ordu));
	ktsv = new Vector<double>(Math::CreateKnots(Ordv));
	ktsw = new Vector<double>(Math::CreateKnots(Ordw));
	ksetu = new KnotSet(*ktsu,ordu,2*ordu);
	ksetv = new KnotSet(*ktsv,ordv,2*ordv);
	ksetw = new KnotSet(*ktsw,ordw,2*ordw);
	leftLimitU = (*ktsu)[ordu-1];
	rightLimitU = (*ktsu)[ordu];
	leftLimitV = (*ktsv)[ordv-1];
	rightLimitV = (*ktsv)[ordv-1];
	leftLimitW = (*ktsw)[ordw-1];
	rightLimitW = (*ktsw)[ordw-1];
}

// constructor builds a BezVol from a Vector of control points, knots
// an order and number of control points
template<class T>
BezVol<T>::BezVol(const Matrix3D<T>& Cpts, int Ordu, int Ordv, int Ordw, double Lu, double Ru, double Lv, double Rv, double Lw, double Rw) :
			ordu(Ordu), ordv(Ordv), ordw(Ordw), cpts(new Matrix3D<T>(Cpts))
{
	Vector<double> Limitsu(2), Limitsv(2), Limitsw(2);
	Limitsu[0] = Lu;
	Limitsu[1] = Ru;
	Limitsv[0] = Lv;
	Limitsv[1] = Rv;
	Limitsw[0] = Lw;
	Limitsw[1] = Rw;
	ktsu = new Vector<double>(Math::CreateKnots(1, Ordu, Limitsu));
	ktsv = new Vector<double>(Math::CreateKnots(1, Ordv, Limitsv));
	ktsw = new Vector<double>(Math::CreateKnots(1, Ordw, Limitsw));
	ksetu = new KnotSet(*ktsu,ordu,2*ordu);
	ksetv = new KnotSet(*ktsv,ordv,2*ordv);
	ksetw = new KnotSet(*ktsw,ordw,2*ordw);
	leftLimitU = (*ktsu)[ordu-1];
	rightLimitU = (*ktsu)[ordu];
	leftLimitV = (*ktsv)[ordv-1];
	rightLimitV = (*ktsv)[ordv];
	leftLimitW = (*ktsw)[ordw-1];
	rightLimitW = (*ktsw)[ordw];
}


// builds a BezCurv from an order
template<class T>
BezVol<T>::BezVol(int Ordu, int Ordv, int Ordw) : ordu(Ordu), ordv(Ordv), ordw(Ordw), cpts(new Matrix3D<T>(Ordu,Ordv,Ordw))
{ 
	// assign knots 0, 1
	ktsu = new Vector<double>(Math::CreateKnots(ordu));
	ktsv = new Vector<double>(Math::CreateKnots(ordv));
	ktsw = new Vector<double>(Math::CreateKnots(ordw));
	ksetu = new KnotSet(*ktsu,ordu,2*ordu);
	ksetv = new KnotSet(*ktsv,ordv,2*ordv);
	ksetw = new KnotSet(*ktsw,ordw,2*ordw);
	leftLimitU = (*ktsu)[ordu-1];
	rightLimitU = (*ktsu)[ordu];
	leftLimitV = (*ktsv)[ordv-1];
	rightLimitV = (*ktsv)[ordv];
	leftLimitW = (*ktsw)[ordw-1];
	rightLimitW = (*ktsw)[ordw];
}

// KNOT CREATION



// ACCESS FUNCTIONS

// get the order of the BezVol
template<class T>
inline int BezVol<T>::GetOrdU() const { return ordu; }

// get the order of the BezVol
template<class T>
inline int BezVol<T>::GetOrdV() const { return ordv; }

// get the order of the BezVol
template<class T>
inline int BezVol<T>::GetOrdW() const { return ordw; }

// get the number of control points
template<class T>
inline int BezVol<T>::GetNumU() const { return ordu; }

// get the number of control points
template<class T>
inline int BezVol<T>::GetNumV() const { return ordv; }

// get the number of control points
template<class T>
inline int BezVol<T>::GetNumW() const { return ordw; }

// get the control points
template<class T>
inline Matrix3D<T> BezVol<T>::GetCPoints() const { return *cpts; }

// get the knot vector
template<class T>
inline Vector<double> BezVol<T>::GetKnotsU() const { return *ktsu; }

// get the knot vector
template<class T>
inline Vector<double> BezVol<T>::GetKnotsV() const { return *ktsv; }

// get the knot vector
template<class T>
inline Vector<double> BezVol<T>::GetKnotsW() const { return *ktsw; }

// get the knot vector
template<class T>
inline KnotSet BezVol<T>::GetKnotSetU() const { return *ksetu; }

// get the knot vector
template<class T>
inline KnotSet BezVol<T>::GetKnotSetV() const { return *ksetv; }

// get the knot vector
template<class T>
inline KnotSet BezVol<T>::GetKnotSetW() const { return *ksetw; }

// get the knot vector
template<class T>
double BezVol<T>::GetLeftLimitU() const { return (*ktsu)[0]; }


// get the knot vector
template<class T>
double BezVol<T>::GetLeftLimitV() const { return (*ktsv)[0]; }

// get the knot vector
template<class T>
double BezVol<T>::GetLeftLimitW() const { return (*ktsw)[0]; }


// get the knot vector
template<class T>
double BezVol<T>::GetRightLimitU() const { return (*ktsu)[ordu]; }


// get the knot vector
template<class T>
double BezVol<T>::GetRightLimitV() const { return (*ktsv)[ordv]; }

// get the knot vector
template<class T>
double BezVol<T>::GetRightLimitW() const { return (*ktsw)[ordw]; }

// ADDITION & SUBTRACTION

template<class T>
BezVol<T> BezVol<T>::Add(const BezVol<T>& b) const
{
	BezVol<T> c(*this), d(b);
	
	if (ordu > b.GetOrdU()) d = b.ElevateU(ordu-b.GetOrdU())
	else if (b.GetOrdU() > ordu) c = ElevateU(b.GetOrdU()-ordu);
	if (ordv > b.GetOrdV()) d = d.ElevateV(ordv-b.GetOrdV());
	else if (b.GetOrdV() > ordv) c = c.ElevateV(b.GetOrdV()-ordv);
	if (ordw > b.GetOrdW()) d = d.ElevateW(ordw-b.GetOrdW());
	else if (b.GetOrdW() > ordw) c = c.ElevateW(b.GetOrdW()-ordw);

	Matrix3D<T> temp(b.GetOrdU(),b.GetOrdV(),b.GetOrdW());

	for (int i=0; i<b.GetOrdU(); i++)
		for (int j=0; j<b.GetOrdV(); j++)
			for (int k=0; k<b.GetOrdW(); k++)
				temp[k][i][j]=c.GetCPoints()[k][i][j]+d.GetCPoints()[k][i][j];

	return BezVol<T>(temp,d.GetKnotsU(),d.GetKnotsV(),d.GetKnotsW(),d.GetOrdU(),d.GetOrdV(),d.GetOrdW());
}

template<class T>
BezVol<T> BezVol<T>::Subtract(const BezVol<T>& b) const
{
	BezVol<T> c(*this), d(b);
	
	if (ordu > b.GetOrdU()) d = b.ElevateU(ordu-b.GetOrdU());
	else if (b.GetOrdU() > ordu) c = ElevateU(b.GetOrdU()-ordu);
	if (ordv > b.GetOrdV()) d = d.ElevateV(ordv-b.GetOrdV());
	else if (b.GetOrdV() > ordv) c = c.ElevateV(b.GetOrdV()-ordv);
	if (ordw > b.GetOrdW()) d = d.ElevateW(ordw-b.GetOrdW());
	else if (b.GetOrdW() > ordw) c = c.ElevateW(b.GetOrdW()-ordw);

	Matrix3D<T> temp(b.GetOrdU(),b.GetOrdV(),b.GetOrdW());

	for (int i=0; i<b.GetOrdU(); i++)
		for (int j=0; j<b.GetOrdV(); j++)
			for (int k=0; k<b.GetOrdW(); k++)
				temp[k][i][j]=c.GetCPoints()[k][i][j]-d.GetCPoints()[k][i][j];

	return BezVol<T>(temp,d.GetKnotsU(),d.GetKnotsV(),d.GetKnotsW(),d.GetOrdU(),d.GetOrdV(),d.GetOrdW());
}


// EVALUATORS

// evaluate the BezVol at the point x using de Boor algorithm
template<class T>
T BezVol<T>::operator()(double u, double v, double w) const
{
  	Vector<T> v2(ordw);

	for (int k=0; k<ordw; k++) v2[k]= BezSurf<T>((*cpts).GetUV(k),*ktsu,*ktsv,ordu,ordv)(u,v);

    return BezCurv<T>(v2,*ktsw,ordw)(w);
}


// evaluate the BezVol using matrix method 
template<class T>
T BezVol<T>::Eval(double u, double v, double w) const
{
	return operator()(u,v,w);
}

// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
T BezVol<T>::Eval1(double u, double v, double w) const
{
	Vector<double> uvec = ksetu->CreateVectorInterp(u);
	Vector<double> vvec = ksetv->CreateVectorInterp(v);
	Vector<double> wvec = ksetw->CreateVectorInterp(w);

	return Math::mult0(Math::mult3(Math::mult7(uvec,*cpts),vvec),wvec);
}


// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Matrix3D<T> BezVol<T>::Eval2(double u, double v, double w) const
{
	Matrix<double> mu(ksetu->CreateVectorInterp(u));
	Matrix<double> mv(ksetv->CreateVectorInterp(v));
	Matrix<double> mw(ksetw->CreateVectorInterp(w));
	//return Math::mult5(Math::mult5(Math::mult6(mu,*cpts),Math::transpose(mv)),Math::transpose(mw));
	return Math::mult9(mw,Math::mult5(Math::mult6(mu,*cpts),Math::transpose(mv)));
}

// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Matrix3D<T> BezVol<T>::EvalDeriv(int levu, int levv, int levw, double u, double v, double w) const
{
	Matrix<double> mu(ksetu->CreateVectorInterpDeriv(levu,u));
	Matrix<double> mv(ksetv->CreateVectorInterpDeriv(levv,v));
	Matrix<double> mw(ksetw->CreateVectorInterpDeriv(levw,w));

	Matrix3D<T> m1 = Math::mult6(Math::mult1(mu,ksetu->CreateMatrixDeriv(levu)),*cpts);
	Matrix3D<T> m2 = Math::mult5(m1,Math::transpose(Math::mult1(mv,ksetv->CreateMatrixDeriv(levv))));
	return Math::mult9(Math::mult1(mw,ksetw->CreateMatrixDeriv(levw)),m2);
}


/*
// evaluate the BezVol using EXPLICIT matrix method 
template<class T>
Matrix3D<T> BezVol<T>::Eval3(double u, double v, double w) const
{
	Matrix<double> mu(ksetu->CreateVectorInterp(u));
	Matrix<double> mv(ksetv->CreateVectorInterp(v));
	Matrix<double> mw(ksetw->CreateVectorInterp(w));
	//return Math::mult5(Math::mult5(Math::mult6(mu,*cpts),Math::transpose(mv)),Math::transpose(mw));
	return Math::mult5(Math::mult5(Math::mult6(mu,*cpts),Math::transpose(mw)),Math::transpose(mv));
}
*/


// DEGREE ELEVATION



// elevate the degree of the BezVol by level
template<class T>
Matrix3D<T> BezVol<T>::ElevateCPoints(int levu, int levv, int levw) const
{
	Matrix3D<T> mat1 = ElevateCPointsW(levw);
	KnotSet kset1 = (*ksetw).CreateKnotSetElevate(levw);

	Matrix3D<T> mat2 = BezVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw+1).ElevateCPointsV(levv);
	KnotSet kset2 = (*ksetv).CreateKnotSetElevate(levv);

	return BezVol<T>(mat2,*ktsu,kset2.GetKnots(),kset1.GetKnots(),ordu,ordv+levv,ordw+levw).ElevateCPointsU(levu);

}


// DERIVATIVES

template<class T>
T BezVol<T>::operator() (int valu, int valv, int valw, double u, double v, double w) const
{
	return Derive(valu,valv,valw,u,v,w);
}


// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
BezVol<T> BezVol<T>::Derive(int levu, int levv, int levw) const
{
	return (DeriveU(levu).DeriveV(levv)).DeriveW(levw);
}


// evaluate the derivative of the BezVol of order deriv at a point
// x. Computes the derivative as a BezVol and then evaluates this at x
template<class T>
T BezVol<T>::Derive(int levu, int levv, int levw, double u, double v, double w) const
{
	return (DeriveU(levu).DeriveV(levv)).DeriveW(levw)(u,v,w);
}


// INTEGRATION

// integrate the BezVol between the limits x1 and x2. Computes
// the indefinite integral as a BezVol and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T BezVol<T>::Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	BezVol<T> temp = IntegrateUVW();
	
	return temp(u2,v2,w2)-temp(u2,v2,w1)-temp(u2,v1,w2)-temp(u1,v2,w2)+
		temp(u2,v1,w1)+temp(u1,v2,w1)+temp(u1,v1,w2)-temp(u1,v1,w1);
}



// compute the indefinite integral of the BezVol as a BezVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> BezVol<T>::IntegrateCPoints() const
{
	Matrix3D<T> mat1 = IntegrateCPointsW();
	KnotSet kset1 = (*ksetw).CreateKnotSetIntegrate();

	Matrix3D<T> mat2 = BezVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw+1).IntegrateCPointsV();
	KnotSet kset2 = (*ksetv).CreateKnotSetIntegrate();

	return BezVol<T>(mat2,*ktsu,kset2.GetKnots(),kset1.GetKnots(),ordu,ordv+1,ordw+1).IntegrateCPointsU();
}	

/*
// PRODUCT

// compute the product of the BezVol with another BezVol and 
// represent the result as a new BezVol 
template<class T>  
BezVol<T> BezVol<T>::Product(const BezVol<T>& b) const
{
	// build the knots sets
	// build the knotset objects
	
	// form the product knot sets
	KnotSet produ=(*ksetu).CreateKnotSetProduct(b.GetKnotSetU());
	KnotSet prodv=(*ksetv).CreateKnotSetProduct(b.GetKnotSetV());
	KnotSet prodw=(*ksetw).CreateKnotSetProduct(b.GetKnotSetW());

	// return the new Surface
	return BezVol<T>(ProductCPoints(b), produ.GetKnots(), prodv.GetKnots(), prodw.GetKnots(), produ.GetOrd(), prodv.GetOrd(), prodw.GetOrd());
}   


// compute the product of the BezVol with another BezVol and 
// represent the result as a new BezVol 
template<class T>  
Matrix3D<T> BezVol<T>::ProductCPoints(const BezVol<T>& b) const
{
	int ordbu = b.GetOrdU();
	int ordbv = b.GetOrdV();
	int ordbw = b.GetOrdW();

	Matrix3D<T> temp(ordbu+ordu-1,ordbv+ordv-1,ordbw+ordw-1);
	
	T sum =0.0;
	for (int k=0; k<ordbu+ordu-1; k++) {
		for (int l=0; l<ordbv+ordv-1; l++) {
			for (int m=0; m<ordbw+ordw-1; m++) {
				for (int s=Math::max1(0,m-ordbw+1); s<=Math::min1(ordw-1,m); s++) {
					for (int i=Math::max1(0,k-ordbu+1); i<=Math::min1(ordu-1,k); i++) {
						for (int j=Math::max1(0,l-ordbv+1); j<=Math::min1(ordv-1,l); j++) {
							sum = sum + (*cpts)[s][i][j]*b.GetCPoints()[m-s][k-i][l-j]*(Math::Combin(ordu-1,i)*Math::Combin(ordbu-1,k-i)*Math::Combin(ordv-1,j)*Math::Combin(ordbv-1,l-j)*Math::Combin(ordw-1,s)*Math::Combin(ordbw-1,m-s)/(Math::Combin(ordu+ordbu-2,k)*Math::Combin(ordv+ordbv-2,l)*Math::Combin(ordw+ordbw-2,m)));
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
*/
// SUBDIVISION

// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BspVol<T> BezVol<T>::Subdivide(int levu, int levv, int levw) const
{
	return SubdivideU(levu).SubdivideV(levv).SubdivideW(levw);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPoints(int levu, int levv, int levw) const
{
	KnotSet kset1 = (*ksetu).CreateKnotSetSubdivide(levu);
	KnotSet kset2 = (*ksetv).CreateKnotSetSubdivide(levv);
	Matrix3D<T> mat = BspVol(SubdivideCPointsU(levu), kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu,ordv,ordw).SubdivideCPointsV(levv);
	return BspVol(mat, kset1.GetKnots(), kset2.GetKnots(), *ktsw, ordu,ordv,ordw,kset1.GetNum()-ordu,kset2.GetNum()-ordv,ordw).SubdivideCPointsW(levw);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPoints(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<T> mat1 = SubdivideCPointsW();
	KnotSet kset1 = (*ksetw).CreateKnotSetSubdivide(w1,w2);
	
	Matrix3D<T> mat2 = BezVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw).SubdivideCPointsV(v1,v2);
	KnotSet kset2 = (*ksetv).CreateKnotSetSubdivide(v1,v2);
		
	return BezVol<T>(mat2,*ktsu,kset2.GetKnots(),kset1.GetKnots(),ordu,ordv,ordw).SubdivideCPointsU(u1,u2);
}

// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BezVol<T> BezVol<T>::Subdivide(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	return (SubdivideU(u1,u2).SubdivideV(v1,v2)).SubdivideW(w1,w2);
}

// INSERT KNOTS


 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspVol<T> BezVol<T>::InsertKnotU(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// calculate the new knot set
	KnotSet kset = (*ksetu).CreateKnotSetInsert(t,mult,n);

	// create and return the new BspSurf
	return BspVol<T>(InsertKnotCPointsU(t,mult,n),kset.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,kset.GetNum()-ordu,ordv,ordw);
}
          
 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Matrix3D<T> BezVol<T>::InsertKnotCPointsU(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// create a local copy of the BspSurf
	BspVol temp = BezVol(*cpts,*ktsu,*ktsv,*ktsw,ordu,ordv,ordw).InsertKnotU(t[0],mult[0]);
    
	for (int i=1; i<n; i++) temp = temp.InsertKnotU(t[i],mult[i]);	
    return temp.GetCPoints();
}


// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspVol<T> BezVol<T>::InsertKnotV(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// calculate the new knot set
	KnotSet kset = (*ksetv).CreateKnotSetInsert(t,mult,n);

	// create and return the new BspSurf
	return BspVol<T>(InsertKnotCPointsV(t,mult,n),*ktsu,kset.GetKnots(),*ktsw,ordu,ordv,ordw,ordu,kset.GetNum()-ordv,ordw);
}
          
 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Matrix3D<T> BezVol<T>::InsertKnotCPointsV(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// create a local copy of the BspSurf
	BspVol temp = BezVol(cpts,*ktsu,*ktsv,*ktsw,ordu,ordv,ordw).InsertKnotV(t[0],mult[0]);
    
	for (int i=1; i<n; i++) temp = temp.InsertKnotV(t[i],mult[i]);	
    return temp.GetCPoints();
}



 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
BspVol<T> BezVol<T>::InsertKnotW(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// calculate the new knot set
	KnotSet kset = (*ksetw).CreateKnotSetInsert(t,mult,n);

	// create and return the new BspSurf
	return BspVol<T>(InsertKnotCPointsW(t,mult,n),*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw,ordu,ordv,kset.GetNum()-ordv);
}
          
 
// inserts a vector of knots into the BspCurv according to a vector
// of multiplicities giving the level of each knot to be inserted
template<class T>
Matrix3D<T> BezVol<T>::InsertKnotCPointsW(const Vector<double>& t,const Vector<int>& mult,int n) const
{
	// create a local copy of the BspSurf
	BspVol temp = BezVol(*cpts,*ktsu,*ktsv,*ktsw,ordu,ordv,ordw).InsertKnotW(t[0],mult[0]);
    
	for (int i=1; i<n; i++) temp = temp.InsertKnotW(t[i],mult[i]);	
    return temp.GetCPoints();
}

// CONVERSION

template<class T>
PolyVol<T> BezVol<T>::ConvertPolyVol() const
{
	Matrix3D<T> m =Math::mult6(Math::ComputePolyCurvMatrix(ordu),*cpts);
	Matrix3D<T> m1 = Math::mult5(m,Math::ComputePolyCurvMatrixTranspose(ordv));
		
	return PolyVol<T>(Math::mult9(Math::ComputePolyCurvMatrix(ordw),m1),ordu,ordv,ordw,(*ktsu)[ordu-1],(*ktsu)[ordu],(*ktsv)[ordv-1],(*ktsv)[ordv],(*ktsw)[ordw-1],(*ktsw)[ordw]);
}


template<class T>
BspVol<T> BezVol<T>::ConvertBspVol() const
{
	return BspVol<T>(*cpts, *ktsu, *ktsv, *ktsw, ordu, ordv, ordw, ordu, ordv, ordw);
}

// READ and WRITE

template <class T>
void BezVol<T>::write(std::ostream& os) 
{
	os << "Bezier Volume\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "order in w is " << ordw << "\n";
	os << "number of control points in u is " << ordu << "\n";
	os << "number of control points in v is " << ordv << "\n";
	os << "number of control points in w is " << ordw << "\n";
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
void BezVol<T>::read(std::istream& is)
{
	int Ordu, Ordv, Ordw;
	std::cout << "order of Bezier Volume in u and v and w";
	is >> Ordu >> Ordv >> Ordw;
	Matrix3D<T> Cpts(Ordu,Ordv,Ordw);
	
	std::cout << "\ninput control points\n";
	is >> Cpts;
	*this = BezVol<T>(Cpts,Ordu,Ordv,Ordw);
} 


template <class T>
void BezVol<T>::writefile(std::ofstream& ofs) 
{
	ofs << "Bezier Volume\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "order in w is " << ordw << "\n";
	ofs << "number of control points in u is " << ordu << "\n";
	ofs << "number of control points in v is " << ordv << "\n";
	ofs << "number of control points in w is " << ordw << "\n";
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are \n";
	ofs << *ktsv;
	
	ofs << "\nknots in w are \n";
	ofs << *ktsw;
	ofs << "\ncontrol points\n";
	ofs << *cpts;
}


template <class T>
void BezVol<T>::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw;
	
	ifs >> Ordu >> Ordv >> Ordw;
	
	Matrix3D<T> Cpts(Ordu,Ordv,Ordw);

	ifs >> Cpts;
	
	*this = BezVol<T>(Cpts,Ordu,Ordv,Ordw);
} 



// PRIVATE FUNCTIONS


// elevate the degree of the BezVol by level
template<class T>
Matrix3D<T> BezVol<T>::ElevateCPointsU(int level) const
{
	Matrix3D<T> ncpts(ordu+level,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetU(j,k),*ktsu,ordu).ElevateCPoints(level);
			// extract control points
			for (int i=0; i<ordu+level; i++) ncpts[k][i][j]=temp[i];
		}

	// return the new points
	return ncpts;
}


// elevate the degree of the BezVol by level
template<class T>
BezVol<T> BezVol<T>::ElevateU(int level) const
{
	// create the knot set
	KnotSet kset = (*ksetu).CreateKnotSetElevate(level);

	// return the new BezVol
	return BezVol<T>(ElevateCPointsU(level), kset.GetKnots(), *ktsv, *ktsw, ordu+level, ordv, ordw);	
}


// elevate the degree of the BezVol by level
template<class T>
Matrix3D<T> BezVol<T>::ElevateCPointsV(int level) const
{
	Matrix3D<T> ncpts(ordu,ordv+level,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetV(i,k),*ktsv,ordv).ElevateCPoints(level);
			// extract control points
			for (int j=0; j<ordv+level; j++) ncpts[k][i][j]=temp[j];
		}
	// return the new points
	return ncpts;
}



// elevate the degree of the BezVol by level
template<class T>
BezVol<T> BezVol<T>::ElevateV(int level) const
{
	// create the knot set
	KnotSet kset = (*ksetv).CreateKnotSetElevate(level);

	// return the new BezVol
	return BezVol<T>(ElevateCPointsV(level), *ktsu, kset.GetKnots(), *ktsw, ordu, ordv+level, ordw);
}


// elevate the degree of the BezVol by level
template<class T>
Matrix3D<T> BezVol<T>::ElevateCPointsW(int level) const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw+level);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetW(i,j),*ktsw,ordw).ElevateCPoints(level);
			// extract control points
			for (int k=0; k<ordw+level; k++) ncpts[k][i][j]=temp[k];
		}
	// return the new points
	return ncpts;
}


// elevate the degree of the BezVol by level
template<class T>
BezVol<T> BezVol<T>::ElevateW(int level) const
{
	// create the knot set
	KnotSet kset = (*ksetw).CreateKnotSetElevate(level);

	// return the new BezVol
	return BezVol<T>(ElevateCPointsW(level), *ktsu, *ktsv, kset.GetKnots(), ordu, ordv, ordw+level);
}


// elevate the degree of the BezVol by level
template<class T>
BezVol<T> BezVol<T>::ElevateUV(int levu, int levv) const
{
	return ElevateU(levu).ElevateV(levv);
}

// elevate the degree of the BezVol by level
template<class T>
BezVol<T> BezVol<T>::ElevateUW(int levu, int levw) const
{
	return ElevateU(levu).ElevateW(levw);
}

// elevate the degree of the BezVol by level
template<class T>
BezVol<T> BezVol<T>::ElevateVW(int levv, int levw) const
{
	return ElevateV(levv).ElevateW(levw);
}

// elevate the degree of the BezVol by level
template<class T>
BezVol<T> BezVol<T>::ElevateUVW(int levu, int levv, int levw) const
{
	return ElevateU(levu).ElevateV(levv).ElevateW(levw);
}


// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
Matrix3D<T> BezVol<T>::DeriveCPointsU(int level) const
{   
	Matrix3D<T> ncpts(ordu-level,ordv,ordw);

	// only need to compute knot vector once
	// returns control points only
	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetU(j,k),*ktsu,ordu).DeriveCPoints(level);
			// extract control points
			for (int i=0; i<ordu-level; i++) ncpts[k][i][j]=temp[i];
		}
	// return the new points
	return ncpts;
}

// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
BezVol<T> BezVol<T>::DeriveU(int level) const
{   
	// create knot set
	KnotSet kset = (*ksetu).CreateKnotSetDeriv(level);

	// return the new BezVol
	return BezVol<T>(DeriveCPointsU(level), kset.GetKnots(), *ktsv, *ktsw, ordu-level, ordv, ordw);
}



// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
Matrix3D<T> BezVol<T>::DeriveCPointsV(int level) const
{
	Matrix3D<T> ncpts(ordu,ordv+level,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetV(i,k),*ktsv,ordv).DeriveCPoints(level);
			// extract control points
			for (int j=0; j<ordv+level; j++) ncpts[k][i][j]=temp[j];
		}
	// return the new BezVol
	return ncpts;
}



// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
BezVol<T> BezVol<T>::DeriveV(int level) const
{
	// create the knot set
	KnotSet kset = (*ksetv).CreateKnotSetDeriv(level);

	// return the new BezVol
	return BezVol<T>(DeriveCPointsV(level), *ktsu, kset.GetKnots(), *ktsw, ordu, ordv+level, ordw);
}


// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
Matrix3D<T> BezVol<T>::DeriveCPointsW(int level) const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw+level);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetW(i,j),*ktsw,ordw).DeriveCPoints(level);
			// extract control points
			for (int k=0; k<ordw+level; k++) ncpts[k][i][j]=temp[k];
		}
	// return the new points
	return ncpts;
}



// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
BezVol<T> BezVol<T>::DeriveW(int level) const
{
	// create the knot set
	KnotSet kset = (*ksetw).CreateKnotSetDeriv(level);

	// return the new BezVol
	return BezVol<T>(DeriveCPointsW(level), *ktsu, *ktsv, kset.GetKnots(), ordu, ordv, ordw+level);
}

// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
BezVol<T> BezVol<T>::DeriveUV(int levu, int levv) const
{
	return DeriveU(levu).DeriveV(levv);
}

// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
BezVol<T> BezVol<T>::DeriveUW(int levu, int levw) const
{
	return DeriveU(levu).DeriveW(levw);
}


// compute the derivative of the BezVol of order deriv and
// represent the result as another BezVol
template<class T>
BezVol<T> BezVol<T>::DeriveVW(int levv, int levw) const
{
	return DeriveV(levv).DeriveW(levw);
}


// evaluate the derivative of the BezVol of order deriv at a point
// x. Computes the derivative as a BezVol and then evaluates this at x
template<class T>
T BezVol<T>::DeriveU(int level, double u, double v, double w) const
{
	return DeriveU(level)(u,v,w);
}


// evaluate the derivative of the BezVol of order deriv at a point
// x. Computes the derivative as a BezVol and then evaluates this at x
template<class T>
T BezVol<T>::DeriveV(int level, double u, double v, double w) const
{
	return DeriveV(level)(u,v,w);
}


// evaluate the derivative of the BezVol of order deriv at a point
// x. Computes the derivative as a BezVol and then evaluates this at x
template<class T>
T BezVol<T>::DeriveW(int level, double u, double v, double w) const
{
	return DeriveW(level)(u,v,w);
}


// evaluate the derivative of the BezVol of order deriv at a point
// x. Computes the derivative as a BezVol and then evaluates this at x
template<class T>
T BezVol<T>::DeriveUV(int levu, int levv, double u, double v, double w) const
{
	return (DeriveU(levu).DeriveV(levv))(u,v,w);
}


// evaluate the derivative of the BezVol of order deriv at a point
// x. Computes the derivative as a BezVol and then evaluates this at x
template<class T>
T BezVol<T>::DeriveUW(int levu, int levw, double u, double v, double w) const
{
	return (DeriveU(levu).DeriveW(levw))(u,v,w);
}


// evaluate the derivative of the BezVol of order deriv at a point
// x. Computes the derivative as a BezVol and then evaluates this at x
template<class T>
T BezVol<T>::DeriveVW(int levv, int levw, double u, double v, double w) const
{
	return (DeriveV(levv).DeriveW(levw))(u,v,w);
}




// integrate the BezVol between the limits x1 and x2. Computes
// the indefinite integral as a BezVol and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T BezVol<T>::IntegrateUVW1(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Vector<T> w(ordw);

	for (int k=0; k<ordw; k++)
		w[k] = BezSurf<T>((*cpts).GetUV(k),*ktsu,*ktsv,ordu,ordv).Integrate(u1,u2,v1,v2);
	return BezCurv<T>(w,*ktsw,ordw).Integrate(w1,w2);
}

// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
Matrix3D<T> BezVol<T>::IntegrateCPointsU() const
{
	Matrix3D<T> ncpts(ordu+1,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = BezCurv<T>((*cpts).GetU(j,k),*ktsu,ordu).IntegrateCPoints();
			// extract control points
			for (int i=0; i<ordu+1; i++) ncpts[k][i][j]=temp[i];
		}	
	// return new points
	return ncpts;
}	



// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
BezVol<T> BezVol<T>::IntegrateU() const
{
	// create knot sets
	KnotSet ksetu = (*ksetu).CreateKnotSetIntegrate();

	// return new BezVol
	return BezVol<T>(IntegrateCPointsU(),ksetu.GetKnots(),*ktsv,*ktsw,ordu+1,ordv,ordw);
}	



// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
Matrix3D<T> BezVol<T>::IntegrateCPointsV() const
{
	Matrix<T> ncpts(ordu,ordv+1,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> v = BezCurv<T>((*cpts).GetV(i,k),*ktsv,ordv).IntegrateCPoints();
			// extract control points
			for (int j=0; j<ordv+1; j++) ncpts[k][i][j]=temp[j];
		}	
	return ncpts;
}	


// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
BezVol<T> BezVol<T>::IntegrateV() const
{
	// create knot sets
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();

	return BezVol<T>(IntegrateCPointsV(),*ktsu,kset.GetKnots(),*ktsw,ordu,ordv+1,ordw);	
}	


// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
Matrix3D<T> BezVol<T>::IntegrateCPointsW() const
{
	Matrix<T> ncpts(ordu,ordv,ordw+1);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> v = BezCurv<T>((*cpts).GetW(i,j),*ktsw,ordw).IntegrateCPoints();
			// extract control points
			for (int k=0; k<ordw+1; k++) ncpts[k][i][j]=temp[k];
		}	
	return ncpts;
}	


// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
BezVol<T> BezVol<T>::IntegrateW() const
{
	// create knot sets
	KnotSet kset = (*ksetw).CreateKnotSetIntegrate();

	// return new BezVol
	return BezVol<T>(IntegrateCPointsW(),*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw+1);	
}	


// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
BezVol<T> BezVol<T>::IntegrateUV() const
{
	return IntegrateU().IntegrateV();
}	

// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
BezVol<T> BezVol<T>::IntegrateUW() const
{
	return IntegrateU().IntegrateW();
}	


// compute the indefinite integral of the BezVol and represent
// it as a BezVol of one higher degree
template<class T>
BezVol<T> BezVol<T>::IntegrateVW() const
{
	return IntegrateV().IntegrateW();
}	


	
// compute the indefinite integral of the BezVol as a BezVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> BezVol<T>::IntegrateCPointsUV() const
{
	Matrix3D<T> mat = IntegrateCPointsV();
	KnotSet kset = (*ksetv).CreateKnotSetIntegrate();

	return BezVol<T>(mat,*ktsu,kset.GetKnots(),*ktsw,ordu,ordv+1,ordw).IntegrateCPointsU();
}	




// compute the indefinite integral of the BezVol as a BezVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> BezVol<T>::IntegrateCPointsUW() const
{
	Matrix3D<T> mat = IntegrateCPointsW();
	KnotSet kset = (*ksetw).CreateKnotSetIntegrate();

	return BezVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw+1).IntegrateCPointsU();
}	


// compute the indefinite integral of the BezVol as a BezVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> BezVol<T>::IntegrateCPointsVW() const
{
	Matrix3D<T> mat = IntegrateCPointsW();
	KnotSet kset = (*ksetw).CreateKnotSetIntegrate();

	return BezVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw+1).IntegrateCPointsV();
}	

// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BspVol<T> BezVol<T>::SubdivideU(int level) const
{
	// create knot set
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(level);
	// return the new BezVol
	return BspVol<T>(SubdivideCPointsU(level), kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu,ordv,ordw);
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsU(int level) const
{
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(level);
	Matrix3D<T> ncpts(kset.GetNum()-ordu,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetU(j,k),*ktsu,ordu).SubdivideCPoints(level);
			// extract control points
			for (int i=0; i<kset.GetNum()-ordu; i++) ncpts[k][i][j]=temp[i];
		}
	// return the new points
	return ncpts;
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BspVol<T> BezVol<T>::SubdivideV(int level) const
{
	// create knot set
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(level);
	// return the new BezVol
	return BspVol<T>(SubdivideCPointsV(level), *ktsu, kset.GetKnots(), *ktsw, ordu, ordv, ordw, ordu, kset.GetNum()-ordv,ordw);
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsV(int level) const
{
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(level);
	Matrix3D<T> ncpts(ordu,kset.GetNum()-ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetV(i,k),*ktsv,ordv).SubdivideCPoints(level);
			// extract control points
			for (int j=0; j<kset.GetNum()-ordv; j++) ncpts[k][i][j]=temp[j];
		}
	// return the new points
	return ncpts;
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BspVol<T> BezVol<T>::SubdivideW(int level) const
{
	// create knot set
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(level);
	// return the new BezVol
	return BspVol<T>(SubdivideCPointsV(level), *ktsu, *ktsv, kset.GetKnots(), ordu, ordv, ordw, ordu, ordv, kset.GetNum()-ordw);
}

// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BspVol<T> BezVol<T>::SubdivideUV(int levu, int levv) const
{
	return SubdivideU(levu).SubdivideV(levv);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BspVol<T> BezVol<T>::SubdivideUW(int levu, int levw) const
{
	return SubdivideU(levu).SubdivideW(levw);
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BspVol<T> BezVol<T>::SubdivideVW(int levv, int levw) const
{
	return SubdivideV(levv).SubdivideW(levw);
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsUV(int levu, int levv) const
{
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(levu);
	return BspVol(SubdivideCPointsU(levu), kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu,ordv,ordw).SubdivideCPointsV(levv);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsVW(int levv, int levw) const
{
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(levv);
	return BspVol(SubdivideCPointsV(levv), *ktsu, kset.GetKnots(), *ktsw, ordu, ordv, ordw, ordu, kset.GetNum()-ordv,ordw).SubdivideCPointsW(levw);
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsUW(int levu, int levw) const
{
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(levu);
	return BspVol(SubdivideCPointsU(levu), kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw, kset.GetNum()-ordu,ordv,ordw).SubdivideCPointsW(levw);
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsW(int level) const
{
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(level);
	Matrix3D<T> ncpts(ordu,ordv,kset.GetNum()-ordw);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetW(i,j),*ktsw,ordw).SubdivideCPoints(level);
			// extract control points
			for (int k=0; k<kset.GetNum()-ordw; k++) ncpts[k][i][j]=temp[k];
		}
	// return the new points
	return ncpts;
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsU(double u1, double u2) const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetU(j,k),*ktsu,ordu).SubdivideCPoints(u1,u2);
			// extract control points
			for (int i=0; i<ordu; i++) ncpts[k][i][j]=temp[i];
		}
	// return the new points
	return ncpts;
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BezVol<T> BezVol<T>::SubdivideU(double u1, double u2) const
{
	// create knot set
	KnotSet kset = (*ksetu).CreateKnotSetSubdivide(u1,u2);
	// return the new BezVol
	return BezVol<T>(SubdivideCPointsU(), kset.GetKnots(), *ktsv, *ktsw, ordu, ordv, ordw);
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVolace after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsV(double v1, double v2) const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetV(i,k),*ktsv,ordv).SubdivideCPoints(v1,v2);
			// extract control points
			for (int j=0; j<ordv; j++) ncpts[k][i][j]=temp[j];
		}
	// return the new points
	return ncpts;
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVolace after the knot refinement
template<class T>
BezVol<T> BezVol<T>::SubdivideV(double v1, double v2) const
{
	// create the knot set
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);
	// return the new BezVol
	return BezVol<T>(SubdivideCPointsV(), *ktsu, kset.GetKnots(), *ktsw, ordu, ordv, ordw);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVolace after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsW(double w1, double w2) const
{
	Matrix3D<T> ncpts(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetW(i,j),*ktsw,ordw).SubdivideCPoints(w1,w2);
			// extract control points
			for (int k=0; k<ordw; k++) ncpts[k][i][j]=temp[k];
		}
	// return the new points
	return ncpts;
}




// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVolace after the knot refinement
template<class T>
BezVol<T> BezVol<T>::SubdivideW(double w1, double w2) const
{
	// create the knot set
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(w1,w2);
	// return the new BezVol
	return BezVol<T>(SubdivideCPointsW(), *ktsu, *ktsv, temp.GetKnots(), ordu, ordv, ordw);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsUV(double u1, double u2, double v1, double v2) const
{
	Matrix3D<T> mat = SubdivideCPointsV();
	KnotSet kset = (*ksetv).CreateKnotSetSubdivide(v1,v2);

	return BezVol<T>(mat,*ktsu,kset.GetKnots(),*ktsw,ordu,ordv,ordw).SubdivideCPointsU(u1,u2);
}



// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsUW(double u1, double u2, double w1, double w2) const
{
	Matrix3D<T> mat = SubdivideCPointsW();
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(w1,w2);

	return BezVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw).SubdivideCPointsU(u1,u2);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
Matrix3D<T> BezVol<T>::SubdivideCPointsVW(double v1, double v2, double w1, double w2) const
{
	Matrix3D<T> mat = SubdivideCPointsW();
	KnotSet kset = (*ksetw).CreateKnotSetSubdivide(w1,w2);

	return BezVol<T>(mat,*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw).SubdivideCPointsV(v1,v2);
}

// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BezVol<T> BezVol<T>::SubdivideUV(double u1, double u2, double v1, double v2) const
{
	return SubdivideU(u1,u2).SubdivideV(v1,v2);
}


// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BezVol<T> BezVol<T>::SubdivideUW(double u1, double u2, double w1, double w2) const
{
	return SubdivideU(u1,u2).SubdivideW(w1,w2);
}

// subdivide the BVolace upto the order given by level. Method
// inserts level knots between each distinct knot in the original
// knot set. Returns the BVol after the knot refinement
template<class T>
BezVol<T> BezVol<T>::SubdivideVW(double v1, double v2, double w1, double w2) const
{
	return SubdivideV(v1,v2).SubdivideW(w1,w2);
}


// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
Matrix3D<T> BezVol<T>::InsertKnotCPointsU(double x) const
{	
	Matrix3D<T> ncpts(ordu+1,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>(cpts.GetU(j,k),ktsu,ordu).InsertKnotCPoints(x);
			// extract control points
			for (int i=0; i<ordu+1; i++) ncpts[k][i][j]=temp[i];
		}
	// return the new points
	return ncpts;
}
     
// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
BspVol<T> BezVol<T>::InsertKnotU(double x) const
{       
	// calculate the new knot set
	KnotSet kset = (*ksetu).CreateKnotSetInsert(x);

	// create and return the new BspSurf
	return BspVol<T>(InsertKnotCPointsU(x),kset.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,ordu+1,ordv,ordw);
}              


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
Matrix3D<T> BezVol<T>::InsertKnotCPointsU(double x, int level) const
{	
	Matrix3D<T> ncpts(ordu+level,ordv,ordw);

	for (int j=0; j<ordv; j++) 
		for (int k=0; k<ordw; k++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetU(j,k),*ktsu,ordu).InsertKnotCPoints(x,level);
			// extract control points
			for (int i=0; i<ordu+level; i++) ncpts[k][i][j]=temp[i];
		}
	// return the new points
	return ncpts;
}
	


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
BspVol<T> BezVol<T>::InsertKnotU(double x, int level) const
{	
	// calculate the new knot set
	KnotSet kset = (*ksetu).CreateKnotSetInsert(x,level);

	// create and return the new BspSurf
	return BspVol<T>(InsertKnotCPointsU(x,level),kset.GetKnots(),*ktsv,*ktsw,ordu,ordv,ordw,ordu+level,ordv,ordw);
}


// inserts a knot x into the BspCurv  and returns the new BspCurv
template<class T>
BspVol<T> BezVol<T>::InsertKnotW(double x) const
{       
	// calculate the new knot set
	KnotSet kset = (*ksetw).CreateKnotSetInsert(x);

	// create and return the new BspSurf
	return BspVol<T>(InsertKnotCPointsW(x),*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw,ordu,ordv,ordv+1);
}              


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
Matrix3D<T> BezVol<T>::InsertKnotCPointsW(double x, int level) const
{	
	Matrix3D<T> ncpts(ordu,ordv,ordw+level);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			Vector<T> temp = BezCurv<T>((*cpts).GetW(i,j),*ktsw,ordw).InsertKnotCPoints(x,level);
			// extract control points
			for (int k=0; k<ordw+level; k++) ncpts[k][i][j]=temp[k];
		}
	// return the new points
	return ncpts;
}


// inserts a knot x 'level' times into the BspCurv and returns
// the new BspCurv
template<class T>
BspVol<T> BezVol<T>::InsertKnotW(double x, int level) const
{	
	// calculate the new knot set
	KnotSet kset = (*ksetw).CreateKnotSetInsert(x,level);

	// create and return the new BspSurf
	return BspVol<T>(InsertKnotCPointsW(x,level),*ktsu,*ktsv,kset.GetKnots(),ordu,ordv,ordw,ordu,ordv,ordw+level);
}



#endif

