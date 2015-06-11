#ifndef COMPOLYSURF
#define COMPOLYSURF


#include "mathfunctions.h"
#include "PolySurf.h"
#include "CompBezCurv.h"
#include "CompBezSurf.h"
#include "BspSurf.h"

// can be of different orders, IGES spec

template<class T>
class CompPolySurf : public Surf<T> {
private:
	// data
	int numu; // number of segments in u
	int numv; // number of segments in v
	int ordu; // max ord in u
	int ordv; // max ord in v
	Ptr<Vector<int> > ordsu; // orders in u
	Ptr<Vector<int> > ordsv;	// orders in v
	Ptr<Matrix<T> > coeffs;	// coeffs
	Ptr<Vector<double> > limitsu; // ranges in u
	Ptr<Vector<double> > limitsv;	// ranges in v
	Ptr<std::set<double> > limitsetu;
	Ptr<std::set<double> > limitsetv;


	// private functions
	CompPolySurf<T> DeriveU(int level) const; 
	CompPolySurf<T> DeriveV(int level) const; 
	CompPolySurf<T> ElevateU(int levu) const;
	CompPolySurf<T> ElevateV(int levv) const;
	CompPolySurf<T> IntegrateU() const;
	CompPolySurf<T> IntegrateV() const;
	Matrix<T> IntegrateCPointsU() const;
	Matrix<T> IntegrateCPointsV() const;
	CompPolySurf<T> Product2(const CompPolySurf<T>& c) const;
	template<class T1>
	CompPolySurf<T> MakeCompatable(const CompPolySurf<T1>& b) const
	{
		return (ConvertBspSurf().MakeBreakCompatable(b.ConvertBspSurf())).ConvertCompPolySurf();
	}
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class CompPolySurf<double>");
		else {
			std::string s(typeid(T).name()), s1("class CompPolySurf<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	CompPolySurf();
	CompPolySurf(int Numu, int Numv, int Ordu, int Ordv, const Matrix<T>& Coeffs, const Vector<double>& Limitsu, const Vector<double>& Limitsv);
	CompPolySurf(int Numu, int Numv, int Ordu, int Ordv, const Matrix<T>& Coeffs);
	CompPolySurf(int Numu, int Numv, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Matrix<T>& Coeffs, const Vector<double>& Limitsu, const Vector<double>& Limitsv);
	CompPolySurf(int Numu, int Numv, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Matrix<T>& Coeffs);

	// access functions
	int GetMaxOrdU() const;
	int GetMaxOrdV() const;
	Vector<double> GetLimitsU() const;
	Vector<double> GetLimitsV() const;
	double GetLeftLimitU(int indexu, int indexv) const;
	double GetRightLimitU(int indexu, int indexv) const;
	double GetLeftLimitV(int indexu, int indexv) const;
	double GetRightLimitV(int indexu, int indexv) const;
	double GetLeftLimitU() const;
	double GetRightLimitU() const;
	double GetLeftLimitV() const;
	double GetRightLimitV() const;
	Matrix<T> GetCoeffsPoly(int indexu, int indexv) const;
	int GetNumU() const;
	int GetNumV() const;
	int FindSegmentU(double x) const;
	int FindSegmentV(double x) const;
	Matrix<T> GetCoeffs() const;
	Vector<int> GetOrdsU() const;
	Vector<int> GetOrdsV() const;
	int GetOrdU() const;
	int GetOrdV() const;
	PolySurf<T> GetPatch(int indexu, int indexv) const;
	PolySurf<T> GetPatchOriginal(int indexu, int indexv) const;
	
	
	// evaluators
	virtual T operator()(double u, double v) const;
	T operator()(double u, double v, int indexu, int indexv) const;
	virtual T operator()(int, int, double u, double v) const;
	T Eval(double u, double v) const;
	T Eval(double u, double v, int indexu, int indexv) const;
	Matrix<T> ComputePoints(int n, int m) const;

	// conversion
	Matrix<T> Convert() const;
	CompBezSurf<T> ConvertCompBezSurf() const;
	BspSurf<T> ConvertBspSurf() const;
	
	// derivatives
	CompPolySurf<T> Derive(int levu, int levv) const;
	virtual T Derive(int levu, int levv, double u, double v) const;

	// degree elevation
	CompPolySurf<T> Elevate(int levu, int levv) const;

	// integration
	CompPolySurf<T> Integrate() const;
	Matrix<T> IntegrateCPoints() const;
	T Integrate(double u1, double u2, double v1, double v2) const;
	
	// product
	template<class T1>
	CompPolySurf<T> Product(const CompPolySurf<T1>& b) const
	{
		CompPolySurf<T> d(*this), e(b);

		if (numu != p.GetNumU() || numv != p.GetNumV()) {
			d = MakeCompatable(b);
			e = b.MakeCompatable(d);
		}
		// find number of new control points
		int sum1=d.GetNumU()*(ordu+p.GetOrdU()-1);
	
		int sum2=d.GetNumV()*(ordv+p.GetOrdV()-1);
	
		Matrix<T> ncoeffs(sum1,sum2);

		int count1;
		int count2;
		int rowstart;
		int colstart;
		// find product patch by patch
		for (int i=0; i<d.GetNumU(); i++) {
			rowstart = i*(ordu+p.GetOrdU()-1);
			count1=rowstart;
			for (int j=0; j<d.GetNumV(); j++) {
				PolySurf<T> p1 = d.GetPatch(i+1,j+1);
				PolySurf<T> p2 = e.GetPatch(i+1,j+1);
				Matrix<T> temp = p1.ProductCPoints(p2);
				colstart = j*(ordv+p.GetOrdV()-1);
				count2 = colstart;
				for (int k=0; k<ordu+p.GetOrdU()-1; k++) {
					for (int l=0; l<ordv+p.GetOrdV()-1; l++) {
						ncoeffs[count1][count2] = temp[k][l];
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
		Vector<int> Ordsu(d.GetNumU());
		Vector<int> Ordsv(d.GetNumV());

		for (int i=0; i<d.GetNumU(); i++) Ordsu[i]=(*ordsu)[i]+p.GetOrdsU()[i]-1;
		for (int i=0; i<d.GetNumV(); i++) Ordsv[i]=(*ordsv)[i]+p.GetOrdsV()[i]-1;
		// form new Ordsu and Ordsv
		// create and return the surface
		return CompPolySurf<T>(d.GetNumU(),d.GetNumV(),Ordsu,Ordsv,ncoeffs,*limitsu,*limitsv);
	}
	
	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

// CONSTRUCTORS

// default constructor
template<class T>
CompPolySurf<T>::CompPolySurf() : numu(0), numv(0), ordu(0), ordv(0), ordsu(), ordsv(), coeffs(), limitsu(),limitsv()
{
}

// constructor builds a CompPolySurf from a Vector of control points,
// an order and number of control points and limits in u and v
template<class T>
CompPolySurf<T>::CompPolySurf(int Numu, int Numv, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Matrix<T>& Coeffs, 
							  const Vector<double>& Limitsu, const Vector<double>& Limitsv) :
	numu(Numu), numv(Numv), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), coeffs(new Matrix<T>(Coeffs)), limitsu(new Vector<double>(Limitsu)), limitsv(new Vector<double>(Limitsv)),
	limitsetu(new std::set<double>((*limitsu).begin(),(*limitsu).end())),limitsetv(new std::set<double>((*limitsv).begin(),(*limitsv).end()))
{
	// find max order in u
	ordu = (*ordsu)[0];
	for (int i=1; i<numu; i++) if (ordu < (*ordsu)[i]) ordu=(*ordsu)[i];

	// find max order in v
	ordv = (*ordsv)[0];
	for (int j=1; j<numv; j++) if (ordv < (*ordsv)[j]) ordv=(*ordsv)[j];
	coeffs = new Matrix<T>(Convert());
}

// constructor builds a CompPolySurf from a Vector of control points,
// an order and number of control points in u and v
template<class T>
CompPolySurf<T>::CompPolySurf(int Numu, int Numv, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Matrix<T>& Coeffs) :
numu(Numu), numv(Numv), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), coeffs(new Matrix<T>(Coeffs)), limitsu(new Vector<double>(numu+1)), limitsv(new Vector<double>(numv+1))
{
	// find max order in u
	ordu = (*ordsu)[0];
	for (int i=1; i<numu; i++) if (ordu < (*ordsu)[i]) ordu=(*ordsu)[i];

	// find max order in v
	ordv = (*ordsv)[0];
	for (int j=1; j<numv; j++) if (ordv < (*ordsv)[j]) ordv=(*ordsv)[j];
	
	// limits are increasing integers in u and v
	for (int i=0; i<=numu; i++) (*limitsu)[i]=(double)i;
	for (int j=0; j<=numv; j++) (*limitsv)[j]=(double)j;
	coeffs = new Matrix<T>(Convert());
	limitsetu = new std::set<double>((*limitsu).begin(),(*limitsu).end());
	limitsetv = new std::set<double>((*limitsv).begin(),(*limitsv).end());

}


// constructor builds a CompPolySurf from a Vector of control points,
// an order and number of control points and limits in u and v 
template<class T>
CompPolySurf<T>::CompPolySurf(int Numu, int Numv, int Ordu, int Ordv, const Matrix<T>& Coeffs, const Vector<double>& Limitsu, const Vector<double>& Limitsv) :
numu(Numu), numv(Numv), ordu(Ordu), ordv(Ordv), ordsu(new Vector<int>(numu)), ordsv(new Vector<int>(numv)), coeffs(new Matrix<T>(Coeffs)), limitsu(new Vector<double>(Limitsu)), limitsv(new Vector<double>(Limitsv)),
limitsetu(new std::set<double>((*limitsu).begin(),(*limitsu).end())),limitsetv(new std::set<double>((*limitsv).begin(),(*limitsv).end()))
{
	// assign orders in u
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	// assign orders in v
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;
}


// constructor builds a CompPolySurf from a Vector of control points,
// an order and number of control points in u and v 
template<class T>
CompPolySurf<T>::CompPolySurf(int Numu, int Numv, int Ordu, int Ordv, const Matrix<T>& Coeffs) :
numu(Numu), numv(Numv), ordu(Ordu), ordv(Ordv), ordsu(new Vector<int>(numu)), 
ordsv(new Vector<int>(numv)), coeffs(new Matrix<T>(Coeffs)), 
limitsu(new Vector<double>(numu+1)), limitsv(new Vector<double>(numv+1))
{
	// assign orders in u
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	// assign orders in v
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;

	// find limits in u
	for (int i=0; i<=numu; i++) (*limitsu)[i]=(double)i;
	// find limits in v
	for (int j=0; j<=numv; j++) (*limitsv)[j]=(double)j;
	limitsetu = new std::set<double>((*limitsu).begin(),(*limitsu).end());
	limitsetv = new std::set<double>((*limitsv).begin(),(*limitsv).end());

}


// ACCESS FUNCTIONS


// get the order of the CompPolySurf
template<class T>
inline int CompPolySurf<T>::GetOrdU() const { return ordu; }

// get the order of the CompPolySurf
template<class T>
inline int CompPolySurf<T>::GetOrdV() const { return ordv; }


// get the orders in u
template<class T>
inline Vector<int> CompPolySurf<T>::GetOrdsU() const { return *ordsu; }


// get the orders in v
template<class T>
inline Vector<int> CompPolySurf<T>::GetOrdsV() const { return *ordsv; }


// get the number of segments in u
template<class T>
inline int CompPolySurf<T>::GetNumU() const { return numu; }

// get the number of segments in v
template<class T>
inline int CompPolySurf<T>::GetNumV() const { return numv; }


// get the coefficients
template<class T>
inline Matrix<T> CompPolySurf<T>::GetCoeffs() const { return *coeffs; }


// get the limit set in u
template<class T>
inline Vector<double> CompPolySurf<T>::GetLimitsU() const { return *limitsu; }

// get the limit set in v
template<class T>
inline Vector<double> CompPolySurf<T>::GetLimitsV() const { return *limitsv; }


template<class T>
double CompPolySurf<T>::GetLeftLimitU() const 
{ 
	return (*limitsu)[0];
}


template<class T>
double CompPolySurf<T>::GetRightLimitU() const 
{ 
	return (*limitsu)[numu];
}

template<class T>
double CompPolySurf<T>::GetLeftLimitV() const 
{ 
	return (*limitsv)[0];
}


template<class T>
double CompPolySurf<T>::GetRightLimitV() const 
{ 
	return (*limitsv)[numv];
}



template<class T>
int CompPolySurf<T>::FindSegmentU(double x) const
{
	std::set<double>::iterator s1 = (*limitsetu).begin();
	std::set<double>::iterator s2 = (*limitsetu).upper_bound(x);
	if (s2 == (*limitsetu).end()) return (*limitsetu).size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);
	return count;		
}


template<class T>
int CompPolySurf<T>::FindSegmentV(double x) const
{
	std::set<double>::iterator s1 = (*limitsetv).begin();
	std::set<double>::iterator s2 = (*limitsetv).upper_bound(x);
	if (s2 == (*limitsetv).end()) return (*limitsetu).size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);
	return count;
}


// get the i,j patch of the composite surface
template<class T>
PolySurf<T> CompPolySurf<T>::GetPatch(int indi, int indj) const
{
	// find start index for control points in u and v
	int sum1=indi*ordu;
	
	int sum2=indj*ordv;
	
	// extract control points
	Matrix<T> ncoeffs(ordu,ordv);
	for (int k=0; k<ordu; k++) 
		for (int l=0; l<ordv; l++)
			ncoeffs[k][l] = (*coeffs)[sum1-ordu+k][sum2-ordv+l];

	// create and return the PolySurf
	return PolySurf<T>(ncoeffs, ordu, ordv, (*limitsu)[indi-1], (*limitsu)[indi], (*limitsv)[indj-1], (*limitsv)[indj]);
}



// get the i,j patch of the composite surface
template<class T>
PolySurf<T> CompPolySurf<T>::GetPatchOriginal(int indi, int indj) const
{
	// order in u and v
	int ordU = (*ordsu)[indi-1];
	int ordV = (*ordsv)[indj-1];

	// find start index for control points in u and v
	int sum1=0;
	for (int i=0; i<indi; i++) sum1 = sum1+(*ordsu)[i];

	int sum2=0;
	for (int j=0; j<indj; j++) sum2 = sum2+(*ordsv)[j];
	
	// extract control points
	Matrix<T> ncoeffs(ordu,ordv);
	for (int k=0; k<ordu; k++) 
		for (int l=0; l<ordv; l++)
			ncoeffs[k][l] = (*coeffs)[sum1-ordu+k][sum2-ordv+l];

	// create and return the PolySurf
	return PolySurf<T>(ncoeffs, ordU, ordV, (*limitsu)[indi-1], (*limitsu)[indi], (*limitsv)[indj-1], (*limitsv)[indj]);
}



// CONVERSIONS

// convert to BspSurf
template<class T>
BspSurf<T> CompPolySurf<T>::ConvertBspSurf() const
{
	return ConvertCompBezSurf().ConvertBspSurf();
}

// convert to composite Bezier form
template<class T>
CompBezSurf<T> CompPolySurf<T>::ConvertCompBezSurf() const
{
	// take each segment and convert to BezSurf form

	// find the sum of the ords in u
	int sum1=numu*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	Matrix<T> ncpts(sum1,sum2);
	int count1;
	int count2;
	int rowstart;
	int colstart;
	// convert each patch
	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> b = GetPatch(i+1,j+1).ConvertBezSurfCPoints();
			colstart = j*ordv;
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv; l++) {
					ncpts[count1][count2]=b[k][l];
					count2++;
				}
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	// create and return the composite Bezier surface
	return CompBezSurf<T>(numu, numv, *ordsu, *ordsv, ncpts, *limitsu, *limitsv);
}


// convert each patch to max degree in u and v
template<class T>
Matrix<T> CompPolySurf<T>::Convert() const
{
	// elevate the degree of each segment to max
	// build array of control points
	int count1=0;
	int count2=0;
	int rowstart;
	int colstart;

	Matrix<T> ncoeffs(numu*ordu,numv*ordv);
	// elevate each patch as necessary
	for (int i=0; i<numu; i++) {
		rowstart=i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {

			
			Matrix<T> mat = GetPatchOriginal(i+1,j+1).Elevate(ordu-(*ordsu)[i],ordv-(*ordsv)[j]).GetCoeffs();
			colstart=j*ordv;
			count2=colstart;
			for (int k=0; k<ordu; k++) { 
				for (int l=0; l<ordv; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	// create and return the CompPolySurf
	return ncoeffs;
}


// EVALUATORS

// evaluate the CompPolySurf at the point x 
template<class T>
T CompPolySurf<T>::operator()(double u, double v) const
{
	// Find segment in u and v
	int indexu = FindSegmentU(u);
	int indexv = FindSegmentV(v);

	// extract and evaluate segment
	return GetPatch(indexu,indexv)(u,v);
}

// evaluate the segment index at the point x 
template<class T>
T CompPolySurf<T>::operator()(double u, double v, int indexu, int indexv) const
{
	// extract and evaluate segment
	return GetPatch(indexu,indexv)(u,v);
}


// evaluate the CompPolySurf using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompPolySurf<T>::Eval(double u, double v) const
{
	// Find segment
	int indu = FindSegment(u);
	int indv = FindSegment(v);

	// evaluate segemnt
	return GetPatch(indu,indv).Eval(u,v);
}

// evaluate the CompPolySurf using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompPolySurf<T>::Eval(double u, double v, int indexu, int indexv) const
{
	// extract and evaluate segment
	return GetPatch(indexu,indexv).Eval(u,v);
}

// DEGREE ELEVATION


// elevate the degree of the CompPolySurf by level
template<class T>
CompPolySurf<T> CompPolySurf<T>::Elevate(int levu, int levv) const
{
	// elevate in u then v
	return ElevateU(levu).ElevateV(levv);
}


// DERIVATIVES


// evaluate the derivative of the CompPolySurf of order deriv at a point
// x. Computes the derivative as a CompPolySurf and then evaluates this at x
template<class T>
CompPolySurf<T> CompPolySurf<T>::Derive(int levu, int levv) const
{
	// derive in u then v
	return DeriveU(levu).DeriveV(levv);
}

template<class T>
T CompPolySurf<T>::operator() (int valu, int valv, double u, double v) const
{
	return Derive(valu,valv,u,v);
}



// evaluate the derivative of the CompPolySurf of order deriv at a point
// x. Computes the derivative as a CompPolySurf and then evaluates this at x
template<class T>
T CompPolySurf<T>::Derive(int levu, int levv, double u, double v) const
{
	return Derive(levu,levv)(u,v);
}


// INTEGRATION

// TO DO
// integrate the CompPolySurf between the limits x1 and x2. Computes
// the indefinite integral as a CompPolySurf and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T CompPolySurf<T>::Integrate(double u1, double u2, double v1, double v2) const
{
	// create the indefinite integral;
	CompPolySurf<T> p = Integrate(); 

	// evaluate and subtract
	return (p(u2,v2)-p(u1,v2)-p(u2,v1)+p(u1,v1));
}



// evaluate the derivative of the CompPolySurf of order deriv at a point
// x. Computes the derivative as a CompPolySurf and then evaluates this at x
template<class T>
CompPolySurf<T> CompPolySurf<T>::Integrate() const
{
	// integrate in u then v
	return IntegrateU().IntegrateV();
}


/*
template<class T>
CompPolySurf<T> CompPolySurf<T>::MakeCompatable(const CompPolySurf<T>& b) const
{
	return (ConvertBspSurf().MakeBreakCompatable(b.ConvertBspSurf())).ConvertCompPolySurf();
}		
*/
/*
// PRODUCT

// compute the product of the CompPolySurf with another CompPolySurf and 
// represent the result as a new CompPolySurf 
// must be compatable in the sense of having the same number of patches
// in u and v
template<class T>  
CompPolySurf<T> CompPolySurf<T>::Product(const CompPolySurf<T>& p) const
{
	CompPolySurf<T> d(*this), e(b);

	if (numu != p.GetNumU() || numv != p.GetNumV()) {
		d = MakeCompatable(b);
		e = b.MakeCompatable(d);
	}
	// find number of new control points
	int sum1=d.GetNumU()*(ordu+p.GetOrdU()-1);
	
	int sum2=d.GetNumV()*(ordv+p.GetOrdV()-1);
	
	Matrix<T> ncoeffs(sum1,sum2);

	int count1;
	int count2;
	int rowstart;
	int colstart;
	// find product patch by patch
	for (int i=0; i<d.GetNumU(); i++) {
		rowstart = i*(ordu+p.GetOrdU()-1);
		count1=rowstart;
		for (int j=0; j<d.GetNumV(); j++) {
			PolySurf<T> p1 = d.GetPatch(i+1,j+1);
			PolySurf<T> p2 = e.GetPatch(i+1,j+1);
			Matrix<T> temp = p1.ProductCPoints(p2);
			colstart = j*(ordv+p.GetOrdV()-1);
			count2 = colstart;
			for (int k=0; k<ordu+p.GetOrdU()-1; k++) {
				for (int l=0; l<ordv+p.GetOrdV()-1; l++) {
					ncoeffs[count1][count2] = temp[k][l];
					count2++;
				}
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	Vector<int> Ordsu(d.GetNumU());
	Vector<int> Ordsv(d.GetNumV());

	for (int i=0; i<d.GetNumU(); i++) Ordsu[i]=(*ordsu)[i]+p.GetOrdsU()[i]-1;
	for (int i=0; i<d.GetNumV(); i++) Ordsv[i]=(*ordsv)[i]+p.GetOrdsV()[i]-1;
	// form new Ordsu and Ordsv
	// create and return the surface
	return CompPolySurf<T>(d.GetNumU(),d.GetNumV(),Ordsu,Ordsv,ncoeffs,*limitsu,*limitsv);
}
*/
// READ and WRITE

template <class T>
void CompPolySurf<T>::write(std::ostream& os) 
{
	os << "CompPoly Surface\n";
	os << "number of segments in u and v\n";
	os << numu << " " << numv;
	os << "orders in u";
	os << *ordsu;
	os << "orders in v";
	os << *ordsv;
	os << "limits in u are\n";
	os << *limitsu;
	os << "\nlimits in v are\n";
	os << *limitsv;
	os << "\ncoefficients are\n";
	os << *coeffs;
}

template <class T>
void CompPolySurf<T>::read(std::istream& is)
{
	int Numu, Numv;
	std::cout << "number of sgements in u and v\n";
	is >> Numu >> Numv;
	Vector<int> Ordsu(Numu), Ordsv(Numv);
	Vector<double> Limitsu(Numu+1), Limitsv(Numv+1);
	std::cout << "orders of Poly Surfaces in u and v";
	is >> Ordsu >> Ordsv;
	std::cout << "limits in u and v\n";
	is >> Limitsu >> Limitsv;
	int sum1=0;
	int sum2=0;
	for (int i=0; i<Numu; i++) sum1+=(*ordsu)[i];
	for (int j=0; j<Numv; j++) sum2+=(*ordsv)[j];
	Matrix<T> Coeffs(sum1,sum2);
	std::cout << "input coefficients\n";
	is >> Coeffs;
	*this = CompPolySurf<T>(Numu,Numv,Ordsu,Ordsv,Coeffs,Limitsu,Limitsv);
} 


template <class T>
void CompPolySurf<T>::writefile(std::ofstream& ofs) 
{
	ofs << "CompPoly Surface\n";
	ofs << "number of segments in u and v\n";
	ofs << numu << " " << numv;
	ofs << "orders in u";
	ofs << *ordsu;
	ofs << "orders in v";
	ofs << *ordsv;
	ofs << "limits in u are\n";
	ofs << *limitsu;
	ofs << "\nlimits in v are\n";
	ofs << *limitsv;
	
	ofs << "\ncoefficients are\n";
	ofs << *coeffs;
}

template <class T>
void CompPolySurf<T>::readfile(std::ifstream& ifs)
{
	int Numu, Numv;
	ifs >> Numu >> Numv;
	Vector<int> Ordsu(Numu), Ordsv(Numv);
	Vector<double> Limitsu(Numu+1), Limitsv(Numv+1);
	ifs >> Ordsu >> Ordsv;
	ifs >> Limitsu >> Limitsv;
	int sum1=0;
	int sum2=0;
	for (int i=0; i<Numu; i++) sum1+=(*ordsu)[i];
	for (int j=0; j<Numv; j++) sum2+=(*ordsv)[j];
	Matrix<T> Coeffs(sum1,sum2);
	ifs >> Coeffs;
	*this = CompPolySurf<T>(Numu,Numv,Ordsu,Ordsv,Coeffs,Limitsu,Limitsv);
} 



// PRIVATE FUNCTIONS

// elevate the degree of the CompPolySurf by level
template<class T>
CompPolySurf<T> CompPolySurf<T>::ElevateU(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*(ordu+level);
	
	int sum2 = numv*ordv;
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncoeffs(sum1,sum2);

	// elevate each patch in u
	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu+level);
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat  = GetPatch(i+1,j+1).ElevateCPointsU(level);
			colstart=j*ordv;
			count2=colstart;
			for (int k=0; k<ordu+level; k++) {
				for (int l=0; l<ordv; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	// construct array of ord
	// create and return the CompPolySurf
	Vector<int> Ordsu(numu);

	for (int i=0; i<numu; i++) Ordsu[i] = (*ordsu)[i]+level;

	return CompPolySurf<T>(numu, numv, Ordsu, *ordsv, ncoeffs, *limitsu, *limitsv);
}


// elevate the degree of the CompPolySurf by level
template<class T>
CompPolySurf<T> CompPolySurf<T>::ElevateV(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv+level);
	
	int count1, count2;
	int colstart, rowstart;

	Matrix<T> ncoeffs(sum1,sum2);
	// elevate each patch in v
	for (int i=0; i<numu; i++) {
		rowstart=i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).ElevateCPointsV(level);
			colstart=j*(ordv+level);
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv+level; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	Vector<int> Ordsv(numv);

	for (int i=0; i<numv; i++) Ordsv[i] = (*ordsv)[i]+level;
	
	// create and return the CompPolySurf
	return CompPolySurf<T>(numu, numv, *ordsu, Ordsv, ncoeffs, *limitsu, *limitsv);
}

// compute the derivative of the CompPolySurf of order deriv and
// represent the result as another CompPolySurf
template<class T>
CompPolySurf<T> CompPolySurf<T>::DeriveU(int level) const
{   
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*(ordu-level);
	int sum2=numv*ordv;
	
	int count1, count2;
	int rowstart, colstart;
	Matrix<T> ncoeffs(sum1,sum2);
	// derive each patch in u
	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu-level);
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).DeriveCPointsU(level);
			colstart=j*ordv;
			count2=colstart;
			for (int k=0; k<ordu-level; k++) {
				for (int l=0; l<ordv; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	
	Vector<int> Ordsu(numu);

	for (int i=0; i<numu; i++) Ordsu[i] = (*ordsu)[i]-level;

	// create and return the new CompPolySurf
	return CompPolySurf<T>(numu, numv, Ordsu, *ordsv, ncoeffs, *limitsu, *limitsv);
}

// elevate the degree of the CompPolySurf by level
template<class T>
CompPolySurf<T> CompPolySurf<T>::DeriveV(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv-level);
	
	int count1, count2;
	int rowstart, colstart;
	Matrix<T> ncoeffs(sum1,sum2);

	// derive each patch in v
	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).DeriveCPointsV(level);
			colstart=j*(ordv-level);
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv-level; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	
	Vector<int> Ordsv(numv);

	for (int i=0; i<numu; i++) Ordsv[i] = (*ordsv)[i]-level;

	// create and return new ComPolySurf
	return CompPolySurf<T>(numu, numv, *ordsu, Ordsv, ncoeffs, *limitsu, *limitsv);
}

// compute the indefinite integral of the CompPolySurf and represent
// it as a CompPolySurf of one higher degree
template<class T>
CompPolySurf<T> CompPolySurf<T>::IntegrateU() const
{
	// find number of new control points
	int sum1=numu*(ordu+1);
	
	int sum2=numv*ordv;
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncoeffs(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart=i*(ordu+1);
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).IntegrateCPointsU();
			colstart=j*ordv;
			count2=colstart;
			for (int k=0; k<ordu+1; k++) {
				for (int l=0; l<ordv; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	
	Vector<int> Ordsu(numu);

	for (int i=0; i<numu; i++) Ordsu[i] = (*ordsu)[i]+1;

	// create and return CompPolySurf
	return CompPolySurf<T>(numu, numv, Ordsu, *ordsv, ncoeffs, *limitsu, *limitsv);
}	

// integrate the surface in v
template<class T>
CompPolySurf<T> CompPolySurf<T>::IntegrateV() const
{
	// find number of new control points
	int sum1=ordu*numu;
	
	int sum2=numv*(ordv+1);
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncoeffs(sum1,sum2);
	// integrate each patch in v
	for (int i=0; i<numu; i++) {
		rowstart=i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).IntegrateCPointsV();
			colstart=j*(ordv+1);
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv+1; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	
	Vector<int> Ordsv(numv);

	for (int i=0; i<numv; i++) Ordsv[i] = (*ordsv)[i]+1;

	// create and return the new surface
	return CompPolySurf<T>(numu, numv, *ordsu, Ordsv, ncoeffs, *limitsu, *limitsv);
}

// compute the indefinite integral of the CompPolySurf as a CompPolySurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> CompPolySurf<T>::IntegrateCPointsU() const
{
	// find number of new control points
	int sum1=numu*(ordu+1);
	int sum2=numv*ordv;
	
	int count1;
	int count2;
	int rowstart;
	int colstart;
	
	Matrix<T> ncoeffs(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart=i*(ordu+1);
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).IntegrateCPointsU();
			colstart=j*ordv;
			count2=colstart;
			for (int k=0; k<ordu+1; k++) {
				for (int l=0; l<ordv; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	return ncoeffs;
}	


// compute the indefinite integral of the CompPolySurf as a CompPolySurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> CompPolySurf<T>::IntegrateCPointsV() const
{

	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv+1);
	
	int count1=0;
	int count2=0;

	Matrix<T> ncoeffs(sum1,sum2);

	// integrate each patch
	for (int i=0; i<numu; i++) {
		rowstart=i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).IntegrateCPointsV();
			colstart=j*(ordv+1);
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv+1; l++) {
					ncoeffs[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	return ncoeffs;
}	

// find product of two CompPolySurfs assuming differing numbers of
// paches in u and v
template<class T>  
CompPolySurf<T> CompPolySurf<T>::Product2(const CompPolySurf<T>& p) const
{
	// convert to CompBezSurf and multiply

	CompBezSurf<T> b1 = ConvertCompBezSurf();
	CompBezSurf<T> b2 = b.ConvertCompBezSurf();
	CompBezSurf<T> prod = b1.Product(b2);
	// return composite poly surf
	return (prod.ConvertCompPolySurf());
}


#endif










