#ifndef COMPBEZSURF
#define COMPBEZSURF

#include "mathfunctions.h"
#include "BezSurf.h"
#include "CompPolySurf.h"
//#include "BspSurf.h"


// can be of different orders, IGES spec

template<class T>
class CompBezSurf : public Surf<T>  {
private:
	// data
	int numu; // number of segments in u
	int numv; // number of segments in v
	int ordu; // order in u
	int ordv; // order in v
	Ptr<Vector<int> > ordsu;  // array of orders in u
	Ptr<Vector<int> > ordsv;	// array of orders in v
	Ptr<Matrix<T> > cpts;		// control points
	Ptr<Vector<double> > ktsu;	// knots in u
	Ptr<Vector<double> > ktsv;	// knots in v
	Ptr<KnotSet> ksetu;
	Ptr<KnotSet> ksetv;
	Ptr<std::set<double> > limitsetu;
	Ptr<std::set<double> > limitsetv;
	
	// private functions
	CompBezSurf<T> DeriveU(int levu) const; 
	CompBezSurf<T> DeriveV(int levv) const;
	CompBezSurf<T> ElevateU(int levu) const;
	CompBezSurf<T> ElevateV(int levv) const;
	CompBezSurf<T> IntegrateU() const;
	CompBezSurf<T> IntegrateV() const;
	Matrix<T> IntegrateCPointsU() const;
	Matrix<T> IntegrateCPointsV() const;
	CompBezSurf<T> Product2(const CompBezSurf<T>& c) const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class CompBezSurf<double>");
		else {
			std::string s(typeid(T).name()), s1("class CompBezSurf<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	CompBezSurf();
	CompBezSurf(int Numu, int Numv, int Ordu, int Ordv, const Matrix<T>& Cpts, 
		const Vector<double>& Ktsu, const Vector<double>& Ktsv);
	CompBezSurf(int Numu, int Numv, int Ordu, int Ordv, const Matrix<T>& Cpts);
	CompBezSurf(int Numu, int Numv, const Vector<int>& Ordsu, 
		const Vector<int>& Ordsv, const Matrix<T>& Cpts);
	CompBezSurf(int Numu, int Numv, const Vector<int>& Ordsu, 
		const Vector<int>& Ordsv, const Vector<double>& Limitsu, const Vector<double>& Limitsv,
		const Matrix<T>& Cpts);
	
	CompBezSurf(int Numu, int Numv, int Ordu, int Ordv, const Vector<double>& Limitsu, const Vector<double>& Limitsv, const Matrix<T>& Cpts);
	CompBezSurf(const Matrix<BezSurf<T> >& mat, int Numu, int Numv, const Vector<double>& Limitsu, const Vector<double>& Limitsv);
	// access functions
	int GetOrdU() const;
	int GetOrdV() const;
	int GetNumU() const;
	int GetNumV() const;
	int GetNumCPointsU() const;
	int GetNumCPointsV() const;
	int FindSegmentU(double x) const;
	int FindSegmentV(double x) const;
	Matrix<T> GetCPoints() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	Vector<int> GetOrdsU() const;
	Vector<int> GetOrdsV() const;
	BezSurf<T> GetPatch(int i, int j) const;
	BezSurf<T> GetPatchOriginal(int i, int j) const;
	Vector<double> GetLimitsU() const;
	Vector<double> GetLimitsV() const;
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
	
	
	// evaluators
	virtual T operator()(double u, double v) const;
	virtual T operator()(int, int, double u, double v) const;
	T Eval(double u, double v) const;

	// derivatives
	CompBezSurf<T> Derive(int levu, int levv) const;
	virtual T Derive(int levu, int levv, double u, double v) const;

	// degree elevation
	CompBezSurf<T> Elevate(int levu, int levv) const;

	// integration
	CompBezSurf<T> Integrate() const;
	T Integrate(double u1, double u2, double v1, double v2) const;
	Matrix<T> IntegrateCPoints() const;

	// product

	template<class T1>
	CompBezSurf<T> MakeCompatable(const CompBezSurf<T1>& b) const
	{
		return (ConvertBspSurf().MakeBreakCompatable(b.ConvertBspSurf())).ConvertCompBezSurf();
	}
	template<class T1>
	CompBezSurf<T> Product(const CompBezSurf<T1>& b) const
	{
		CompBezSurf<T> d(*this);
		CompBezSurf<T1> e(b);

		if (numu != b.GetNumU() || numv != b.GetNumV()) {
			// needs to add distinct knots from b not in current object
			// into current object (independent of order)
			d = MakeCompatable(b);
			e = b.MakeCompatable(d);
		}

		// find number of new control points
	
		int ordU = ordu+b.GetOrdU()-1;
		int ordV = ordv+b.GetOrdV()-1;

		int sum1= d.GetNumU()*ordU;
		int sum2= d.GetNumV()*ordV;
	
		Matrix<T> ncpts(sum1,sum2);

		int count1;
		int count2;
		int rowstart;
		int colstart;
		for (int i=0; i<d.GetNumU(); i++) {
			rowstart=i*ordU;
			count1=rowstart;
			for (int j=0; j<d.GetNumV(); j++) {
				BezSurf<T> bez1 = d.GetPatch(i+1,j+1);
				BezSurf<T1> bez2 = e.GetPatch(i+1,j+1);
				Matrix<T> temp = bez1.ProductCPoints(bez2);
				colstart = j*ordV;
				count2=colstart;
				for (int k=0; k<ordU; k++) {
					for (int l=0; l<ordV; l++) {
						ncpts[count1][count2] = temp[k][l];
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	
		Vector<int> Ordsu(d.GetNumU()), Ordsv(d.GetNumV());
	//	for (int i=0; i<d.GetNumU(); i++) Ordsu[i] = (*ordsu)[i]+b.GetOrdsU()[i]-1;
	//	for (int i=0; i<d.GetNumV(); i++) Ordsv[i] = (*ordsv)[i]+b.GetOrdsV()[i]-1;
		for (int i=0; i<d.GetNumU(); i++) Ordsu[i] = d.GetOrdsU()[i]+e.GetOrdsU()[i]-1;
		for (int i=0; i<d.GetNumV(); i++) Ordsv[i] = d.GetOrdsV()[i]+e.GetOrdsV()[i]-1;
		// create and return the CompBezSurf
		return CompBezSurf<T>(d.GetNumU(),d.GetNumV(),Ordsu,Ordsv,GetLimitsU(),GetLimitsV(), ncpts);
	}

	// conversion
	Matrix<T> Convert() const;
	BspSurf<T> ConvertBspSurf() const;
	BspSurf<T> ConvertBspSurfRemoveKnots() const;
	CompPolySurf<T> ConvertCompPolySurf() const;

	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};


// CONSTRUCTORS

// default constructor
template<class T>
CompBezSurf<T>::CompBezSurf() : numu(0), numv(0), ordu(0), ordv(0), ordsu(), ordsv(), cpts(), ktsu(), ktsv()
{
}


// constructor builds a CompBezSurf from a matrix of control points
// an order and number of control points in u and v
template<class T>
CompBezSurf<T>::CompBezSurf(int Numu, int Numv, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Matrix<T>& Cpts) :
numu(Numu), numv(Numv), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), cpts(new Matrix<T>(Cpts))
{
	// find max order in u
	ordu = (*ordsu)[0];
	for (int i=1; i<numu; i++) if (ordu < (*ordsu)[i]) ordu=(*ordsu)[i];
	// find max order in v
	ordv = (*ordsv)[0];
	for (int j=1; j<numv; j++) if (ordv < (*ordsv)[j]) ordv=(*ordsv)[j];



	// create knots in u and v
	ktsu = new Vector<double>(Math::CreateKnots(Numu,ordu));
	ktsv = new Vector<double>(Math::CreateKnots(Numv,ordv));
	ksetu = new KnotSet(*ktsu,ordu,(ordu+1)*numu);
	ksetv = new KnotSet(*ktsv,ordv,(ordv+1)*numv);
	cpts = new Matrix<T>(Convert());	
	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsu).end());
	limitsetv = new std::set<double>((*ktsv).begin(),(*ktsv).end());
}


// constructor builds a CompBezSurf from a matrix of control points
// an order and number of control points in u and v
template<class T>
CompBezSurf<T>::CompBezSurf(int Numu, int Numv, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Vector<double>& Limitsu, const Vector<double>& Limitsv,
					const Matrix<T>& Cpts) :
numu(Numu), numv(Numv), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), cpts(new Matrix<T>(Cpts))
{
	// find max order in u
	ordu = (*ordsu)[0];
	for (int i=1; i<numu; i++) if (ordu < (*ordsu)[i]) ordu=(*ordsu)[i];
	// find max order in v
	ordv = (*ordsv)[0];
	for (int j=1; j<numv; j++) if (ordv < (*ordsv)[j]) ordv=(*ordsv)[j];

	// create knots in u and v
	ktsu = new Vector<double>(Math::CreateKnots(Numu,ordu,Limitsu));
	ktsv = new Vector<double>(Math::CreateKnots(Numv,ordv,Limitsv));
	ksetu = new KnotSet(*ktsu,ordu,(ordu+1)*numu);
	ksetv = new KnotSet(*ktsv,ordv,(ordv+1)*numv);
	cpts = new Matrix<T>(Convert());
	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsu).end());
	limitsetv = new std::set<double>((*ktsv).begin(),(*ktsv).end());
}

template<class T>
CompBezSurf<T>::CompBezSurf(const Matrix<BezSurf<T> >& mat, int Numu, int Numv, const Vector<double>& Limitsu, const Vector<double>& Limitsv) :
numu(Numu), numv(Numv), ordsu(new Vector<int>(Numu)), ordsv(new Vector<int>(Numv))
{
	int sum1 = 0;
	int sum2 = 0;
	for (int i=0; i<numu; i++) (*ordsu)[i] = mat[i][0].GetOrdU();
	for (int j=0; j<numv; j++) (*ordsv)[j] = mat[0][j].GetOrdV();

	ordu = (*ordsu)[0];
	for (int i=1; i<numu; i++) if (ordu < (*ordsu)[i]) ordu=(*ordsu)[i];
	// find max order in v
	ordv = (*ordsv)[0];
	for (int j=1; j<numv; j++) if (ordv < (*ordsv)[j]) ordv=(*ordsv)[j];

	for (int i=0; i<numu; i++) sum1 = sum1+ (*ordsu)[i];
	for (int j=0; j<numv; j++) sum2 = sum2+ (*ordsv)[j];

	cpts = new Matrix<T>(sum1,sum2);
	int count1 = 0;
	int count2 = 0;
	int rowstart, colstart;
	for (int i=0; i<numu; i++) {
		rowstart = i*(*ordsu)[i];
		count1 = rowstart;
		for (int j=0; j<numv; j++) {
			colstart = j*(*ordsv)[j];
			count2 = colstart;
			for (int k=0; k<(*ordsu)[i]; k++) {
				for (int l=0; l<(*ordsv)[j]; l++) {
					(*cpts)[count1][count2] = mat[i][j].GetCPoints()[k][l];
					count2++;
				}
				count2 = colstart;
				count1++;
			}
			count1 = rowstart;
		}
	}
	
	// create knots in u and v
	ktsu = new Vector<double>(Math::CreateKnots(Numu,ordu,Limitsu));
	ktsv = new Vector<double>(Math::CreateKnots(Numv,ordv,Limitsv));
	ksetu = new KnotSet(*ktsu,ordu,(ordu+1)*numu);
	ksetv = new KnotSet(*ktsv,ordv,(ordv+1)*numv);
	cpts = new Matrix<T>(Convert());
	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsu).end());
	limitsetv = new std::set<double>((*ktsv).begin(),(*ktsv).end());
}

// constructor builds a CompBezSurf from a matrix of control points,
// knot vectors an order and a number of control points in u and v
template<class T>
CompBezSurf<T>::CompBezSurf(int Numu, int Numv, int Ordu, int Ordv, const Matrix<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv) :
numu(Numu), numv(Numv), ordu(Ordu), ordv(Ordv), ordsu(new Vector<int>(numu)), ordsv(new Vector<int>(numv)), cpts(new Matrix<T>(Cpts)), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)),
ksetu(new KnotSet(*ktsu,ordu,ordu*(numu+1))), ksetv(new KnotSet(*ktsv,ordv,ordv*(numv+1))), 
limitsetu(new std::set<double>(Ktsu.begin(),Ktsu.end())), limitsetv(new std::set<double>(Ktsv.begin(),Ktsv.end()))
{
	// assign order arrays in u and v 
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;
}

// constructor builds a CompBezSurf from a matrix of control points
// orders in u and v
template<class T>
CompBezSurf<T>::CompBezSurf(int Numu, int Numv, int Ordu, int Ordv, const Matrix<T>& Cpts) :
numu(Numu), numv(Numv), ordu(Ordu), ordv(Ordv), ordsu(new Vector<int>(numu)), ordsv(new Vector<int>(numv)), cpts(new Vector<T>(Cpts))
{
	// set orders in u
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	// set orders in v
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;
	// create knots in u and v
	ktsu = new Vector<double>(Math::CreateKnots(numu,ordu));
	ktsv = new Vector<double>(Math::CreateKnots(numv,ordv));

	ksetu = new KnotSet((*ktsu,ordu,(numu+1)*ordu));
	ksetv = new KnotSet((*ktsv,ordv,(numv+1)*ordv));

	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsv).end());
	limitsetv = new std::set<double>((*ktsu).begin(),(*ktsv).end());
}

template<class T>
CompBezSurf<T>::CompBezSurf(int Numu, int Numv, int Ordu, int Ordv, const Vector<double>& Limitsu, const Vector<double>& Limitsv, const Matrix<T>& Cpts) :
numu(Numu), numv(Numv), ordu(Ordu), ordv(Ordv), ordsu(new Vector<int>(numu)), ordsv(new Vector<int>(numv)), cpts(new Vector<T>(Cpts))
{
	// set orders in u
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	// set orders in v
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;
	// create knots in u and v
	ktsu = new Vector<double>(Math::CreateKnots(numu,ordu,Limitsu));
	ktsv = new Vector<double>(Math::CreateKnots(numv,ordv,Limitsv));

	ksetu = new KnotSet((*ktsu,ordu,(numu+1)*ordu));
	ksetv = new KnotSet((*ktsv,ordv,(numv+1)*ordv));

	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsv).end());
	limitsetv = new std::set<double>((*ktsu).begin(),(*ktsv).end());
}




// ACCESS FUNCTIONS
template<class T>
int CompBezSurf<T>::GetNumCPointsU() const
{
	return (numu+1)*ordu;
}

template<class T>
int CompBezSurf<T>::GetNumCPointsV() const
{
	return (numv+1)*ordv;
}



// get the order of the CompBezSurf
template<class T>
inline int CompBezSurf<T>::GetOrdU() const { return ordu; }

// get the order of the CompBezSurf
template<class T>
inline int CompBezSurf<T>::GetOrdV() const { return ordv; }

// get the orders in u
template<class T>
inline Vector<int> CompBezSurf<T>::GetOrdsU() const { return *ordsu; }

// get the orders in v
template<class T>
inline Vector<int> CompBezSurf<T>::GetOrdsV() const { return *ordsv; }

// get the number of patches in u
template<class T>
inline int CompBezSurf<T>::GetNumU() const { return numu; }

// get the number of patches in v
template<class T>
inline int CompBezSurf<T>::GetNumV() const { return numv; }


// get the control points
template<class T>
inline Matrix<T> CompBezSurf<T>::GetCPoints() const { return *cpts; }

// get the knot vector in u
template<class T>
inline Vector<double> CompBezSurf<T>::GetKnotsU() const { return *ktsu; }

// get the knot vector in v
template<class T>
inline Vector<double> CompBezSurf<T>::GetKnotsV() const { return *ktsv; }

// get the knot vector in v
template<class T>
inline Vector<double> CompBezSurf<T>::GetLimitsU() const 
{ 
	Vector<double> limits((*limitsetu).size());
	std::copy((*limitsetu).begin(),(*limitsetu).end(),limits.begin());

	return limits;
}

// get the knot vector in v
template<class T>
inline Vector<double> CompBezSurf<T>::GetLimitsV() const 
{ 
	Vector<double> limits((*limitsetv).size());
	std::copy((*limitsetv).begin(),(*limitsetv).end(),limits.begin());

	return limits;
}


template<class T>
double CompBezSurf<T>::GetLeftLimitU() const 
{ 
	return (*ktsu)[0];
}


template<class T>
double CompBezSurf<T>::GetRightLimitU() const 
{ 
	return (*ktsu)[numu*ordu];
}

template<class T>
double CompBezSurf<T>::GetLeftLimitV() const 
{ 
	return (*ktsv)[0];
}


template<class T>
double CompBezSurf<T>::GetRightLimitV() const 
{ 
	return (*ktsv)[numv*ordv];
}



template<class T>
int CompBezSurf<T>::FindSegmentU(double x) const
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
int CompBezSurf<T>::FindSegmentV(double x) const
{
	std::set<double>::iterator s1 = (*limitsetv).begin();
	std::set<double>::iterator s2 = (*limitsetv).upper_bound(x);
	if (s2 == (*limitsetv).end()) return (*limitsetv).size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);
	return count;
}


// get the i,jth patch of the surface and return as a BezSurf
template<class T>
BezSurf<T> CompBezSurf<T>::GetPatchOriginal(int indu, int indv) const
{
	int ordU = (*ordsu)[indu-1];
	int ordV = (*ordsv)[indv-1];

	// find start index for control points in u and v
	int sum1=0;
	for (int i=0; i<indu; i++) sum1 = sum1+(*ordsu)[i];

	int sum2=0;
	for (int j=0; j<indv; j++) sum2 = sum2+(*ordsv)[j];
	

	// extract control points
	Matrix<T> ncpts(ordU,ordV);

	for (int k=0; k<ordU; k++) 
		for (int l=0; l<ordV; l++)
			ncpts[k][l] = (*cpts)[sum1-ordU+k][sum2-ordV+l];

	// construct knots
	
	// create and return the BezSurf
	return BezSurf<T>(ncpts, ordU, ordV, (*ksetu).GetDistinctKnots()[indu-1],
		(*ksetu).GetDistinctKnots()[indu],(*ksetv).GetDistinctKnots()[indv-1],(*ksetv).GetDistinctKnots()[indv]);
}


// get the i,jth patch of the surface and return as a BezSurf
template<class T>
BezSurf<T> CompBezSurf<T>::GetPatch(int indu, int indv) const
{
	// find start index for control points in u and v
	int sum1= indu*ordu;
	
	int sum2= indv*ordv;
	

	// extract control points
	Matrix<T> ncpts(ordu,ordv);
	for (int k=0; k<ordu; k++) 
		for (int l=0; l<ordv; l++)
			ncpts[k][l] = (*cpts)[sum1-ordu+k][sum2-ordv+l];

	// create and return the BezSurf
	return BezSurf<T>(ncpts, ordu, ordv, (*ksetu).GetDistinctKnots()[indu-1],
		(*ksetu).GetDistinctKnots()[indu],(*ksetv).GetDistinctKnots()[indv-1],(*ksetv).GetDistinctKnots()[indv]);
}




// CONVERSION

// convert to BspSurf form
template<class T>
BspSurf<T> CompBezSurf<T>::ConvertBspSurf() const
{
	return BspSurf<T>(*cpts,*ktsu,*ktsv,ordu,ordv,ordu*numu,ordv*numv);
}

// convert to BspSurf removing all possible knots in u and v
template<class T>
BspSurf<T> CompBezSurf<T>::ConvertBspSurfRemoveKnots() const
{
	return BspSurf<T>(*cpts,*ktsu,*ktsv,ordu,ordv,ordu*numu,ordv*numv).RemovePossibleKnotsUV();
}


// convert each BezSurf patch to the max degree in u and v
template<class T>
Matrix<T> CompBezSurf<T>::Convert() const
{
	// elevate the degree of each patch to max in u and v
	// extract BezSurf patch
	// build array of control points

	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncpts(numu*ordu,numv*ordv);
	// elevate each patch to max order in u and v
	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatchOriginal(i+1,j+1).ElevateCPointsUV(ordu-(*ordsu)[i],ordv-(*ordsv)[j]);	
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv; l++) {
					ncpts[count1][count2]=mat[k][l];
					count2++;
				}
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}

	// create and return the new CompBezSurf
	return ncpts;
}


// convert to a CompPolySurf
template<class T>
CompPolySurf<T> CompBezSurf<T>::ConvertCompPolySurf() const
{
	// take each segment and convert to BezSurf form

	// find the sum of the ords
	int sum1= numu*ordu;
	
	int sum2= numv*ordv;
	
	Matrix<T> ncoeffs(sum1,sum2);
	int count1;
	int count2;
	int rowstart;
	int colstart;
	// extract each patch and convert
	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) {
			PolySurf<T> p = GetPatch(i+1,j+1).ConvertPolySurf();
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv; l++) {
					ncoeffs[count1][count2]=p.GetCoeffs()[k][l];
					count2++;
				}
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	// create and return the CompBezSurf
	return CompPolySurf<T>(numu, numv, *ordsu, *ordsv, ncoeffs, GetLimitsU(), GetLimitsV());
}


// EVALUATORS

// evaluate the CompBezSurf at the point x using de Boor algorithm
template<class T>
T CompBezSurf<T>::operator()(double u, double v) const
{
	// find the segment in u and v
	int i = FindSegmentU(u);
	int j = FindSegmentV(v);
	// get the patch and evaluate
	return GetPatch(i,j)(u,v);
}


// evaluate the CompBezSurf using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompBezSurf<T>::Eval(double u, double v) const
{
	// find the segment in u and v
	int i = FindSegmentU(u);
	int j = FindSegmentV(v);
	// get the patch and evaluate
	return GetPatch(i,j).Eval(u,v);
}


// DEGREE ELEVATION


// elevate the degree of the CompBezSurf by level
template<class T>
CompBezSurf<T> CompBezSurf<T>::Elevate(int levu, int levv) const
{
	// elevate in U and V
	return ElevateU(levu).ElevateV(levv);
}


// DERIVATIVES


// evaluate the derivative of the CompBezSurf of order deriv at a point
// x. Computes the derivative as a CompBezSurf and then evaluates this at x
template<class T>
CompBezSurf<T> CompBezSurf<T>::Derive(int levu, int levv) const
{
	// derive u then v
	return DeriveU(levu).DeriveV(levv);
}

template<class T>
T CompBezSurf<T>::operator() (int valu, int valv, double u, double v) const
{
	return Derive(valu,valv,u,v);
}


// evaluate the derivative of the CompBezSurf of order deriv at a point
// x. Computes the derivative as a CompBezSurf and then evaluates this at x
template<class T>
T CompBezSurf<T>::Derive(int levu, int levv, double u, double v) const
{
	// derive u then v and evaluate
	return Derive(levu,levv)(u,v);
}


// INTEGRATION

// integrate the CompBezSurf between the limits x1 and x2. Computes
// the indefinite integral as a CompBezSurf and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T CompBezSurf<T>::Integrate(double u1, double u2, double v1, double v2) const
{
	// create the indefinite integral;
	CompBezSurf<T> b = IntegrateUV(); 

	// evaluate and subtract
	return (b(u2,v2)-b(u1,v2)-b(u2,v1)+b(u1,v1));
}



// evaluate the derivative of the CompBezSurf of order deriv at a point
// x. Computes the derivative as a CompBezSurf and then evaluates this at x
template<class T>
CompBezSurf<T> CompBezSurf<T>::Integrate() const
{
	return IntegrateU().IntegrateV();
}



// Integrates with respect to u and v and returns the CompBezSurf
template<class T>
Matrix<T> CompBezSurf<T>::IntegrateCPoints() const
{
	Matrix<T> mat = IntegrateCPointsV();
	
	// create and return CompBezSurf
	return CompBezSurf<T>(numu, numv, ordu, ordv+1, mat).IntegrateCPointsU();
}

/*
template<class T>
CompBezSurf<T> CompBezSurf<T>::MakeCompatable(const CompBezSurf<T>& b) const
{
	return (ConvertBspSurf().MakeBreakCompatable(b.ConvertBspSurf())).ConvertCompBezSurf();
}		
*/

/*
// PRODUCT

// compute the product of the CompBezSurf with another CompBezSurf and 
// represent the result as a new CompBezSurf 
// must be compatible, i.e have the same number of patches in u and v
template<class T>  
CompBezSurf<T> CompBezSurf<T>::Product(const CompBezSurf<T>& b) const
{

	CompBezSurf<T> d(*this), e(b);

	if (numu != b.GetNumU() || numv != b.GetNumV()) {
		// needs to add distinct knots from b not in current object
		// into current object (independent of order)
		d = MakeCompatable(b);
		e = b.MakeCompatable(d);
	}

	// find number of new control points
	
	int ordU = ordu+b.GetOrdU()-1;
	int ordV = ordv+b.GetOrdV()-1;

	int sum1= d.GetNumU()*ordU;
	int sum2= d.GetNumV()*ordV;
	
	Matrix<T> ncpts(sum1,sum2);

	int count1;
	int count2;
	int rowstart;
	int colstart;
	for (int i=0; i<d.GetNumU(); i++) {
		rowstart=i*ordU;
		count1=rowstart;
		for (int j=0; j<d.GetNumV(); j++) {
			BezSurf<T> bez1 = d.GetPatch(i+1,j+1);
			BezSurf<T> bez2 = e.GetPatch(i+1,j+1);
			Matrix<T> temp = bez1.ProductCPoints(bez2);
			colstart = j*ordV;
			count2=colstart;
			for (int k=0; k<ordU; k++) {
				for (int l=0; l<ordV; l++) {
					ncpts[count1][count2] = temp[k][l];
					count2++;
				}
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	
	Vector<int> Ordsu(d.GetNumU()), Ordsv(d.GetNumV());
	for (int i=0; i<d.GetNumU(); i++) Ordsu[i] = (*ordsu)[i]+b.GetOrdsU()[i]-1;
	for (int i=0; i<d.GetNumV(); i++) Ordsv[i] = (*ordsv)[i]+b.GetOrdsV()[i]-1;
	// create and return the CompBezSurf
	return CompBezSurf<T>(d.GetNumU(),d.GetNumV(),Ordsu,Ordsv,GetLimitsU(),GetLimitsV(), ncpts);
}
*/

// READ and WRITE

template <class T>
void CompBezSurf<T>::write(std::ostream& os) 
{
	os << "Comp Bezier Surface\n";
	os << "number of segments in u and v\n";
	os << numu << " " << numv;
	os << "orders in u";
	os << *ordsu;
	os << "orders in v";
	os << *ordsv;
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are\n";
	os << *ktsv;
	os << "\ncontrol points are\n";
	os << *cpts;
}

template <class T>
void CompBezSurf<T>::read(std::istream& is)
{
	int Numu, Numv;
	std::cout << "number of sgements in u and v\n";
	is >> Numu >> Numv;
	Vector<int> Ordsu(Numu), Ordsv(Numv);
	std::cout << "orders of Bezier Surfaces in u and v";
	is >> Ordsu >> Ordsv;
	int sum1=0;
	int sum2=0;
	for (int i=0; i<Numu; i++) sum1+=(*ordsu)[i];
	for (int j=0; j<Numv; j++) sum2+=(*ordsv)[j];
	Matrix<T> Cpts(sum1,sum2);
	std::cout << "input control points";
	is >> Cpts;
	*this = CompBezSurf<T>(Numu,Numv,Ordsu,Ordsv,Cpts);
} 


template <class T>
void CompBezSurf<T>::writefile(std::ofstream& ofs)
{
	ofs << "Comp Bezier Surface\n";
	ofs << "number of segments in u and v\n";
	ofs << numu << " " << numv;
	ofs << "orders in u";
	ofs << *ordsu;
	ofs << "orders in v";
	ofs << *ordsv;
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are\n";
	ofs << *ktsv;
	
	ofs << "\ncontrol points are\n";
	ofs << *cpts;
}

template <class T>
void CompBezSurf<T>::readfile(std::ifstream& ifs)
{
	int Numu, Numv;
	ifs >> Numu >> Numv;
	Vector<int> Ordsu(Numu), Ordsv(Numv);
	ifs >> Ordsu >> Ordsv;
	int sum1=0;
	int sum2=0;
	for (int i=0; i<Numu; i++) sum1+=(*ordsu)[i];
	for (int j=0; j<Numv; j++) sum2+=(*ordsv)[j];
	Matrix<T> Cpts(sum1,sum2);
	ifs >> Cpts;
	*this = CompBezSurf<T>(Numu,Numv,Ordsu,Ordsv,Cpts);
} 



// elevate the degree of the CompBezSurf by level
template<class T>
CompBezSurf<T> CompBezSurf<T>::ElevateU(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*(ordu+level);
	
	int sum2=numv*ordv;
	
	int count1;
	int count2;
	int rowstart;
	int colstart;
	Matrix<T> ncpts(sum1,sum2);

	// extract and elevate each patch
	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu+level);
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).ElevateCPointsU(level);
			colstart = j*ordv;
			count2=colstart;
			for (int k=0; k<ordu+level; k++) {
				for (int l=0; l<ordv; l++) {
					ncpts[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}

	Vector<int> Ordsu(numu);
	for (int i=0; i<numu; i++) Ordsu[i] = (*ordsu)[i]+level;
	
	// create and return the CompBezSurf
	return CompBezSurf<T>(numu, numv, Ordsu, *ordsv, GetLimitsU(), GetLimitsV(),ncpts);
}


// elevate the degree of the CompBezSurf by level
template<class T>
CompBezSurf<T> CompBezSurf<T>::ElevateV(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv+level);
	
	int count1;
	int count2;
	int rowstart;
	int colstart;
	Matrix<T> ncpts(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart=i*ordu;	
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).ElevateCPointsV(level);
			colstart = j*(ordv+level);
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv+level; l++) {
					ncpts[count1][count2]=mat[k][l];
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
	
	// create and return the CompBezSurf
	return CompBezSurf<T>(numu, numv, *ordsu, Ordsv, GetLimitsU(), GetLimitsV(),ncpts);
}
// compute the derivative of the CompBezSurf of order deriv and
// represent the result as another CompBezSurf
template<class T>
CompBezSurf<T> CompBezSurf<T>::DeriveU(int level) const
{   
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*(ordu-level);
	
	int sum2=numv*ordv;
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncpts(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu-level);
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).DeriveCPointsU(level);
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<ordu-level; k++) {
				for (int l=0; l<ordv; l++) {
					ncpts[count1][count2]=mat[k][l];
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
	
	// create and return the CompBezSurf
	return CompBezSurf<T>(numu, numv, Ordsu, *ordsv, GetLimitsU(), GetLimitsV(), ncpts);
}

// elevate the degree of the CompBezSurf by level
template<class T>
CompBezSurf<T> CompBezSurf<T>::DeriveV(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv-level);
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncpts(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart=i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).DeriveCPointsV(level);
			colstart = j*(ordv-level);
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv-level; l++) {
					ncpts[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	Vector<int> Ordsv(numv);
	for (int i=0; i<numv; i++) Ordsv[i] = (*ordsv)[i]-level;
	
	// create and return BezSurf
	return CompBezSurf<T>(numu, numv, *ordsu, Ordsv, GetLimitsU(), GetLimitsV(), ncpts);
}

// compute the indefinite integral of the CompBezSurf and represent
// it as a CompBezSurf of one higher degree
template<class T>
CompBezSurf<T> CompBezSurf<T>::IntegrateU() const
{
	// find number of new control points
	int sum1=numu*(ordu+1);
	
	int sum2=numv*ordv;
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncpts(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu+1);
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).IntegrateCPointsU();
			colstart = j*ordv;
			count2=colstart;
			for (int k=0; k<ordu+1; k++) {
				for (int l=0; l<ordv; l++) {
					ncpts[count1][count2]=mat[k][l];
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
	
	// create and return the CompBezSurf
	return CompBezSurf<T>(numu, numv, Ordsu, *ordsv, GetLimitsU(), GetLimitsV(), ncpts);
}	

// elevate the degree of the CompBezSurf by level
template<class T>
CompBezSurf<T> CompBezSurf<T>::IntegrateV() const
{
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv+1);
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncpts(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).IntegrateCPointsV();
			colstart=j*(ordv+1);
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv+1; l++) {
					ncpts[count1][count2]=mat[k][l];
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
	
	// create and return CompBezSurf
	return CompBezSurf<T>(numu, numv, *ordsu, Ordsv, GetLimitsU(), GetLimitsV(), ncpts);
}

// compute the indefinite integral of the CompBezSurf as a CompBezSurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> CompBezSurf<T>::IntegrateCPointsU() const
{
	// find number of new control points
	int sum1=numu*(ordu+1);
	
	int sum2=numv*ordv;
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncpts(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart=i*(ordu+1);
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).IntegrateCPointsU();
			colstart=j*ordv;
			count2=colstart;
			for (int k=0; k<ordu+1; k++) {
				for (int l=0; l<ordv; l++) {
					ncpts[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;
				count1++;
			}
			count1=rowstart;
		}
	}
	return ncpts;
}	


// compute the indefinite integral of the CompBezSurf as a CompBezSurf
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix<T> CompBezSurf<T>::IntegrateCPointsV() const
{

	// find number of new control points
	int sum1=ordu*numu;
	
	int sum2=numv*(ordv+1);
	
	int count1;
	int count2;
	int rowstart;
	int colstart;

	Matrix<T> ncpts(sum1,sum2);

	for (int i=0; i<numu; i++) {
		rowstart=i*ordu;
		count1=rowstart;
		for (int j=0; j<numv; j++) {
			Matrix<T> mat = GetPatch(i+1,j+1).IntegrateCPointsV();
			colstart=j*(ordv+1);
			count2=colstart;
			for (int k=0; k<ordu; k++) {
				for (int l=0; l<ordv+1; l++) {
					ncpts[count1][count2]=mat[k][l];
					count2++;
				}	
				count2=colstart;;
				count1++;
			}
			count1=rowstart;
		}
	}
	return ncpts;
}	

// compute the product of the CompBezSurf with another CompBezSurf and 
// represent the result as a new CompBezSurf 
// need not be compatable
template<class T>  
CompBezSurf<T> CompBezSurf<T>::Product2(const CompBezSurf<T>& b) const
{
	// convert to CompBspSurf and multiply

	BspSurf<T> b1 = ConvertBspSurf();
	BspSurf<T> b2 = b.ConvertBspSurf();
	BspSurf<T> prod = b1.Product(b2);
	// return composite poly surf
	return (prod.ConvertCompBezSurf());
}


#endif