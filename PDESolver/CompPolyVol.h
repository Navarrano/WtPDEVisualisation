#ifndef COMPOLYVOL
#define COMPOLYVOL

#include "PolyVol.h"
#include "CompBezCurv.h"
#include "CompBezVol.h"
#include "BspVol.h"
#include "mathfunctions.h"

template<class T>
class CompBezVol;

// can be of different orders, IGES spec

template<class T>
class CompPolyVol : public Vol<T> {
private:
	// data
	int numu; // number of segments in u
	int numv; // number of segments in v
	int numw;
	int ordu; // max ord in u
	int ordv; // max ord in v
	int ordw;
	Ptr<Vector<int> > ordsu; // orders in u
	Ptr<Vector<int> > ordsv;	// orders in v
	Ptr<Vector<int> > ordsw;
	Ptr<Matrix3D<T> > coeffs;	// coeffs
	Ptr<Vector<double> > limitsu; // ranges in u
	Ptr<Vector<double> > limitsv;	// ranges in v
	Ptr<Vector<double> > limitsw;
	Ptr<std::set<double> > limitsetu;
	Ptr<std::set<double> > limitsetv;
	Ptr<std::set<double> > limitsetw;



	// private functions
	CompPolyVol<T> DeriveU(int level) const; 
	CompPolyVol<T> DeriveV(int level) const; 
	CompPolyVol<T> DeriveW(int level) const; 
	CompPolyVol<T> DeriveUV(int levu, int levv) const;
	CompPolyVol<T> DeriveUW(int levu, int levw) const;
	CompPolyVol<T> DeriveVW(int levv, int levw) const;
	CompPolyVol<T> ElevateU(int levu) const;
	CompPolyVol<T> ElevateV(int levv) const;
	CompPolyVol<T> ElevateW(int levw) const;
	CompPolyVol<T> ElevateUV(int levu, int levw) const;
	CompPolyVol<T> ElevateUW(int levu, int levw) const;
	CompPolyVol<T> ElevateVW(int levv, int levw) const;
	CompPolyVol<T> IntegrateU() const;
	CompPolyVol<T> IntegrateV() const;
	CompPolyVol<T> IntegrateW() const;
	CompPolyVol<T> IntegrateUV() const;
	CompPolyVol<T> IntegrateUW() const;
	CompPolyVol<T> IntegrateVW() const;
	Matrix3D<T> IntegrateCPointsU() const;
	Matrix3D<T> IntegrateCPointsV() const;
	Matrix3D<T> IntegrateCPointsW() const;
	Matrix3D<T> IntegrateCPointsUV() const;
	Matrix3D<T> IntegrateCPointsUW() const;
	Matrix3D<T> IntegrateCPointsVW() const;
	CompPolyVol<T> Product2(const CompPolyVol<T>& c) const;
	template<class T1>
	CompPolyVol<T> MakeCompatable(const CompPolyVol<T1>& b) const
	{
		return (ConvertBspVol().MakeBreakCompatable(b.ConvertBspVol())).ConvertCompPolyVol();
	}
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class CompPolyVol<double>");
		else {
			std::string s(typeid(T).name()), s1("class CompPolyVol<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	CompPolyVol();
	CompPolyVol(int Numu, int Numv, int Numw, int Ordu, int Ordv, int Ordw, const Matrix3D<T>& Coeffs, const Vector<double>& Limitsu, 
		const Vector<double>& Limitsv, const Vector<double>& Limitsw);
	CompPolyVol(int Numu, int Numv, int Numw, int Ordu, int Ordv, int Ordw, const Matrix3D<T>& Coeffs);
	CompPolyVol(int Numu, int Numv, int Numw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, 
		const Vector<int>& Ordsw, const Matrix3D<T>& Coeffs, const Vector<double>& Limitsu, const Vector<double>& Limitsv, const Vector<double>& Limitsw);
	CompPolyVol(int Numu, int Numv, int Numw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Vector<int>& Ordsw, const Matrix3D<T>& Coeffs);


	// access functions
	int GetOrdU() const;
	int GetOrdV() const;
	int GetOrdW() const;
	Vector<double> GetLimitsU() const;
	Vector<double> GetLimitsV() const;
	Vector<double> GetLimitsW() const;
	int FindSegmentU(double x) const;
	int FindSegmentV(double x) const;
	int FindSegmentW(double x) const;
	double GetLeftLimitU(int indexu, int indexv, int indexw) const;
	double GetRightLimitU(int indexu, int indexv, int indexw) const;
	double GetLeftLimitV(int indexu, int indexv, int indexw) const;
	double GetRightLimitV(int indexu, int indexv, int indexw) const;
	double GetLeftLimitW(int indexu, int indexv, int indexw) const;
	double GetRightLimitW(int indexu, int indexv, int indexw) const;

	double GetLeftLimitU() const;
	double GetRightLimitU() const;
	double GetLeftLimitV() const;
	double GetRightLimitV() const;
	double GetLeftLimitW() const;
	double GetRightLimitW() const;
	int GetNumU() const;
	int GetNumV() const;
	int GetNumW() const;
	Matrix3D<T> GetCoeffs() const;
	Vector<int> GetOrdsU() const;
	Vector<int> GetOrdsV() const;
	Vector<int> GetOrdsW() const;
	Matrix3D<T> GetCoeffsPoly(int indexu, int indexv, int indexw) const;
	PolyVol<T> GetVolume(int indexu, int indexv, int indexw) const;
	PolyVol<T> GetVolumeOriginal(int i, int j, int k) const;

	// evaluators
	virtual T operator()(double u, double v, double w) const;
	virtual T operator()(int levu, int levv, int levw, double u, double v, double w) const;
	T operator()(double u, double v, double w, int indexu, int indexv, int indexw) const;
	T Eval(double u, double v, double w) const;
	T Eval(double u, double v, double w, int indexu, int indexv, int indexw) const;
	Matrix3D<T> ComputePoints(int,int,int) const;

	// conversion
	Matrix3D<T> Convert() const;
	CompBezVol<T> ConvertCompBezVol() const;
	BspVol<T> ConvertBspVol() const;
	
	// derivatives
	CompPolyVol<T> Derive(int levu, int levv, int levw) const;
	virtual T Derive(int levu, int levv, int levw, double u, double v, double w) const;

	// degree elevation
	CompPolyVol<T> Elevate(int levu, int levv, int levw) const;

	// integration
	Matrix3D<T> IntegrateCPoints() const;
	CompPolyVol<T> Integrate() const;
	T Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const;

	// product
	template<class T1>
	CompPolyVol<T> Product(const CompPolyVol<T1>& b) const
	{
		CompPolyVol<T> d(*this), e(b);

		if (numu != p.GetNumU() || numv != p.GetNumV() || numw != p.GetNumW()) {
			d = MakeCompatable(b);
			e = b.MakeCompatable(d);
		}


		// find number of new control points
		int sum1=d.GetNumU()*(ordu+p.GetOrdU()-1);	
		int sum2=d.GetNumV()*(ordv+p.GetOrdV()-1);
		int sum3=d.GetNumW()*(ordw+p.GetOrdW()-1);
	
		Matrix3D<T> ncoeffs(sum1,sum2,sum3);

		int count1;
		int count2;
		int count3;
		int rowstart;
		int colstart;
		int polstart;

		// find product patch by patch
		for (int i=0; i<d.GetNumU(); i++) {
			rowstart = i*(ordu+p.GetOrdU()-1);
			count1=rowstart;
			for (int j=0; j<d.GetNumV(); j++) {
				colstart=j*(ordv+p.GetOrdV()-1);
				count2=colstart;
				for (int k=0; k<d.GetNumW(); k++) {
					polstart=k*(ordw+p.GetOrdW()-1);
					count3=polstart;
					PolyVol<T> p1 = d.GetVolume(i+1,j+1,k+1);
					PolyVol<T> p2 = e.GetVolume(i+1,j+1,k+1);
					Matrix3D<T> temp = p1.ProductCPoints(p2);
					for (int l=0; l<ordu+p.GetOrdU()-1; l++) {
						for (int m=0; m<ordv+p.GetOrdV()-1; m++) {
							for (int q=0; q<ordw+p.GetOrdW()-1; q++) {
								ncoeffs[count3][count1][count2] = temp[q][l][m];
								count3++;
							}
							count3=polstart;
							count2++;
						}
						count2=colstart;
						count1++;
					}
					count1=rowstart;
				}
			}
		}

		Vector<int> Ordsu(d.GetNumU());
		for (int i=0; i<d.GetNumU(); i++) Ordsu[i]=(*ordsu)[i]+p.GetOrdsU()[i]-1;


		Vector<int> Ordsv(d.GetNumV());
		for (int i=0; i<d.GetNumV(); i++) Ordsv[i]=(*ordsv)[i]+p.GetOrdsV()[i]-1;
	
		Vector<int> Ordsw(d.GetNumW());
		for (int i=0; i<d.GetNumW(); i++) Ordsw[i]=(*ordsw)[i]+p.GetOrdsW()[i]-1;
	
		// create and return the Vol
		return CompPolyVol<T>(d.GetNumU(),d.GetNumV(),d.GetNumW(),Ordsu,Ordsv,Ordsw,ncoeffs,*limitsu,*limitsv,*limitsw);
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
CompPolyVol<T>::CompPolyVol() : numu(0), numv(0), numw(0), ordu(0), ordv(0), ordw(0), ordsu(), ordsv(), ordsw(), coeffs(), limitsu(),limitsv(), limitsw()
{
}

// constructor builds a CompPolyVol from a Vector of control points,
// an order and number of control points and limits in u and v
template<class T>
CompPolyVol<T>::CompPolyVol(int Numu, int Numv, int Numw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, 
				const Vector<int>& Ordsw, const Matrix3D<T>& Coeffs, 
							  const Vector<double>& Limitsu, const Vector<double>& Limitsv, const Vector<double>& Limitsw) :
	numu(Numu), numv(Numv), numw(Numw), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), ordsw(new Vector<int>(Ordsw)), 
		coeffs(new Matrix3D<T>(Coeffs)), limitsu(new Vector<double>(Limitsu)), 
		limitsv(new Vector<double>(Limitsv)), 
		limitsw(new Vector<double>(Limitsw)),
		limitsetu(new std::set<double>((*limitsu).begin(),(*limitsu).end())),
		limitsetv(new std::set<double>((*limitsv).begin(),(*limitsv).end())),
		limitsetw(new std::set<double>((*limitsw).begin(),(*limitsw).end()))
{
	// find max order in u
	ordu = (*ordsu)[0];
	for (int i=1; i<numu; i++) if (ordu < (*ordsu)[i]) ordu=(*ordsu)[i];

	// find max order in v
	ordv = (*ordsv)[0];
	for (int j=1; j<numv; j++) if (ordv < (*ordsv)[j]) ordv=(*ordsv)[j];

	// find max order in w
	ordw = (*ordsw)[0];
	for (int k=1; k<numw; k++) if (ordw < (*ordsw)[k]) ordw=(*ordsw)[k];
	coeffs = new Matrix3D<T>(Convert());
}

// constructor builds a CompPolyVol from a Vector of control points,
// an order and number of control points in u and v
template<class T>
CompPolyVol<T>::CompPolyVol(int Numu, int Numv, int Numw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Vector<int>& Ordsw,
		const Matrix3D<T>& Coeffs) :
numu(Numu), numv(Numv), numw(Numw), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), ordsw(new Vector<int>(Ordsw)), coeffs(new Matrix3D<T>(Coeffs)), 
limitsu(new Vector<double>(numu+1)), limitsv(new Vector<double>(numv+1)), 
limitsw(new Vector<double>(numw+1))
{
	// find max order in u
	ordu = (*ordsu)[0];
	for (int i=1; i<numu; i++) if (ordu < (*ordsu)[i]) ordu=(*ordsu)[i];

	// find max order in v
	ordv = (*ordsv)[0];
	for (int j=1; j<numv; j++) if (ordv < (*ordsv)[j]) ordv=(*ordsv)[j];

	// find max order in w
	ordw = (*ordsw)[0];
	for (int k=1; k<numw; k++) if (ordw < (*ordsw)[k]) ordw=(*ordsw)[k];

	// limits are increasing integers in u and v
	for (int i=0; i<=numu; i++) (*limitsu)[i]=(double)i;
	for (int j=0; j<=numv; j++) (*limitsv)[j]=(double)j;
	for (int k=0; k<=numw; k++) (*limitsw)[k]=(double)k;
	limitsetu = new std::set<double>((*limitsu).begin(),(*limitsu).end());
	limitsetv = new std::set<double>((*limitsv).begin(),(*limitsv).end());
	limitsetw = new std::set<double>((*limitsw).begin(),(*limitsw).end());
	coeffs = new Matrix3D<T>(Convert());
}


// constructor builds a CompPolyVol from a Vector of control points,
// an order and number of control points and limits in u and v 
template<class T>
CompPolyVol<T>::CompPolyVol(int Numu, int Numv, int Numw, int Ordu, int Ordv, int Ordw, 
const Matrix3D<T>& Coeffs, const Vector<double>& Limitsu, const Vector<double>& Limitsv, const Vector<double>& Limitsw) :
numu(Numu), numv(Numv), numw(Numw), ordu(Ordu), ordv(Ordv), ordw(Ordw), 
ordsu(new Vector<int>(Numu)), ordsv(new Vector<int>(Numv)), 
ordsw(new Vector<int>(Numw)), coeffs(new Matrix3D<T>(Coeffs)), 
limitsu(new Vector<double>(Limitsu)), limitsv(new Vector<double>(Limitsv)), 
limitsw(new Vector<double>(Limitsw)),
limitsetu(new std::set<double>((*limitsu).begin(),(*limitsu).end())),
limitsetv(new std::set<double>((*limitsv).begin(),(*limitsv).end())),
limitsetw(new std::set<double>((*limitsw).begin(),(*limitsw).end()))
{
	// assign orders in u
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	// assign orders in v
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;

	// assign orders in w
	for (int k=0; k<numw; k++) (*ordsw)[k]=ordw;
}


// constructor builds a CompPolyVol from a Vector of control points,
// an order and number of control points in u and v 
template<class T>
CompPolyVol<T>::CompPolyVol(int Numu, int Numv, int Numw, int Ordu, int Ordv, int Ordw, const Matrix3D<T>& Coeffs) :
numu(Numu), numv(Numv), numw(Numw), ordu(Ordu), ordv(Ordv), ordw(Ordw), 
ordsu(new Vector<int>(Numu)), ordsv(new Vector<int>(Numv)), ordsw(new Vector<int>(Numw)), 
coeffs(new Matrix3D<T>(Coeffs)), 
limitsu(new Vector<double>(Numu+1)), 
limitsv(new Vector<double>(Numv+1)),
limitsw(new Vector<double>(Numw+1))
{
	// assign orders in u
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	// assign orders in v
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;

	// assign orders in w
	for (int k=0; k<numw; k++) (*ordsw)[k]=ordw;

	// find limits in u
	for (int i=0; i<=numu; i++) (*limitsu)[i]=(double)i;
	// find limits in v
	for (int j=0; j<=numv; j++) (*limitsv)[j]=(double)j;

	// find limits in w
	for (int k=0; k<=numw; k++) (*limitsw)[k]=(double)k;

	limitsetu = new std::set<double>((*limitsu).begin(),(*limitsu).end());
	limitsetv = new std::set<double>((*limitsv).begin(),(*limitsv).end());
	limitsetw = new std::set<double>((*limitsw).begin(),(*limitsw).end());
}


// ACCESS FUNCTIONS



// get the order of the CompPolyVol
template<class T>
inline int CompPolyVol<T>::GetOrdU() const { return ordu; }

// get the order of the CompPolyVol
template<class T>
inline int CompPolyVol<T>::GetOrdV() const { return ordv; }

// get the order of the CompPolyVol
template<class T>
inline int CompPolyVol<T>::GetOrdW() const { return ordw; }


// get the orders in u
template<class T>
inline Vector<int> CompPolyVol<T>::GetOrdsU() const { return *ordsu; }


// get the orders in v
template<class T>
inline Vector<int> CompPolyVol<T>::GetOrdsV() const { return *ordsv; }


// get the orders in w
template<class T>
inline Vector<int> CompPolyVol<T>::GetOrdsW() const { return *ordsw; }

// get the number of segments in u
template<class T>
inline int CompPolyVol<T>::GetNumU() const { return numu; }

// get the number of segments in v
template<class T>
inline int CompPolyVol<T>::GetNumV() const { return numv; }


// get the number of segments in v
template<class T>
inline int CompPolyVol<T>::GetNumW() const { return numw; }

// get the coefficients
template<class T>
inline Matrix3D<T> CompPolyVol<T>::GetCoeffs() const { return *coeffs; }


// get the limit set in u
template<class T>
inline Vector<double> CompPolyVol<T>::GetLimitsU() const { return *limitsu; }

// get the limit set in v
template<class T>
inline Vector<double> CompPolyVol<T>::GetLimitsV() const { return *limitsv; }


// get the limit set in v
template<class T>
inline Vector<double> CompPolyVol<T>::GetLimitsW() const { return *limitsw; }

template<class T>
double CompPolyVol<T>::GetLeftLimitU() const 
{ 
	return (*limitsu)[0];
}


template<class T>
double CompPolyVol<T>::GetRightLimitU() const 
{ 
	return (*limitsu)[numu];
}

template<class T>
double CompPolyVol<T>::GetLeftLimitV() const 
{ 
	return (*limitsv)[0];
}


template<class T>
double CompPolyVol<T>::GetRightLimitV() const 
{ 
	return (*limitsv)[numv];
}

template<class T>
double CompPolyVol<T>::GetLeftLimitW() const 
{ 
	return (*limitsw)[0];
}


template<class T>
double CompPolyVol<T>::GetRightLimitW() const 
{ 
	return (*limitsw)[numw];
}


template<class T>
int CompPolyVol<T>::FindSegmentU(double x) const
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
int CompPolyVol<T>::FindSegmentV(double x) const
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


template<class T>
int CompPolyVol<T>::FindSegmentW(double x) const
{
	std::set<double>::iterator s1 = (*limitsetw).begin();
	std::set<double>::iterator s2 = (*limitsetw).upper_bound(x);
	if (s2 == (*limitsetw).end()) return (*limitsetw).size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);
	return count;
}

// get the i,j patch of the composite Volace
template<class T>
PolyVol<T> CompPolyVol<T>::GetVolumeOriginal(int i, int j, int k) const
{
	// order in u and v
	int ordU = (*ordsu)[i-1];
	int ordV = (*ordsv)[j-1];
	int ordW = (*ordsw)[k-1];

	// find start index for control points in u and v
	int sum1=0;
	for (int l=0; l<i; l++) sum1 = sum1+(*ordsu)[l];

	int sum2=0;
	for (int l=0; l<j; l++) sum2 = sum2+(*ordsv)[l];
	
	int sum3=0;
	for (int l=0; l<k; l++) sum3 = sum3+(*ordsw)[l];

	// extract control points
	Matrix3D<T> ncoeffs(ordU,ordV,ordW);
	for (int l=0; l<ordU; l++) 
		for (int p=0; p<ordV; p++)
			for (int q=0; q<ordW; q++) 
				ncoeffs[q][l][p] = (*coeffs)[sum3-ordW+q][sum1-ordU+l][sum2-ordV+p];

	// create and return the PolyVol
	return PolyVol<T>(ncoeffs, ordU, ordV, ordW, (*limitsu)[i-1], (*limitsu)[i], (*limitsv)[j-1], (*limitsv)[j], (*limitsw)[k-1], (*limitsw)[k]);
}


// get the i,j patch of the composite Volace
template<class T>
PolyVol<T> CompPolyVol<T>::GetVolume(int i, int j, int k) const
{
	// find start index for control points in u and v
	int sum1=i*ordu;
	
	int sum2=j*ordv;
	
	int sum3=k*ordw;
	
	// extract control points
	Matrix3D<T> ncoeffs(ordu,ordv,ordw);
	for (int l=0; l<ordu; l++) 
		for (int p=0; p<ordv; p++)
			for (int q=0; q<ordw; q++) 
				ncoeffs[q][l][p] = (*coeffs)[sum3-ordw+q][sum1-ordu+l][sum2-ordv+p];

	// create and return the PolyVol
	return PolyVol<T>(ncoeffs, ordu, ordv, ordw, (*limitsu)[i-1], (*limitsu)[i], (*limitsv)[j-1], (*limitsv)[j], (*limitsw)[k-1], (*limitsw)[k]);
}


// CONVERSIONS

// convert to BspVol
template<class T>
BspVol<T> CompPolyVol<T>::ConvertBspVol() const
{
	return ConvertCompBezVol().ConvertBspVol();
}

// convert to composite Bezier form
template<class T>
CompBezVol<T> CompPolyVol<T>::ConvertCompBezVol() const
{
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw*ordw
	

	Matrix3D<T> ncpts(sum1,sum2,sum3);
	
	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> temp = GetVolume(i+1,j+1,k+1).ConvertBezVolCPoints();
				cout << "temp" << temp;
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=temp[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}
	// create and return the composite Bezier Volace
	return CompBezVol<T>(numu, numv, numw, *ordsu, *ordsv, *ordsw, ncpts, *limitsu, *limitsv, *limitsw);
}


// convert each patch to max degree in u and v
template<class T>
Matrix3D<T> CompPolyVol<T>::Convert() const
{
	// elevate the degree of each segment to max
	// build array of control points
	
	Matrix3D<T> ncoeffs(numu*ordu,numv*ordv,numw*ordw);
	
	// elevate each patch as necessary
	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*(*ordsu)[i];
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*(*ordsv)[j];
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolumeOriginal(i+1,j+1,k+1).ElevateUVW(ordu-(*ordsu)[i],ordv-(*ordsv)[j],ordw-(*ordsw)[k]).GetCoeffs();
				int polstart = k*(*ordsw)[k];
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}
	// create and return the composite Poly Vol
	return ncoeffs;
}

// EVALUATORS

// evaluate the CompPolyVol at the point x 
template<class T>
T CompPolyVol<T>::operator()(double u, double v, double w) const
{
	// Find segment in u and v
	int indexu = FindSegmentU(u);
	int indexv = FindSegmentV(v);
	int indexw = FindSegmentW(w);

	// extract and evaluate segment
	return GetVolume(indexu,indexv,indexw)(u,v,w);
}

// evaluate the segment index at the point x 
template<class T>
T CompPolyVol<T>::operator()(double u, double v, double w, int indexu, int indexv, int indexw) const
{
	// extract and evaluate segment
	return GetVolume(indexu,indexv,indexw)(u,v,w);
}


// evaluate the CompPolyVol using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompPolyVol<T>::Eval(double u, double v, double w) const
{
	// Find segment
	int indu = FindSegment(u);
	int indv = FindSegment(v);
	int indw = FindSegment(w);

	// evaluate segemnt
	return GetVolume(indu,indv,indw).Eval(u,v,w);
}

// evaluate the CompPolyVol using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompPolyVol<T>::Eval(double u, double v, double w, int indexu, int indexv, int indexw) const
{
	// extract and evaluate segment
	return GetVolume(indexu,indexv,indexw).Eval(u,v,w);
}

// DEGREE ELEVATION

// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::Elevate(int levu, int levv, int levw) const
{
	// elevate in u then v then w
	return (ElevateU(levu).ElevateV(levv)).ElevateW(levw);
}


// DERIVATIVES

// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
CompPolyVol<T> CompPolyVol<T>::Derive(int levu, int levv, int levw) const
{
	// derive in u then v then w
	return (DeriveU(levu).DeriveV(levv)).DeriveW(levw);
}


template<class T>
T CompPolyVol<T>::operator() (int valu, int valv, int valw, double u, double v, double w) const
{
	return Derive(valu,valv,valw,u,v,w);
}



// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
T CompPolyVol<T>::Derive(int levu, int levv, int levw, double u, double v, double w) const
{
	return Derive(levu,levv,levw)(u,v,w);
}


// INTEGRATION

// TO DO
// integrate the CompPolyVol between the limits x1 and x2. Computes
// the indefinite integral as a CompPolyVol and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T CompPolyVol<T>::Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	// create the indefinite integral;
	CompPolyVol<T> p = Integrate(); 
	
	return p(u2,v2,w2)-p(u2,v2,w1)-p(u2,v1,w2)-p(u1,v2,w2)+
		p(u2,v1,w1)+p(u1,v2,w1)+p(u1,v1,w2)-p(u1,v1,w1);
}


// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
CompPolyVol<T> CompPolyVol<T>::Integrate() const
{
	// integrate in u then v
	return (IntegrateU().IntegrateV()).IntegrateW();
}


// READ and WRITE

template <class T>
void CompPolyVol<T>::write(std::ostream& os)
{
	os << "CompPoly Volume\n";
	os << "number of segments in u and v and w\n";
	os << numu << " " << numv << " " << numw;
	os << "orders in u";
	os << *ordsu;
	os << "orders in v";
	os << *ordsv;
	os << "orders in w";
	os << *ordsw;
	os << "limits in u are\n";
	os << *limitsu;
	os << "\nlimits in v are\n";
	os << *limitsv;
	os << "\nlimits in w are\n";
	os << *limitsw;
	os << "\ncoefficients are\n";
	os << *coeffs;
}

template <class T>
void CompPolyVol<T>::read(std::istream& is)
{
	int Numu, Numv, Numw;
	std::cout << "number of sgements in u and v and w\n";
	is >> Numu >> Numv >> Numw;
	Vector<int> Ordsu(Numu), Ordsv(Numv), Ordsw(Numw);
	Vector<double> Limitsu(Numu+1), Limitsv(Numv+1), Limitsw(Numw+1);
	std::cout << "orders of Poly Volumes in u and v and w";
	is >> Ordsu >> Ordsv >> Ordsw;
	std::cout << "limits in u and v and w\n";
	is >> Limitsu >> Limitsv >> Limitsw;
	int sum1=0;
	int sum2=0;
	int sum3=0;
	for (int i=0; i<Numu; i++) sum1+=(*ordsu)[i];
	for (int j=0; j<Numv; j++) sum2+=(*ordsv)[j];
	for (int k=0; k<Numw; k++) sum3+=(*ordsw)[k];
	Matrix3D<T> Coeffs(sum1,sum2,sum3);
	std::cout << "input coefficients\n";
	is >> Coeffs;
	*this = CompPolyVol<T>(Numu,Numv,Numw,Ordsu,Ordsv,Ordsw,Coeffs,Limitsu,Limitsv,Limitsw);
} 


template <class T>
void CompPolyVol<T>::writefile(std::ofstream& ofs) 
{
	ofs << "CompPoly Volume\n";
	ofs << "number of segments in u and v and w\n";
	ofs << numu << " " << numv << " " << numw;
	ofs << "orders in u";
	ofs << *ordsu;
	ofs << "orders in v";
	ofs << *ordsv;
	ofs << "orders in w";
	ofs << *ordsw;
	ofs << "limits in u are\n";
	ofs << *limitsu;
	ofs << "\nlimits in v are\n";
	ofs << *limitsv;
	ofs << "\nlimits in w are\n";
	ofs << *limitsw;
	ofs << "\ncoefficients are\n";
	ofs << *coeffs;
}

template <class T>
void CompPolyVol<T>::readfile(std::ifstream& ifs)
{
	int Numu, Numv, Numw;
	ifs >> Numu >> Numv >> Numw;
	Vector<int> Ordsu(Numu), Ordsv(Numv), Ordsw(Numw);
	Vector<double> Limitsu(Numu+1), Limitsv(Numv+1), Limitsw(Numw+1);
	ifs >> Ordsu >> Ordsv >> Ordsw;
	ifs >> Limitsu >> Limitsv >> Limitsw;
	int sum1=0;
	int sum2=0;
	int sum3=0;
	for (int i=0; i<Numu; i++) sum1+=(*ordsu)[i];
	for (int j=0; j<Numv; j++) sum2+=(*ordsv)[j];
	for (int k=0; k<Numw; k++) sum3+=(*ordsw)[k];
	Matrix3D<T> Coeffs(sum1,sum2,sum3);
	ifs >> Coeffs;
	*this = CompPolyVol<T>(Numu,Numv,Numw,Ordsu,Ordsv,Ordsw,Coeffs,Limitsu,Limitsv,Limitsw);
} 

// PRIVATE FUNCTIONS

// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::ElevateU(int level) const
{
	if (level <=0) return *this;
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*(ordu+level);
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw*ordw;
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*(ordu+level);
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).ElevateCPointsU(level);
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu+level; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsu(numu);
	for (int i=0; i<numu; i++) (*ordsu)[i]=(*ordsu)[i]+level;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, Ordsu, *ordsv, *ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);
}


// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::ElevateV(int level) const
{
	if (level <0) return *this;
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*(ordv+level);
	
	// find the sum of the ords in w
	int sum3=numw*ordw;
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*(ordv+level);
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).ElevateCPointsV(level);
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv+level; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsv(numv);
	for (int j=0; j<numv; j++) Ordsv[j]=(*ordsv)[j]+level;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, *ordsu, Ordsv, *ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);

}


// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::ElevateW(int level) const
{
	if (level <=0) return *this;
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw*(ordw+level);
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).ElevateCPointsW(level);
				int polstart = k*(ordw+level);
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw+level; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsw(numw);
	for (int k=0; k<numw; k++) Ordsw[k]=(*ordsw)[k]+level;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, *ordsu, *ordsv, Ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);

}


// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::ElevateUV(int levu, int levv) const
{
	// elevate in u then v
	return ElevateU(levu).ElevateV(levv);
}

// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::ElevateUW(int levu, int levw) const
{
	// elevate in u then w
	return ElevateU(levu).ElevateW(levw);
}


// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::ElevateVW(int levv, int levw) const
{
	// elevate in v then w
	return ElevateV(levv).ElevateW(levw);
}

// compute the derivative of the CompPolyVol of order deriv and
// represent the result as another CompPolyVol
template<class T>
CompPolyVol<T> CompPolyVol<T>::DeriveU(int level) const
{
	if (level <=0) return *this;
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*(ordu-level);
	
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw*ordw;

	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*(ordu-level);
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).DeriveCPointsU(level);
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu-level; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsu(numu);
	for (int i=0; i<numu; i++) Ordsu[i]=(*ordsu)[i]-level;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, Ordsu, *ordsv, *ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);
   
}

// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::DeriveV(int level) const
{
	if (level <=0) return *this;
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*(ordv-level);
	
	// find the sum of the ords in w
	int sum3=numw*ordw;
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*(ordv-level);
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).DeriveCPointsV(level);
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv-level; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsv(numv);
	for (int j=0; j<numv; j++) Ordsv[j]=(*ordsv)[j]-level;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, *ordsu, Ordsv, *ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);
}


// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::DeriveW(int level) const
{
	if (level <=0) return *this;
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*ordu;
		// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw*(ordw-level);
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);
	
	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).DeriveCPointsW(level);
				int polstart = k*(ordw-level);
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw-level; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsw(numw);
	for (int k=0; k<numw; k++) Ordsw[k]=(*ordsw)[k]-level;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, *ordsu, *ordsv, Ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);
}



// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
CompPolyVol<T> CompPolyVol<T>::DeriveUV(int levu, int levv) const
{
	// derive in u then v
	return DeriveU(levu).DeriveV(levv);
}



// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
CompPolyVol<T> CompPolyVol<T>::DeriveUW(int levu, int levw) const
{
	// derive in u then w
	return DeriveU(levu).DeriveW(levw);
}

// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
CompPolyVol<T> CompPolyVol<T>::DeriveVW(int levv, int levw) const
{
	// derive in u then v
	return DeriveV(levv).DeriveW(levw);
}


// compute the indefinite integral of the CompPolyVol and represent
// it as a CompPolyVol of one higher degree
template<class T>
CompPolyVol<T> CompPolyVol<T>::IntegrateU() const
{
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*(ordu+1);
	
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw*ordw;
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*(ordu+1);
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).IntegrateCPointsU();
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu+1; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsu(numu);
	for (int i=0; i<numu; i++) Ordsu[i]=(*ordsu)[i]+1;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, Ordsu, *ordsv, *ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);
}	

// integrate the Volace in v
template<class T>
CompPolyVol<T> CompPolyVol<T>::IntegrateV() const
{
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*(ordv+1);
	
	// find the sum of the ords in w
	int sum3=numw*ordw;
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*(ordv+1);
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).DeriveCPointsV(level);
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv+1; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsv(numv);
	for (int j=0; j<numv; j++) Ordsv[j]=(*ordsv)[j]+1;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, *ordsu, Ordsv, *ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);
}


// elevate the degree of the CompPolyVol by level
template<class T>
CompPolyVol<T> CompPolyVol<T>::IntegrateW() const
{
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw*(ordw+1);
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).IntegrateCPointsW();
				int polstart = k*(ordw+1);
				int count3=polstart;
				for (int l=0; l<(*ordsu)[i]; l++) {
					for (int r=0; r<(*ordsv)[j]; r++) {
						for (int q=0; q<(*ordsw)[k]+1; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}

	// construct array of ord
	Vector<int> Ordsw(numw);
	for (int k=0; k<numw; k++) Ordsw[k]=(*ordsw)[k]+1;

	// create and return the CompPolyVol
	return CompPolyVol<T>(numu, numv, numw, *ordsu, *ordsv, Ordsw, ncoeffs, *limitsu, *limitsv, *limitsw);

}


// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
CompPolyVol<T> CompPolyVol<T>::IntegrateUV() const
{
	// integrate in u then v
	return IntegrateU().IntegrateV();
}


// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
CompPolyVol<T> CompPolyVol<T>::IntegrateUW() const
{
	// integrate in u then v
	return IntegrateU().IntegrateW();
}

// evaluate the derivative of the CompPolyVol of order deriv at a point
// x. Computes the derivative as a CompPolyVol and then evaluates this at x
template<class T>
CompPolyVol<T> CompPolyVol<T>::IntegrateVW() const
{
	// integrate in u then v
	return IntegrateV().IntegrateW();
}


// compute the indefinite integral of the CompPolyVol as a CompPolyVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> CompPolyVol<T>::IntegrateCPointsU() const
{
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*(ordu+1);
	
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw*ordw;
	

	Matrix3D<T> ncoeffs(sum1,sum2,sum3);
	
	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*(ordu+1);
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).IntegrateCPointsU();
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu+1; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}
	return ncoeffs;
}	


// compute the indefinite integral of the CompPolyVol as a CompPolyVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> CompPolyVol<T>::IntegrateCPointsV() const
{
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=numu*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*(ordv+1);
	
	// find the sum of the ords in w
	int sum3=numw*ordw;
	
	Matrix3D<T> ncoeffs(sum1,sum2,sum3);
	
	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*(ordv+1);
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).IntegrateCPointsV();
				int polstart = k*ordw;
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv+1; r++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}
	return ncoeffs;
}


// compute the indefinite integral of the CompPolyVol as a CompPolyVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> CompPolyVol<T>::IntegrateCPointsW() const
{
	// take each segment and convert to BezVol form

	// find the sum of the ords in u
	int sum1=num*ordu;
	
	// find the sum of the ords in v
	int sum2=numv*ordv;
	
	// find the sum of the ords in w
	int sum3=numw(*ordw+1);
	


	Matrix3D<T> ncoeffs(sum1,sum2,sum3);

	// convert each patch
	for (int i=0; i<numu; i++) {
		int rowstart = i*ordu;
		int count1=rowstart;
		for (int j=0; j<numv; j++) {
			int colstart = j*ordv;
			int count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> mat = GetVolume(i+1,j+1,k+1).IntegrateCPointsW();
				int polstart = k*(ordw+1);
				int count3=polstart;
				for (int l=0; l<ordu; l++) {
					for (int r=0; r<ordv; r++) {
						for (int q=0; q<ordw+1; q++) {
							ncoeffs[count3][count1][count2]=mat[q][l][r];
							count3++;
						}
						count3=polstart;
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
	}
	return ncoeffs;
}


// find product of two CompPolyVols assuming differing numbers of
// paches in u and v
template<class T>  
CompPolyVol<T> CompPolyVol<T>::Product2(const CompPolyVol<T>& p) const
{
	// convert to CompBezVol and multiply

	CompBezVol<T> b1 = ConvertCompBezVol();
	CompBezVol<T> b2 = p.ConvertCompBezVol();
	CompBezVol<T> prod = b1.Product(b2);
	
	// return composite poly Vol
	return (prod.ConvertCompPolyVol());
}


#endif










