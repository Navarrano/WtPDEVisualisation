#ifndef COMPBEZVOL
#define COMPBEZVOL

#include "BezVol.h"
#include "CompPolyVol.h"
#include "BspVol.h"
#include "mathfunctions.h"

template<class T>
class CompPolyVol;


// can be of different orders, IGES spec

template<class T>
class CompBezVol : public Vol<T> {
private:
	// data
	int numu; // number of segments in u
	int numv; // number of segments in v
	int numw; // number of sgements in w
	int ordu; // order in u
	int ordv; // order in v
	int ordw; // order in w
	Ptr<Vector<int> > ordsu;  // array of orders in u
	Ptr<Vector<int> > ordsv;	// array of orders in v
	Ptr<Vector<int> > ordsw;	// array of orders in w
	Ptr<Matrix3D<T> > cpts;		// control points
	Ptr<Vector<double> > ktsu;	// knots in u
	Ptr<Vector<double> > ktsv;	// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<KnotSet> ksetu;
	Ptr<KnotSet> ksetv;
	Ptr<KnotSet> ksetw;
	Ptr<std::set<double> > limitsetu;
	Ptr<std::set<double> > limitsetv;
	Ptr<std::set<double> > limitsetw;
	

	// private functions
	CompBezVol<T> ElevateU(int levu) const;
	CompBezVol<T> ElevateV(int levv) const;
	CompBezVol<T> ElevateW(int levw) const;
	CompBezVol<T> ElevateUV(int levu, int levv) const;
	CompBezVol<T> ElevateUW(int levu, int levw) const;
	CompBezVol<T> ElevateVW(int levv, int levw) const;
	CompBezVol<T> ElevateUVW(int levu, int levv, int levw) const;
	CompBezVol<T> DeriveU(int levu) const; 
	CompBezVol<T> DeriveV(int levv) const;
	CompBezVol<T> DeriveW(int levw) const;
	CompBezVol<T> DeriveUV(int levu, int levv) const;
	CompBezVol<T> DeriveUW(int levu, int levw) const;
	CompBezVol<T> DeriveVW(int levv, int levw) const;
	T DeriveU(int levu, double u, double v, double w) const;
	T DeriveV(int levv, double u, double v, double w) const;
	T DeriveW(int levw, double u, double v, double w) const;
	T DeriveUV(int levu, int levv, double u, double v, double w) const;
	T DeriveUW(int levu, int levw, double u, double v, double w) const;
	T DeriveVW(int levv, int levw, double u, double v, double w) const;
	CompBezVol<T> IntegrateU() const;
	CompBezVol<T> IntegrateV() const;
	CompBezVol<T> IntegrateW() const;
	CompBezVol<T> IntegrateUV() const;
	CompBezVol<T> IntegrateUW() const;
	CompBezVol<T> IntegrateVW() const;
	
	Matrix3D<T> IntegrateCPointsU() const;
	Matrix3D<T> IntegrateCPointsV() const;
	Matrix3D<T> IntegrateCPointsW() const;
	Matrix3D<T> IntegrateCPointsUV() const;
	Matrix3D<T> IntegrateCPointsUW() const;
	Matrix3D<T> IntegrateCPointsVW() const;
	CompBezVol<T> Product2(const CompBezVol<T>& c) const;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class CompBezVol<double>");
		else {
			std::string s(typeid(T).name()), s1("class CompBezVol<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:

	// constructors
	CompBezVol();
	CompBezVol(int Numu, int Numv, int Numw, int Ordu, int Ordv, int Ordw, const Matrix3D<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw);
	CompBezVol(int Numu, int Numv, int Numw, int Ordu, int Ordv, int Ordw, const Matrix3D<T>& Cpts);
	CompBezVol(int Numu, int Numv, int Numw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Vector<int>& Ordsw, const Matrix3D<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw);
	CompBezVol(int Numu, int Numv, int Numw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Vector<int>& Ordsw, const Matrix3D<T>& Cpts);
	CompBezVol(int Numu, int Numv, int Numw, const Vector<double>& Limitsu, const Vector<double>& Limitsv, const Vector<double>& Limitsw,
	const Vector<int>& Ordsu, const Vector<int>& Ordsv, const Vector<int>& Ordsw, const Matrix3D<T>& Cpts);
	CompBezVol(const Matrix3D<BezVol<T> >& mat, int Numu, int Numv, int Numw, const Vector<double>& Limitsu, const Vector<double>& Limitsv, const Vector<double>& Limitsw);

	// access functions
	int GetOrdU() const;
	int GetOrdV() const;
	int GetOrdW() const;
	int GetNumU() const;
	int GetNumV() const;
	int GetNumW() const;
	int GetNumCPointsU() const;
	int GetNumCPointsV() const;
	int GetNumCPointsW() const;
	int FindSegmentW(double x) const;
	int FindSegmentV(double x) const;
	int FindSegmentU(double x) const;
	Matrix3D<T> GetCPoints() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	Vector<double> GetKnotsW() const;
	Vector<double> GetLimitsU() const;
	Vector<double> GetLimitsV() const;
	Vector<double> GetLimitsW() const;

	virtual double GetLeftLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetLeftLimitW() const;

	virtual double GetRightLimitU() const;
	virtual double GetRightLimitV() const;
	virtual double GetRightLimitW() const;

	Vector<int> GetOrdsU() const;
	Vector<int> GetOrdsV() const;
	Vector<int> GetOrdsW() const;
	BezVol<T> GetVolume(int i, int j, int k) const;
	BezVol<T> GetVolumeOriginal(int i, int j, int k) const;

	// conversion
	Matrix3D<T> Convert() const;
	BspVol<T> ConvertBspVol() const;
	BspVol<T> ConvertBspVolRemoveKnots() const;
	CompPolyVol<T> ConvertCompPolyVol() const;
	
	// evaluators
	virtual T operator()(double u, double v, double w) const;
	virtual T operator()(int, int, int, double u, double v, double w) const;
	T Eval(double u, double v, double w) const;
	Matrix3D<T> ComputePoints(int, int , int) const;

	// degree elevation
	CompBezVol<T> Elevate(int levu, int levv, int levw) const;

	// derivatives
	CompBezVol<T> Derive(int levu, int levv, int levw) const;
	virtual T Derive(int levu, int levv, int levw, double u, double v, double w) const;

	// integration
	CompBezVol<T> Integrate() const;
	T Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<T> IntegrateCPoints() const;

	// product

		template<class T1>
	CompBezVol<T> MakeCompatable(const CompBezVol<T1>& b) const
	{
		return (ConvertBspVol().MakeBreakCompatable(b.ConvertBspVol())).ConvertCompBezVol();
	}

	template<class T1>
	CompBezVol<T> Product(const CompBezVol<T1>& b) const
	{
		CompBezVol<T> d(*this);
		CompBezVol<T1> e(b);

		if (numu != b.GetNumU() || numv != b.GetNumV() || numw != b.GetNumW()) {
			d = MakeCompatable(b);
			e = b.MakeCompatable(d);
		}


		// find number of new control points
		int sum1=d.GetNumU()*(ordu+b.GetOrdU()-1);
	
		int sum2=d.GetNumV()*(ordv+b.GetOrdV()-1);
	
		int sum3=d.GetNumW()*(ordw+b.GetOrdW()-1);
	
		Matrix3D<T> ncpts(sum1,sum2,sum3);

		int count1;
		int count2;
		int count3;

		int rowstart;
		int colstart;
		int polstart;

		for (int i=0; i<d.GetNumU(); i++) {
			rowstart=i*(ordu+b.GetOrdU()-1);
			count1=rowstart;
			for (int j=0; j<d.GetNumV(); j++) {
				colstart=j*(ordv+b.GetOrdV()-1); 
				count2=colstart;
				for (int k=0; k<d.GetNumW(); k++) {
					BezVol<T> bez1 = GetVolume(i+1,j+1,k+1);
					BezVol<T1> bez2 = b.GetVolume(i+1,j+1,k+1);
					Matrix3D<T> temp = bez1.ProductCPoints(bez2);
					polstart = k*(ordw+b.GetOrdW()-1);
					count3=polstart;
					for (int l=0; l<ordu+b.GetOrdU()-1; l++) {
						for (int p=0; p<ordv+b.GetOrdV()-1; p++)  {
							for (int q=0; q<ordw+b.GetOrdW()-1; q++) {
								ncpts[count3][count1][count2] = temp[q][l][p];
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
		// form new Ordsu and Ordsv
		Vector<int> Ordsu(d.GetNumU());
		Vector<int> Ordsv(d.GetNumV());
		Vector<int> Ordsw(d.GetNumW());

		for (int i=0; i<d.GetNumU(); i++) Ordsu[i] = (*ordsu)[i]+b.GetOrdsU()[i]-1;
		for (int j=0; j<d.GetNumV(); j++) Ordsv[j] = (*ordsv)[j]+b.GetOrdsV()[j]-1;
		for (int k=0; k<d.GetNumW(); k++) Ordsw[k] = (*ordsw)[k]+b.GetOrdsW()[k]-1;

	
		// create and return the CompBezVol
		return CompBezVol<T>(d.GetNumU(),d.GetNumV(),d.GetNumW(),GetLimitsU(),GetLimitsV(),GetLimitsW(),Ordsu,Ordsv,Ordsw,ncpts);
	}
	
	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
};

// CONSTRUCTORS

// default constructor
template<class T>
CompBezVol<T>::CompBezVol() : numu(0), numv(0), numw(0), ordu(0), ordv(0), ordw(0), 
ordsu(), ordsv(), ordsw(), cpts(), ktsu(), ktsv(), ktsw()
{
}

// constructor builds a CompBezVol from a matrix of control points, knots
// an order and number of control points in u and v
template<class T>
CompBezVol<T>::CompBezVol(int Numu, int Numv, int Numw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, 
							const Vector<int>& Ordsw, const Matrix3D<T>& Cpts, 
							const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw) :
numu(Numu), numv(Numv), numw(Numw), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), 
ordsw(new Vector<int>(Ordsw)), cpts(new Matrix3D<T>(Cpts)), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)), ktsw(new Vector<double>(Ktsw)),
ksetu(new KnotSet(*ktsu,ordu,ordu*(numu+1))), ksetv(new KnotSet(*ktsv,ordv,ordv*(numv+1))), ksetw(new KnotSet(*ktsw,ordw,ordw*(numw+1))),
limitsetu(new std::set<double>(Ktsu.begin(),Ktsu.end())), limitsetv(new std::set<double>(Ktsv.begin(),Ktsv.end())),
limitsetw(new std::set<double>(Ktsw.begin(),Ktsw.end()))
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
	cpts = new Matrix3D<T>(Convert());
}

// constructor builds a CompBezVol from a matrix of control points
// an order and number of control points in u and v
template<class T>
CompBezVol<T>::CompBezVol(int Numu, int Numv, int Numw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, 
		const Vector<int>& Ordsw, const Matrix3D<T>& Cpts) :
	numu(Numu), numv(Numv), numw(Numw), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), ordsw(new Vector<int>(Ordsw)), 
		cpts(new Matrix3D<T>(Cpts))
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
	// create knots in u and v, w

	ktsu = new Vector<double>(Math::CreateKnots(Numu,ordu));
	ktsv = new Vector<double>(Math::CreateKnots(Numv,ordv));
	ktsw = new Vector<double>(Math::CreateKnots(Numw,ordw));
	ksetu = new KnotSet(*ktsu,ordu,(ordu+1)*numu);
	ksetv = new KnotSet(*ktsv,ordv,(ordv+1)*numv);
	ksetw = new KnotSet(*ktsw,ordw,(ordw+1)*numw);
	cpts = new Matrix3D<T>(Convert());	
	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsu).end());
	limitsetv = new std::set<double>((*ktsv).begin(),(*ktsv).end());
	limitsetw = new std::set<double>((*ktsw).begin(),(*ktsw).end());
}


template<class T>
CompBezVol<T>::CompBezVol(int Numu, int Numv, int Numw, const Vector<double>& Limitsu, const Vector<double>& Limitsv,
		const Vector<double>& Limitsw, const Vector<int>& Ordsu, const Vector<int>& Ordsv, 
		const Vector<int>& Ordsw, const Matrix3D<T>& Cpts) :
	numu(Numu), numv(Numv), numw(Numw), ordsu(new Vector<int>(Ordsu)), ordsv(new Vector<int>(Ordsv)), ordsw(new Vector<int>(Ordsw)), 
		cpts(new Matrix3D<T>(Cpts))
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
	// create knots in u and v, w

	ktsu = new Vector<double>(Math::CreateKnots(Numu,ordu,Limitsu));
	ktsv = new Vector<double>(Math::CreateKnots(Numv,ordv,Limitsv));
	ktsw = new Vector<double>(Math::CreateKnots(Numw,ordw,Limitsw));
	ksetu = new KnotSet(*ktsu,ordu,(ordu+1)*numu);
	ksetv = new KnotSet(*ktsv,ordv,(ordv+1)*numv);
	ksetw = new KnotSet(*ktsw,ordw,(ordw+1)*numw);
	cpts = new Matrix3D<T>(Convert());	
	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsu).end());
	limitsetv = new std::set<double>((*ktsv).begin(),(*ktsv).end());
	limitsetw = new std::set<double>((*ktsw).begin(),(*ktsw).end());
}



	
// constructor builds a CompBezVol from a matrix of control points,
// knot vectors an order and a number of control points in u and v
template<class T>
CompBezVol<T>::CompBezVol(int Numu, int Numv, int Numw, int Ordu, int Ordv, int Ordw, 
	const Matrix3D<T>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw) :
numu(Numu), numv(Numv), numw(Numw), ordu(Ordu), ordv(Ordv), ordw(Ordw), ordsu(new Vector<int>(numu)), ordsv(new Vector<int>(numv)), 
ordsw(new Vector<int>(numw)), cpts(new Matrix3D<T>(Cpts)), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)), ktsw(new Vector<double>(Ktsw))
{
	// assign order arrays in u and v 
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;
	for (int k=0; k<numw; k++) (*ordsw)[k]=ordw;
	ksetu = new KnotSet(*ktsu,ordu,(ordu+1)*numu);
	ksetv = new KnotSet(*ktsv,ordv,(ordv+1)*numv);
	ksetw = new KnotSet(*ktsw,ordw,(ordw+1)*numw);
	
	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsu).end());
	limitsetv = new std::set<double>((*ktsv).begin(),(*ktsv).end());
	limitsetw = new std::set<double>((*ktsw).begin(),(*ktsw).end());
}

// constructor builds a CompBezVol from a matrix of control points
// orders in u and v
template<class T>
CompBezVol<T>::CompBezVol(int Numu, int Numv, int Numw, int Ordu, int Ordv, int Ordw, const Matrix3D<T>& Cpts) :
numu(Numu), numv(Numv), numw(Numw), ordu(Ordu), ordv(Ordv), ordw(Ordw), 
ordsu(new Vector<int>(Numu)), ordsv(new Vector<int>(Numv)), ordsw(new Vector<int>(numw)), 
cpts(new Matrix3D<T>(Cpts))
{
	// set orders in u
	for (int i=0; i<numu; i++) (*ordsu)[i]=ordu;
	// set orders in v
	for (int j=0; j<numv; j++) (*ordsv)[j]=ordv;
	// create knots in u and v
	for (int k=0; k<numw; k++) (*ordsw)[k]=ordw;
	// create knots in u and v, w
	ktsu = new Vector<double>(Math::CreateKnots(Numu,ordu));
	ktsv = new Vector<double>(Math::CreateKnots(Numv,ordv));
	ktsw = new Vector<double>(Math::CreateKnots(Numw,ordw));
	ksetu = new KnotSet(*ktsu,ordu,(ordu+1)*numu);
	ksetv = new KnotSet(*ktsv,ordv,(ordv+1)*numv);
	ksetw = new KnotSet(*ktsw,ordw,(ordw+1)*numw);
	
	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsu).end());
	limitsetv = new std::set<double>((*ktsv).begin(),(*ktsv).end());
	limitsetw = new std::set<double>((*ktsw).begin(),(*ktsw).end());
}


template<class T>
CompBezVol<T>::CompBezVol(const Matrix3D<BezVol<T> >& mat, int Numu, int Numv, int Numw, const Vector<double>& Limitsu, const Vector<double>& Limitsv, const Vector<double>& Limitsw) :
numu(Numu), numv(Numv), numw(Numw), ordsu(new Vector<int>(Numu)), ordsv(new Vector<int>(Numv)), ordsw(new Vector<int>(Numw))
{
//	std::ofstream ofs("bez.dat");

//	ofs << mat;
//	ofs.close();
	int sum1 = 0;
	int sum2 = 0;
	int sum3 = 0;
	for (int i=0; i<numu; i++) (*ordsu)[i] = mat[0][i][0].GetOrdU();
	for (int j=0; j<numv; j++) (*ordsv)[j] = mat[0][0][j].GetOrdV();
	for (int k=0; k<numw; k++) (*ordsw)[k] = mat[k][0][0].GetOrdW();

	ordu = (*ordsu)[0];
	for (int i=1; i<numu; i++) if (ordu < (*ordsu)[i]) ordu=(*ordsu)[i];
	// find max order in v
	ordv = (*ordsv)[0];
	for (int j=1; j<numv; j++) if (ordv < (*ordsv)[j]) ordv=(*ordsv)[j];

	// find max order in w
	ordw = (*ordsw)[0];
	for (int k=1; k<numw; k++) if (ordw < (*ordsw)[k]) ordw=(*ordsw)[k];

	for (int i=0; i<numu; i++) sum1 = sum1+ (*ordsu)[i];
	for (int j=0; j<numv; j++) sum2 = sum2+ (*ordsv)[j];
	for (int k=0; k<numw; k++) sum3 = sum3+ (*ordsw)[k];

	cpts = new Matrix3D<T>(sum1,sum2,sum3);
	int count1 = 0;
	int count2 = 0;
	int count3 = 0;
	int rowstart, colstart, polstart;
	for (int i=0; i<numu; i++) {
		rowstart = i*(*ordsu)[i];
		count1 = rowstart;
		for (int j=0; j<numv; j++) {
			colstart = j*(*ordsv)[j];
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				polstart = k*(*ordsw)[k];
				count3 = polstart;

				for (int l=0; l<(*ordsu)[i]; l++) {
					for (int m=0; m<(*ordsv)[j]; m++) {
						for (int n=0; n<(*ordsw)[k]; n++) {
							(*cpts)[count3][count1][count2] = mat[k][i][j].GetCPoints()[n][l][m];
							count3++;
						}
						count3 = polstart;
						count2++;
					}
					count2 = colstart;
					count1++;
				}
				count1 = rowstart;
			}
		}
	}
	
	// create knots in u and v
	ktsu = new Vector<double>(Math::CreateKnots(Numu,ordu,Limitsu));
	ktsv = new Vector<double>(Math::CreateKnots(Numv,ordv,Limitsv));
	ktsw = new Vector<double>(Math::CreateKnots(Numw,ordw,Limitsw));
	ksetu = new KnotSet(*ktsu,ordu,(ordu+1)*numu);
	ksetv = new KnotSet(*ktsv,ordv,(ordv+1)*numv);
	ksetw = new KnotSet(*ktsw,ordw,(ordw+1)*numw);
	cpts = new Matrix3D<T>(Convert());
	limitsetu = new std::set<double>((*ktsu).begin(),(*ktsu).end());
	limitsetv = new std::set<double>((*ktsv).begin(),(*ktsv).end());
	limitsetw = new std::set<double>((*ktsw).begin(),(*ktsw).end());
}


// ACESS FUNCTIONS

template<class T>
int CompBezVol<T>::GetNumCPointsU() const
{
	return numu*ordu;
}

template<class T>
int CompBezVol<T>::GetNumCPointsV() const
{
	return numv*ordv;
}

template<class T>
int CompBezVol<T>::GetNumCPointsW() const
{
	return numw*ordw;
}


// get the order of the CompBezVol
template<class T>
inline int CompBezVol<T>::GetOrdU() const { return ordu; }

// get the order of the CompBezVol
template<class T>
inline int CompBezVol<T>::GetOrdV() const { return ordv; }

// get the order of the CompBezVol
template<class T>
inline int CompBezVol<T>::GetOrdW() const { return ordw; }


// get the orders in u
template<class T>
inline Vector<int> CompBezVol<T>::GetOrdsU() const { return *ordsu; }

// get the orders in v
template<class T>
inline Vector<int> CompBezVol<T>::GetOrdsV() const { return *ordsv; }

// get the orders in v
template<class T>
inline Vector<int> CompBezVol<T>::GetOrdsW() const { return *ordsw; }

// get the number of patches in u
template<class T>
inline int CompBezVol<T>::GetNumU() const { return numu; }

// get the number of patches in v
template<class T>
inline int CompBezVol<T>::GetNumV() const { return numv; }

// get the number of patches in v
template<class T>
inline int CompBezVol<T>::GetNumW() const { return numw; }


// get the control points
template<class T>
inline Matrix3D<T> CompBezVol<T>::GetCPoints() const { return *cpts; }

// get the knot vector in u
template<class T>
inline Vector<double> CompBezVol<T>::GetKnotsU() const { return *ktsu; }


// get the knot vector in v
template<class T>
inline Vector<double> CompBezVol<T>::GetKnotsV() const { return *ktsv; }

// get the knot vector in v
template<class T>
inline Vector<double> CompBezVol<T>::GetKnotsW() const { return *ktsw; }


// get the knot vector in v
template<class T>
inline Vector<double> CompBezVol<T>::GetLimitsU() const { return (*ksetu).GetDistinctKnots(); }

// get the knot vector in v
template<class T>
inline Vector<double> CompBezVol<T>::GetLimitsV() const { return (*ksetv).GetDistinctKnots(); }

// get the knot vector in v
template<class T>
inline Vector<double> CompBezVol<T>::GetLimitsW() const { return (*ksetw).GetDistinctKnots(); }

// get the knot vector
template<class T>
double CompBezVol<T>::GetLeftLimitU() const { return (*ktsu)[0]; }


// get the knot vector
template<class T>
double CompBezVol<T>::GetLeftLimitV() const { return (*ktsv)[0]; }

// get the knot vector
template<class T>
double CompBezVol<T>::GetLeftLimitW() const { return (*ktsw)[0]; }


// get the knot vector
template<class T>
double CompBezVol<T>::GetRightLimitU() const { return (*ktsu)[numu*ordu]; }


// get the knot vector
template<class T>
double CompBezVol<T>::GetRightLimitV() const { return (*ktsv)[numv*ordv]; }

// get the knot vector
template<class T>
double CompBezVol<T>::GetRightLimitW() const { return (*ktsw)[numw*ordw]; }


template<class T>
int CompBezVol<T>::FindSegmentU(double x) const
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
int CompBezVol<T>::FindSegmentV(double x) const
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
int CompBezVol<T>::FindSegmentW(double x) const
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

// get the i,jth patch of the Volace and return as a BezVol
template<class T>
BezVol<T> CompBezVol<T>::GetVolumeOriginal(int indu, int indv, int indw) const
{
	// order
	int ordU = (*ordsu)[indu-1];
	int ordV = (*ordsv)[indv-1];
	int ordW = (*ordsw)[indw-1];

	// find start index for control points in u and v
	int sum1=0;
	for (int l=0; l<indu; l++) sum1 = sum1+(*ordsu)[l];

	int sum2=0;
	for (int l=0; l<indv; l++) sum2 = sum2+(*ordsv)[l];
	
	int sum3=0;
	for (int l=0; l<indw; l++) sum3 = sum3+(*ordsw)[l];

	// extract control points
	Matrix3D<T> ncpts(ordU,ordV,ordW);
	for (int l=0; l<ordU; l++) 
		for (int p=0; p<ordV; p++)
			for (int q=0; q<ordW; q++)
				ncpts[q][l][p] = (*cpts)[sum3-ordW+q][sum1-ordU+l][sum2-ordV+p];

	
	// create and return the BezVol
	return BezVol<T>(ncpts, ordU, ordV, ordW,(*ksetu).GetDistinctKnots()[indu-1],
		(*ksetu).GetDistinctKnots()[indu],(*ksetv).GetDistinctKnots()[indv-1],(*ksetv).GetDistinctKnots()[indv],
		(*ksetw).GetDistinctKnots()[indw-1],(*ksetw).GetDistinctKnots()[indw]);
}

// get the i,jth patch of the Volace and return as a BezVol
template<class T>
BezVol<T> CompBezVol<T>::GetVolume(int indu, int indv, int indw) const
{
	// extract control points
	int sum1 = indu*ordu;
	int sum2 = indv*ordv;
	int sum3 = indw*ordw;

	Matrix3D<T> ncpts(ordu,ordv,ordw);
	for (int l=0; l<ordu; l++) 
		for (int p=0; p<ordv; p++)
			for (int q=0; q<ordw; q++)
				ncpts[q][l][p] = (*cpts)[sum3-ordw+q][sum1-ordu+l][sum2-ordv+p];

	// construct knots
	
	// create and return the BezVol
	return BezVol<T>(ncpts, ordu, ordv, ordw,(*ksetu).GetDistinctKnots()[indu-1],
		(*ksetu).GetDistinctKnots()[indu],(*ksetv).GetDistinctKnots()[indv-1],(*ksetv).GetDistinctKnots()[indv],
		(*ksetw).GetDistinctKnots()[indw-1],(*ksetw).GetDistinctKnots()[indw]);
}


// CONVERSION

// convert to BspVol form
template<class T>
BspVol<T> CompBezVol<T>::ConvertBspVol() const
{
	// create and return the bspVol
	return BspVol<T>(*cpts,*ktsu,*ktsv,*ktsw,ordu,ordv,ordw,ordu*numu,ordv*numv,ordw*numw);
}

// convert to BspVol removing all possible knots in u and v
template<class T>
BspVol<T> CompBezVol<T>::ConvertBspVolRemoveKnots() const
{
	// create and return BspVol removing all possible knots in u and v
	return BspVol<T>(*cpts,*ktsu,*ktsv,*ktsw,ordu,ordv,ordw,ordu*numu,ordv*numv,ordw*numw).RemovePossibleKnotsUVW();
}


// convert each BezVol patch to the max degree in u and v
template<class T>
Matrix3D<T> CompBezVol<T>::Convert() const
{
	// elevate the degree of each patch to max in u and v
	// extract BezVol patch
	// build array of control points
	int count1;
	int count2;
	int count3;
	int rowstart;
	int colstart;
	int polstart;


	Matrix3D<T> ncpts(numu*ordu,numv*ordv,numw*ordw);
	// elevate each patch to max order in u and v, w
	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				BezVol<T> b = GetVolumeOriginal(i+1,j+1,k+1).ElevateUVW(ordu-(*ordsu)[i],ordv-(*ordsv)[j],ordw-(*ordsw)[k]);	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=b.GetCPoints()[q][l][p];
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

	// create and return the new CompBezVol
	return ncpts;
}


// convert to a CompPolyVol
template<class T>
CompPolyVol<T> CompBezVol<T>::ConvertCompPolyVol() const
{
	// take each segment and convert to BezVol form

	// find the sum of the ords
	int sum1=numu*ordu;
	int sum2=numv*ordv;
	int sum3=numw*ordw;

	Matrix3D<T> ncoeffs(sum1,sum2,sum3);
	int count1;
	int count2;
	int count3;
	int rowstart;
	int colstart;
	int polstart;

	// elevate each patch to max order in u and v, w
	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				PolyVol<T> p1 = GetVolume(i+1,j+1,k+1).ConvertPolyVol();	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw; q++) {
							ncoeffs[count3][count1][count2]=p1.GetCoeffs()[q][l][p];
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

	// create and return the new CompBezVol
	return CompPolyVol<T>(numu, numv, numw, ,GetLimitsU(),GetLimitsV(),GetLimitsW(), *ordsu, *ordsv, *ordsw, ncoeffs);
}


// EVALUATORS

// evaluate the CompBezVol at the point x using de Boor algorithm
template<class T>
T CompBezVol<T>::operator()(double u, double v, double w) const
{
	// find the segment in u and v
	int i = FindSegmentU(u);
	int j = FindSegmentV(v);
	int k = FindSegmentW(w);

	// get the patch and evaluate
	return GetVolume(i,j,k)(u,v,w);
}


// evaluate the CompBezVol using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompBezVol<T>::Eval(double u, double v, double w) const
{
	// find the segment in u and v
	int i = FindSegmentU(u);
	int j = FindSegmentV(v);
	int k = FindSegmentW(w);

	// get the patch and evaluate
	return GetVolume(i,j,k).Eval(u,v,w);
}


// DEGREE ELEVATION


// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::Elevate(int levu, int levv, int levw) const
{
	// elevate in U and V
	return (ElevateU(levu).ElevateV(levv)).ElevateW(levw);
}

// DERIVATIVES


// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
CompBezVol<T> CompBezVol<T>::Derive(int levu, int levv, int levw) const
{
	// derive u then v
	return (DeriveU(levu).DeriveV(levv)).DeriveW(levw);
}

// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
T CompBezVol<T>::Derive(int levu, int levv, int levw, double u, double v, double w) const
{
	// derive u then v
	return (DeriveU(levu).DeriveV(levv)).DeriveW(levw)(u,v,w);
}

template<class T>
T CompBezVol<T>::operator() (int valu, int valv, int valw, double u, double v, double w) const
{
	return Derive(valu,valv,valw,u,v,w);
}



// INTEGRATE

// integrate the CompBezVol between the limits x1 and x2. Computes
// the indefinite integral as a CompBezVol and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T CompBezVol<T>::Integrate(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	// create the indefinite integral;
	return ConvertBspVol().Integrate(u1,u2,v1,v2,w1,w2);
}


// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
CompBezVol<T> CompBezVol<T>::Integrate() const
{
	return IntegrateU().IntegrateV().IntegrateW();
}

/*
template<class T>
CompBezVol<T> CompBezVol<T>::MakeCompatable(const CompBezVol<T>& b) const
{
	return (ConvertBspVol().MakeBreakCompatable(b.ConvertBspVol())).ConvertCompBezVol();
}		
*/

/*
// PRODUCT

// compute the product of the CompBezVol with another CompBezVol and 
// represent the result as a new CompBezVol 
// must be compatible, i.e have the same number of patches in u and v
template<class T>  
CompBezVol<T> CompBezVol<T>::Product(const CompBezVol<T>& b) const
{
	CompBezVol<T> d(*this), e(b);

	if (numu != b.GetNumU() || numv != b.GetNumV() || numw != b.GetNumW()) {
		d = MakeCompatable(b);
		e = b.MakeCompatable(d);
	}


	// find number of new control points
	int sum1=d.GetNumU()*(ordu+b.GetOrdU()-1);
	
	int sum2=d.GetNumV()*(ordv+b.GetOrdV()-1);
	
	int sum3=d.GetNumW()*(ordw+b.GetOrdW()-1);
	
	Matrix3D<T> ncpts(sum1,sum2,sum3);

	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	for (int i=0; i<d.GetNumU(); i++) {
		rowstart=i*(ordu+b.GetOrdU()-1);
		count1=rowstart;
		for (int j=0; j<d.GetNumV(); j++) {
			colstart=j*(ordv+b.GetOrdV()-1); 
			count2=colstart;
			for (int k=0; k<d.GetNumW(); k++) {
				BezVol<T> bez1 = GetVolume(i+1,j+1,k+1);
				BezVol<T> bez2 = b.GetVolume(i+1,j+1,k+1);
				Matrix3D<T> temp = bez1.ProductCPoints(bez2);
				polstart = k*(ordw+b.GetOrdW()-1);
				count3=polstart;
				for (int l=0; l<ordu+b.GetOrdU()-1; l++) {
					for (int p=0; p<ordv+b.GetOrdV()-1; p++)  {
						for (int q=0; q<ordw+b.GetOrdW()-1; q++) {
							ncpts[count3][count1][count2] = temp[q][l][p];
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
	// form new Ordsu and Ordsv
	Vector<int> Ordsu(d.GetNumU());
	Vector<int> Ordsv(d.GetNumV());
	Vector<int> Ordsw(d.GetNumW());

	for (int i=0; i<d.GetNumU(); i++) Ordsu[i] = (*ordsu)[i]+b.GetOrdsU()[i]-1;
	for (int j=0; j<d.GetNumV(); j++) Ordsv[j] = (*ordsv)[j]+b.GetOrdsV()[j]-1;
	for (int k=0; k<d.GetNumW(); k++) Ordsw[k] = (*ordsw)[k]+b.GetOrdsW()[k]-1;

	
	// create and return the CompBezVol
	return CompBezVol<T>(d.GetNumU(),d.GetNumV(),d.GetNumW(),GetLimitsU(),GetLimitsV(),GetLimitsW(),Ordsu,Ordsv,Ordsw,ncpts);
}
*/

// READ & WRITE

template <class T>
void CompBezVol<T>::write(std::ostream& os) 
{
	os << "Comp Bezier Volume\n";
	os << "number of segments in u, v, w\n";
	os << numu << " " << numv << " " << numw;
	os << "orders in u\n";
	os << ordsu;
	os << "\norders in v\n";
	os << ordsv;
	os << "\norders in w\n";
	os << ordsw;
	os << "knots in u are\n";
	os << ktsu;
	os << "\nknots in v are \n";
	os << ktsv;
	
	os << "\nknots in w are \n";
	os << ktsw;
	
	os << "\ncontrol points are\n";
	os << cpts;
}

template <class T>
void CompBezVol<T>::read(std::istream& is)
{
	int Numu, Numv, Numw;
	std::cout << "number of segments in u,v,w\n";
	is >> Numu >> Numv >> Numw;
	Vector<int> Ordsu(Numu), Ordsv(Numv), Ordsw(Numw);
	std::cout << "order of Bezier Volume in u and v and w";
	is >> Ordsu >> Ordsv >> Ordsw;
	int sum1=0;
	int sum2=0;
	int sum3=0;
	for (int i=0; i<Numu; i++) sum1+=(*ordsu)[i];
	for (int j=0; j<Numv; j++) sum2+=(*ordsv)[j];
	for (int k=0; k<Numw; k++) sum3+=(*ordsw)[k];
	Matrix3D<T> Cpts(sum1,sum2,sum3);
	std::cout << "\ninput control points\n";
	is >> Cpts;
	*this = CompBezVol<T>(Numu,Numv,Numw,Ordsu,Ordsv,Ordsw,Cpts);
} 


template <class T>
void CompBezVol<T>::writefile(std::ofstream& ofs) 
{
	ofs << "Comp Bezier Volume\n";
	ofs << "number of segments in u, v, w\n";
//	ofs << numu << " " << numv << " " << numw;
	ofs << "orders in u\n";
	ofs << ordsu;
	ofs << "\norders in v\n";
	ofs << ordsv;
	ofs << "\norders in w\n";
	ofs << ordsw;
	ofs << "knots in u are\n";
	ofs << ktsu;
	ofs << "\nknots in v are \n";
	ofs << ktsv;
	
	ofs << "\nknots in w are \n";
	ofs << ktsw;
	ofs << "\ncontrol points\n";
	ofs << cpts;
}


template <class T>
void CompBezVol<T>::readfile(std::ifstream& ifs)
{
	
	int Numu, Numv, Numw;
	ifs >> Numu >> Numv >> Numw;
	
	Vector<int> Ordsu(Numu), Ordsv(Numv), Ordsw(Numw);
	ifs >> Ordsu >> Ordsv >> Ordsw;
	int sum1=0;
	int sum2=0;
	int sum3=0;
	for (int i=0; i<Numu; i++) sum1+=(*ordsu)[i];
	for (int j=0; j<Numv; j++) sum2+=(*ordsv)[j];
	for (int k=0; k<Numw; k++) sum3+=(*ordsw)[k];
	
	Matrix3D<T> Cpts(sum1,sum2,sum3);

	ifs >> Cpts;
	
	*this = CompBezVol<T>(Numu,Numv,Numw,Ordsu,Ordsv,Ordsw,Cpts);
} 


// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::ElevateU(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*(ordu+level);
	
	int sum2=numv*ordv
	
	int sum3=numw*ordw;

	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu+level);
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).ElevateCPointsU(level);	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu+level; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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
	for (int i=0; i<numu; i++) Ordsu[i]=(*ordsu)[i]+level;

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, , GetLimitsU(), GetLimitsV(), GetLimitsW(), Ordsu, *ordsv, *ordsw, ncpts);
}


// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::ElevateV(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv+level);
	
	int sum3=numw*ordw;
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*(ordv+level);
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).ElevateCPointsU(level);	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv+level; p++)  {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, GetLimitsU(), GetLimitsV(), GetLimitsW(), *ordsu, Ordsv, *ordsw, ncpts);
}


// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::ElevateW(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*ordv;
	
	int sum3=numw*(ordw+level);
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).ElevateCPointsW(level);	
				polstart = k*(ordw+level);
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw+level; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, GetLimitsU(), GetLimitsV(), GetLimitsW(), *ordsu, *ordsv, Ordsw, ncpts);
}


// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::ElevateUV(int levu, int levv) const
{
	// elevate in U and V
	return ElevateU(levu).ElevateV(levv);
}

// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::ElevateUW(int levu, int levw) const
{
	// elevate in U and V
	return ElevateU(levu).ElevateW(levw);
}


// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::ElevateVW(int levv, int levw) const
{
	// elevate in U and V
	return ElevateV(levv).ElevateW(levw);
}

// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::ElevateUVW(int levu, int levv, int levw) const
{
	// elevate in U and V
	return ElevateU(levu).ElevateV(levv).ElevateW(levw);
}

// compute the derivative of the CompBezVol of order deriv and
// represent the result as another CompBezVol
template<class T>
CompBezVol<T> CompBezVol<T>::DeriveU(int level) const
{   
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*(ordu-level);
	
	int sum2=numv*ordv;
	
	int sum3=numw*ordw;
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu-level);
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).DeriveCPointsU(level);	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu-level; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, GetLimitsU(), GetLimitsV(), GetLimitsW(),Ordsu, *ordsv, *ordsw, ncpts);
}



// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::DeriveV(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv-level);
	
	int sum3=numw*ordw;
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*(ordv-level);
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).DeriveCPointsV(level);	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv-level; p++) {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, GetLimitsU(), GetLimitsV(), GetLimitsW(), *ordsu, Ordsv, *ordsw, ncpts);
}

// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::DeriveW(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*ordv;
	
	int sum3=numw*(ordw-level);
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).DeriveCPointsW(level);	
				polstart = k*(ordw-level);
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw-level; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, GetLimitsU(), GetLimitsV(), GetLimitsW(), *ordsu, *ordsv, Ordsw, ncpts);
}

// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
CompBezVol<T> CompBezVol<T>::DeriveUV(int levu, int levv) const
{
	// derive u then v
	return DeriveU(levu).DeriveV(levv);
}

// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
CompBezVol<T> CompBezVol<T>::DeriveUW(int levu, int levw) const
{
	// derive u then v
	return DeriveU(levu).DeriveW(levw);
}


// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
CompBezVol<T> CompBezVol<T>::DeriveVW(int levv, int levw) const
{
	// derive u then v
	return DeriveV(levv).DeriveW(levw);
}


// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
T CompBezVol<T>::DeriveUV(int levu, int levv, double u, double v, double w) const
{
	// derive u then v and evaluate
	return DeriveUV(levu,levv)(u,v,w);
}


// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
T CompBezVol<T>::DeriveUW(int levu, int levw, double u, double v, double w) const
{
	// derive u then v and evaluate
	return DeriveUW(levu,levw)(u,v,w);
}


// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
T CompBezVol<T>::DeriveVW(int levv, int levw, double u, double v, double w) const
{
	// derive u then v and evaluate
	return DeriveVW(levv,levw)(u,v,w);
}

// compute the indefinite integral of the CompBezVol and represent
// it as a CompBezVol of one higher degree
template<class T>
CompBezVol<T> CompBezVol<T>::IntegrateU() const
{
	// find number of new control points
	int sum1=numu*(ordu+1);
	
	int sum2=numv*ordv;
	
	int sum3=numw*ordw;
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu+1);
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).IntegrateCPointsU();	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu+1; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, Ordsu, GetLimitsU(), GetLimitsV(), GetLimitsW(),*ordsv, *ordsw, ncpts);
}	

// elevate the degree of the CompBezVol by level
template<class T>
CompBezVol<T> CompBezVol<T>::IntegrateV() const
{
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv+1);
	
	int sum3=numw*ordw;
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*(ordv+1);
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).IntegrateCPointsV();	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv+1; p++) {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, GetLimitsU(), GetLimitsV(), GetLimitsW(), *ordsu, Ordsv, *ordsw, ncpts);
}


template<class T>
CompBezVol<T> CompBezVol<T>::IntegrateW() const
{
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*ordv;
	
	int sum3=numw*(ordw+1);
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).IntegrateCPointsW();	
				polstart = k*(ordw+1);
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw+1; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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
	for (int k=0; k<numw; k++) Ordsw[k]=(*ordsw)[k]+1;

	// create and return the CompBezVol
	return CompBezVol<T>(numu, numv, numw, GetLimitsU(), GetLimitsV(), GetLimitsW(), *ordsu, *ordsv, Ordsw, ncpts);
}


// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
CompBezVol<T> CompBezVol<T>::IntegrateUV() const
{
	return IntegrateU().IntegrateV();
}

// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
CompBezVol<T> CompBezVol<T>::IntegrateUW() const
{
	return IntegrateU().IntegrateW();
}

// evaluate the derivative of the CompBezVol of order deriv at a point
// x. Computes the derivative as a CompBezVol and then evaluates this at x
template<class T>
CompBezVol<T> CompBezVol<T>::IntegrateVW() const
{
	return IntegrateV().IntegrateW();
}


// compute the indefinite integral of the CompBezVol as a CompBezVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> CompBezVol<T>::IntegrateCPointsU() const
{
	// find number of new control points
	int sum1=numu*(ordu+1);
	
	int sum2=numv*ordv;
	
	int sum3=numw*ordw;
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*(ordu+1);
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).IntegrateCPointsU();	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu+1; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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
	return ncpts;
}	


// compute the indefinite integral of the CompBezVol as a CompBezVol
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Matrix3D<T> CompBezVol<T>::IntegrateCPointsV() const
{

	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*(ordv+1);
	
	int sum3=numw*ordw;
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*(ordv+1);
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).IntegrateCPointsV();	
				polstart = k*ordw;
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (p=0; p<ordv+1; p++) {
						for (q=0; q<ordw; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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
	return ncpts;
}	


template<class T>
Matrix3D<T> CompBezVol<T>::IntegrateCPointsW() const
{
	// find number of new control points
	int sum1=numu*ordu;
	
	int sum2=numv*ordv;
	
	int sum3=numw*(ordw+1);
	
	
	int count1;
	int count2;
	int count3;

	int rowstart;
	int colstart;
	int polstart;

	Matrix3D<T> ncpts(sum1,sum2,sum3);

	for (int i=0; i<numu; i++) {
		rowstart = i*ordu;
		count1 = rowstart;
		for (int j=0; j<numv; j++) { 
			colstart = j*ordv;
			count2 = colstart;
			for (int k=0; k<numw; k++) {
				Matrix3D<T> v = GetVolume(i+1,j+1,k+1).IntegrateCPointsW();	
				polstart = k*(ordw+1);
				count3 = polstart;
				for (int l=0; l<ordu; l++) {
					for (int p=0; p<ordv; p++) {
						for (int q=0; q<ordw+1; q++) {
							ncpts[count3][count1][count2]=v[q][l][p];
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
	return ncpts;
}

// compute the product of the CompBezVol with another CompBezVol and 
// represent the result as a new CompBezVol 
// need not be compatable
template<class T>  
CompBezVol<T> CompBezVol<T>::Product2(const CompBezVol<T>& b) const
{
	// convert to CompBspVol and multiply

	BspVol<T> b1 = ConvertBspVol();
	BspVol<T> b2 = b.ConvertBspVol();
	BspVol<T> prod = b1.Product(b2);
	// return composite poly Vol
	return (prod.ConvertCompBezVol());
}


#endif