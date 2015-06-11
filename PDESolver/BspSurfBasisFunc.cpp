
#include "bspsurf.h"
//#include "compbezsurf.h"

// CONSTRUCTORS

// default constructor
BspSurfBasisFunc::BspSurfBasisFunc() : ordu(0), ordv(0), ktsu(), ktsv() {}

// constructor taking orders amd knot sets in u and v
BspSurfBasisFunc::BspSurfBasisFunc(int Ordu, int Ordv, const Vector<double>& Ktsu,const Vector<double>& Ktsv) : 
	ordu(Ordu), ordv(Ordv), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)), b(new BspSurf<double>(CreateBspSurf())) { }


// ACCESS FUNCTIONS

// compute the dimension of the u knot set
int BspSurfBasisFunc::ComputeDimU() const
{
	// create KnotSet object 
	KnotSet kset(*ktsu,ordu,ordu+1);
	int num = kset.GetNumDistinct();
	int num_cond=0;
	// find the dimension
	for (int i=1; i<num-1; i++) num_cond = num_cond + ordu - kset.GetMult()[i];
	return ((num-1)*ordu-num_cond);
}

// compute the dimension of the v knot set
int BspSurfBasisFunc::ComputeDimV() const
{
	// create KnotSet object
	KnotSet kset(*ktsv,ordv,ordv+1);
	int num = kset.GetNumDistinct();
	int num_cond=0;
	// find the dimension
	for (int i=1; i<num-1; i++) num_cond = num_cond + ordv - kset.GetMult()[i];
	return ((num-1)*ordv-num_cond);
}

// evaluate the basis function at u,v
BspSurf<double> BspSurfBasisFunc::GetBspSurf() const
{    
	return *b;
}

// create the BspSurf representation of the basis function
BspSurf<double> BspSurfBasisFunc::CreateBspSurf() const
{
	// compute dimension in u and v
	int dimu = ComputeDimU();
	int dimv = ComputeDimV();
	
	Matrix<double> cpts(dimu,dimv,0.0);

	KnotSet ksetu(*ktsu,ordu,ordu+1);
	KnotSet ksetv(*ktsv,ordv,ordv+1);
	
	int numu = ksetu.GetNumDistinct();
	int numv = ksetv.GetNumDistinct();

	// create distinct knot and multiplicity vectors 
	Vector<double> dts1(ksetu.GetDistinctKnots());
	Vector<int> mult1(ksetu.GetMult());
	
	Vector<double> dts2(ksetv.GetDistinctKnots());
	Vector<int> mult2(ksetv.GetMult());

	// assign multiplicity ordu,ordv at the two ends
	mult1[0]=ordu;
	mult1[numu-1]=ordu;

	mult2[0]=ordv;
	mult2[numv-1]=ordv;

	// compute offsets
	int start1=ksetu.GetMult()[0];
	int start2=ksetv.GetMult()[0];

	KnotSet kset1(dts1,mult1,ordu,dimu);
	KnotSet kset2(dts2,mult2,ordv,dimv);

	// assign 1.0 to appropriate control point
	cpts[ordu-start1][ordv-start2]=1.0;
	
	//  create and return the BspSurf
	return BspSurf<double>(cpts,kset1.GetKnots(),kset2.GetKnots(),ordu,ordv,dimu,dimv);
}

Vector<double> BspSurfBasisFunc::GetKnotsU() const
{
	return *ktsu;
}


Vector<double> BspSurfBasisFunc::GetKnotsV() const
{
	return *ktsv;
}


// EVALUATORS

// evaluate the basis function at u,v
double BspSurfBasisFunc::Eval(double u, double v) const
{    
	// create the BspSurf for the basis function and evaluate it at u,v
	if (u < (*ktsu)[0] || u > (*ktsu)[ordu]) return 0;
	if (v < (*ktsv)[0] || v > (*ktsv)[ordv]) return 0;
	return (*b).Eval(u,v);
}

// evaluate the basis function at u,v
double BspSurfBasisFunc::operator()(double u, double v) const
{    
	// create the BspSurf for the basis function and evaluate it at u,v
	if (u < (*ktsu)[0] || u > (*ktsu)[ordu]) return 0;
	if (v < (*ktsv)[0] || v > (*ktsv)[ordv]) return 0;
	return (*b)(u,v);
}


// evaluate the basis function at u,v
double BspSurfBasisFunc::operator()(int valu, int valv, double u, double v) const
{    
	// create the BspSurf for the basis function and evaluate it at u,v
	return (*b).Derive(valu, valv, u,v);	
}


double BspSurfBasisFunc::Derive(int m, int n, double u, double v) const
{
	return (*b).Derive(m,n,u,v);
}

double BspSurfBasisFunc::GetLeftLimitU() const
{
	return (*ktsu)[0];
}

double BspSurfBasisFunc::GetRightLimitU() const
{
	return (*ktsu)[ordu];
}

double BspSurfBasisFunc::GetLeftLimitV() const
{
	return (*ktsv)[0];
}

double BspSurfBasisFunc::GetRightLimitV() const
{
	return (*ktsv)[ordv];
}

// READ and WRITE

void BspSurfBasisFunc::write(std::ostream& os) 
{
	os << "Bspline Surface basis function\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	
	os << "knots in u are\n";
	(*ktsu).write(os);
	os << "\nknots in v are \n";
	(*ktsv).write(os);
	
	os << "\nBspSurf representation\n";
	(*b).write(os);
}


void BspSurfBasisFunc::read(std::istream& is)
{
	int Ordu, Ordv;
	std::cout << "order of Bspline Surface basis function in u and v";
	is >> Ordu >> Ordv;
	Vector<double> Ktsu(Ordu+1), Ktsv(Ordv+1);
	std::cout << "knots in u";
	is >> Ktsu;
	std::cout << "knots in v";
	is >> Ktsv;
	*this = BspSurfBasisFunc(Ordu,Ordv,Ktsu,Ktsv);
} 



void BspSurfBasisFunc::writefile(std::ofstream& ofs) 
{
	ofs << "Bspline Surface basis function\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	
	ofs << "knots in u are\n";
	(*ktsu).writefile(ofs);
	ofs << "\nknots in v are\n";
	(*ktsv).writefile(ofs);
	
	ofs << "\nBspSurf representation\n";
	(*b).writefile(ofs);
}


void BspSurfBasisFunc::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv;
	ifs >> Ordu >> Ordv;
	Vector<double> Ktsu(Ordu+1), Ktsv(Ordv+1);
	ifs >> Ktsu;
	ifs >> Ktsv;
	*this = BspSurfBasisFunc(Ordu,Ordv,Ktsu,Ktsv);
} 


// CONSTRUCTORS

// default constructor
BspSurfBasisFuncSet::BspSurfBasisFuncSet() : ordu(0), ordv(0), numu(0), numv(0), ktsu(), ktsv(), b() {}


// constructor taking order in u, order in v, knots in u and v
BspSurfBasisFuncSet::BspSurfBasisFuncSet(int Ordu, int Ordv, int Numu, int Numv, const Vector<double>& Ktsu, const Vector<double>& Ktsv) : 
	ordu(Ordu), ordv(Ordv), numu(Numu), numv(Numv), ktsu(new Vector<double>(Ktsu)), 
		ktsv(new Vector<double>(Ktsv)), b(new Matrix<BspSurfBasisFunc>(Numu-Ordu,Numv-Ordv,BspSurfBasisFunc()))
{
	Vector<double> v1(ordu+1), v2(ordv+1);
	for (int i=0; i<numu-ordu; i++) {
		for (int j=0; j<numv-ordv; j++) {
			for (int k=0; k<=ordu; k++) v1[k] = (*ktsu)[i+k];
			for (int l=0; l<=ordv; l++) v2[l] = (*ktsv)[j+l];
			(*b)[i][j] = BspSurfBasisFunc(ordu,ordv,v1,v2);
		}
	}
}



// EVALUATORS
// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
Matrix<double> BspSurfBasisFuncSet::Eval(double u, double v) const
{    
	Matrix<double> mat(numu-ordu,numv-ordv);

	for (int i=0; i<numu-ordu; i++) 
		for (int j=0; j<numv-ordv; j++) 
			mat[i][j] = (*b)[i][j].Eval(u,v);

	return mat;
}


// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
Matrix<double> BspSurfBasisFuncSet::operator()(double u, double v) const
{    
	Matrix<double> mat(numu-ordu,numv-ordv);

	for (int i=0; i<numu-ordu; i++) 
		for (int j=0; j<numv-ordv; j++) 
			mat[i][j] = (*b)[i][j](u,v);//.Eval(u,v);

	return mat;
}



Matrix<double> BspSurfBasisFuncSet::CreateMatrixKnotAveragesU() const
{
	return BspCurvBasisFuncSet(*ktsu,ordu,numu).CreateMatrixKnotAverages();
}


Matrix<double> BspSurfBasisFuncSet::CreateMatrixKnotAveragesV() const
{
	return BspCurvBasisFuncSet(*ktsv,ordv,numv).CreateMatrixKnotAverages();
}




Matrix<double> BspSurfBasisFuncSet::CreateIntegralNewU(const BspSurf<double>& s, const BspSurf<double>& norm, double u1, double u2, double v1, double v2) const
{
	Matrix<double> res(numu-ordu,numv-ordv,0.0);
	
	for (int j=0; j<numu-ordu; j++) 
		for (int k=0; k<numv-ordv; k++) {
			for (int i=1; i<=s.GetOrdU(); i++) {
				BspSurf<double> surf = (*b)[j][k].GetBspSurf().IntegrateU(i);
				double y1=surf.GetLeftLimitU();
				double y2=surf.GetRightLimitU();
				BspSurf<double> s1 = s.DeriveU(i-1).Subdivide(surf.GetLeftLimitU(),surf.GetRightLimitU(),surf.GetLeftLimitV(),surf.GetRightLimitV());
				BspCurv<double> c2 = surf.GetIsoparametricU(y2);
				BspCurv<double> c3 = s1.GetIsoparametricU(y2);
				BspCurv<double> c4 = surf.GetIsoparametricU(y1);
				BspCurv<double> c5 = s1.GetIsoparametricU(y1);
				BspCurv<double> prod1 = c2.Product(c3);
				BspCurv<double> prod2 = c4.Product(c5);
				BspCurv<double> sub = prod1.Subtract(prod2);
				res[j][k] = res[j][k] + pow(-1.0,double(i-1))*sub.Integrate(v1,v2);
		}
	}


	return res;
}


Matrix<double> BspSurfBasisFuncSet::CreateIntegralNewV(const BspSurf<double>& s, const BspSurf<double>& norm, double u1, double u2, double v1, double v2) const
{
	Matrix<double> res(numu-ordu,numv-ordv,0.0);

	for (int j=0; j<numu-ordu; j++) 
		for (int k=0; k<numv-ordv; k++) {
			for (int i=1; i<=s.GetOrdV(); i++) {
				BspSurf<double> surf = (*b)[j][k].GetBspSurf().IntegrateV(i);
				double y1=surf.GetLeftLimitV();
				double y2=surf.GetRightLimitV();
				BspSurf<double> s1 = s.DeriveV(i-1).Subdivide(surf.GetLeftLimitU(),surf.GetRightLimitU(),surf.GetLeftLimitV(),surf.GetRightLimitV());
				BspCurv<double> c2 = surf.GetIsoparametricV(y2);
				BspCurv<double> c3 = s1.GetIsoparametricV(y2);
				BspCurv<double> c4 = surf.GetIsoparametricV(y1);
				BspCurv<double> c5 = s1.GetIsoparametricV(y1);
				BspCurv<double> prod1 = c2.Product(c3);
				BspCurv<double> prod2 = c4.Product(c5);
				BspCurv<double> sub = prod1.Subtract(prod2);
				res[j][k] = res[j][k] + pow(-1.0,double(i-1))*sub.Integrate(u1,u2);
		}
	}


	return res;
}

Matrix<double> BspSurfBasisFuncSet::CreateIntegral(double u1, double u2, double v1, double v2) const
{
	Matrix<double> res(numu-ordu,numv-ordv);

	for (int i=0; i<numu-ordu; i++) 
		for (int j=0; j<numv-ordv; j++)
			res[i][j] = (*b)[i][j].GetBspSurf().Integrate(u1,u2,v1,v2);

	return res;
}






Matrix<double> BspSurfBasisFuncSet::ComputeUBasisMatrix() const
{

	
	Matrix<double> mat(numu-ordu,numu-ordu);
	Vector<double> u(KnotSet((*ktsu),ordu,numu).ComputeKnotSetAver());


	BspCurvBasisFuncSet b((*ktsu),ordu,numu);

	for (int i=0; i<numu-ordu; i++) 
			mat[i] = b.ComputeVecBasis(u[i]);

	return mat;
}

Matrix<double> BspSurfBasisFuncSet::ComputeVBasisMatrix() const
{


	Matrix<double> mat(numv-ordv,numv-ordv);
	Vector<double> v(KnotSet((*ktsv),ordv,numv).ComputeKnotSetAver());

	BspCurvBasisFuncSet b((*ktsv),ordv,numv);

	for (int i=0; i<numv-ordv; i++)
			mat[i] = b.ComputeVecBasis(v[i]);

	return mat;
}



// READ AND WRITE

void BspSurfBasisFuncSet::write(std::ostream& os) 
{
	os << "Bspline Surface basis function set\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are\n";
	os << *ktsv;
	os << "\nBspSurf representation\n";
	for (int i=0; i<numu-ordu; i++)
		for (int j=0; j<numv-ordv; j++) (*b)[i][j].write(os);
}


void BspSurfBasisFuncSet::read(std::istream& is)
{
	int Ordu, Ordv, Numu, Numv;
	std::cout << "order of Bspline Surface  basis function in u and v";
	is >> Ordu >> Ordv;
	std::cout << "number of knots in u and v\n";
	is >> Numu >> Numv;
	Vector<double> Ktsu(Numu);
	Vector<double> Ktsv(Numv);
	is >> Ktsu;
	is >> Ktsv;
	*this = BspSurfBasisFuncSet(Ordu,Ordv,Numu,Numv,Ktsu,Ktsv);
} 



void BspSurfBasisFuncSet::writefile(std::ofstream& ofs) 
{
	ofs << "Bspline Surface basis function\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are\n";
	ofs << *ktsv;
	
	ofs << "\nBspSurf representation is\n";
	for (int i=0; i<numu-ordu; i++)
		for (int j=0; j<numv-ordv; j++) (*b)[i][j].writefile(ofs);
}


void BspSurfBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Numu, Numv;
	ifs >> Ordu >> Ordv;
	ifs >> Numu >> Numv;
	Vector<double> Ktsu(Numu), Ktsv(Numv);
	ifs >> Ktsu; 
	ifs >> Ktsv;
	*this = BspSurfBasisFuncSet(Ordu,Ordv,Numu,Numv,Ktsu,Ktsv);
} 

