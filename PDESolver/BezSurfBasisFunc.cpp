
#include "bezsurf.h"

// CONSTRUCTORS

// default constructor
BezSurfBasisFunc::BezSurfBasisFunc() : ordu(0), ordv(0), ktsu(), ktsv() {}

// constructor taking order in u, order in v, knots in u and v
BezSurfBasisFunc::BezSurfBasisFunc(int Ordu, int Ordv, const Vector<double>& Ktsu, const Vector<double>& Ktsv) : 
	ordu(Ordu), ordv(Ordv), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)), 
		b(new BezSurf<double>(CreateBezSurf())) { }


// EVALUATORS


// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
double BezSurfBasisFunc::Eval(double u, double v) const
{    
	return (*b).Eval(u,v);
}


// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
double BezSurfBasisFunc::operator()(double u, double v) const
{    
	return (*b)(u,v);
}


// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
double BezSurfBasisFunc::operator()(int valu, int valv, double u, double v) const
{    
	return (*b).Derive(valu,valv,u,v);
}


// ACCESS FUNCTIONS

// evaluate the basis function at u,v
BezSurf<double> BezSurfBasisFunc::GetBezSurf() const
{    
	return (*b);
}


int BezSurfBasisFunc::GetOrdU() const
{
	return ordu;
}

int BezSurfBasisFunc::GetOrdV() const
{
	return ordv;
}

Vector<double> BezSurfBasisFunc::GetKnotsU() const
{
	return *ktsu;
}


Vector<double> BezSurfBasisFunc::GetKnotsV() const
{
	return *ktsv;
}

// compute the dimension of the u knot set
int BezSurfBasisFunc::ComputeDimU() const
{
	// create KnotSet object 
	return ordu;
}

// compute the dimension of the v knot set
int BezSurfBasisFunc::ComputeDimV() const
{
	return ordv;
}


// creates the BezSurf representation of the basis function
BezSurf<double> BezSurfBasisFunc::CreateBezSurf() const
{

	// matrix for control points
	Matrix<double> cpts(ordu,ordv,0.0);

	// create KnotSet objects
	KnotSet ksetu(*ktsu,ordu,ordu+1);
	KnotSet ksetv(*ktsv,ordv,ordv+1);
	
	int numu = ksetu.GetNumDistinct();
	int numv = ksetv.GetNumDistinct();

	// create Vectors for distinct knots and multiplicities
	Vector<double> dts1(ksetu.GetDistinctKnots());
	Vector<int> mult1(ksetu.GetMult());
	
	Vector<double> dts2(ksetv.GetDistinctKnots());
	Vector<int> mult2(ksetv.GetMult());

	// assign ordu, ordv multiplicities at the ends
	mult1[0]=ordu;
	mult1[numu-1]=ordu;

	mult2[0]=ordv;
	mult2[numv-1]=ordv;

	// find offsets
	int start1=ksetu.GetMult()[0];
	int start2=ksetv.GetMult()[0];

	// create KnotSet objects for surface
	KnotSet kset1(dts1,mult1,ordu,numu);
	KnotSet kset2(dts2,mult2,ordv,numv);

	// assign 1.0 to appropriate control point
	cpts[ordu-start1][ordv-start2]=1.0;
	
	// create and return the BezSurf
	return BezSurf<double>(cpts,kset1.GetKnots(),kset2.GetKnots(),ordu,ordv);
}


double BezSurfBasisFunc::Derive(int m, int n, double u, double v) const
{
	return (*b).Derive(m,n,u,v);
}

double BezSurfBasisFunc::GetLeftLimitU() const
{
	return (*ktsu)[0];
}

double BezSurfBasisFunc::GetRightLimitU() const
{
	return (*ktsu)[ordu];
}

double BezSurfBasisFunc::GetLeftLimitV() const
{
	return (*ktsv)[0];
}

double BezSurfBasisFunc::GetRightLimitV() const
{
	return (*ktsv)[ordv];
}

// READ AND WRITE

void BezSurfBasisFunc::write(std::ostream& os) 
{
	os << "Bezier Surface basis function\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are\n";
	os << *ktsv;
	os << "\nBezSurf representation\n";
	b->write(os);
}


void BezSurfBasisFunc::read(std::istream& is)
{
	int Ordu, Ordv;
	std::cout << "order of Bezier Surface  basis function in u and v";
	is >> Ordu >> Ordv;
	Vector<double> Ktsu(Ordu+1);
	Vector<double> Ktsv(Ordv+1);
	is >> Ktsu;
	is >> Ktsv;
	*this = BezSurfBasisFunc(Ordu,Ordv,Ktsu,Ktsv);
} 



void BezSurfBasisFunc::writefile(std::ofstream& ofs) 
{
	ofs << "Bezier Surface basis function\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are\n";
	ofs << *ktsv;
	
	ofs << "\nBezSurf representation is\n";
	b->writefile(ofs);
}


void BezSurfBasisFunc::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv;
	ifs >> Ordu >> Ordv;
	Vector<double> Ktsu(Ordu+1), Ktsv(Ordv+1);
	ifs >> Ktsu; 
	ifs >> Ktsv;
	*this = BezSurfBasisFunc(Ordu,Ordv,Ktsu,Ktsv);
} 




// CONSTRUCTORS

// default constructor
BezSurfBasisFuncSet::BezSurfBasisFuncSet() : ordu(0), ordv(0), ktsu(), ktsv(), b() {}


// constructor taking order in u, order in v, knots in u and v
BezSurfBasisFuncSet::BezSurfBasisFuncSet(int Ordu, int Ordv, const Vector<double>& Ktsu, const Vector<double>& Ktsv) : 
	ordu(Ordu), ordv(Ordv), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)),
		b(new Matrix<BezSurfBasisFunc>(Ordu,Ordv,BezSurfBasisFunc())) 
{
	Vector<double> v1(ordu+1), v2(ordv+1);
	for (int i=0; i<ordu; i++) {
		for (int k=0; k<=ordu; k++) v1[k] = (*ktsu)[i+k];
		for (int j=0; j<ordv; j++) {
			for (int l=0; l<=ordv; l++) v2[l] = (*ktsv)[j+l];
			(*b)[i][j] = BezSurfBasisFunc(ordu,ordv,v1,v2);
		}
	}
}


int BezSurfBasisFuncSet::GetOrdU() const
{
	return ordu;
}

int BezSurfBasisFuncSet::GetOrdV() const
{
	return ordv;
}

Vector<double> BezSurfBasisFuncSet::GetKnotsU() const
{
	return *ktsu;
}


Vector<double> BezSurfBasisFuncSet::GetKnotsV() const
{
	return *ktsv;
}

BezSurfBasisFunc BezSurfBasisFuncSet::GetBezSurfBasisFunc(int i, int j) const
{
	Vector<double> v1(ordu+1), v2(ordv+1);
	
	for (int k=0; k<=ordu; k++) v1[k] = (*ktsu)[i-1+k];
	for (int l=0; l<=ordv; l++) v2[l] = (*ktsv)[j-1+l];
	return BezSurfBasisFunc(ordu,ordv,v1,v2);
}

// EVALUATORS
// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
Matrix<double> BezSurfBasisFuncSet::Eval(double u, double v) const
{    
	Matrix<double> mat(ordu,ordv);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) 
			mat[i][j] = (*b)[i][j].Eval(u,v);

	return mat;
}


// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
Matrix<double> BezSurfBasisFuncSet::operator()(double u, double v) const
{    
	Matrix<double> mat(ordu,ordv);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) 
			mat[i][j] = (*b)[i][j].Eval(u,v);

	return mat;
}


// READ ADN WRITE

void BezSurfBasisFuncSet::write(std::ostream& os) 
{
	os << "Bezier Surface basis function set\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are\n";
	os << *ktsv;
	os << "\nBezSurf representation\n";
	for (int i=0; i<ordu; i++)
		for (int j=0; j<ordv; j++) ((*b)[i][j]).write(os);
}


void BezSurfBasisFuncSet::read(std::istream& is)
{
	int Ordu, Ordv;
	std::cout << "order of Bezier Surface  basis function in u and v";
	is >> Ordu >> Ordv;
	Vector<double> Ktsu(Ordu+Ordu);
	Vector<double> Ktsv(Ordv+Ordv);
	is >> Ktsu;
	is >> Ktsv;
	*this = BezSurfBasisFuncSet(Ordu,Ordv,Ktsu,Ktsv);
} 



void BezSurfBasisFuncSet::writefile(std::ofstream& ofs) 
{
	ofs << "Bezier Surface basis function\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are\n";
	ofs << *ktsv;
	
	ofs << "\nBezSurf representation is\n";
	for (int i=0; i<ordu; i++)
		for (int j=0; j<ordv; j++) ((*b)[i][j]).writefile(ofs);
}


void BezSurfBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv;
	ifs >> Ordu >> Ordv;
	Vector<double> Ktsu(Ordu+Ordu), Ktsv(Ordv+Ordv);
	ifs >> Ktsu; 
	ifs >> Ktsv;
	*this = BezSurfBasisFuncSet(Ordu,Ordv,Ktsu,Ktsv);
} 

