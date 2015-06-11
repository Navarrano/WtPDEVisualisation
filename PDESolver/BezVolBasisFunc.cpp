
#include "bezvol.h"

// CONSTRUCTORS

// default constructor
BezVolBasisFunc::BezVolBasisFunc() : ordu(0), ordv(0), ordw(0), ktsu(), ktsv(), ktsw() {}

// constructor taking order in u, order in v, knots in u and v
BezVolBasisFunc::BezVolBasisFunc(int Ordu, int Ordv, int Ordw, const Vector<double>& Ktsu, 
	const Vector<double>& Ktsv, const Vector<double>& Ktsw) : 
	ordu(Ordu), ordv(Ordv), ordw(Ordw), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)), ktsw(new Vector<double>(Ktsw)), 
		b(new BezVol<double>(CreateBezVol())) { }


// EVALUATORS

// evaluate the basis function at u, v. Creates the BezVol
// representation and evaluates it
double BezVolBasisFunc::Eval(double u, double v, double w) const
{    
	return (*b).Eval(u,v,w);
}


// evaluate the basis function at u, v. Creates the BezVol
// representation and evaluates it
double BezVolBasisFunc::operator()(double u, double v, double w) const
{    
	return (*b)(u,v,w);
}

// evaluate the basis function at u, v. Creates the BezVol
// representation and evaluates it
double BezVolBasisFunc::operator()(int valu, int valv, int valw, double u, double v, double w) const
{    
	return (*b).Derive(valu,valv,valw,u,v,w);
}


// ACCESS FUNCTIONS

// evaluate the basis function at u,v
BezVol<double> BezVolBasisFunc::GetBezVol() const
{    
	return *b;
}

// compute the dimension of the u knot set
int BezVolBasisFunc::ComputeDimU() const
{
	return ordw;
}

// compute the dimension of the v knot set
int BezVolBasisFunc::ComputeDimV() const
{
	return ordv;
}


// compute the dimension of the v knot set
int BezVolBasisFunc::ComputeDimW() const
{
	return ordw;
}


// creates the BezVol representation of the basis function
BezVol<double> BezVolBasisFunc::CreateBezVol() const
{

	// matrix for control points
	Matrix3D<double> cpts(ordu,ordv,ordw,0.0);

	// create KnotSet objects
	KnotSet ksetu(*ktsu,ordu,ordu+1);
	KnotSet ksetv(*ktsv,ordv,ordv+1);
	KnotSet ksetw(*ktsw,ordw,ordw+1);	

	int numu = ksetu.GetNumDistinct();
	int numv = ksetv.GetNumDistinct();
	int numw = ksetw.GetNumDistinct();

	// create Vectors for distinct knots and multiplicities
	Vector<double> dts1(ksetu.GetDistinctKnots());
	Vector<int> mult1(ksetu.GetMult());
	
	Vector<double> dts2(ksetv.GetDistinctKnots());
	Vector<int> mult2(ksetv.GetMult());

	Vector<double> dts3(ksetw.GetDistinctKnots());
	Vector<int> mult3(ksetw.GetMult());

	// assign ordu, ordv multiplicities at the ends
	mult1[0]=ordu;
	mult1[numu-1]=ordu;

	mult2[0]=ordv;
	mult2[numv-1]=ordv;

	mult3[0]=ordw;
	mult3[numw-1]=ordw;

	// find offsets
	int start1=ksetu.GetMult()[0];
	int start2=ksetv.GetMult()[0];
	int start3=ksetw.GetMult()[0];

	// create KnotSet objects for Volace
	KnotSet kset1(dts1,mult1,ordu,numu);
	KnotSet kset2(dts2,mult2,ordv,numv);
	KnotSet kset3(dts3,mult3,ordw,numw);

	// assign 1.0 to appropriate control point
	cpts[ordw-start3][ordu-start1][ordv-start2]=1.0;
	
	// create and return the BezVol
	return BezVol<double>(cpts,kset1.GetKnots(),kset2.GetKnots(),kset3.GetKnots(),ordu,ordv,ordw);
}


double BezVolBasisFunc::Derive(int m, int n, int p, double u, double v, double w) const
{
	return (*b).Derive(m,n,p,u,v,w);
}

double BezVolBasisFunc::GetLeftLimitU() const
{
	return (*ktsu)[0];
}

double BezVolBasisFunc::GetRightLimitU() const
{
	return (*ktsu)[ordu];
}

double BezVolBasisFunc::GetLeftLimitV() const
{
	return (*ktsv)[0];
}

double BezVolBasisFunc::GetRightLimitV() const
{
	return (*ktsv)[ordv];
}

double BezVolBasisFunc::GetLeftLimitW() const
{
	return (*ktsw)[0];
}

double BezVolBasisFunc::GetRightLimitW() const
{
	return (*ktsw)[ordw];
}


// READ and WRITE

void BezVolBasisFunc::write(std::ostream& os) 
{
	os << "Bezier Basis Volume function\n";
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
	
	os << "\nBezVol representation is\n";
	os << *b;
}


void BezVolBasisFunc::read(std::istream& is)
{
	int Ordu, Ordv, Ordw;
	std::cout << "order of Bezier Volume in u and v and w";
	is >> Ordu >> Ordv >> Ordw;
	Vector<double> Ktsu(Ordu+1), Ktsv(Ordv+1), Ktsw(Ordw+1);
	std::cout << "input knots in u\n";
	is >> Ktsu;
	std::cout << "input knots in v\n";
	is >> Ktsv;
	std::cout << "input knots in w\n";
	is >> Ktsw;
	*this = BezVolBasisFunc(Ordu,Ordv,Ordw,Ktsu,Ktsv,Ktsw);
} 



void BezVolBasisFunc::writefile(std::ofstream& ofs) 
{
	ofs << "Bezier Volume basis function\n";
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
	ofs << "\nBezVol representation is\n";
	(*b).writefile(ofs);
}



void BezVolBasisFunc::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw;

	ifs >> Ordu >> Ordv >> Ordw;
	
	Vector<double> Ktsu(Ordu+1),Ktsv(Ordv+1), Ktsw(Ordw+1);

	ifs >> Ktsu >> Ktsv >> Ktsw;
	
	*this = BezVolBasisFunc(Ordu,Ordv,Ordw,Ktsu,Ktsv,Ktsw);
} 



// CONSTRUCTORS

// default constructor
BezVolBasisFuncSet::BezVolBasisFuncSet() : ordu(0), ordv(0), ordw(0), ktsu(), ktsv(), ktsw(), b() {}


// constructor taking order in u, order in v, knots in u and v
BezVolBasisFuncSet::BezVolBasisFuncSet(int Ordu, int Ordv, int Ordw, 
const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw) : 
	ordu(Ordu), ordv(Ordv), ordw(Ordw), ktsu(new Vector<double>(Ktsu)), 
	ktsv(new Vector<double>(Ktsv)), ktsw(new Vector<double>(Ktsw)),
	b(new Matrix3D<BezVolBasisFunc>(Ordu,Ordv,Ordw,BezVolBasisFunc()))
{
	Vector<double> v1(ordu+1), v2(ordv+1), v3(ordw+1);
	for (int i=0; i<ordu; i++) {
		for (int l=0; l<=ordu; l++) v1[l] = (*ktsu)[i+l];
		for (int j=0; j<ordv; j++) {
			for (int m=0; m<=ordv; m++) v2[m] = (*ktsv)[j+m];
			for (int k=0; k<ordw; k++) {
				for (int n=0; n<=ordw; n++) v3[n] = (*ktsw)[k+n];
				(*b)[k][i][j] = BezVolBasisFunc(ordu,ordv,ordw,v1,v2,v3);
			}
		}
	}
}


// EVALUATORS
// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
Matrix3D<double> BezVolBasisFuncSet::Eval(double u, double v, double w) const
{    
	Matrix3D<double> mat(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) 
			for (int k=0; k<ordw; k++)
				mat[k][i][j] = (*b)[k][i][j].Eval(u,v,w);

	return mat;
}


// evaluate the basis function at u, v. Creates the BezSurf
// representation and evaluates it
Matrix3D<double> BezVolBasisFuncSet::operator()(double u, double v, double w) const
{    
	Matrix3D<double> mat(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) 
			for (int k=0; k<ordw; k++)
				mat[k][i][j] = (*b)[k][i][j].Eval(u,v,w);

	return mat;
}


// READ ADN WRITE

void BezVolBasisFuncSet::write(std::ostream& os) 
{
	os << "Bezier Volume basis function set\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "order in w is " << ordw << "\n";
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are\n";
	os << *ktsv;
	os << "\nknots in w are\n";
	os << *ktsw;
	os << "\nBezVol representation\n";
	for (int i=0; i<ordu; i++)
		for (int j=0; j<ordv; j++) 
			for (int k=0; k<ordw; k++) os << (*b)[k][i][j];
}


void BezVolBasisFuncSet::read(std::istream& is)
{
	int Ordu, Ordv, Ordw;
	std::cout << "order of Bezier Volume basis function in u,v,w";
	is >> Ordu >> Ordv >> Ordw;
	Vector<double> Ktsu(Ordu+Ordu);
	Vector<double> Ktsv(Ordv+Ordv);
	Vector<double> Ktsw(Ordw+Ordw);
	is >> Ktsu;
	is >> Ktsv;
	is >> Ktsw;
	*this = BezVolBasisFuncSet(Ordu,Ordv,Ordw,Ktsu,Ktsv,Ktsw);
} 



void BezVolBasisFuncSet::writefile(std::ofstream& ofs) 
{
	ofs << "Bezier Volumes basis function set\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "order in w is " << ordw << "\n";
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are\n";
	ofs << *ktsv;
	ofs << "\nknots in w are\n";
	ofs << *ktsw;
	ofs << "\nBezVol representation is\n";
	for (int i=0; i<ordu; i++)
		for (int j=0; j<ordv; j++) 
			for (int k=0; k<ordw; k++) (*b)[k][i][j].writefile(ofs);
}


void BezVolBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw;
	ifs >> Ordu >> Ordv >> Ordw;
	Vector<double> Ktsu(Ordu+Ordu), Ktsv(Ordv+Ordv), Ktsw(Ordw+Ordw);
	ifs >> Ktsu; 
	ifs >> Ktsv;
	ifs >> Ktsw;
	*this = BezVolBasisFuncSet(Ordu,Ordv,Ordw,Ktsu,Ktsv,Ktsw);
} 



