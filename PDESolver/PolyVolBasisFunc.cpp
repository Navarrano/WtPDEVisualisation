
#include "polyvol.h"


// CONSTRUCTORS

// default constructor
PolyVolBasisFunc::PolyVolBasisFunc() : ordu(0), ordv(0), ordw(0), indexu(0), indexv(0), indexw(0), leftlimitu(0), leftlimitv(0), rightlimitu(0), rightlimitv(0), 
leftlimitw(0), rightlimitw(0) {}

// alternate constructor
PolyVolBasisFunc::PolyVolBasisFunc(int Ordu, int Ordv, int Ordw, int Indexu, int Indexv, int Indexw, double LeftLimitU, double RightLimitU, double LeftLimitV, double RightLimitV, 
						   double LeftLimitW, double RightLimitW) : 
	ordu(Ordu), ordv(Ordv), ordw(Ordw), indexu(Indexu), indexv(Indexv), indexw(Indexw),
		leftlimitu(LeftLimitU), leftlimitv(LeftLimitV), 
		rightlimitu(RightLimitU), rightlimitv(RightLimitV), 
		leftlimitw(LeftLimitW), rightlimitw(RightLimitW), b(new PolyVol<double>(CreatePolyVol())) { }


// EVALUATORS

// evaluate the basis function
double PolyVolBasisFunc::Eval(double u, double v, double w) const
{    
	return (*b).Eval(u,v,w);
}

// evaluate the basis function
double PolyVolBasisFunc::operator()(double u, double v, double w) const
{    
	return (*b)(u,v,w);
}

double PolyVolBasisFunc::operator()(int valu, int valv, int valw, double u, double v, double w) const
{    
	return (*b).Derive(valu,valv,valw,u,v,w);
}



// ACCESS FUNCTIONS

PolyVol<double> PolyVolBasisFunc::GetPolyVol() const
{
	return (*b);
}

// create PolyVol representation of the basis function
PolyVol<double> PolyVolBasisFunc::CreatePolyVol() const
{	
	Matrix3D<double> coeffs(ordu,ordv,ordw,0.0);

	coeffs[indexu-1][indexv-1][indexw-1]=1.0;
	return PolyVol<double>(coeffs,ordu,ordv,ordw,leftlimitu,rightlimitu,leftlimitv,rightlimitv,leftlimitw,rightlimitw);
}

double PolyVolBasisFunc::Derive(int m, int n, int p, double u, double v, double w) const
{
	return (*b).Derive(m,n,p,u,v,w);
}

double PolyVolBasisFunc::GetLeftLimitU() const
{
	return leftlimitu;
}

double PolyVolBasisFunc::GetRightLimitU() const
{
	return rightlimitu;
}

double PolyVolBasisFunc::GetLeftLimitV() const
{
	return leftlimitv;
}

double PolyVolBasisFunc::GetRightLimitV() const
{
	return rightlimitv;
}

double PolyVolBasisFunc::GetLeftLimitW() const
{
	return leftlimitw;
}

double PolyVolBasisFunc::GetRightLimitW() const
{
	return rightlimitw;
}


// READ and WRITE

void PolyVolBasisFunc::write(std::ostream& os) 
{
	os << "Poly Basis Volume function\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "order in w is " << ordw << "\n";
	os << "index is " << indexu << " " << indexv << " " << indexw << "\n";
	os << "left right limit u,v,w\n";
	os << leftlimitu << rightlimitu << leftlimitv << rightlimitv << leftlimitw << rightlimitw;
	
	os << "\nPolyVol representation is\n";
	os << (*b);
}


void PolyVolBasisFunc::read(std::istream& is)
{
	int Ordu, Ordv, Ordw, Indexu, Indexv, Indexw;
	std::cout << "order of Poly Volume in u and v and w";
	is >> Ordu >> Ordv >> Ordw;
	std::cout << "index in u, v and w\n";
	is >> Indexu >> Indexv >> Indexw;
	double leftu, rightu, leftv, rightv, leftw, rightw;
	std::cout << "left right limit u\n";
	is >> leftu >> rightu;
	std::cout << "\nleft right limit v\n";
	is >> leftv >> rightv;
	std::cout << "left right limit w\n";
	is >> leftw >> rightw;
	*this = PolyVolBasisFunc(Ordu,Ordv,Ordw,Indexu,Indexv,Indexw,leftu,rightu,leftv,rightv,leftw,rightw);
} 



void PolyVolBasisFunc::writefile(std::ofstream& ofs) 
{
	ofs << "Poly Volume basis function\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "order in w is " << ordw << "\n";
	ofs << "index in u, v, w " << indexu << " " << indexv << " " << indexw << "\n";
	ofs << "left right limit u\n";
	ofs << leftlimitu << rightlimitu;
	ofs << "\nleft right limit v\n";
	ofs << leftlimitv << rightlimitv;
	ofs << "\nleft right limit w\n";
	ofs << leftlimitw << rightlimitw;
	ofs << "\nPolyVol representation is\n";
	ofs << b;
}



void PolyVolBasisFunc::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw;
	int Indexu, Indexv, Indexw;

	ifs >> Ordu >> Ordv >> Ordw;
	ifs >> Indexu >> Indexv >> Indexw;

	double leftu, rightu, leftv, rightv, leftw, rightw;
	ifs >> leftu >> rightu >> leftv >> rightv >> leftw >> rightw;
	
	*this = PolyVolBasisFunc(Ordu,Ordv,Ordw,Indexu, Indexv, Indexw, leftu,rightu,leftv,rightv,leftw,rightw);
} 



// CONSTRUCTORS
PolyVolBasisFuncSet::PolyVolBasisFuncSet() : ordu(0), ordv(0), ordw(0), b() {}

PolyVolBasisFuncSet::PolyVolBasisFuncSet(int Ordu, int Ordv, int Ordw,double LeftLimitU, double RightLimitU, double LeftLimitV, double RightLimitV,
					double LeftLimitW, double RightLimitW) : ordu(Ordu), ordv(Ordv), ordw(Ordw), 
					leftlimitu(LeftLimitU), rightlimitu(RightLimitU), leftlimitv(LeftLimitV), rightlimitv(RightLimitV),
					leftlimitw(LeftLimitW), rightlimitw(RightLimitW),
					b(new Matrix3D<PolyVolBasisFunc>(Ordu,Ordv,Ordw,PolyVolBasisFunc()))
{
	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++)
			for (int k=0; k<ordw; k++)
				(*b)[k][i][j] = PolyVolBasisFunc(Ordu,Ordv,Ordw, i+1, j+1, k+1,leftlimitu,rightlimitu,leftlimitv,rightlimitv,leftlimitw,rightlimitw);
}


// EVALUATORS
// evaluate the basis Func by converting to BspVol form
Matrix3D<double> PolyVolBasisFuncSet::Eval(double u, double v, double w) const
{    
	Matrix3D<double> mat(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) 
			for (int k=0; k<ordw; k++) mat[k][i][j] = (*b)[k][i][j](u,v,w);
	return mat;
}


// evaluate the basis Func by converting to BspVol form
Matrix3D<double> PolyVolBasisFuncSet::operator()(double u, double v, double w) const
{   
	Matrix3D<double> mat(ordu,ordv,ordw);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) 
			for (int k=0; k<ordw; k++) mat[k][i][j] = (*b)[k][i][j](u,v,w);
	return mat;
}



// READ and WRITE
void PolyVolBasisFuncSet::write(std::ostream& os) 
{
	os << "Polyier Vol basis Func set\n";
	os << "order of basis Funcs in u is " << ordu << "\n";
	os << "order of basis Funcs in v is " << ordv << "\n";
	os << "order of basis Funcs in w is " << ordw << "\n";
	os << "number of basis Func in u  is " << ordu << "\n";
	os << "number of basis Func in v  is " << ordv << "\n";
	os << "number of basis Func in w  is " << ordw << "\n";
	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) 
			for (int k=0; k<ordw; k++) {
				os << "\n" << i << j << k << "th Basis Func" << "\n";
				os << (*b)[k][i][j];
			}
}

void PolyVolBasisFuncSet::read(std::istream& is)
{
	int Ordu, Ordv, Ordw;
	std::cout << "order of Poly Vol";
	is >> Ordu >> Ordv >> Ordw;
	double Lu, Ru, Lv, Rv, Lw, Rw;
	is >> Lu >> Ru >> Lv >> Rv >> Lw >> Rw;
	*this = PolyVolBasisFuncSet(Ordu,Ordv,Ordw,Lu,Ru,Lv,Rv,Lw,Rw);
} 


void PolyVolBasisFuncSet::writefile(std::ofstream& ofs) 
{
	ofs << "Polyier Vol basis Func set\n";
	ofs << "order of basis Funcs in u is " << ordu << "\n";
	ofs << "order of basis Funcs in v is " << ordv << "\n";
	ofs << "order of basis Funcs in w is " << ordv << "\n";
	ofs << "number of basis Func in u  is " << ordu << "\n";
	ofs << "number of basis Func in v  is " << ordv << "\n";
	ofs << "number of basis Func in w  is " << ordw << "\n";
	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) 
			for (int k=0; k<ordw; k++) {
				ofs << "\n" <<  i << j << k << "th Basis Func" << "\n";
				(*b)[k][i][j].writefile(ofs);
	}	
}


void PolyVolBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw;
	ifs >> Ordu >> Ordv >> Ordw;
	double Lu, Ru, Lv, Rv, Lw, Rw;
	ifs >> Lu >> Ru >> Lv >> Rv >> Lw >> Rw;
	*this = PolyVolBasisFuncSet(Ordu,Ordv,Ordw,Lu,Ru,Lv,Rv,Lw,Rw);
} 

