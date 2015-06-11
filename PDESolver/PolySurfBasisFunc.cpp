
#include "polysurf.h"


// CONSTRUCTORS

// default constructor
PolySurfBasisFunc::PolySurfBasisFunc() : ordu(0), ordv(0), indexu(0), indexv(0), leftlimitu(0), leftlimitv(0), rightlimitu(0), rightlimitv(0) {}

// alternate constructor
PolySurfBasisFunc::PolySurfBasisFunc(int Ordu, int Ordv, int Indexu, int Indexv, double LeftLimitU, double RightLimitU, double LeftLimitV, double RightLimitV) : 
	ordu(Ordu), ordv(Ordv), indexu(indexu), indexv(Indexv), leftlimitu(LeftLimitU), leftlimitv(LeftLimitV), rightlimitu(RightLimitU), rightlimitv(RightLimitV),
	b(new PolySurf<double>(CreatePolySurf())) 
{ }


// EVALUATORS

// evaluate the basis function
double PolySurfBasisFunc::Eval(double u, double v) const
{    
	return (*b).Eval(u,v);
}

// evaluate the basis function
double PolySurfBasisFunc::operator()(double u, double v) const
{    
	return (*b)(u,v);
}

double PolySurfBasisFunc::operator()(int valu, int valv, double u, double v) const
{    
	return (*b).Derive(valu,valv,u,v);
}

// ACCESS FUNCTIONS

PolySurf<double> PolySurfBasisFunc::GetPolySurf() const
{
	return (*b);
}

// create PolySurf representation of the basis function
PolySurf<double> PolySurfBasisFunc::CreatePolySurf() const
{	
	Matrix<double> coeffs(ordu,ordv,0.0);

	coeffs[indexu-1][indexv-1]=1.0;
	return PolySurf<double>(coeffs,ordu,ordv,leftlimitu,rightlimitu,leftlimitv,rightlimitv);
}

double PolySurfBasisFunc::Derive(int m, int n, double u, double v) const
{
	return (*b).Derive(m,n,u,v);
}

double PolySurfBasisFunc::GetLeftLimitU() const
{
	return leftlimitu;
}

double PolySurfBasisFunc::GetRightLimitU() const
{
	return leftlimitv;
}

double PolySurfBasisFunc::GetLeftLimitV() const
{
	return rightlimitu;
}

double PolySurfBasisFunc::GetRightLimitV() const
{
	return rightlimitv;
}

// READ and WRITE

void PolySurfBasisFunc::write(std::ostream& os) 
{
	os << "Poly Surface basis function\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "left right limit in u\n";
	os << leftlimitu << " " << rightlimitu << "\n";
	os << "\nleft right limit in v\n";
	os << leftlimitv << " " << rightlimitv << "\n";
	os << "\nPolySurf representation\n";
	os << (*b);
}


void PolySurfBasisFunc::read(std::istream& is)
{
	int Ordu, Ordv;
	int Indexu, Indexv;
	std::cout << "order of Poly Surface  basis function in u and v";
	is >> Ordu >> Ordv;
	is >> Indexu >> Indexv;
	double leftu, leftv, rightu, rightv;
	is >> leftu >> rightu;
	is >> leftv >> rightv;
	*this = PolySurfBasisFunc(Ordu,Ordv,Indexu,Indexv,leftu,rightu,leftv,rightv);
} 



void PolySurfBasisFunc::writefile(std::ofstream& ofs) 
{
	ofs << "Poly Surface basis function\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "limits in u\n";
	ofs << leftlimitu << " "  << rightlimitu;
	ofs << "\nlimits in v\n";
	ofs << leftlimitv << " " << rightlimitv;
	
	ofs << "\nPolySurf representation is\n";
	ofs << b;
}


void PolySurfBasisFunc::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv;
	int Indexu, Indexv;
	double leftu, leftv, rightu, rightv;
	ifs >> Ordu >> Ordv;
	ifs >> Indexu >> Indexv;
	ifs >> leftu >> rightu; 
	ifs >> leftv >> rightv;
	*this = PolySurfBasisFunc(Ordu,Ordv,Indexu, Indexv, leftu,rightu,leftv,rightv);
} 




// CONSTRUCTORS
PolySurfBasisFuncSet::PolySurfBasisFuncSet() : ordu(0), ordv(0), b() {}

PolySurfBasisFuncSet::PolySurfBasisFuncSet(int Ordu, int Ordv,double LeftLimitU, double RightLimitU, double LeftLimitV, double RightLimitV) : 
ordu(Ordu), ordv(Ordv), 
leftlimitu(LeftLimitU), rightlimitu(RightLimitU), leftlimitv(LeftLimitV), rightlimitv(RightLimitV),
b(new Matrix<PolySurfBasisFunc>(Ordu,Ordv,PolySurfBasisFunc()))
{
	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++)
			(*b)[i][j] = PolySurfBasisFunc(Ordu,Ordv, i+1, j+1, leftlimitu,rightlimitu,leftlimitv,rightlimitv);
}


// EVALUATORS
// evaluate the basis Func by converting to BspSurf form
Matrix<double> PolySurfBasisFuncSet::Eval(double u, double v) const
{    
	Matrix<double> mat(ordu,ordv);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) mat[i][j] = (*b)[i][j](u,v);

	return mat;
}


// evaluate the basis Func by converting to BspSurf form
Matrix<double> PolySurfBasisFuncSet::operator()(double u, double v) const
{   
	Matrix<double> mat(ordu,ordv);

	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) mat[i][j] = (*b)[i][j](u,v);

	return mat;
}



// READ and WRITE
void PolySurfBasisFuncSet::write(std::ostream& os) 
{
	os << "Polyier Surf basis Func set\n";
	os << "order of basis Funcs in u is " << ordu << "\n";
	os << "order of basis Funcs in v is " << ordv << "\n";
	os << "number of basis Func in u  is " << ordu << "\n";
	os << "number of basis Func in v  is " << ordv << "\n";
	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			os << "\n" <<  i << j << "th Basis Func" << "\n";
			os << (*b)[i][j];
	}
}

void PolySurfBasisFuncSet::read(std::istream& is)
{
	int Ordu, Ordv;
	double Lu,Lv,Ru,Rv;
	std::cout << "order of Poly Surf";
	is >> Ordu >> Ordv;
	std::cout << "limits in u and v\n";
	is >> Lu >> Lv >> Ru >> Rv;
	*this = PolySurfBasisFuncSet(Ordu,Ordv,Lu,Lv,Ru,Rv);
} 


void PolySurfBasisFuncSet::writefile(std::ofstream& ofs)  
{
	ofs << "Polyier Surf basis Func set\n";
	ofs << "order of basis Funcs in u is " << ordu << "\n";
	ofs << "order of basis Funcs in v is " << ordv << "\n";
	ofs << "number of basis Func in u  is " << ordu << "\n";
	ofs << "number of basis Func in v  is " << ordv << "\n";
	for (int i=0; i<ordu; i++) 
		for (int j=0; j<ordv; j++) {
			ofs << "\n" <<  i << j << "th Basis Func" << "\n";
			(*b)[i][j].writefile(ofs);
	}	
}


void PolySurfBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv;
	double Lu,Lv,Ru,Rv;
	ifs >> Ordu >> Ordv;
	ifs >> Lu >> Lv >> Ru >> Rv;
	*this = PolySurfBasisFuncSet(Ordu,Ordv,Lu,Lv,Ru,Rv);
} 


