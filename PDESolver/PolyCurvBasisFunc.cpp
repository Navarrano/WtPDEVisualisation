
#include "polycurv.h"

// CONSTRUCTORS

// default constructor
PolyCurvBasisFunc::PolyCurvBasisFunc() : ord(0), leftlimit(0), rightlimit(0) {}

// constructor taking an order, left and right limits
PolyCurvBasisFunc::PolyCurvBasisFunc(int Ord, int Index, double LeftLimit, double RightLimit) : 
	ord(Ord), index(Index), leftlimit(LeftLimit), rightlimit(RightLimit), b(new PolyCurv<double>(CreatePolyCurv())) { }


// EVALUATORS

// evaluate basis function at x
double PolyCurvBasisFunc::Eval(double x) const
{    
	return (*b).Eval(x);
}


// evaluate basis function at x
double PolyCurvBasisFunc::operator()(double x) const
{    
	return (*b)(x);
}

double PolyCurvBasisFunc::operator()(int val, double x) const
{    
	return (*b).Derive(val,x);
}


// ACCESS FUNCTIONS

PolyCurv<double> PolyCurvBasisFunc::GetPolyCurv() const
{
	return (*b);
}

// create the PolyCurv representation of basis function
PolyCurv<double> PolyCurvBasisFunc::CreatePolyCurv() const
{	
	Vector<double> coeffs(ord,0.0);

	coeffs[index-1]=1.0;
	return PolyCurv<double>(coeffs,ord,leftlimit,rightlimit);
}

double PolyCurvBasisFunc::GetLeftLimit() const
{
	return leftlimit;
}


double PolyCurvBasisFunc::GetRightLimit() const
{
	return rightlimit;
}

double PolyCurvBasisFunc::Derive(int n, double x) const
{
	return (*b).Derive(n,x);
}



// READ and WRITE

void PolyCurvBasisFunc::write(std::ostream& os) 
{
	os << "PolyCurv Basis Function\n";
	os << "order is " << ord << "\n";
	os << "left and right limits are\n";
	os << leftlimit << " " << rightlimit;
	os << "index is " << index;
	os << "\nPolyCurv representation is\n";
	(*b).write(os);
}
	
	
void PolyCurvBasisFunc::read(std::istream& is)
{
	int Ord;
	std::cout << "order of basis function\n";
	is >> Ord;
	std::cout << "input left and right limit\n";
	double left, right;
	is >> left >> right;
	std::cout << "input index\n";
	is >> index;
	*this = PolyCurvBasisFunc(Ord,index,left,right);
} 



void PolyCurvBasisFunc::writefile(std::ofstream& ofs) 
{
	ofs << "PolyCurv basis function\n";
	ofs << "order is " << ord << "\n";
	ofs << "left and right limit\n";
	ofs << leftlimit << " " << rightlimit;
	ofs << "index is " << index;
	ofs << "\nPolyCurv representation is\n";
	(*b).writefile(ofs);
}


void PolyCurvBasisFunc::readfile(std::ifstream& ifs)
{
	int Ord;
	ifs >> Ord;
	double left, right;
	ifs >> left >> right;
	int index;
	ifs >> index;
	*this = PolyCurvBasisFunc(Ord,index,left,right);
} 




// CONSTRUCTORS
PolyCurvBasisFuncSet::PolyCurvBasisFuncSet() : ord(0), b() {}

PolyCurvBasisFuncSet::PolyCurvBasisFuncSet(int Ord, double LeftLimit, double RightLimit) : 
ord(Ord), leftlimit(LeftLimit), rightlimit(RightLimit),
b(new Vector<PolyCurvBasisFunc>(Ord,PolyCurvBasisFunc()))
{
	for (int i=0; i<ord; i++) 
		(*b)[i] = PolyCurvBasisFunc(Ord,i+1,leftlimit,rightlimit);
}


// EVALUATORS
// evaluate the basis Func by converting to BspCurv form
Vector<double> PolyCurvBasisFuncSet::Eval(double x) const
{    
	Vector<double> v(ord);

	for (int j=0; j<ord; j++) v[j] = (*b)[j](x);
	return v;
}


// evaluate the basis Func by converting to BspCurv form
Vector<double> PolyCurvBasisFuncSet::operator()(double x) const
{   
	Vector<double> v(ord);

	for (int j=0; j<ord; j++) v[j] = (*b)[j].Eval(x);
	return v;
}



// READ and WRITE
void PolyCurvBasisFuncSet::write(std::ostream& os) 
{
	os << "Polyier Curve basis Func set\n";
	os << "order of basis Funcs is " << ord << "\n";
	os << "number of basis Funcs is " << ord << "\n";
	for (int i=0; i<ord; i++) {
		os <<  i << "th Basis Func" << "\n";
		(*b)[i].write(os);
	}
}

void PolyCurvBasisFuncSet::read(std::istream& is)
{
	int Ord;
	double L, R;
	std::cout << "order of Bspline curve";
	is >> Ord;
	std::cout << "limits\n";
	is >> L >> R;
	*this = PolyCurvBasisFuncSet(Ord,L,R);
} 


void PolyCurvBasisFuncSet::writefile(std::ofstream& ofs)  
{
	ofs << "Polyier Curve basis Func set\n";
	ofs << "order of basis Funcs is " << ord << "\n";
	ofs << "number of basis Funcs is " << ord << "\n";
	for (int i=0; i<ord; i++) {
		ofs << i << "th Basis Func" << "\n";
		(*b)[i].writefile(ofs);
	}
}


void PolyCurvBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ord;
	double L,R;
	ifs >> Ord;
	ifs >> L >> R;
	*this = PolyCurvBasisFuncSet(Ord,L,R);
} 


