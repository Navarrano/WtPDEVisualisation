
#include "bezcurv.h"


 // CONSTRUCTORS

// default constructor
BezCurvBasisFunc::BezCurvBasisFunc() : ord(0), kts() {}

// constructor building a BezCurv basis function from 
// a Vector of knots and an order
BezCurvBasisFunc::BezCurvBasisFunc(const Vector<double>& Kts, int Ord) : 
	ord(Ord), kts(new Vector<double>(Kts)), b(new BezCurv<double>(CreateBezCurv())) 
{ 
}

// EVALUATORS

// Basis function evaluator
// Creates BezCurv representation of basis function and evaluates it
double BezCurvBasisFunc::Eval(double x) const
{    
	return (*b).Eval(x);
}

// Basis function evaluator
// Creates BezCurv representation of basis function and evaluates it
double BezCurvBasisFunc::operator()(double x) const
{    
	return (*b)(x);
}


// Basis function evaluator
// Creates BezCurv representation of basis function and evaluates it
double BezCurvBasisFunc::operator()(int val, double x) const
{    
	return (*b).Derive(val,x);
}


// ACCESS FUNCTIONS

BezCurv<double> BezCurvBasisFunc::GetBezCurv() const 
{
	return *b;
}


// compute the dimension of the basis function knot set
int BezCurvBasisFunc::ComputeDim() const
{
	return ord;
}



double BezCurvBasisFunc::GetLeftLimit() const
{
	return (*kts)[0];
}


double BezCurvBasisFunc::GetRightLimit() const
{
	return (*kts)[ord];
}

Vector<double> BezCurvBasisFunc::GetKnots() const
{
	return *kts;
}

int BezCurvBasisFunc::GetOrd() const
{
	return ord;
}

double BezCurvBasisFunc::Derive(int n, double x) const
{
	return (*b).Derive(n,x);
}



// create the BezCurv representation of the basis function
BezCurv<double> BezCurvBasisFunc::CreateBezCurv() const
{	
	Vector<double> cpts(ord,0.0);

	// create the knot set
	KnotSet kset(*kts,ord,ord+1);
	
	int num = kset.GetNumDistinct();
	
	// create Vectors of distinct knots and multiplicities
	Vector<double> dts(kset.GetDistinctKnots());
	Vector<int> mult(kset.GetMult());
	
	// set multiplicity to be ord at two ends
	mult[0]=ord;
	mult[num-1]=ord;

	// find offset
	int start=kset.GetMult()[0];
	// create knot set for curve
	KnotSet kset1(dts,mult,ord,num);
	// assign 1.0 to appropriate control point (rest=0.0)
	cpts[ord-start]=1.0;
	// create and return the BezCurv
	return BezCurv<double>(cpts,kset1.GetKnots(),ord);
}


// READ AND WRITE

void BezCurvBasisFunc::write(std::ostream& os) 
{
	os << "Bezier Curve Basis Function\n";
	os << "order is " << ord << "\n";
	os << "knots are\n";
	(*kts).write(os);	
	os << "\nBezCurv representation is\n";
	(*b).write(os);
}


void BezCurvBasisFunc::read(std::istream& is)
{
	int Ord;
	std::cout << "order of basis function\n";
	is >> Ord;
	Vector<double> Kts(Ord+1);
	std::cout << "input knots\n";
	is >> Kts;
	*this = BezCurvBasisFunc(Kts,Ord);
} 



void BezCurvBasisFunc::writefile(std::ofstream& ofs) 
{
	ofs << "Bezier basis function\n";
	ofs << "order is " << ord << "\n";
	ofs << "knots are\n";
	(*kts).writefile(ofs);
	ofs << "\nBezCurv representation is\n";
	(*b).writefile(ofs);
}


void BezCurvBasisFunc::readfile(std::ifstream &ifs)
{
	int Ord;
	ifs >> Ord;
	Vector<double> Kts(Ord+1);
	ifs >> Kts;
	*this = BezCurvBasisFunc(Kts,Ord);
} 







// CONSTRUCTORS
BezCurvBasisFuncSet::BezCurvBasisFuncSet() : ord(0), kts(), b() {}


BezCurvBasisFuncSet::BezCurvBasisFuncSet(const Vector<double>& Kts, int Ord) : ord(Ord), kts(new Vector<double>(Kts)), b(new Vector<BezCurvBasisFunc>(Ord,BezCurvBasisFunc()))
{
	Vector<double> v(Ord+1);
	for (int i=0; i<ord; i++) {
		for (int j=0; j<=ord; j++) v[j] = (*kts)[i+j];
		(*b)[i] = BezCurvBasisFunc(v,Ord);
	}
}


// EVALUATORS
// evaluate the basis Func by converting to BspCurv form
Vector<double> BezCurvBasisFuncSet::Eval(double x) const
{    
	Matrix<double> v(ord,ord,0.0);
	Vector<double> dp(ord,0.0), dm(ord,0.0);

//	int ind = kset.Find_index(x);
	double xnew = (x-(*kts)[ord-1])/((*kts)[ord]-(*kts)[ord-1]);

	v[0][0] = 1.0;
	for (int j=0; j<ord-1; j++) {
		dp[j] = 1.0-xnew;
		dm[j] = xnew;
		for (int i=0; i<=j; i++) {
			double m = v[i][j]/(dp[i]+dm[j-i]);
			v[i][j+1] = v[i][j+1]+dp[i]*m;
			v[i+1][j+1]=dm[j-i]*m;
		}
	}
	Vector<double> res(ord);
	for (int i=0; i<ord; i++) res[i] = v[i][ord-1];
	return res;
}


Vector<double> BezCurvBasisFuncSet::GetKnots() const
{
	return *kts;
}

int BezCurvBasisFuncSet::GetOrd() const
{
	return ord;
}
	
BezCurvBasisFunc BezCurvBasisFuncSet::GetBasisFunc(int i) const
{
	Vector<double> v(ord+1);

	for (int j=0; j<=ord; j++) v[j] = (*kts)[i-1+j]; 
	
	return BezCurvBasisFunc(v,ord);
}



// evaluate the basis Func by converting to BspCurv form
Vector<double> BezCurvBasisFuncSet::operator()(double x) const
{   
	Vector<double> v(ord);

	for (int j=0; j<ord; j++) v[j] = (*b)[j].Eval(x);
	return v;
}



// READ and WRITE
void BezCurvBasisFuncSet::write(std::ostream& os) 
{
	os << "Bezier Curve basis Func set\n";
	os << "order of basis Funcs is " << ord << "\n";
	os << "number of basis Funcs is " << ord << "\n";
	os << "knot set is\n";
	(*kts).write(os);
	for (int i=0; i<ord; i++) {
		os <<  i << "th Basis Func" << "\n";
		(*b)[i].write(os);
	}
}

void BezCurvBasisFuncSet::read(std::istream& is)
{
	int Ord;
	std::cout << "order of Bspline curve";
	is >> Ord;
	Vector<double> Kts(Ord+Ord);
	std::cout << "input knots";
	is >> Kts;
	*this = BezCurvBasisFuncSet(Kts,Ord);
} 


void BezCurvBasisFuncSet::writefile(std::ofstream& ofs) 
{
	ofs << "Bezier Curve basis Func set\n";
	ofs << "order of basis Funcs is " << ord << "\n";
	ofs << "number of basis Funcs is " << ord << "\n";
	ofs << "knot set is\n";
	(*kts).writefile(ofs);
	for (int i=0; i<ord; i++) {
		ofs <<  i << "th Basis Func" << "\n";
		(*b)[i].writefile(ofs);
	}
}


void BezCurvBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ord;
	ifs >> Ord;
	Vector<double> Kts(Ord+Ord);
	ifs >> Kts;
	*this = BezCurvBasisFuncSet(Kts,Ord);
} 

