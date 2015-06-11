
#include "bspcurv.h"
#include "mathfunctions.h"
#include "compbezcurv.h"

// CONSTRUCTORS

// default constructor
BspCurvBasisFunc::BspCurvBasisFunc() : ord(0), kts(), b() {}

// constructor taking an order and a Vector of knots
BspCurvBasisFunc::BspCurvBasisFunc(const Vector<double>& Kts, int Ord) : 
	ord(Ord), kts(new Vector<double>(Kts)), b(new BspCurv<double>(CreateBspCurv())) 
{

}


// ACCESS FUCNTIONS

// compute the dimension of the basis Func knot std::set
int BspCurvBasisFunc::ComputeDim() const
{
	// what about the case when num = 1?
	// create KnotSet object
	KnotSet kset(*kts,ord,ord+1);
	int num = kset.GetNumDistinct();
	int num_cond=0;
	// compute dimension
	for (int i=1; i<num-1; i++) num_cond = num_cond + ord - kset.GetMult()[i];
	return ((num-1)*ord-num_cond);
}

BspCurv<double> BspCurvBasisFunc::GetBspCurv() const
{
	return *b;
}


Vector<double> BspCurvBasisFunc::GetKnots() const
{
	return *kts;
}

int BspCurvBasisFunc::GetOrd() const
{
	return ord;
}

// create the BspCurv representation of the basis Func
BspCurv<double> BspCurvBasisFunc::CreateBspCurv() const
{
	int dim = ComputeDim();	
	
	Vector<double> cpts(dim,0.0);

	// create KnotSet object
	KnotSet kset(*kts,ord,ord+1);
	
	int num = kset.GetNumDistinct();
	
	// create temporary Vectors storing multiplicities and distinct knots
	Vector<double> dts(kset.GetDistinctKnots());
	Vector<int> mult(kset.GetMult());
	
	// std::set multiplicities at the end equal to ord
	mult[0]=ord;
	mult[num-1]=ord;

	// compute offstd::set
	int start=kset.GetMult()[0];
	// create KnotSet object for curve
	KnotSet kset1(dts,mult,ord,num);
	// assign 1.0 to approriate control point
	cpts[ord-start]=1.0;
	// create and return curve
	return BspCurv<double>(cpts,kset1.GetKnots(),ord,dim);
}



// EVALUATORS

// evaluate the basis Func by converting to BspCurv form
double BspCurvBasisFunc::Eval(double x) const
{    
	if (x < (*kts)[0] || x > (*kts)[ord]) return 0;
	else return (*b).Eval(x);
}

// evaluate the basis Func by converting to BspCurv form
double BspCurvBasisFunc::operator()(double x) const
{    
	if (x < (*kts)[0] || x > (*kts)[ord]) return 0;
	return (*b)(x);
}

// evaluate the basis Func by converting to BspCurv form
double BspCurvBasisFunc::operator()(int val, double x) const
{    
	return (*b).Derive(val,x);
}

double BspCurvBasisFunc::CreateIntegral(double x1, double x2) const
{
	return GetBspCurv().Integrate(x1,x2);
}


double BspCurvBasisFunc::GetLeftLimit() const
{
	return (*kts)[0];
}


double BspCurvBasisFunc::GetRightLimit() const
{
	return (*kts)[ord];
}

double BspCurvBasisFunc::Derive(int n, double x) const
{
	return (*b).Derive(n,x);
}



// READ and WRITE
void BspCurvBasisFunc::write(std::ostream& os)
{
	os << "Bspline Curve basis Func\n";
	os << "order is " << ord << "\n";
	os << "number of control points is " << ord << "\n";
	os << "knots are\n";
	os << *kts;
	os << "\nBspCurv representation\n";
	(*b).write(os);
}

void BspCurvBasisFunc::read(std::istream& is)
{
	int Ord;
	std::cout << "order of Bspline curve";
	is >> Ord;
	Vector<double> Kts(Ord+1);
	std::cout << "input knots";
	is >> Kts;
	*this = BspCurvBasisFunc(Kts,Ord);
} 


void BspCurvBasisFunc::writefile(std::ofstream& ofs) 
{
	ofs << "Bspline Curve\n";
	ofs << "order is " << ord << "\n";
	ofs << "knots are\n";
	ofs << *kts;
	ofs << "\nBspCurv representation\n";
	(*b).writefile(ofs);
}


void BspCurvBasisFunc::readfile(std::ifstream& ifs)
{
	int Ord;
	ifs >> Ord;
	Vector<double> Kts(Ord+1);
	ifs >> Kts;
	*this = BspCurvBasisFunc(Kts,Ord);
} 



// CONSTRUCTORS
BspCurvBasisFuncSet::BspCurvBasisFuncSet() : ord(0), num(0), kts(), b()  {}

BspCurvBasisFuncSet::BspCurvBasisFuncSet(const Vector<double>& Kts, int Ord, int Num) : kts(new Vector<double>(Kts)), ord(Ord), num(Num), 
b(new Vector<BspCurvBasisFunc>(Num-Ord,BspCurvBasisFunc())) 
{
	Vector<double> v(Ord+1);
	for (int i=0; i<Num-Ord; i++) {
		for (int j=0; j<=Ord; j++) v[j] = (*kts)[i+j]; 
		(*b)[i] = BspCurvBasisFunc(v,Ord);
	}
}

BspCurvBasisFuncSet::BspCurvBasisFuncSet(const KnotSet& kset) : kts(new Vector<double>(kset.GetKnots())), ord(kset.GetOrd()), num(kset.GetNum()),b(new Vector<BspCurvBasisFunc>(num-ord,BspCurvBasisFunc())) 
{
	Vector<double> v(ord+1);
	for (int i=0; i<num-ord; i++) {
		for (int j=0; j<=ord; j++) v[j] = (*kts)[i+j]; 
		(*b)[i] = BspCurvBasisFunc(v,ord);
	}

}


Vector<double> BspCurvBasisFuncSet::Eval(double x) const
{
	Vector<double> v(num-ord);

	for (int i=0; i<num-ord; i++) v[i] = (*b)[i].Eval(x);

	return v;
}


// EVALUATORS
// evaluate the basis Func by converting to BspCurv form
Vector<double> BspCurvBasisFuncSet::EvalNonZero(double x) const
{    
	Matrix<double> v(ord,ord,0.0);
	Vector<double> dp(ord,0.0), dm(ord,0.0);

	int ind = KnotSet(*kts,ord,num).Find_index(x);
	

	v[0][0] = 1.0;
	for (int j=0; j<ord-1; j++) {
		dp[j] = (*kts)[ind+j+1]-x;
		dm[j] = x - (*kts)[ind+1-j-1];
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



Vector<double> BspCurvBasisFuncSet::GetKnots() const
{
	return *kts;
}

int BspCurvBasisFuncSet::GetOrd() const
{
	return ord;
}
	
BspCurvBasisFunc BspCurvBasisFuncSet::GetBasisFunc(int i) const
{
	Vector<double> v(ord+1);

	for (int j=0; j<=ord; j++) v[j] = (*kts)[i-1+j]; 
	
	return BspCurvBasisFunc(v,ord);
}

Vector<double> BspCurvBasisFuncSet::CreateIntegralNew(const BspCurv<double>& c, double x1, double x2) const
{
	Vector<double> res(num-ord,0.0);

	for (int j=0; j<num-ord; j++) {
		for (int i=1; i<=c.GetOrd(); i++) {
			double y1=x1;
			double y2=x2;
	
			BspCurv<double> curve = (*b)[j].GetBspCurv().Integrate(i);
	
			KnotSet kset = curve.GetKnotSet();
			if (x1 < kset.GetKnots()[curve.GetOrd()-1]) y1 = kset.GetKnots()[curve.GetOrd()-1];
			if (x2 > kset.GetKnots()[curve.GetNum()]) y2 = kset.GetKnots()[curve.GetNum()];
			BspCurv<double> c1 = c.Derive(i-1);
		
			res[j] = res[j] + pow(-1.0,(double)(i-1))*curve(y2)*c1(curve.GetRightLimit())-curve(y1)*c1(curve.GetLeftLimit());
		}
	}


	return res;
}





Vector<double> BspCurvBasisFuncSet::CreateIntegralProduct(const Curve<double>& c, double x1, double x2) const
{
	Vector<double> res(num-ord);


	for (int i=0; i<num-ord; i++) {
		double y1=x1;
		double y2=x2;
		BspCurv<double> curve = (*b)[i].GetBspCurv().Integrate();
		KnotSet kset = curve.GetKnotSet();
		if (x1 < kset.GetKnots()[curve.GetOrd()-1]) y1 = kset.GetKnots()[curve.GetOrd()-1];
		if (x2 > kset.GetKnots()[curve.GetNum()]) y2 = kset.GetKnots()[curve.GetNum()];
		res[i] = curve(y2)*c(curve.GetRightLimit())-curve(y1)*c(curve.GetLeftLimit());
	}

	return res;
}

Vector<double> BspCurvBasisFuncSet::CreateIntegralIntegralProduct(const Curve<double>& c, double x1, double x2) const
{
	Vector<double> res(num-ord);


	for (int i=0; i<num-ord; i++) {
		double y1=x1;
		double y2=x2;
		BspCurv<double> curve = (*b)[i].GetBspCurv().Integrate().Integrate();
		KnotSet kset = curve.GetKnotSet();
		if (x1 < kset.GetKnots()[curve.GetOrd()-1]) y1 = kset.GetKnots()[curve.GetOrd()-1];
		if (x2 > kset.GetKnots()[curve.GetNum()]) y2 = kset.GetKnots()[curve.GetNum()];
		res[i] = curve(y2)*c(curve.GetRightLimit())-curve(y1)*c(curve.GetLeftLimit());
	}

	return res;
}

Vector<double> BspCurvBasisFuncSet::CreateIntegral(double x1, double x2) const
{
	Vector<double> res(num-ord);

	for (int i=0; i<num-ord; i++) 
		res[i] = (*b)[i].GetBspCurv().Integrate(x1,x2);

	return res;
}


Vector<double> BspCurvBasisFuncSet::CreateIntegralIntegral(double x1, double x2) const
{
	Vector<double> res(num-ord);

	for (int i=0; i<num-ord; i++) 
		res[i] = (*b)[i].GetBspCurv().Integrate().Integrate(x1,x2);

	return res;
}

Vector<double> BspCurvBasisFuncSet::CreateIntegralIntegralIntegral(double x1, double x2) const
{
	Vector<double> res(num-ord);

	for (int i=0; i<num-ord; i++) 
		res[i] = (*b)[i].GetBspCurv().Integrate().Integrate().Integrate(x1,x2);

	return res;
}


// evaluate the basis Func by converting to BspCurv form
Matrix<double> BspCurvBasisFuncSet::ComputeNMatrix(double x) const
{    
	Matrix<double> mat(ord,ord,0.0);
	Vector<double> dp(ord,0.0), dm(ord,0.0);

	int ind = KnotSet(*kts,ord,num).Find_index(x);


	mat[0][0] = 1.0;
	for (int j=0; j<ord-1; j++) {
		dp[j] = (*kts)[ind+j+1]-x;
		dm[j] = x - (*kts)[ind-j];
		for (int i=0; i<=j; i++) {
			double m = mat[i][j]/(dp[i]+dm[j-i]);
			mat[i][j+1] = mat[i][j+1]+dp[i]*m;
			mat[i+1][j+1]=dm[j-i]*m;
		}
	}
	return mat;
}


// evaluate the basis Func by converting to BspCurv form
Vector<double> BspCurvBasisFuncSet::operator()(double x) const
{   
	Vector<double> v(num-ord);
	for (int j=0; j<num-ord; j++) v[j] = (*b)[j](x); 
	return v;
}

Matrix<double> BspCurvBasisFuncSet::CreateMatrixKnotAverages() const
{
	Matrix<double> mat(num-ord,num-ord);
	Vector<double> knotav = KnotSet(*kts,ord,num).ComputeKnotSetAver();

	for (int i=0; i<num-ord; i++)
		for (int j=0; j<num-ord; j++) mat[i][j] = (*b)[j](knotav[i]);

	return mat;
}

Matrix<double> BspCurvBasisFuncSet::CreateMatrixKnotAverages(Vector<double>& vec) const
{
	Matrix<double> mat(vec.GetNum(),num-ord);
	//Vector<double> knotav = KnotSet(*kts,ord,num).ComputeKnotSetAver();

	for (int i=0; i<vec.GetNum(); i++)
		for (int j=0; j<num-ord; j++) mat[i][j] = (*b)[j](vec[i]);

	return mat;
}



Matrix<double> BspCurvBasisFuncSet::CreateMatrixMinimisation(int lev, double x1, double x2) const
{
	Matrix<double> mat = CreateMatrixIntegral(lev, x1, x2);

	Matrix<double> deriv = KnotSet(*kts,ord,num).CreateMatrixDeriv(lev);
	
	Matrix<double> derivt = Math::transpose(deriv);

	return Math::mult1(Math::mult1(derivt,mat),deriv);
}


Matrix<double> BspCurvBasisFuncSet::CreateMatrixMinimisation(int lev) const
{
	Matrix<double> mat = CreateMatrixIntegral(lev);

	Matrix<double> deriv = KnotSet(*kts,ord,num).CreateMatrixDeriv(lev);
	
	Matrix<double> derivt = Math::transpose(deriv);

	return Math::mult1(Math::mult1(derivt,mat),deriv);
}

Matrix<double> BspCurvBasisFuncSet::CreateMatrixMinimisation(int lev1, int lev2) const
{
	Matrix<double> mat = CreateMatrixIntegral(lev1,lev2);

	Matrix<double> deriv1 = KnotSet(*kts,ord,num).CreateMatrixDeriv(lev1);
	Matrix<double> deriv2 = KnotSet(*kts,ord,num).CreateMatrixDeriv(lev2);
	
	Matrix<double> deriv1t = Math::transpose(deriv1);

	return Math::mult1(Math::mult1(deriv1t,mat),deriv2);
}

Matrix<double> BspCurvBasisFuncSet::CreateMatrixMinimisation(int lev1, int lev2, double x1, double x2) const
{
	Matrix<double> mat = CreateMatrixIntegral(lev1,lev2,x1,x2);

	Matrix<double> deriv1 = KnotSet(*kts,ord,num).CreateMatrixDeriv(lev1);
	Matrix<double> deriv2 = KnotSet(*kts,ord,num).CreateMatrixDeriv(lev2);
	
	Matrix<double> deriv1t = Math::transpose(deriv1);

	return Math::mult1(Math::mult1(deriv1t,mat),deriv2);
}

Matrix<double> BspCurvBasisFuncSet::CreateMatrixIntegral(int lev1, int lev2, double x1, double x2) const
{
	KnotSet kset1 = KnotSet(*kts,ord,num).CreateKnotSetDeriv(lev1);
	KnotSet kset2 = KnotSet(*kts,ord,num).CreateKnotSetDeriv(lev2);
	
	Matrix<double> mat(kset1.GetNum()-(ord-lev1),kset2.GetNum()-(ord-lev2));
	
	BspCurvBasisFuncSet basis1(kset1.GetKnots(),ord-lev1,kset1.GetNum());
	BspCurvBasisFuncSet basis2(kset2.GetKnots(),ord-lev2,kset2.GetNum());

	// find the first basis function for which the last distinct
	// knot is greater than x1

	int ind1=-1;
	do {
		ind1++;
	} while ((*(basis1.b))[ind1].GetKnots()[ord-lev1] <= x1);


	// find the last basis function for which the first distinct
	// knot is less than x2

	int ind2=-1;
	do {
		ind2++;
	} while (ind2 < kset1.GetNum()-ord+lev1 && (*(basis1.b))[ind2].GetKnots()[0] < x2);


	int ind3=-1;
	do {
		ind3++;
	} while ((*(basis2.b))[ind3].GetKnots()[ord-lev2] <= x1);


	// find the last basis function for which the first distinct
	// knot is less than x2

	int ind4=-1;
	do {
		ind4++;
	} while (ind4 < kset2.GetNum()-ord+lev2 && (*(basis2.b))[ind2].GetKnots()[0] < x2);


	Matrix<double> mat1(kset1.GetNum()-ord+lev1,kset2.GetNum()-ord+lev2,0.0);

	int i1, i2;
	if (ind1 < ind3) i1 = ind1; else i1 = ind3;
	if (ind2 > ind4) i2 = ind2; else i2 = ind4;

	for (int i=i1; i<=i2-1; i++)
		for (int j=i1; j<=i2-1; j++) {
		
			
			// create the two std::sets representing the two knot std::sets
			Vector<double> temp1((*(basis1.b))[i].GetKnots());
			Vector<double> temp2((*(basis2.b))[j].GetKnots());
			std::set<double> s1(temp1.begin(),temp1.end());
			std::set<double> s2(temp2.begin(),temp2.end());

			if (*(--s2.end()) > *(s1.begin())) {
				// form the intersection
				std::set<double> s3;
				std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(s3,s3.begin()));
			
				// if there is an intersection
				if (s3.size() > 1) {
					Vector<double> v(s3.size());
					std::set<double>::iterator s = s3.begin();

					// copy the elements into a vector
					for (unsigned int k=0; k<s3.size(); k++) v[k] = *s++;
				
					// create the compbezcurvs
					Vector<BezCurv<double> > vec1(s3.size()-1);
					Vector<BezCurv<double> > vec2(s3.size()-1);

					BspCurv<double> b1((*(basis1.b))[i].GetBspCurv()), b2((*(basis2.b))[j].GetBspCurv());
				
					// find the segments of intersection
					for (unsigned int k=0; k<s3.size()-1; k++) {
						int segb1 = b1.GetKnotSet().Find_segment(v[k]);
						int segb2 = b2.GetKnotSet().Find_segment(v[k]);
						
						vec1[k] = b1.GetSegment(segb1);
						vec2[k] = b2.GetSegment(segb2);
					}
				
					CompBezCurv<double> cb1(vec1,s3.size()-1,v);
					CompBezCurv<double> cb2(vec2,s3.size()-1,v);
					CompBezCurv<double> prod = cb1.Product(cb2);
				
					mat1[i][j] = prod.ConvertBspCurv().Integrate(x1,x2);
				}
			}
		}
	

	return mat1;
}




Matrix<double> BspCurvBasisFuncSet::CreateMatrixIntegral(int lev, double x1, double x2) const
{
	

	KnotSet kset = KnotSet(*kts,ord,num).CreateKnotSetDeriv(lev);
	
	Matrix<double> mat(kset.GetNum()-(ord-lev),kset.GetNum()-(ord-lev));
	
	BspCurvBasisFuncSet basis(kset.GetKnots(),ord-lev,kset.GetNum());

	// find the first basis function for which the last distinct
	// knot is greater than x1

	int ind1=-1;
	do {
		ind1++;
	} while ((*(basis.b))[ind1].GetKnots()[ord-lev] <= x1);


	// find the last basis function for which the first distinct
	// knot is less than x2

	int ind2=-1;
	do {
		ind2++;
	} while (ind2 < kset.GetNum()-ord+lev && (*(basis.b))[ind2].GetKnots()[0] < x2);


	Matrix<double> mat1(kset.GetNum()-ord+lev,kset.GetNum()-ord+lev,0.0);

	for (int i=ind1; i<=ind2-1; i++)
		for (int j=ind1; j<=i; j++) {
		
			
			// create the two std::sets representing the two knot std::sets
			Vector<double> temp1((*(basis.b))[i].GetKnots());
			Vector<double> temp2((*(basis.b))[j].GetKnots());
			std::set<double> s1(temp1.begin(),temp1.end());
			std::set<double> s2(temp2.begin(),temp2.end());

			if (*(--s2.end()) > *(s1.begin())) {
				// form the intersection
				std::set<double> s3;
				std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(s3,s3.begin()));
			
				// if there is an intersection
				if (s3.size() > 1) {
					Vector<double> v(s3.size());
					std::set<double>::iterator s = s3.begin();

					// copy the elements into a vector
					for (unsigned int k=0; k<s3.size(); k++) v[k] = *s++;
				
					// create the compbezcurvs
					Vector<BezCurv<double> > vec1(s3.size()-1);
					Vector<BezCurv<double> > vec2(s3.size()-1);

					BspCurv<double> b1((*(basis.b))[i].GetBspCurv()), b2((*(basis.b))[j].GetBspCurv());
				
					// find the segments of intersection
					for (unsigned int k=0; k<s3.size()-1; k++) {
						int segb1 = b1.GetKnotSet().Find_segment(v[k]);
						int segb2 = b2.GetKnotSet().Find_segment(v[k]);
						
						vec1[k] = b1.GetSegment(segb1);
						vec2[k] = b2.GetSegment(segb2);
					}
				
					CompBezCurv<double> cb1(vec1,s3.size()-1,v);
					CompBezCurv<double> cb2(vec2,s3.size()-1,v);
					CompBezCurv<double> prod = cb1.Product(cb2);
				
					mat1[i][j] = prod.ConvertBspCurv().Integrate(x1,x2);
				}
			}
		}
	
	for (int i=ind1; i<=ind2-2; i++)
		for (int j=i+1; j<=ind2-1; j++) mat1[i][j] = mat1[j][i];
	
	return mat1;
}

Matrix<double> BspCurvBasisFuncSet::CreateMatrixIntegral(int lev1, int lev2) const
{
	KnotSet kset1 = KnotSet(*kts,ord,num).CreateKnotSetDeriv(lev1);
	KnotSet kset2 = KnotSet(*kts,ord,num).CreateKnotSetDeriv(lev2);
	
	Matrix<double> mat(kset1.GetNum()-(ord-lev1),kset2.GetNum()-(ord-lev2));
	
	BspCurvBasisFuncSet basis1(kset1.GetKnots(),ord-lev1,kset1.GetNum());
	BspCurvBasisFuncSet basis2(kset2.GetKnots(),ord-lev2,kset2.GetNum());

	// find the first basis function for which the last distinct
	// knot is greater than x1



	Matrix<double> mat1(kset1.GetNum()-ord+lev1,kset2.GetNum()-ord+lev2,0.0);



	for (int i=0; i<kset1.GetNum()-ord+lev1; i++)
		for (int j=0; j<kset2.GetNum()-ord+lev2; j++) {
		
			
			// create the two std::sets representing the two knot std::sets
			Vector<double> temp1((*(basis1.b))[i].GetKnots());
			Vector<double> temp2((*(basis2.b))[j].GetKnots());
			std::set<double> s1(temp1.begin(),temp1.end());
			std::set<double> s2(temp2.begin(),temp2.end());

			if (*(--s2.end()) > *(s1.begin())) {
				// form the intersection
				std::set<double> s3;
				std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(s3,s3.begin()));
			
				// if there is an intersection
				if (s3.size() > 1) {
					Vector<double> v(s3.size());
					std::set<double>::iterator s = s3.begin();

					// copy the elements into a vector
					for (unsigned int k=0; k<s3.size(); k++) v[k] = *s++;
				
					// create the compbezcurvs
					Vector<BezCurv<double> > vec1(s3.size()-1);
					Vector<BezCurv<double> > vec2(s3.size()-1);

					BspCurv<double> b1((*(basis1.b))[i].GetBspCurv()), b2((*(basis2.b))[j].GetBspCurv());
				
					// find the segments of intersection
					for (unsigned int k=0; k<s3.size()-1; k++) {
						int segb1 = b1.GetKnotSet().Find_segment(v[k]);
						int segb2 = b2.GetKnotSet().Find_segment(v[k]);
						
						vec1[k] = b1.GetSegment(segb1);
						vec2[k] = b2.GetSegment(segb2);
					}
				
					CompBezCurv<double> cb1(vec1,s3.size()-1,v);
					CompBezCurv<double> cb2(vec2,s3.size()-1,v);
					CompBezCurv<double> prod = cb1.Product(cb2);
				
					mat1[i][j] = prod.ConvertBspCurv().Integrate((*kts)[ord-1],(*kts)[num-ord]);
				}
			}
		}
	
	
	
	return mat1;
}


Matrix<double> BspCurvBasisFuncSet::CreateMatrixIntegral(int lev) const
{

	KnotSet kset = KnotSet(*kts,ord,num).CreateKnotSetDeriv(lev);
	
	Matrix<double> mat(kset.GetNum()-(ord-lev),kset.GetNum()-(ord-lev));
	
	BspCurvBasisFuncSet basis(kset.GetKnots(),ord-lev,kset.GetNum());

	Matrix<double> mat1(kset.GetNum()-ord+lev,kset.GetNum()-ord+lev,0.0);
	for (int i=0; i<kset.GetNum()-ord+lev; i++)
		for (int j=0; j<=i; j++) {
			// create the two std::sets representing the two knot std::sets
			
			Vector<double> temp1((*(basis.b))[i].GetKnots());
			Vector<double> temp2((*(basis.b))[j].GetKnots());
			std::set<double> s1(temp1.begin(),temp1.end());
			std::set<double> s2(temp2.begin(),temp2.end());
			
			// if there is an intersection
			if (*(--s2.end()) > *(s1.begin())) {
				// form the intersection
				std::set<double> s3;
				std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(s3,s3.begin()));
			
				// if there is an intersection
				if (s3.size() > 1) {
					Vector<double> v(s3.size());
					std::set<double>::iterator s = s3.begin();

					// copy the elements into a vector
					for (unsigned int k=0; k<s3.size(); k++) v[k] = *s++;

					// create the compbezcurvs
					Vector<BezCurv<double> > vec1(s3.size()-1);
					Vector<BezCurv<double> > vec2(s3.size()-1);

					BspCurv<double> b1((*(basis.b))[i].GetBspCurv()), b2((*(basis.b))[j].GetBspCurv());
					// find the segments of intersection
					for (unsigned int k=0; k<s3.size()-1; k++) {
						int segb1 = b1.GetKnotSet().Find_segment(v[k]);
						int segb2 = b2.GetKnotSet().Find_segment(v[k]);
						vec1[k] = b1.GetSegment(segb1);
						vec2[k] = b2.GetSegment(segb2);
					}
				
					CompBezCurv<double> cb1(vec1,s3.size()-1,v);
					CompBezCurv<double> cb2(vec2,s3.size()-1,v);
					CompBezCurv<double> prod = cb1.Product(cb2);
	
					mat1[i][j] = prod.ConvertBspCurv().Integrate((*kts)[ord-1],(*kts)[num-ord]);
				}
			}
		}
	for (int i=0; i<kset.GetNum()-ord+lev-1; i++)
		for (int j=i+1; j<kset.GetNum()-ord+lev; j++) mat1[i][j] = mat1[j][i];
	
	return mat1;
}

Matrix<double> BspCurvBasisFuncSet::ComputeBasisMatrix() const
{
	Matrix<double> mat(num-ord,num-ord);
	Vector<double> u(KnotSet((*kts),ord,num).ComputeKnotSetAver());

	for (int i=0; i<num-ord; i++)
			mat[i] = ComputeVecBasis(u[i]);

	return mat;
}


Vector<double> BspCurvBasisFuncSet::ComputeVecBasis(double x1) const
{
	Vector<double> v(num-ord);

	for (int i=0; i<num-ord; i++) v[i] = (*b)[i](x1);


	return v;
}


// READ and WRITE
void BspCurvBasisFuncSet::write(std::ostream& os)
{
	os << "Bspline Curve basis Func std::set\n";
	os << "order of basis Funcs is " << ord << "\n";
	os << "number of basis Funcs is " << num-ord << "\n";
	os << "knot std::set is\n";
	os << *kts;
	for (int i=0; i<num-ord; i++) {
		os <<  i << "th Basis Func" << "\n";
		os << (*b)[i];
	}
}

void BspCurvBasisFuncSet::read(std::istream& is)
{
	int Ord;
	std::cout << "order of Bspline curve\n";
	is >> Ord;
	int Num;
	std::cout << "number of knots\n";
	is >> Num;

	Vector<double> Kts(Num);
	std::cout << "input knots\n";
	is >> Kts;
	*this = BspCurvBasisFuncSet(Kts,Ord,Num);
} 


void BspCurvBasisFuncSet::writefile(std::ofstream& ofs) 
{
	ofs << "Bspline Curve basis Func std::set\n";
	ofs << "order of basis Funcs is " << ord << "\n";
	ofs << "number of basis Funcs is " << num-ord << "\n";
	ofs << "knot std::set is\n";
	ofs << *kts;
	for (int i=0; i<num-ord; i++) {
		ofs << i << "th Basis Func" << "\n";
		(*b)[i].writefile(ofs);
	}
}


void BspCurvBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ord, Num;
	ifs >> Ord;
	ifs >> Num;
	Vector<double> Kts(Num);
	ifs >> Kts;
	*this = BspCurvBasisFuncSet(Kts,Ord,Num);
} 


