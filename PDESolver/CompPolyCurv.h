#ifndef COMPOLYCURV
#define COMPOLYCURV


// can be of different orders, IGES spec

#include "vector.h"
#include "bspcurv.h"
#include "compbezcurv.h"


template<class T>
class CompPolyCurv : public Curve<T> {
private:
	// data
	int num; // number of segments
	int ord; // max order of the segments
	Ptr<Vector<int> > ords; // orders of segments
	Ptr<Vector<T> > coeffs; // coefficients
	Ptr<Vector<double> > limits; // range of the segments
	Ptr<std::set<double> > limitset;

	// private functions
	CompPolyCurv<T> Product2(const CompPolyCurv<T>& c) const;
	template<class T1>
	CompPolyCurv<T> MakeCompatable(const CompPolyCurv<T1>& b) const
	{
		return (ConvertBspCurv().MakeBreakCompatable(b.ConvertBspCurv())).ConvertCompPolyCurv();
	}
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class CompPolyCurv<double>");
		else {
			std::string s(typeid(T).name()), s1("class CompPolyCurv<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	CompPolyCurv();
	CompPolyCurv(int Num, int Ord, const Vector<T>& Coeffs, const Vector<double>& Limits);
	CompPolyCurv(int Num, int Ord, const Vector<T>& Coeffs);
	CompPolyCurv(int Num, const Vector<int>& Ords, const Vector<T>& Coeffs, const Vector<double>& Limits);
	CompPolyCurv(int Num, const Vector<int>& Ords, const Vector<T>& Coeffs);
	CompPolyCurv(const Vector<PolyCurv<T> >& vec, int nseg, const Vector<double>& Limits);
	CompPolyCurv(const Vector<PolyCurv<T> >& vec, int nseg);
	
	// access functions
	int GetOrd() const;
	Vector<double> GetLimits() const;
	virtual double GetLeftLimit(int index) const;
	virtual double GetRightLimit(int index) const;
	Vector<T> GetCoeffsPoly(int index) const;
	int FindSegment(double x) const;
	double GetLeftLimit() const;
	double GetRightLimit() const;
	int GetNum() const;
	Vector<T> GetCoeffs() const;
	Vector<int> GetOrds() const;
	PolyCurv<T> GetSegment(int seg) const;
	PolyCurv<T> GetSegmentOriginal(int seg) const;

	// evaluators
	virtual T operator()(double x) const;
	T operator()(double x, int index) const;
	virtual T operator()(int, double x) const;
	T Eval(double x) const;
	T Eval(double x, int index) const;
	Vector<T> ComputePoints(int n) const;

	// conversion
	Vector<T> Convert() const;
	CompBezCurv<T> ConvertCompBezCurv() const;
	BspCurv<T> ConvertBspCurv() const;

	// derivatives
	CompPolyCurv<T> Derive(int level) const; 
	virtual T Derive(int level, double x) const;

	// degree elevation
	CompPolyCurv<T> Elevate(int level) const;

	// integration
	CompPolyCurv<T> Integrate() const;
	T Integrate(double x1, double x2) const;
	Vector<T> IntegrateCPoints() const;

	// product
	template<class T1>
	CompPolyCurv<T> Product(const CompPolyCurv<T1>& b) const
	{
		CompPolyCurv<T> d(*this), e(b);

		if (num != b.GetNum()) {
			d = MakeCompatable(b);
			e = b.MakeCompatable(d);
		}

		// create array for new coeffs
		Vector<T> ncoeffs(d.GetNum()*(ord+b.GetOrd()-1));
		int count=0;

		// multiply corresponding segments
		for (int i=0; i<d.GetNum(); i++) {
			PolyCurv<T> p1 = d.GetSegment(i+1);
			PolyCurv<T> p2 = e.GetSegment(i+1);
			Vector<T> temp = p1.ProductCPoints(p2);
			for (int j=0; j<ord+b.GetOrd()-1; j++) {
				ncoeffs[count] = temp[j];
				count++;
			}
		}
		Vector<int> Ords(d.GetNum());
		for (int i=0; i<d.GetNum(); i++) Ords[i] = (*ords)[i]+b.GetOrds()[i]-1;
	
		// return composite poly
		return CompPolyCurv<T>(d.GetNum(),Ords,ncoeffs,*limits);
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
CompPolyCurv<T>::CompPolyCurv() : num(0), ord(0), ords(), coeffs(), limits() 
{
}

template<class T>
CompPolyCurv<T>::CompPolyCurv(const Vector<PolyCurv<T> >& vec, int nseg, const Vector<double>& Limits) : num(nseg), ords(new Vector<int>(nseg)), limits(new Vector<double>(Limits))
{
	int sum = 0;
	for (int i=0; i<num; i++) (*ords)[i] = vec[i].GetOrd();
	for (int i=0; i<num; i++) sum = sum+ (*ords)[i];
	coeffs = new Vector<T>(sum);

	int count = 0;
	for (int i=0; i<num; i++) 
		for (int j=0; j<(*ords)[i]; j++) {
			(*coeffs)[count] = vec[i].GetCoeffs()[j];
			count++;
		}

	// find max order
	ord = (*ords)[0];
	for (int i=1; i<num; i++)
		if (ord < (*ords)[i]) ord=(*ords)[i];
	coeffs = new Vector<T>(Convert());
	limitset = new std::set<double>((*limits).begin(),(*limits).end());
}


template<class T>
CompPolyCurv<T>::CompPolyCurv(const Vector<PolyCurv<T> >& vec, int nseg) : num(nseg), ords(new Vector<int>(nseg)), limits(new Vector<double>(nseg+1))
{
	int sum = 0;
	for (int i=0; i<num; i++) (*ords)[i] = vec[i].GetOrd();
	for (int i=0; i<=num; i++) (*limits)[i]=(double)i;
	for (int i=0; i<num; i++) sum = sum+ (*ords)[i];
	coeffs = new Vector<T>(sum);
	int count = 0;
	for (int i=0; i<num; i++) 
		for (int j=0; j<(*ords)[i]; j++) {
			(*coeffs)[count] = vec[i].GetCoeffs()[j];
			count++;
		}
	// find max order
	ord = (*ords)[0];
	for (int i=1; i<num; i++) if (ord < (*ords)[i]) ord=(*ords)[i]; 
	coeffs = neww Vector<T>(Convert());
	limitset = new std::set<double>((*limits).begin(),(*limits).end());
}


// constructor builds a CompPolyCurv from a Vector of control points, 
// an order and number of control points and limits
template<class T>
CompPolyCurv<T>::CompPolyCurv(int Num, const Vector<int>& Ords, const Vector<T>& Coeffs, const Vector<double>& Limits) :
num(Num), ords(new Vector<int>(Ords)), coeffs(new Vector<T>(Coeffs)), limits(new Vector<double>(Limits))
{
	// find max order
	ord = (*ords)[0];
	for (int i=1; i<num; i++) 
		if (ord < (*ords)[i]) ord=(*ords)[i]; 
	coeffs = new Vector<T>(Convert());
	limitset = new std::set<double>((*limits).begin(),(*limits).end());
}


// constructor builds a CompPolyCurv from a Vector of control points, 
// an order and number of control points
template<class T>
CompPolyCurv<T>::CompPolyCurv(int Num, const Vector<int>& Ords, const Vector<T>& Coeffs) :
num(Num), ords(new Vector<int>(Ords)), coeffs(new Vector<T>(Coeffs)), limits(new Vector<double>(Num+1))
{
	// find max order
	ord = (*ords)[0];
	for (int i=1; i<num; i++) 
		if (ord < (*ords)[i]) ord=(*ords)[i];

	// limits are increasing integers
	for (int i=0; i<num+1; i++) (*limits)[i] = (double)i;
	coeffs = new Vector<T>(Convert());
	limitset = new std::set<double>((*limits).begin(),(*limits).end());
}


// constructor building a CompOloyCurv from coefficients, limits an order and 
// number of segments
template<class T>
CompPolyCurv<T>::CompPolyCurv(int Num, int Ord, const Vector<T>& Coeffs, const Vector<double>& Limits) :
num(Num), ord(Ord), ords(new Vector<int>(Num)), coeffs(new Vector<T>(Coeffs)), limits(new Vector<double>(Limits))
{
	// assign the orders
	for (int i=0; i<num; i++) (*ords)[i]=ord;
	limitset = new std::set<double>((*limits).begin(),(*limits).end());
}


// constructor building CompPolyCurv from coefficients, an order and
// number of sgements
template<class T>
CompPolyCurv<T>::CompPolyCurv(int Num, int Ord, const Vector<T>& Coeffs) :
num(Num), ord(Ord), ords(new Vector<int>(Num)), coeffs(new Vector<T>(Coeffs)), limits(new Vector<double>(Num+1))
{
	// assign orders
	for (int i=0; i<num; i++) (*ords)[i]=ord;
	//assign limits
	for (int i=0; i<=num; i++) (*limits)[i]=(double)i;
	limitset = new std::set<double>((*limits).begin(),(*limits).end());
}


// ACCESS FUNCTIONS


// get the ith segment of the CompPolyCurv
template<class T>
PolyCurv<T> CompPolyCurv<T>::GetSegmentOriginal(int seg) const
{
	// extract control points
	//range();if (i <=0 || i > num) throw();
	int Ord = (*ords)[seg-1];
	Vector<T> ncoeffs(Ord);
	
	int sum = 0;
	for (int i=0; i<seg; i++) sum=sum+(*ords)[i];
	sum=sum-Ord;

	for (int k=0; k<Ord; k++) ncoeffs[k] = (*coeffs)[sum+k];

	// construct Poly
	return PolyCurv<T>(ncoeffs,ord,(*limits)[seg-1],(*limits)[seg]);
}



// get the ith segment of the CompPolyCurv
template<class T>
PolyCurv<T> CompPolyCurv<T>::GetSegment(int seg) const
{
	// extract control points
	//range();if (i <=0 || i > num) throw();
	Vector<T> ncoeffs(ord);
	int sum = (seg-1)*ord;
	for (int k=0; k<ord; k++) ncoeffs[k] = (*coeffs)[sum+k];

	// construct Poly
	return PolyCurv<T>(ncoeffs,ord,(*limits)[seg-1],(*limits)[seg]);
}


template<class T>
double CompPolyCurv<T>::GetLeftLimit() const
{
	return (*limits)[0];
}

template<class T>
double CompPolyCurv<T>::GetRightLimit() const
{
	return (*limits)[num];
}

// get the order of the CompPolyCurv
template<class T>
inline int CompPolyCurv<T>::GetOrd() const { return ord; }

// get the orders
template<class T>
inline Vector<int> CompPolyCurv<T>::GetOrds() const { return *ords; }

// get the number of segments
template<class T>
inline int CompPolyCurv<T>::GetNum() const { return num; }

// get the coefficients
template<class T>
inline Vector<T> CompPolyCurv<T>::GetCoeffs() const { return *coeffs; }

// get the limit set
template<class T>
inline Vector<double> CompPolyCurv<T>::GetLimits() const { return *limits; }

// get the left limit of index segment
template<class T>
inline double CompPolyCurv<T>::GetLeftLimit(int index) const
{
	return (*limits)[index-1];
}

// get the right limit of index segment
template<class T>
inline double CompPolyCurv<T>::GetRightLimit(int index) const
{
	return (*limits)[index];
}





template<class T>
int CompPolyCurv<T>::FindSegment(double x) const
{
	// range(x);
	std::set<double>::iterator s1 = (*limitset).begin();
	std::set<double>::iterator s2 = (*limitset).upper_bound(x);
	if (s2 == (*limitset).end()) return (*limitset).size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);
	return count;
}




// CONVERSIONS

// convert to CompBezCurv
template<class T>
CompBezCurv<T> CompPolyCurv<T>::ConvertCompBezCurv() const
{
	// take each segment and convert to Bezier form
	Vector<T> ncpts(num*ord);
	int count=0;
	
	// convert each segment
	for (int i=0; i<num; i++) {
		Vector<T> v = GetSegment(i+1).ConvertBezCurvCPoints();
		for (int j=0; j<ord; j++) {
			ncpts[count]=v[j];
			count++;
		}
	}

	return CompBezCurv<T>(num, *ords, ncpts, *limits);
}


// convert to BspCurv
template<class T>
BspCurv<T> CompPolyCurv<T>::ConvertBspCurv() const
{
	return ConvertCompBezCurv().ConvertBspCurv();
}

// convert all segments of curve to a common degree
template<class T>
Vector<T> CompPolyCurv<T>::Convert() const
{
	// elevate the degree of each segment to max
	// build array of control points

	int count=0;
	Vector<T> ncoeffs(num*ord);
	// convert each segment by degree elevation
	for (int i=0; i<num; i++) {
		Vector<T> v = GetSegmentOriginal(i+1).ElevateCPoints(ord-(*ords)[i]);
		for (int j=0; j<ord; j++) {
			ncoeffs[count]=v[j];
			count++;
		}
	}

	// create and return the new CompPolyCurv
	return ncoeffs;
}


// EVALUATORS

// evaluate the CompPolyCurv at the point x 
template<class T>
T CompPolyCurv<T>::operator()(double x) const // throw
{
	// Find segment
	int index = FindSegment(x);

	// extract and evaluate segment
	return GetSegment(index)(x);
}

// evaluate the segment index at the point x 
template<class T>
T CompPolyCurv<T>::operator()(double x, int index) const // throw
{
	// extract and evaluate segment
	return GetSegment(index)(x);
}


// evaluate the CompPolyCurv using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompPolyCurv<T>::Eval(double x) const // throw
{
	// Find segment
	int ind = FindSegment(x);
	// extract segment and evaluate
	return GetSegment(ind).Eval(x);
}

// evaluate the CompPolyCurv using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompPolyCurv<T>::Eval(double x, int index) const // throw
{
	// extract and evaluate segment
	return GetSegment(index).Eval(x);
}


// DEGREE ELEVATION

// elevate the degree of the CompPolyCurv by level
template<class T>
CompPolyCurv<T> CompPolyCurv<T>::Elevate(int level) const
{
	if (level <=0) return *this;
	// find number of new control points
	
	int count=0;
	Vector<T> ncoeffs(num*(ord+level));
	// elevate each segment
	for (int i=0; i<num; i++) {
		Vector<T> v = GetSegment(i+1).ElevateCPoints(level);
		for (int j=0; j<ord+level; j++) {
			ncoeffs[count]=v[j];
			count++;
		}
	}

	Vector<int> Ords(num);
	for (int i=0; i<num; i++) Ords[i] = (*ords)[i]+level;
	// create and return the CompPolyCurv
	return CompPolyCurv<T>(num, Ords, ncoeffs, *limits);
}

// DERIVATIVES

// compute the derivative of the CompPolyCurv of order deriv and
// represent the result as another CompPolyCurv
template<class T>
CompPolyCurv<T> CompPolyCurv<T>::Derive(int level) const
{   
	if (level <=0) return *this;
	// if (level >= ord) return 0 CompPolyCurv
	Vector<T> ncoeffs(num*(ord-level));

	// derive each segment
	int count=0;
	for (int i=0; i<num; i++) {
		Vector<T> v = GetSegment(i+1).DeriveCPoints(level);
		for (int j=0; j<ord-level; j++) {
			ncoeffs[count]=v[j];
			count++;
		}
	}

	Vector<int> Ords(num);
	for (int i=0; i<num; i++) Ords[i] = (*ords)[i]-level;
	
	// create and return the derived curve	
	return CompPolyCurv<T>(num, Ords, ncoeffs, *limits);
}

template<class T>
T CompPolyCurv<T>::operator()(int lev, double x) const
{
	return Derive(lev)(x);
}

// evaluate the derivative of the CompPolyCurv of order deriv at a point
// x. Computes the derivative as a CompPolyCurv and then evaluates this at x
template<class T>
T CompPolyCurv<T>::Derive(int level, double x) const
{
	return Derive(level)(x);
}


// INTEGRATION

// integrate the CompPolyCurv between the limits x1 and x2. Computes
// the indefinite integral as a CompPolyCurv and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T CompPolyCurv<T>::Integrate(double x1, double x2) const
{
	if (x2 < (*limits)[0] || x1 > (*limits)[num]) return 0;

	if (x1 < (*limits)[0]) x1 = (*limits)[0];
	if (x2 > (*limits)[num]) x2 = (*limits)[num];
	
	// create the indefinite integral;
	CompPolyCurv<T> intCurve = Integrate(); 

	// evaluate and subtract
	return (intCurve(x2) - intCurve(x1));
}


// compute the indefinite integral of the CompPolyCurv and represent
// it as a CompPolyCurv of one higher degree
template<class T>
CompPolyCurv<T> CompPolyCurv<T>::Integrate() const
{
	// find number of new control points
	Vector<T> ncoeffs(num*(ord+1));

	// integrate each segment
	int count=0;
	for (int i=0; i<num; i++) {
		Vector<T> v = GetSegment(i+1).IntegrateCPoints();
		for (int j=0; j<ord+1; j++) {
			ncoeffs[count]=v[j];
			count++;
		}
	}

	Vector<int> Ords(num);
	for (int i=0; i<num; i++) Ords[i] = (*ords)[i]+1;
	
	// create and return the curve
	return CompPolyCurv<T>(num, Ords, ncpts, *limits);
}	


// compute the indefinite integral of the CompPolyCurv as a CompPolyCurv
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Vector<T> CompPolyCurv<T>::IntegrateCPoints() const
{
	// find number of new control points
	
	Vector<T> ncoeffs(num*(ord+1));

	// integrate each segment
	int count=0;
	for (int i=0; i<num; i++) {
		Vector<T> v = GetSegment(i+1).IntegrateCPoints();
		for (int j=0; j<ord+1; j++) {
			ncoeffs[count]=v[j];
			count++;
		}
	}
	// return the array of control points
	return ncoeffs;
}	

/*
template<class T>
CompPolyCurv<T> CompPolyCurv<T>::MakeCompatable(const CompPolyCurv<T>& b) const
{
	return (ConvertBspCurv().MakeBreakCompatable(b.ConvertBspCurv())).ConvertCompPolyCurv();
}		
*/
/*
// PRODUCT

// compute the product of the CompPolyCurv with another CompPolyCurv and 
// represent the result as a new CompPolyCurv 
// assumes both curves have the same number of segments
template<class T>  
CompPolyCurv<T> CompPolyCurv<T>::Product(const CompPolyCurv<T>& b) const
{
	CompPolyCurv<T> d(*this), e(b);

	if (num != b.GetNum()) {
		d = MakeCompatable(b);
		e = b.MakeCompatable(d);
	}

	// create array for new coeffs
	Vector<T> ncoeffs(d.GetNum()*(ord+b.GetOrd()-1));
	int count=0;

	// multiply corresponding segments
	for (int i=0; i<d.GetNum(); i++) {
		PolyCurv<T> p1 = d.GetSegment(i+1);
		PolyCurv<T> p2 = e.GetSegment(i+1);
		Vector<T> temp = p1.ProductCPoints(p2);
		for (int j=0; j<ord+b.GetOrd()-1; j++) {
			ncoeffs[count] = temp[j];
			count++;
		}
	}
	Vector<int> Ords(d.GetNum());
	for (int i=0; i<d.GetNum(); i++) Ords[i] = (*ords)[i]+b.GetOrds()[i]-1;
	
	// return composite poly
	return CompPolyCurv<T>(d.GetNum(),Ords,ncoeffs,*limits);
}
*/

// READ and WRITE

template <class T>
void CompPolyCurv<T>::write(std::ostream& os) 
{
	os << "Comp Poly Curve\n";
	os << "number of segments\n";
	os << num;
	os << "\norders are\n";
	os << *ords;
	os << "\nlimits are\n";
	os << *limits;
	os << "\ncoefficients are\n";
	os << *coeffs;
}
	
	
template <class T>
void CompPolyCurv<T>::read(std::istream& is)
{
	std::cout << "number of segments\n";
	int Num;
	is >> Num;
	Vector<int> Ords(Num);
	std::cout << "orders of comp poly curve\n";
	Vector<double> Limits(Num+1);
	is >> Ords;
	std::cout << "limits\n";
	is >> Limits;
	int sum =0;
	for (int i=0; i<Num; i++) sum+=Ords[i];
	Vector<T> Coeffs(sum);
	std::cout << "input coefficients\n";
	is >> Coeffs;
	*this = CompPolyCurv<T>(Num,Ords,Coeffs,Limits);
} 


template <class T>
void CompPolyCurv<T>::writefile(std::ofstream& ofs)
{
	ofs << "CompPolyCurve\n";
	ofs << "number of segments are\n";
	ofs << num;
	ofs << "orders are\n";
	ofs << *ords;
	ofs << "limits are\n";
	ofs << *limits;
	ofs << "\ncoefficients are\n";
	ofs << *coeffs;
}

template <class T>
void CompPolyCurv<T>::readfile(std::ifstream& ifs)
{
	int Num;
	ifs >> Num;
	Vector<int> Ords(Num);
	ifs >> Ords;
	Vector<double> Limits(Num+1);
	ifs >> Limits;
	int sum =0;
	for (int i=0; i<Num; i++) sum+=Ords[i];
	Vector<T> Coeffs(sum);
	ifs >> Coeffs;
	*this = CompPolyCurv<T>(Num,Ords,Coeffs,Limits);
} 


// compute the product of the CompPolyCurv with another CompPolyCurv and 
// represent the result as a new CompPolyCurv 
// assumes both curves have different numbers of segments
template<class T>  
CompPolyCurv<T> CompPolyCurv<T>::Product2(const CompPolyCurv<T>& b) const
{
	// segment numbers are different
	// convert to CompBezCurv and multiply
	CompBezCurv<T> b1 = ConvertCompBezCurv();
	CompBezCurv<T> b2 = b.ConvertCompBezCurv();
	CompBezCurv<T> prod = b1.Product(b2);
	// return composite poly
	return (prod.ConvertCompPolyCurv());
}


#endif











