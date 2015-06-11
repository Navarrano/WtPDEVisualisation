#ifndef COMPBEZCURV
#define COMPBEZCURV


#include "bezcurv.h"
#include "bspcurv.h"
#include "comppolycurv.h"


// can be of different orders, IGES spec

template<class T>
class CompBezCurv : public Curve<T> {
private:
	// data
	int num; // number of segments
	int ord; // maximum order
	Ptr<Vector<int> > ords; // vector of orders
	Ptr<Vector<T> > cpts;	// control points
	Ptr<Vector<double> > kts; // knots
	Ptr<KnotSet> kset; // just added
	Ptr<std::set<double> > limitset;
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class CompBezCurv<double>");
		else {
			std::string s(typeid(T).name()), s1("class CompBezCurv<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:
	// constructors
	CompBezCurv();
	CompBezCurv(int Num, int Ord, const Vector<T>& Cpts, const Vector<double>& Kts);
	CompBezCurv(int Num, int Ord, const Vector<T>& Cpts);
	CompBezCurv(int Num, const Vector<int>& Ords, const Vector<T>& Cpts, const Vector<double>& Kts);
	CompBezCurv(int Num, const Vector<int>& Ords, const Vector<T>& Cpts);
	CompBezCurv(int Num, const Vector<double>& limits, int Ord, const Vector<T>& Cpts);
	CompBezCurv(int Num, const Vector<double>& limits, const Vector<int>& Ords, const Vector<T>& Cpts);
	CompBezCurv(const Vector<BezCurv<T> >& vec, int nseg);
	CompBezCurv(const Vector<BezCurv<T> >& vec, int nseg, const Vector<double>& limits);
	
	// access functions
	int GetNum() const;
	int GetNumCPoints() const;
	Vector<T> GetCPoints() const;
	Vector<double> GetKnots() const;
	KnotSet GetKnotSet() const;
	int GetOrd() const;
	int FindSegment(double x) const;
	Vector<int> GetOrds() const;
	Vector<T> Convert() const;
	BezCurv<T> GetSegment(int i) const;
	BezCurv<T> GetSegmentOriginal(int i) const;
	Vector<double> GetLimits() const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;

	// evaluators
	virtual T operator()(double x) const;
	virtual T operator()(int, double x) const;
	T Eval(double x) const;
	Vector<T> ComputePoints(int n) const;

	// derivative
	CompBezCurv<T> Derive(int level) const; 
	virtual T Derive(int level, double x) const;

	// make compatable with another CompBezCurv
	template<class T1>
	CompBezCurv<T> MakeCompatable(const CompBezCurv<T1>& b) const
	{
		return (ConvertBspCurv().MakeBreakCompatable(b.ConvertBspCurv())).ConvertCompBezCurv();
	}

	// degree elevation
	CompBezCurv<T> Elevate(int level) const;

	// integration
	CompBezCurv<T> Integrate() const;
	T Integrate(double x1, double x2) const;
	Vector<T> IntegrateCPoints() const;

	// product
	template<class T1>
	CompBezCurv<T> Product(const CompBezCurv<T1>& b) const
	{
		// test whether num segments is the same
		CompBezCurv<T> d(*this);
		CompBezCurv<T1> e(b);

		if (!(*kset).IsSameDistinctKnotSet(b.GetKnotSet())) {
			// needs to add distinct knots from b not in current object
			// into current object (independent of order)
		
			d = MakeCompatable(b);
			e = b.MakeCompatable(d);
		}

		Vector<T> ncpts((d.GetOrd()+e.GetOrd()-1)*d.GetNum());

		// extract and multiply coresponding segments
		int count=0;
	
		for (int i=0; i<d.GetNum(); i++) {
			BezCurv<T> bez1 = d.GetSegment(i+1);
			BezCurv<T1> bez2 = e.GetSegment(i+1);
			Vector<T> v = bez1.ProductCPoints(bez2);
			for (int j=0; j<d.GetOrd()+e.GetOrd()-1; j++) {
				ncpts[count] = v[j];
				count++;
			}
		}

		Vector<int> Ords(d.GetNum());
		for (int i=0; i<d.GetNum(); i++) Ords[i]=(*ords)[i]+b.GetOrds()[i]-1;
	//	for (int i=0; i<d.GetNum(); i++) Ords[i]=d.GetOrds()[i]+e.GetOrds()[i]-1;
	
		return CompBezCurv<T>(d.GetNum(),GetLimits(),Ords,ncpts);
	}

	template<class T1>
	CompBezCurv<T> Product2(const CompBezCurv<T1>& c) const
	{
		// segment numbers are different
		// convert to CompBezCurv and multiply
		BspCurv<T> b1 = ConvertBspCurv();
		BspCurv<T> b2 = b.ConvertBspCurv();
		BspCurv<T> prod = b1.Product(b2);
		// return composite poly
		return (prod.ConvertCompBezCurv());
	}

	// conversion
	BspCurv<T> ConvertBspCurv() const;
	CompPolyCurv<T> ConvertCompPolyCurv() const;
	
	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};

// CONSTRUCTORS

// default constructor
template<class T>
CompBezCurv<T>::CompBezCurv() : num(0), ord(0), ords(), cpts(), kts(), kset()
{
}

// constructor builds a CompBezCurv from a Vector of control points, knots
// an order and number of control points
template<class T>
CompBezCurv<T>::CompBezCurv(int Num, const Vector<double>& Limits, const Vector<int>& Ords, const Vector<T>& Cpts) :
num(Num), ords(new Vector<int>(Ords)), cpts(new Vector<T>(Cpts))
{
	// find maximum order
	ord = (*ords)[0];
	for (int i=1; i<num; i++)  if (ord < (*ords)[i]) ord=(*ords)[i];

	kts = new Vector<double>(Math::CreateKnots(num,ord,Limits));
	kset = new KnotSet(*kts,ord,ord*(num+1));
	cpts = new Vector<T>(Convert());
	kset = new KnotSet(*kts,ord,(num+1)*ord);
	limitset = new std::set<double>((*kts).begin(),(*kts).end());
}


// constructor builds a CompBezCurv from a Vector of control points, knots
// an order and number of control points
template<class T>
CompBezCurv<T>::CompBezCurv(int Num, const Vector<int>& Ords, const Vector<T>& Cpts) :
num(Num), ords(new Vector<int>(Ords)), cpts(new Vector<T>(Cpts))
{
	// find maximum order
	ord = (*ords)[0];
	for (int i=1; i<num; i++)  if (ord < (*ords)[i]) ord=(*ords)[i];

	kts = new Vector<double>(Math::CreateKnots(num,ord));
	kset = new KnotSet(*kts,ord,ord*(num+1));
	cpts = new Vector<T>(Convert());
	limitset = new std::set<double>((*kts).begin(),(*kts).end());
}


// constructor builds a CompBezCurv from a Vector of control points, knots
// an order and number of control points
template<class T>
CompBezCurv<T>::CompBezCurv(int Num, const Vector<int>& Ords, const Vector<T>& Cpts, const Vector<double>& Kts) :
num(Num), ords(new Vector<int>(Ords)), cpts(new Vector<T>(Cpts)), kts(new Vector<double>(Kts)),
kset(new KnotSet(*kts,ord,ord*(num+1)))
{
	// find maximum order
	ord = (*ords)[0];
	for (int i=1; i<num; i++)  if (ord < (*ords)[i]) ord=(*ords)[i];
	cpts = new Vector<T>(Convert());
	limitset = new std::set<double>((*kts).begin(),(*kts).end());
}




// constructor builds a CompBezCurv from control points, knots, order and number
// of segments
template<class T>
CompBezCurv<T>::CompBezCurv(int Num, int Ord, const Vector<T>& Cpts, const Vector<double>& Kts) :
num(Num), ord(Ord), ords(new Vector<int>(Num)), cpts(new Vector<T>(Cpts)), 
kts(new Vector<double>(Kts)), kset(new KnotSet(*kts,ord,(num+1)*ord)),
limitset(new std::set<double>(Kts.begin(),Kts.end()))
{
	// assign vector of orders	
	for (int i=0; i<num; i++) (*ords)[i]=ord;
}

// constructor builds a CompBezCurv from control points an order and
// a number of segments
template<class T>
CompBezCurv<T>::CompBezCurv(int Num, int Ord, const Vector<T>& Cpts) :
num(Num), ord(Ord), ords(new Vector<int>(Num)), cpts(new Vector<T>(Cpts)), 	kts(new Vector<double>(Math::CreateKnots(num,ord))),
kset(new KnotSet(*kts,ord,(num+1)*ord)), limitset(new std::set<double>((*kts).begin(),(*kts).end()))
{
	// set orders
	for (int i=0; i<num; i++) (*ords)[i]=ord;
	// create knots
}


template<class T>
CompBezCurv<T>::CompBezCurv(int Num, const Vector<double>& Limits, int Ord, const Vector<T>& Cpts) :
num(Num), ord(Ord), ords(new Vector<int>(Num)), cpts(new Vector<T>(Cpts)), kts(new Vector<double>(Math::CreateKnots(num,ord,Limits))),
kset(new KnotSet(*kts,ord,ord*(num+1))), limitset(new std::set<double>((*kts).begin(),(*kts).end()))
{
	for (int i=0; i<=num; i++) (*ords)[i]=ord;
}

template<class T>
CompBezCurv<T>::CompBezCurv(const Vector<BezCurv<T> >& vec, int nseg) : num(nseg), ords(new Vector<int>(nseg))
{
	int sum = 0;
	for (int i=0; i<num; i++) (*ords)[i] = vec[i].GetOrd();

	for (int i=0; i<num; i++) sum = sum+ (*ords)[i];
	cpts = new Vector<T>(sum);
	
	int count = 0;
	for (int i=0; i<num; i++) 
		for (int j=0; j<(*ords)[i]; j++) {
			(*cpts)[count] = vec[i].GetCPoints()[j];
			count++;
		}

	// find max order
	ord = (*ords)[0];
	for (int i=1; i<num; i++) if (ord < (*ords)[i]) ord=(*ords)[i];

	kts = new Vector<double>(Math::CreateKnots(num,ord));
	
	kset = new KnotSet(*kts,ord,(num+1)*ord);
	limitset = new std::set<double>((*kts).begin(),(*kts).end());
	cpts = new Vector<T>(Convert());
}

template<class T>
CompBezCurv<T>::CompBezCurv(const Vector<BezCurv<T> >& vec, int nseg, const Vector<double>& limits) : num(nseg), ords(new Vector<int>(nseg))
{
	int sum = 0;
	for (int i=0; i<num; i++) (*ords)[i] = vec[i].GetOrd();
	for (int i=0; i<num; i++) sum = sum+ (*ords)[i];

	cpts = new Vector<T>(sum);
	int count = 0;
	for (int i=0; i<num; i++) 
		for (int j=0; j<(*ords)[i]; j++) {
			(*cpts)[count] = vec[i].GetCPoints()[j];
			count++;
		}
	
	// find max order
	ord = (*ords)[0];
	for (int i=1; i<num; i++) if (ord < (*ords)[i]) ord=(*ords)[i]; 

	kts = new Vector<double>(Math::CreateKnots(num,ord,limits));
	kset = new KnotSet(*kts,ord,(num+1)*ord);
	limitset = new std::set<double>((*kts).begin(),(*kts).end());
	cpts = new Vector<T>(Convert());
}



// ACCESS FUNCTIONS
template<class T>
int CompBezCurv<T>::GetNumCPoints() const
{
	return num*ord;
}



// get the original ith segment and return as a BezCurv
template<class T>
BezCurv<T> CompBezCurv<T>::GetSegmentOriginal(int ind) const
{
	// find start index for control points
	int sum=0;
	int Ord = (*ords)[ind-1];
	for (int i=0; i<ind; i++) sum = sum + (*ords)[i];
	

	// extract control points
	Vector<T> ncpts(Ord);
	for (int k=0; k<Ord; k++) ncpts[k] = (*cpts)[sum-Ord+k];

	return BezCurv<T>(ncpts,Ord,(*kset).GetDistinctKnots()[ind-1],(*kset).GetDistinctKnots()[ind]);
}


// get the ith segment when all orders are the same and return as a BezCurv
template<class T>
BezCurv<T> CompBezCurv<T>::GetSegment(int ind) const
{
	// find start index for control points
	int sum=ind*ord;
	
	// extract control points
	Vector<T> ncpts(ord);
	for (int k=0; k<ord; k++) ncpts[k] = (*cpts)[sum-ord+k];

	// construct knots
	return BezCurv<T>(ncpts,ord,(*kset).GetDistinctKnots()[ind-1],(*kset).GetDistinctKnots()[ind]);
}


// get the order of the CompBezCurv
template<class T>
inline int CompBezCurv<T>::GetOrd() const { return ord; }

// get the orders
template<class T>
inline Vector<int> CompBezCurv<T>::GetOrds() const { return *ords; }

// get the number of segments
template<class T>
inline int CompBezCurv<T>::GetNum() const { return num; }

// get the control points
template<class T>
inline Vector<T> CompBezCurv<T>::GetCPoints() const { return *cpts; }

// get the knot vector
template<class T>
inline Vector<double> CompBezCurv<T>::GetKnots() const { return *kts; }

template<class T>
inline KnotSet CompBezCurv<T>::GetKnotSet() const { return *kset; }


template<class T>
inline Vector<double> CompBezCurv<T>::GetLimits() const
{
	return (*kset).GetDistinctKnots();
}


template<class T>
inline double CompBezCurv<T>::GetLeftLimit() const
{
	return (*kts)[0];
}

template<class T>
inline double CompBezCurv<T>::GetRightLimit() const
{
	return (*kts)[num*ord];
}



template<class T>
int CompBezCurv<T>::FindSegment(double x) const
{
	std::set<double>::iterator s1 = (*limitset).begin();
	std::set<double>::iterator s2 = (*limitset).upper_bound(x);
	if (s2 == (*limitset).end()) return (*limitset).size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);
	return count;
}


// CONVERSION
// convert to a BspCurv removing all possible knots
template<class T>
BspCurv<T> CompBezCurv<T>::ConvertBspCurv() const
{
	// convert
	return BspCurv<T>(*cpts,*kts,ord,ord*num);
}


// convert all segments to same degree
template<class T>
Vector<T> CompBezCurv<T>::Convert() const
{
	// elevate the degree of each segment to max
	// extract BezCurv segment
	// build array of control points
	
	int count=0;
	Vector<T> ncpts(num*ord);
	
	// elevate and extract control points
	for (int i=0; i<num; i++) {
		Vector<T> b = GetSegmentOriginal(i+1).ElevateCPoints(ord-(*ords)[i]);
		for (int j=0; j<ord; j++) {
			ncpts[count]=b[j];
			count++;
		}
	}

	// create and return the new control points
	return ncpts;
}


// convert to CompPoly form
template<class T>
CompPolyCurv<T> CompBezCurv<T>::ConvertCompPolyCurv() const
{
	// take each segment and convert to PolyCurv form
	// find the sum of the ords
	
	Vector<T> ncoeffs(num*ord);
	int count=0;
	// convert each segment
	for (int i=0; i<num; i++) {
		PolyCurv<T> p = GetSegment(i+1).ConvertPolyCurv().Reparameterise2((*kset).GetDistinctKnots()[i],(*kset).GetDistinctKnots()[i+1]);
		for (int j=0; j<ord; j++) {
			ncoeffs[count]=p.GetCoeffs()[j];
			count++;
		}
	}

	// create and return the CompPolyCurv
	return CompPolyCurv<T>(num, ord, ncoeffs, (*kset).GetDistinctKnots());
}

/*
// MAKE COMPATABLE

template<class T>
CompBezCurv<T> CompBezCurv<T>::MakeCompatable(const CompBezCurv<T>& b) const
{
	return (ConvertBspCurv().MakeBreakCompatable(b.ConvertBspCurv())).ConvertCompBezCurv();
}		
*/

// EVALUATORS
// evaluate the CompBezCurv at the point x using de Boor algorithm
template<class T>
T CompBezCurv<T>::operator()(double x) const
{
	int ind = FindSegment(x); // find the segment
	// extract BezCurv and evaluate
	int sum=ind*ord;

  	return GetSegment(ind)((x-(*kts)[sum-1])/((*kts)[sum]-(*kts)[sum-1]));
}


// evaluate the CompBezCurv using standard method rather than overloaded operator
// identical code to that of overloaded operator
template<class T>
T CompBezCurv<T>::Eval(double x) const
{
	int ind = FindSegment(x); // find the segment
	// extract BezCurv and evaluate
	int sum=ind*ord;

  	return GetSegment(ind).Eval((x-(*kts)[sum-1])/((*kts)[sum]-(*kts)[sum-1]));
}


// DEGREE ELEVATION

// elevate the degree of the CompBezCurv by level
template<class T>
CompBezCurv<T> CompBezCurv<T>::Elevate(int level) const
{
	// find number of new control points
	if (level <= 0) return *this;

	int count=0;

	Vector<T> ncpts(num*(ord+level));
	// elevate each segment
	for (int i=0; i<num; i++) {
		Vector<T> b = GetSegment(i+1).ElevateCPoints(level);
		for (int j=0; j<ord+level; j++) {
			ncpts[count]=b[j];
			count++;
		}
	}

	Vector<int> Ords(num);
	for (int i=0; i<num; i++) Ords[i]=(*ords)[i]+level;
	// create and return the elevated curve
	return CompBezCurv<T>(num, GetLimits(), Ords, ncpts);
}


// DERIVATIVES

// compute the derivative of the CompBezCurv of order deriv and
// represent the result as another CompBezCurv
template<class T>
CompBezCurv<T> CompBezCurv<T>::Derive(int level) const
{   
	// find number of new control points
	if (level <= 0) return *this;

	if (level >= ord) {
		KnotSet kset(*kts,ord,ord*(num+1));
		Vector<double> knots(kset.GetDistinctKnots());
		Vector<T> ncpts(num-1,0.0);
		return CompBezCurv<T>(num,1,ncpts,knots);
	}
	
	
	Vector<T> ncpts(num*(ord-level));

	// derive each segment
	int count=0;
	for (int i=0; i<num; i++) {
		Vector<T> b = GetSegment(i+1).DeriveCPoints(level);
		for (int j=0; j<ord-level; j++) {
			ncpts[count]=b[j];
			count++;
		}
	}

	Vector<int> Ords(num);
	for (int i=0; i<num; i++) Ords[i]=(*ords)[i]-level;
	
	// create and return the derived curve
	return CompBezCurv<T>(num, GetLimits(), Ords, ncpts);
}

template<class T>
T CompBezCurv<T>::operator()(int lev, double x) const
{
	return Derive(lev)(x);
}

// evaluate the derivative of the CompBezCurv of order deriv at a point
// x. Computes the derivative as a CompBezCurv and then evaluates this at x
template<class T>
T CompBezCurv<T>::Derive(int level, double x) const
{
	return Derive(level)(x);
}


// INTEGRATION

// integrate the CompBezCurv between the limits x1 and x2. Computes
// the indefinite integral as a CompBezCurv and then evaluates 
// this at x2 and x1, subtracting the results.
template<class T>
T CompBezCurv<T>::Integrate(double x1, double x2) const
{
	// create the indefinite integral;
	Vector<double> dts = intCurve.GetKnotSet().GetDistinctKnots();
	if (x2 < dts[0] || x1 > dts[limitset.size()-1]) return 0;

	if (x1 < dts[0]) x1 = dts[0];
	if (x2 > dts[limitset.size()-1]) x2 = dts[limitset.size()-1];

	CompBezCurv<T> intCurve = Integrate(); 
	
	// evaluate and subtract
	return (intCurve(x2) - intCurve(x1));
}


// compute the indefinite integral of the CompBezCurv and represent
// it as a CompBezCurv of one higher degree
template<class T>
CompBezCurv<T> CompBezCurv<T>::Integrate() const
{
	// find number of new control points
	
	Vector<T> ncpts(num*(ord+1));

	// integrate each segment
	int count=0;
	for (int i=0; i<num; i++) {
		Vector<T> b = GetSegment(i+1).IntegrateCPoints();
		for (int j=0; j<ord+1; j++) {
			ncpts[count]=b[j];
			count++;
		}
	}

	Vector<int> Ords(num);
	for (int i=0; i<num; i++) Ords[i]=(*ords)[i]+1;
	
	// construct curve
	return CompBezCurv<T>(num,GetLimits(), Ords, ncpts);
}	


// compute the indefinite integral of the CompBezCurv as a CompBezCurv
// and return just the control points of this Curve as a CPoints
// object
template<class T>
Vector<T> CompBezCurv<T>::IntegrateCPoints() const
{
	Vector<T> ncpts(num*(ord+1));

	// integrate each segment
	int count=0;
	for (int i=0; i<num; i++) {
		Vector<T> b = GetSegment(i+1).IntegrateCPoints();
		for (int j=0; j<ord+1; j++) {
			ncpts[count]=b[j];
			count++;
		}
	}
	return ncpts;
}	

/*
// PRODUCT

// compute the product of the CompBezCurv with another CompBezCurv and 
// represent the result as a new CompBezCurv 
// assumes both curves have the same number of segments and same order
template<class T>  
CompBezCurv<T> CompBezCurv<T>::Product(const CompBezCurv<T>& b) const
{
	// test whether num segments is the same

	CompBezCurv<T> d(*this), e(b);


	if (!(*kset).IsSameDistinctKnotSet(b.GetKnotSet())) {
		// needs to add distinct knots from b not in current object
		// into current object (independent of order)
		d = MakeCompatable(b);
		e = b.MakeCompatable(d);
	}

	Vector<T> ncpts((d.GetOrd()+e.GetOrd()-1)*d.GetNum());

	
	// extract and multiply coresponding segments
	int count=0;
	
	for (int i=0; i<d.GetNum(); i++) {
		BezCurv<T> bez1 = d.GetSegment(i+1);
		BezCurv<T> bez2 = e.GetSegment(i+1);
		Vector<T> v = bez1.ProductCPoints(bez2);
		for (int j=0; j<d.GetOrd()+e.GetOrd()-1; j++) {
			ncpts[count] = v[j];
			count++;
		}
	}

	Vector<int> Ords(d.GetNum());
	for (int i=0; i<d.GetNum(); i++) Ords[i]=(*ords)[i]+b.GetOrds()[i]-1;
	
	return CompBezCurv<T>(d.GetNum(),GetLimits(),Ords,ncpts);
}		
*/
// READ and WRITE

template <class T>
void CompBezCurv<T>::write(std::ostream& os) 
{
	os << "Comp Bezier Curve\n";
	os << "number of segments\n";
	os << num;
	os << "\norders are\n";
	os << *ords;
	os << "\nknots are\n";
	os << *kts;
	os << "\ncontrol points are\n";
	os << *cpts;
}
	
	
template <class T>
void CompBezCurv<T>::read(std::istream& is)
{
	std::cout << "number of segments\n";
	int Num;
	is >> Num;
	Vector<int> Ords(Num);
	std::cout << "orders of Bezier curve\n";
	is >> Ords;
	int sum =0;
	for (int i=0; i<Num; i++) sum+=Ords[i];
	Vector<T> Cpts(sum);
	std::cout << "input control points\n";
	is >> Cpts;
	*this = CompBezCurv<T>(Num,Ords,Cpts);
} 


template <class T>
void CompBezCurv<T>::writefile(std::ofstream& ofs)
{
	ofs << "CompBezier Curve\n";
	ofs << "number of segments are\n";
//	ofs << num;
	ofs << "orders are\n";
	ofs << *ords;
	ofs << "knots are\n";
	ofs << *kts;
	ofs << "\ncontrol points are\n";
	ofs << *cpts;
}

template <class T>
void CompBezCurv<T>::readfile(std::ifstream&  ifs)
{
	int Num;
	ifs >> Num;
	Vector<int> Ords(Num);
	ifs >> Ords;
	int sum =0;
	for (int i=0; i<Num; i++) sum+=Ords[i];
	Vector<T> Cpts(sum);
	ifs >> Cpts;
	*this = CompBezCurv<T>(Num,Ords,Cpts);
} 


#endif