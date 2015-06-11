#ifndef CPOINTS
#define CPOINTS

#include "vector.h"
#include "knotset.h"

template<class T>
class CPoints {


	// data
	int num;	// number of control points
	Ptr<Vector<T> > cpts;	// vector of control points
	Ptr<KnotSet> kset; // associated knot set

	// private update functions
	CPoints<T>& UpdateCPointsLeft(); // updates left end
	CPoints<T>& UpdateCPointsRight(); // updates right end
	virtual ObjectID Identity() const 
	{ 
		if (std::string(typeid(T).name()) == std::string("double"))
			return std::string("class CPoints<double>");
		else {
			std::string s(typeid(T).name()), s1("class CPoints<");
			return s1.append(s,0,s.size())+">";
		}
	}
public:

	// constructors
	CPoints();
	CPoints(const Vector<T>& Cpts, int N);
	CPoints(const Vector<T>& Cpts, const KnotSet& Kset, int N);

	// access functions
	Vector<T> GetCPoints() const;
	KnotSet GetKnotSet() const;
	int GetNum() const;
	
	// update functions
	CPoints<T>& UpdateCPoints(); // updates control points for ord
								// multiplicity at ends
	CPoints<T>& UpdateKnotSet(const KnotSet& kset);
	CPoints<T> Removal(int ind1, int ind2) const;

	// modify cpoints
	CPoints<T> CreateCPointsDeriv(int level) const;
	CPoints<T> CreateCPointsInsert(double t) const;
	CPoints<T> CreateCPointsInsert(double t, int level) const;
	CPoints<T> CreateCPointsInsert(const Vector<double>& Kts, int N) const;
	CPoints<T> CreateCPointsInsert(const Vector<double>& Kts, const Vector<int>& Mult, int N) const;
	CPoints<T> CreateCPointsInsert(int level) const;
	CPoints<T> CreateCPointsSubdivide(double t1, double t2) const;
	CPoints<T> CreateCPointsSubdivide(int level) const;
	CPoints<T> CreateCPointsIntegrate() const;
};


// CONSTRUCTORS

// default constructor
template<class T>
CPoints<T>::CPoints() : cpts(), kset() { }


// constructor taking vector of control points and KnotSet object
// and number of control points
template<class T>
CPoints<T>::CPoints(const Vector<T>& Cpts, const KnotSet& Kset, int N) : 
	num(N),cpts(new Vector<T>(Cpts)), kset(new KnotSet(Kset)) 
{
}


// constructor taking vector of control points and number of control points
template<class T>
CPoints<T>::CPoints(const Vector<T>& Cpts, int N) :	
	num(N), cpts(new Vector<T>(Cpts)), kset(new KnotSet()) 
{
	// what about the knot set
}


// get the control points
template<class T>
inline Vector<T> CPoints<T>::GetCPoints()  const { return *cpts; }


// get the KnotSet object
template<class T>
inline KnotSet CPoints<T>::GetKnotSet() const 
{
	return *kset;
}


// get the number of control points
template<class T>
inline int CPoints<T>::GetNum() const
{
	return num;
}


// KNOT and CONTROL POINT REMOVAL

// remove the control point and knot from the curve
// given indices for both
template<class T>
CPoints<T> CPoints<T>::Removal(int ind1, int ind2) const
{
	int ord = (*kset).GetOrd();

	Vector<T> ncpts(num-1);
	Vector<double> kts(num-1+ord);
	for (int i=0; i<ind1; i++) ncpts[i] = (*cpts)[i];
	for (int i=ind1+1; i<num; i++) ncpts[i-1] = (*cpts)[i];
	for (int i=0; i<ind2; i++) kts[i] = (*kset).GetKnots()[i];
	for (int i=ind2+1; i<num+ord; i++) kts[i-1]=(*kset).GetKnots()[i];
	
	return CPoints<T>(ncpts, KnotSet(kts,ord,ord+num-1), num-1);
}


// PRIAVTE FUNCTIONS

// update the knot set, replaces original KnotSet with new one
// used primarily for converting knot set that doesn't have multiplicity
// ord at both ends to such a knot set
template<class T>
inline CPoints<T>& CPoints<T>::UpdateKnotSet(const KnotSet& Kset) 
{
	*kset = Kset;
	return *this;
}

// update the control points to represent a curve that has knot multiplicity
// equal to the order at both ends of the curve
template<class T>
CPoints<T>& CPoints<T>::UpdateCPoints()
{
	// find the multiplicity of the left end knot of index ord-1
	int s1 = (*kset).FindMultiplicity((*kset).GetKnots()[(*kset).GetOrd()-1]);

	// find the multipliicty of the right end knot of index num
	int s2 = (*kset).FindMultiplicity((*kset).GetKnots()[num]);
	
	// update the control points and knot sets
	if (s1 < (*kset).GetOrd()) {
		UpdateCPointsLeft();
		(*kset).UpdateKnotsLeft();
	}
	if (s2 < (*kset).GetOrd()) {
		UpdateCPointsRight();
		(*kset).UpdateKnotsRight();
	}
	// return the updated Cpoint object
	return *this;
}


// updates the cntrol points at the left hand end of the curve
// according to the coincident end knot specification
// Uses the subdivision property of the de Boor algorithm
template<class T>
CPoints<T>& CPoints<T>::UpdateCPointsLeft()
{
	int ord = (*kset).GetOrd();
	Vector<T> temp(ord);
	Vector<T> ntemp(ord-1);
	double lbd;

	int ind = ord;
	
	// extract the relevant control points
	for (int j=0; j<ord; j++) temp[j]=(*cpts)[ind-ord+j];

	//compute the point according to de Boor algorithm 
	for (int i=0; i<ord-1; i++) {
		for (int j=0; j<ord-i-1; j++) {
			lbd = ((*kset).GetKnots()[ord-1] - (*kset).GetKnots()[ind-ord+i+j+1])/((*kset).GetKnots()[ind+j]-(*kset).GetKnots()[ind-ord+i+j+1]);
			temp[j] = lbd*temp[j+1]+(1.0-lbd)*temp[j];
		}
		ntemp[i] = temp[ord-i-2];
	}
	// update the control points
	for (int i=0; i<ord-1; i++) (*cpts)[i] = ntemp[ord-i-2];

	// return the modified CPoint object
	return *this;
}


// updates the cntrol points at the right hand end of the curve
// according to the coincident end knot specification
// Uses the subdivision property of the de Boor algorithm
template<class T>
CPoints<T>& CPoints<T>::UpdateCPointsRight()
{
	int ord = (*kset).GetOrd();
	Vector<T> temp(ord);
	Vector<T> ntemp(ord-1);
	double lbd;

	int ind = num;
  

	for (int j=0; j<ord; j++) temp[j]=(*cpts)[ind-ord+j];

	//compute the point according to de Boor algorithm 
	for (int i=0; i<ord-1; i++) {
		for (int j=0; j<ord-i-1; j++) {
			lbd = ((*kset).GetKnots()[num] - (*kset).GetKnots()[ind-ord+i+j+1])/((*kset).GetKnots()[ind+j]-(*kset).GetKnots()[ind-ord+i+j+1]);
			temp[j] = lbd*temp[j+1]+(1.0-lbd)*temp[j];
		}
		ntemp[i] = temp[0];
	}

	// update control points
	for (int i=0; i<ord-1; i++) (*cpts)[ind-ord+i+1] = ntemp[i];

	// return update CPoint object
	return *this;
}

// MODIFY CPOINTS

// create control points after subdivision of original curve
// upto level 'level'.
template<class T>
CPoints<T> CPoints<T>::CreateCPointsSubdivide(int level) const
{
	double x, x1, x2;
	double step;
	CPoints<T> cpt(*this);

	// find number of distinct knots
	int end = (*kset).GetNumDistinct()-2;

	// between each two distinct knots insert
	// level knots and compute new control points
	for (int i=0; i<=end; i++) {
		x1 = (*kset).GetDistinctKnots()[i];
		x2 = (*kset).GetDistinctKnots()[i+1];
	
		step = (x2-x1)/level;
	
		for (int j=1; j<=level-1; j++) {
			x = x1+j*step;
			cpt = cpt.CreateCPointsInsert(x);
		}
	}
	// return CPoint object with control poinst for subdivided curve
	return cpt;
}


// create control points for derivative Bcurve, represented as a Bcurve
// Returns new Cpoint object with derivative control points
template<class T>
CPoints<T> CPoints<T>::CreateCPointsDeriv(int deriv) const
{
	// set up temporary arrays for calculations 
	if (deriv <= 0) return *this;
	int ord = (*kset).GetOrd();
	int n = (*kset).GetNumDistinct();
	Vector<T> temp(ord);
	Vector<double> dm(ord);
	Vector<T> dd(ord);
	Vector<double> dp(ord);	

	int norder = ord-deriv;

	Vector<int> over(n);
	Vector<int> mult((*kset).GetMult());
	over[0] = 0;
	over[n-1] = 0;
	int k1=1;

	for(int i=1;i<n-1;i++) over[i] = ord-mult[i]-deriv;

	// create control points
	KnotSet nkset = (*kset).CreateKnotSetDeriv(deriv);
	int ns = nkset.GetNum()-norder;
	Vector<T> ncoef(ns);
	double that;
	int ind;

	for(int l=0;l<n-1;l++) {
		that=(*kset).GetDistinctKnots()[l];
		//ind=(*kset).Fndint1(that);
		ind = (*kset).Find_index(that);

		for(int j=0;j<ord;j++) {
			temp[j] = (*cpts)[ind-ord+j+1];
			dp[j] = (*kset).GetKnots()[ind+j+1]-that;
			dm[j] = that-(*kset).GetKnots()[ind-ord+j+1];
		}
	
		// compute the points according to de Boor algorithm 
		for (int i=0; i<deriv; i++) {
			for(int j=i+1;j<ord;j++) 
				dd[j] = (temp[j]-temp[j-1])/(dm[j]+dp[j-i-1]);
			for(int j=1;j<ord;j++)
				temp[j] = dd[j];
		}
		if(over[l]>0) {	k1-=over[l];}

		for(int j=deriv;j<ord;j++) {
			ncoef[k1-1] = temp[j];
			k1++; 
		}
	}

	ind = 1;
	for (int i=1; i<=deriv; i++) ind*=(ord-i);
	for(int i=0;i<ns;i++) ncoef[i] = ncoef[i]*(double)ind;
	// create and return the new CPoint object
	return CPoints<T>(ncoef,nkset,ns); 
}


// fid the new control points after inserting a knot at
// a point x
template<class T>
CPoints<T> CPoints<T>::CreateCPointsInsert(double x) const
{	
	int ord = (*kset).GetOrd();
	Vector<T> newcoef(num+1);

	// find index and multiplicity of knot x in original curve
	//int p= (*kset).Find_index(x); //
	int p = (*kset).Fndint1(x); 
	int s= (*kset).FindMultiplicity(x);
	
	// need to catch the case when x is already a knot of mult = ord


	// calculate the new coefficients qaccording to knot insertion algorithm
	double alpha; 

	for(int i=0;i<=p-ord+s+1;i++)
		newcoef[i] = (*cpts)[i];
	for(int i=p-ord+s+2;i<=p;i++) {
		alpha=(x-(*kset).GetKnots()[i])/((*kset).GetKnots()[i+ord-1]-(*kset).GetKnots()[i]);
		newcoef[i] = alpha*(*cpts)[i] +(1-alpha)*(*cpts)[i-1];
	}
	for(int i=p+1;i<num+1;i++)
		newcoef[i] = (*cpts)[i-1];

	// update knot set
	CPoints<T> cpt(newcoef,num+1);
	KnotSet nkset = (*kset).CreateKnotSetInsert(x);
	return cpt.UpdateKnotSet(nkset);
}


// fid the new control points after inserting a knot at
// a point x 'level' times
template<class T>
CPoints<T> CPoints<T>::CreateCPointsInsert(double x, int level) const
{
	// int	p= (*kset).Find_index(x); 
	int p = (*kset).Fndint1(x); 
	int s= (*kset).FindMultiplicity(x);


	int ord = (*kset).GetOrd();
	if (level <= 0 || s + level > ord) return *this;

	double alpha;
	Vector<T> newcoef(num+level);
	Vector<T> temp(num+level);
	

	for (int i=0; i<num; i++) temp[i] = (*cpts)[i];
	// find new control points according to knot insertion algorithm
	for (int j=0; j<level; j++) {
		for (int i=0; i<=p-ord+s+j+1; i++)
			newcoef[i] = temp[i];
		for (int i=p-ord+s+j+2; i<=p; i++) {
			alpha=(x-(*kset).GetKnots()[i])/((*kset).GetKnots()[i+ord-j-1]-(*kset).GetKnots()[i]);
			newcoef[i] = alpha*temp[i] +(1-alpha)*temp[i-1];
		}
		for(int i=p+1;i<num+j+1;i++)
			newcoef[i] = temp[i-1];
		for (int i=0; i<num+j+1; i++) temp[i] = newcoef[i];
	}
	// create new knot set
	KnotSet nkset = (*kset).CreateKnotSetInsert(x, level);
	// create new CPoint object
	CPoints<T> cpt(newcoef,num+level);
	// assign new knot set
	return cpt.UpdateKnotSet(nkset);
}


template<class T>
CPoints<T> CPoints<T>::CreateCPointsInsert(const Vector<double>& Kts, int N) const
{
	 CPoints<T> cpt(*this);
    
	// perform knot insertion 
	for(int i=0;i<N;i++)
		cpt=cpt.CreateCPointsInsert(Kts[i]);

	return cpt;
}


template<class T>
CPoints<T> CPoints<T>::CreateCPointsInsert(const Vector<double>& Kts, const Vector<int>& Mult, int N) const
{
	CPoints<T> cpt(*this);
    
	// perform knot insertion 
	for(int i=0;i<N;i++)
//		for(int j=0;j<Mult[i];j++) {
		cpt=cpt.CreateCPointsInsert(Kts[i],Mult[i]);
	return cpt;
}

// if multiplicity is ord or ord-1 nothing to do
template<class T>
CPoints<T> CPoints<T>::CreateCPointsSubdivide(double x1, double x2) const
{
	int ord = (*kset).GetOrd();

	// find correct position to insset knot x1
	int s1= (*kset).FindMultiplicity(x1);

	// find correct position to insert knot x2
	int s2 = (*kset).FindMultiplicity(x2);

	
	CPoints<T> cpt1 = CreateCPointsInsert(x1,ord-s1);
	CPoints<T> cpt2 = cpt1.CreateCPointsInsert(x2,ord-s2);
	
	//extract the correct section of curve
	int ind1 = cpt2.GetKnotSet().Find_index(x1);
	int ind2 = cpt2.GetKnotSet().Find_index(x2);
	
	// detect if x2 is the end knot
	std::multiset<double> mset = cpt2.GetKnotSet().GetKnotSet();
	std::multiset<double>::iterator p = mset.begin();
	for (unsigned int i=0; i<mset.size()-1; i++) p++;
	if (x2 == *p) ind2 = ind2 + ord;
	
	Vector<double> newcoef(ind2-ind1);
	
	for (int i=ind1-ord+1; i<=ind2-ord; i++) newcoef[i-ind1+ord-1] = cpt2.GetCPoints()[i];
	return CPoints<T>(newcoef,ind2-ind1);
}




// create CPoints object representing the indefinite integral of the curve
template<class T>
CPoints<T> CPoints<T>::CreateCPointsIntegrate() const
{
	T higher = 0.0;
	int k=(*kset).GetOrd();

	Vector<T> ncpts(num+1);
	ncpts[0] = higher;
	for(int i=0; i<num; i++) {
		higher = 0.0;
		for (int j=0; j<=i; j++)  {
			higher = higher + ((*cpts)[j])*((*kset).GetKnots()[j+k]-(*kset).GetKnots()[j])/(double)k;
		}
		ncpts[i+1]=higher;
	}

	//ncpts[num]=0.0;
	return CPoints<T>(ncpts,num+1);
}



#endif