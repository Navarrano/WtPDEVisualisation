#include "knotset.h"

// CONSTRUCTORS

// default constructor
KnotSet::KnotSet() : ord(0), num(0), n(0), kts(), dts(), mult() { }

// constructor taking an integer Num the number of knots
KnotSet::KnotSet(int Num) : 
	ord(0), num(Num), n(Num), kts(new Vector<double>(Num)), dts(new Vector<double>(Num)), 
		mult(new Vector<int>(Num)), mset(new std::multiset<double>())
{ 
	// create knots somehow
}

// constructor taking an ord and a num of knots
KnotSet::KnotSet(int Ord, int Num) : 
	ord(Ord), num(Num), n(Num), kts(new Vector<double>(Num)), dts(new Vector<double>(Num)), 
		mult(new Vector<int>(Num)) , mset(new std::multiset<double>())
{
	// create knots somehow
}


// constructor taking vector of knots, ord and num
KnotSet::KnotSet(const Vector<double>& Kts, int Ord, int Num) : 
	ord(Ord), num(Num), kts(new Vector<double>(Kts)), mset(new std::multiset<double>((*kts).begin(),(*kts).end())),
		dts(new Vector<double>(ComputeDistinctKnots())), 
		mult(new Vector<int>(ComputeMultiplicities()))
{
	n = (*dts).size();
}

// constructor taking distinct knots and multiplicities, order and num of distinct knots
KnotSet::KnotSet(const Vector<double>& Dts, const Vector<int>& Mult, int Ord, int N) :	
	ord(Ord), n(N), dts(new Vector<double>(Dts)), mult(new Vector<int>(Mult))
{
	kts = new Vector<double>(ComputeKnots());
	mset = new std::multiset<double>((*kts).begin(),(*kts).end());
	num = (*kts).size();
}

// PRIVATE CONSTRUCTOR


// constructor taking vector of knots, distinct knots and multplicities and num of distinct knots
KnotSet::KnotSet(const Vector<double>& Kts, const Vector<double>& Dts, const Vector<int>& Mult, int N) :
ord(0), num(N), n(N), kts(new Vector<double>(Kts)), mset(new std::multiset<double>((*kts).begin(),(*kts).end())), 
dts(new Vector<double>(Dts)), mult(new Vector<int>(Mult))
{
}


// ACCESS FUNCTIONS

// get multiplicity vector
inline Vector<int> KnotSet::GetMult() const 
{
	return *mult;
}

// get num distinct knots
inline int KnotSet::GetNumDistinct() const 
{
	return n;
}

// get the knots
inline Vector<double> KnotSet::GetKnots() const
{
	return *kts;
}

// get the distinct knts()
inline Vector<double> KnotSet::GetDistinctKnots() const
{
	return *dts;
}

inline std::multiset<double> KnotSet::GetKnotSet() const
{
	return *mset;
}

// get the order
inline int KnotSet::GetOrd() const 
{
	return ord;
}

// get the total number of knots
inline int KnotSet::GetNum() const 
{
	return num;
}


bool KnotSet::Contains(const KnotSet& kset) const
{
	return KnotIntersection(kset).IsSameKnotSet(kset);
}

KnotSet KnotSet::KnotIntersection(const KnotSet& kset) const
{
	std::multiset<double> set1(kset.GetKnotSet());
	std::multiset<double> set2;

	std::set_intersection((*mset).begin(),(*mset).end(),set1.begin(),set1.end(),std::inserter(set2,set2.begin()));
	if (set2.size() > 0) {
		Vector<double> nkts(set2.size());
		std::copy(set2.begin(),set2.end(),nkts.begin());
		return KnotSet(nkts,ord,set2.size());
	} else return KnotSet();
}
	

// PRIVATE FUNCTIONS
KnotSet KnotSet::KnotDifference(const KnotSet& kset) const
{
	std::multiset<double> set1(kset.GetKnotSet());
	std::multiset<double> set2;

	std::set_difference((*mset).begin(),(*mset).end(),set1.begin(),set1.end(),std::inserter(set2,set2.begin()));
	if (set2.size() > 0) {
		Vector<double> nkts(set2.size());
		std::copy(set2.begin(),set2.end(),nkts.begin());
		return KnotSet(nkts,ord,set2.size());
	} else return KnotSet();
	//return KnotSet(nkset.GetDistinctKnots(),nkset.GetDistinctKnots(),nkset.GetMult(),set3.size());
}


KnotSet KnotSet::SortKnotSet() const
{
	Vector<double> nkts((*mset).size());
	std::copy((*mset).begin(),(*mset).end(),nkts.begin());

	KnotSet kset(nkts,ord,(*mset).size());

	return KnotSet(kset.GetDistinctKnots(),kset.GetDistinctKnots(),kset.GetMult(),kset.GetNumDistinct());
}
		
		
int KnotSet::FindMultiplicity(double x) const
{
	return (*mset).count(x);
}


int KnotSet::Fndint1(double x) const
{
	std::multiset<double>::iterator r2 = (*mset).lower_bound(x);
	std::multiset<double>::iterator r1 = (*mset).begin();
	
	if (r2 == r1) return ord-1; // added
	int i = 0;
	do { i++; } while (r1++ != r2);

	return i-2;
}


int KnotSet::Find_index(double x) const
{
	std::multiset<double>::iterator r3 = (*mset).upper_bound(x);
	std::multiset<double>::iterator r2 = (*mset).lower_bound(x);
	std::multiset<double>::iterator r1 = (*mset).begin();
	if (r2 == r1 || r2 == (*mset).end()) return ord-1;
	else if (r3 == (*mset).end()) return num-ord-1;
	else {
		int i = 0;
		do { i++; } while (r1++ != r3);
		return i-2;
	}
}


int KnotSet::Find_segment(double x) const
{
	std::set<double> s((*dts).begin(),(*dts).end());

	std::set<double>::iterator s1 = s.begin();
	std::set<double>::iterator s2 = s.upper_bound(x);
	if (s2 == s.end()) return s.size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);
	return count;
}

int KnotSet::ComputeParameterAverUpper(double x) const
{
	Vector<double> v = ComputeKnotSetAver();

	std::set<double> s(v.begin(),v.end());

	std::set<double>::iterator s1 = s.begin();
	std::set<double>::iterator s2 = s.upper_bound(x);
	if (s2 == s.end()) return s.size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);

	return count;
}

int KnotSet::ComputeParameterAverLower(double x) const
{
	Vector<double> v = ComputeKnotSetAver();

	std::set<double> s(v.begin(),v.end());

	std::set<double>::iterator s1 = s.begin();
	std::set<double>::iterator s2 = s.upper_bound(x);
	if (s2 == s.end()) return s.size()-1;
	int count=-1;
	do {
		count++;
	} while (s1++ != s2);

	return count-1;
}

Matrix<double> KnotSet::ComputeLeastSquaresMatrix(const Vector<double>& tau, int m) const
{
	int dim = ord + n - 2;
	Matrix<double> vmat(m,dim);
	for (int i=0; i<m; i++) {
		Vector<double> bval(CreateVectorInterp(tau[i]));
		vmat.InsertRow(bval,i);
	}
	
	return vmat;
}


Vector<double> KnotSet::ComputeDistinctKnots() const
{	
	std::set<double> set1((*kts).begin(),(*kts).end());
	Vector<double> v(set1.size());
	std::copy(set1.begin(),set1.end(),v.begin());
	return v;
}


Vector<int> KnotSet::ComputeMultiplicities() const
{
	std::set<double> set1((*kts).begin(),(*kts).end());
	std::set<double>::iterator p = set1.begin();
	Vector<int> v(set1.size());
	for (unsigned int i=0; i<set1.size(); i++) v[i] = (*mset).count(*p++);

	return v;
}


Vector<double> KnotSet::ComputeKnots() const
{
	Vector<double> v(FindNumKnots());

	int count=0;
	for (int i=0; i<(*mult).size(); i++)
		for (int j=0; j<(*mult)[i]; j++) {
			v[count] = (*dts)[i];
			count++;
		}
	return v;
}

   
int KnotSet::FindNumKnots() const
{
	int count = 0;
	for (int i=0; i<(*mult).size(); i++) count += (*mult)[i];
	return count;
}


int KnotSet::FindNumDistinctKnots() const
{
	return std::set<double>((*kts).begin(),(*kts).end()).size();
}



// CHECK FUNCTIONS
// checks to see that knot set has multiplicity ord at both ends
bool KnotSet::CheckKnotSet() const
{
	int s1=FindMultiplicity((*kts)[ord-1]);

	int s2=FindMultiplicity((*kts)[num-ord]);

	if (s1 < ord || s2 < ord) return true;
	else return false;
}


// updates the knot set to multiplicity ord at left end
KnotSet& KnotSet::UpdateKnotsLeft() 
{
	Vector<double> newknot(num);

	for (int i=0; i<ord; i++) newknot[i]=(*kts)[ord-1];
	for (int i=ord; i<num; i++) newknot[i]=(*kts)[i];
	*this = KnotSet(newknot,ord,num);
	return *this;
}

// updates knot set to multiplicity ord at right end
KnotSet& KnotSet::UpdateKnotsRight() 
{
	Vector<double> newknot(num);

	for (int i=0; i<num-ord; i++) newknot[i] = (*kts)[i];
	for (int i=num-ord; i<num; i++) newknot[i]=(*kts)[num-ord];
	*this = KnotSet(newknot,ord,num);
	return *this;
}


KnotSet KnotSet::UpdateKnotSet()
{
	// find the multiplicity of the left end knot of index ord-1
	int s1= FindMultiplicity((*kts)[ord-1]);

	// find the multipliicty of the right end knot of index num
	int s2= FindMultiplicity((*kts)[num-ord]);
	
	// update the control points and knot sets
	if (s1 < ord) {
		UpdateKnotsLeft();
	}
	if (s2 < ord) {
		UpdateKnotsRight();
	}
	// return the updated Cpoint object
	return *this;
}



// NORMALISATION FUNCTIONS

// normalise the knot set with with respect to [0,1]
KnotSet KnotSet::Normalise() const
{
	Vector<double> nkts(num);
	for (int i=0; i<num; i++) nkts[i] = (*kts)[i]/(*kts)[num-ord];
	return KnotSet(nkts,ord,num);
}


KnotSet KnotSet::Normalise(double a, double b) const
{
	Vector<double> nkts(*kts);

	// create new knots
	double scale = (b-a)/((*kts)[num-ord]-(*kts)[ord-1]);
	for (int i=ord; i<num; i++) nkts[i] = ((*kts)[i]-(*kts)[ord-1])*scale+(*kts)[ord-1];

	double dist = (*kts)[ord-1]-a;
	for (int i=0; i<num; i++) nkts[i] = nkts[i]-dist;
	
	return KnotSet(nkts,ord,num);	
}


// normalise the knot set with respect to the knot set kset
KnotSet KnotSet::Normalise(const KnotSet& kset) const
{
	Vector<double> nkts(*kts);

	// create new knots
	int num1 = kset.GetNum();
	int ord1 = kset.GetOrd();

	double scale = (kset.GetKnots()[num1-ord1]-kset.GetKnots()[ord1-1])/((*kts)[num-ord]-(*kts)[ord-1]);
	for (int i=ord; i<num; i++) nkts[i] = ((*kts)[i]-(*kts)[ord-1])*scale+(*kts)[ord-1];


	double dist = (*kts)[ord-1]-kset.GetKnots()[ord1-1];
	for (int i=0; i<num; i++) nkts[i] = nkts[i]-dist;

	return KnotSet(nkts,ord,num);
}



// TEST FOR EQUALITY
bool KnotSet::IsSameKnotSet(const KnotSet& kset) const
{
	// questionable ???
	KnotSet nkset = Normalise(kset);

	if (nkset.GetKnotSet() == kset.GetKnotSet()) return true;
	else return false;
}

// TEST FOR EQUALITY
bool KnotSet::IsSameDistinctKnotSet(const KnotSet& kset) const
{
	// questionable ???
	KnotSet k1(*this), k2(kset);

//	KnotSet nkset = Normalise(kset);

	Vector<double> ndts1(GetDistinctKnots());
	Vector<double> ndts2(kset.GetDistinctKnots());

	std::set<double> s1(ndts1.begin(),ndts1.end());
	std::set<double> s2(ndts2.begin(),ndts2.end());
	
	if (s1 == s2) {  return true;}
	else { return false;}
}


KnotSet KnotSet::KnotUnion(const KnotSet& kset) const
{
	// questionable ??
//	KnotSet nkset = Normalise(kset);
	KnotSet nkset(*this);
	std::multiset<double> set3;

	// added 2/2/03
	if (kset.GetNumDistinct() > 2) {
		Vector<double> vec(kset.GetNum()-2*kset.GetOrd());
		for (int i=0; i<kset.GetNum()-2*kset.GetOrd(); i++) vec[i] = kset.GetKnots()[i+kset.GetOrd()];
	
	
		std::multiset<double> set1=nkset.GetKnotSet();
		std::multiset<double> set2(vec.begin(),vec.end());
		std::set_union(set1.begin(),set1.end(),set2.begin(),set2.end(),std::inserter(set3,set3.begin()));
	
		Vector<double> nkts(set3.size());
		std::copy(set3.begin(),set3.end(),nkts.begin());
		

		return KnotSet(nkts,ord,set3.size());
	} else return nkset;
}
		


// SPECIFIC KNOT SETS
KnotSet KnotSet::CreateKnotSetCompBezCurv() const
{
	Vector<int> nmult(n);
	Vector<double> ndts(*dts);
	for (int i=0; i<n; i++) nmult[i]=ord;
	return KnotSet(ndts,nmult,ord,n);
}


// create knot set representing subdivision of curve to level
KnotSet KnotSet::CreateKnotSetSubdivide(int level) const
{
	double x, x1, x2;
	double step;
	KnotSet kset(*this);
	
	// insert knots into each interval
	for (int i=0; i<=n-2; i++) {
		x1=(*dts)[i];
		x2=(*dts)[i+1];
		step=(x2-x1)/level;
		for (int j=1; j<=level-1; j++) {
			x=x1+j*step;
			kset=kset.CreateKnotSetInsert(x);
		}
	}
	return kset;
}



// create knot set reprsenting indefinite integral of curve
KnotSet KnotSet::CreateKnotSetIntegrate() const
{
	Vector<double> nkts(num+2);

	for(int i=0;i<ord+1; i++)
		nkts[i]=(*kts)[0];

	for(int i=ord+1; i<num-ord+1; i++)
		nkts[i]=(*kts)[i-1];

	for(int i=num+1-ord; i<num+2; i++)
		nkts[i]=(*kts)[num-1];

	return KnotSet(nkts,ord+1,num+2);
}

// create knot set reprsenting subdivision at x1 x2
// if x1, x2 have multiplicity ord there is nothing to do
KnotSet KnotSet::CreateKnotSetSubdivide(double x1, double x2) const
{
	// find correct position to insert knot x1
	int s1 = FindMultiplicity(x1);
	int s2 = FindMultiplicity(x2);

	KnotSet kset1=CreateKnotSetInsert(x1,ord-s1);
	// find correct position to insert knot x2

	KnotSet kset2=kset1.CreateKnotSetInsert(x2,ord-s2);

	//ofs.close();
	//extract the correct section of curve
	int ind1=kset2.Find_index(x1);
	int ind2=kset2.Find_index(x2);

	// detect if x2 is the end knot
	std::multiset<double>::iterator p = (*mset).begin();
	for (unsigned int i=0; i<(*mset).size()-1; i++) p++;
	if (x2 == *p) ind2 = ind2 + ord;
	
	Vector<double> newknot(ind2-ind1+ord);
	
	for (int i=ind1-ord+1; i<=ind2; i++) newknot[i-ind1+ord-1] = kset2.GetKnots()[i];
	return KnotSet(newknot,ord,ind2-ind1+ord);
}



// create knot set for derivative curve 
KnotSet KnotSet::CreateKnotSetDeriv(int level) const
{
	if (level <= 0) return *this;

	Vector<int> nmult(n);
	Vector<double> nt(n);

	Vector<int> over(n);
	over[0]=0;
	over[n-1]=0;

	for (int i=1;i<n-1;i++) over[i]=ord-(*mult)[i]-level;

	int norder=ord-level;
	nmult[0]=norder;
	nt[0]=(*kts)[0];
	nmult[n-1]=norder;
	nt[n-1]=(*kts)[num-ord];

	for(int i=1;i<n-1;i++) {
		if(over[i]<=0) nmult[i]=norder;
		else nmult[i]=norder-over[i];
		nt[i]=(*dts)[i];
	}

	return KnotSet(nt, nmult, norder, n);
}



// create the knot set for knot insertion 
KnotSet KnotSet::CreateKnotSetInsert(double x) const
{
	// find correct position to insert new knot 
	int p = Fndint1(x);
	int s = FindMultiplicity(x);

	// check to see that multiplicity is <= ord
	if (s > ord-1) return *this;
        
	// calculate new knot set 
	std::vector<double> newknot(num+1);
	std::copy((*kts).begin(),(*kts).end(),newknot.begin());

	std::vector<double>::iterator q = newknot.begin()+p+1;
	newknot.insert(q,x);
	newknot.pop_back();

	return KnotSet(newknot,ord,num+1);
}


// create knot set for knot insertion of x level times
KnotSet KnotSet::CreateKnotSetInsert(double x, int level) const
{
	// find correct position to insert new knot 
	if (level <= 0) return *this;
	int p = Fndint1(x);
	int s = FindMultiplicity(x);

	if (s > ord-level) return *this;
	

	std::vector<double> newknot(num+level);
	std::copy((*kts).begin(),(*kts).end(),newknot.begin());


	std::vector<double>::iterator q = newknot.begin()+p+1;

	newknot.insert(q,level,x);
	// added sun 3rd Feb -level to -(level+1)
	newknot.erase(newknot.end()-level,newknot.end());

	return KnotSet(newknot,ord,num+level);
}


// create knot set for knot insertion of x level times
KnotSet KnotSet::CreateKnotSetInsert(const Vector<double>& Kts, int Num) const
{
	if (Num <= 0) return *this;
	KnotSet kset(*this);
	for (int i=0; i<Num; i++) kset = kset.CreateKnotSetInsert(Kts[i]);
	return kset;
}


KnotSet KnotSet::CreateKnotSetElevate(int level) const
{
	if (level <= 0) return *this;
	Vector<int> nmult(n);
	Vector<double> ndts(*dts);
	for (int i=0; i<n; i++) nmult[i]=(*mult)[i]+level;
	return KnotSet(ndts,nmult,ord+level,n);
}

// create knot set for knot insertion of x level times
KnotSet KnotSet::CreateKnotSetInsert(const Vector<double>& Kts, const Vector<int>& mult, int N) const
{
	if (N <= 0) return *this;
	KnotSet kset(*this);
	for (int i=0; i<N; i++) kset = kset.CreateKnotSetInsert(Kts[i],mult[i]);
	return kset;
}

KnotSet KnotSet::CreateKnotSetRemoval(double knot) const
{
	// find position of knot for removal
	int p = Fndint1(knot);
	int s = FindMultiplicity(knot);

	// compute i, j values 
	int i=p+s-ord+1;
	int j=p;

		// claculate new control point values
	while (j-i > 0)
	{
		i++;
		j--;
	}
	// create CPoints object with new control points
	int num_eqns=ord-s;
	// remove according to whether number of equations is even or odd
	div_t div_result;
	div_result = div(num_eqns,2);
	int ind1, ind2;
	if (div_result.rem == 0) {
		i--;
		j++;
		ind1 = j;
		ind2 = p+s;
	} else {
		ind1 = i;
		ind2 = p+s;
	}

	Vector<double> kts1(num-1);
	for (int i=0; i<ind2; i++) kts1[i] = (*kts)[i];
	for (int i=ind2+1; i<num; i++) kts1[i-1]=(*kts)[i];
	
	return KnotSet(kts1,ord,num-1);
}

KnotSet KnotSet::CreateKnotSetProduct(const KnotSet& kset) const
{
	// normalise knot sets
	KnotSet nkset = Normalise(kset);

	int n1=kset.GetNumDistinct();
	int ord1=kset.GetOrd();
	Vector<SortCollection> v(n+n1);

	int count=0;
	for (int i=0; i<n; i++) {
		v[i].kts = 1;
		v[i].dts = nkset.GetDistinctKnots()[i];
		v[i].mult= nkset.GetMult()[i]; 
	}

	for (int i=0;i<n1; i++) {
		v[n+i].kts=2;
		v[n+i].dts= kset.GetDistinctKnots()[i];
		v[n+i].mult= kset.GetMult()[i];
	}
	
	std::multiset<SortCollection, std::less<SortCollection> > s(v.begin(),v.end());

	Vector<double> dts2(n+n1);
	Vector<int> mult2(n+n1);

	Vector<SortCollection> v1(s.size());
	std::copy(s.begin(),s.end(),v1.begin());
	
	count=0;
	int i=0;
	int temp1, temp2;
	do {
		if (v1[i].dts==v1[i+1].dts) {
			dts2[count]=v1[i].dts;
			if (v1[i].kts==1) {
				temp1=ord1-1+v1[i].mult;
				temp2=ord-1+v1[i+1].mult;
			} else {
				temp1=ord-1+v1[i].mult;
				temp2=ord1-1+v1[i+1].mult;
			}
			if (temp1 > temp2) 
				mult2[count]=temp1;
			else mult2[count]=temp2;
			count++;
			i+=2;
		} else if (v1[i].kts==1) {
			dts2[count]=v1[i].dts;
			mult2[count]=ord1-1+v1[i].mult;
			count++;
			i++;
		} else if (v1[i].kts==2) {
			dts2[count]=v1[i].dts;
			mult2[count]=ord-1+v1[i].mult;
			count++;
			i++;
		}
	} while(i<n+n1-1);

	Vector<double> fdts(count);
	Vector<int> fmult(count);

	std::copy(dts2.begin(),dts2.begin()+count,fdts.begin());
	std::copy(mult2.begin(),mult2.begin()+count,fmult.begin());

	return KnotSet(fdts,fmult,ord+ord1-1,count);
}


// MATRIX FOR DERIVATIVES
	
// create the matrix representing derivative curve control points
// defined in terms of original control points
Matrix<double> KnotSet::CreateMatrixDeriv() const
{
	// set up temporary arrays for calculations 
	Vector<double> dm(ord);
	Vector<double> dp(ord);
	Vector<double> tp1(ord);
	Vector<double> tp2(ord);
	

	int norder = ord-1;

	Vector<int> over(n);
	over[0] = 0;
	over[n-1] = 0;
	
	int k1=1;

	for(int i=1;i<n-1;i++) over[i] = ord-(*mult)[i]-1;

	// create control points
	KnotSet nkset = CreateKnotSetDeriv(1);
	int ns = nkset.GetNum()-norder;

	Matrix<double> d(ns,num-ord,0.0);

	double that;
	int ind;

	for(int l=0;l<n-1;l++) {
		that=(*dts)[l];
		ind=Find_index(that);

		for(int j=0;j<ord;j++) {
			dp[j] = (*kts)[ind+j+1]-that;
			dm[j] = that-(*kts)[ind-ord+j+1];
		}

		// compute the points according to de Boor algorithm 
	
		for(int j=1;j<ord;j++) {
			tp1[j] = 1/(dm[j] + dp[j-1]);
			tp2[j] = -1/(dm[j] + dp[j-1]); 
		}	
		if(over[l]>0) {	k1-=over[l];}

		for(int j=1;j<ord;j++) {
			d[k1-1][k1-1] = tp2[j];
			d[k1-1][k1] = tp1[j]; 
			k1++;
		}
	}
	
	for(int i=0;i<ns;i++)
		for(int j=0;j<num-ord;j++)
			d[i][j] = d[i][j]*(ord-1);

	return d;
}


// create the matrix representing derivative curve control points
// defined in terms of original control points
Vector<double> KnotSet::CreateVectorInterp(double val) const
{

	Matrix<double> v(ord,ord,0.0);
	Vector<double> dp(ord,0.0), dm(ord,0.0);

	int ind = Find_index(val);

	v[0][0] = 1.0;
	for (int j=0; j<ord-1; j++) {
		dp[j] = (*kts)[ind+j+1]-val;
		dm[j] = val - (*kts)[ind+1-j-1];
		for (int i=0; i<=j; i++) {
			double m = v[i][j]/(dp[i]+dm[j-i]);
			v[i][j+1] = v[i][j+1]+dp[i]*m;
			v[i+1][j+1]=dm[j-i]*m;
		}
	}
	Vector<double> res(num-ord,0.0);
	for (int i=0; i<ord; i++) res[i+ind-ord+1] = v[i][ord-1];
	return res;
}


// create the matrix representing derivative curve control points
// defined in terms of original control points
Vector<double> KnotSet::CreateVectorInterpDeriv(int level, double val) const
{
	Matrix<double> m = CreateMatrixDeriv(level);
	KnotSet kset = CreateKnotSetDeriv(level);
	Vector<double> v = kset.CreateVectorInterp(val);
	return Math::mult4(v,m);
}



	
// create the matrix representing derivative curve control points
// defined in terms of original control points
Matrix<double> KnotSet::CreateMatrixDeriv(int level) const
{	
	if (level <= 0) return Math::ComputeIdentityMatrix(num-ord);
	// create matrix for first derivative
	Matrix<double> d = CreateMatrixDeriv();

	// create knot sets for each subsequent derivative up to deriv
	// find new matrix and multiply by previous matrix.
	Vector<double> nkts(*kts);
	KnotSet kset = KnotSet(nkts,ord,num);
	for (int i=1; i<level; i++) {
		kset = kset.CreateKnotSetDeriv(1);
		d = Math::mult1(kset.CreateMatrixDeriv(),d);
	}	
	// return matrix of linear combinations
	return d;
}



// create the matrix representing derivative curve control points
// defined in terms of original control points
Vector<double> KnotSet::ComputeKnotSetAver() const
{	
	// index should be between 1 and num-ord;
	Vector<double> av(num-ord);
	double sum=0.0;
	for (int i=0; i<num-ord; i++) {
		for (int j=i+1; j<=i+ord-1; j++) sum=sum+(*kts)[j];
		av[i] = sum/(double)(ord-1);
		sum=0.0;
	}
	return av;
}

	
// create the matrix representing derivative curve control points
// defined in terms of original control points
KnotSet KnotSet::CreateKnotSetCurveBasis(int index) const
{	
	// index should be between 1 and num-ord;
	Vector<double> knots(ord+1);
	Vector<double> r(*kts);
	std::vector<double> sv = r.convert();

	//Vector<double>::iterator r = (*kts).begin();
	std::vector<double>::iterator vi(sv.begin());
	//std::vector<double> r((*kts).begin());
	std::copy(vi+index-1,vi+ord+index, knots.begin());
	return KnotSet(knots,ord,ord+1);
}



// READ and WRITE
void KnotSet::write(std::ostream& os) 
{
	os << "Knot Set\n";
	os << *kts;
	os << "distinct knots\n";
	os << *dts;
	os << "multiplicities\n";
	os << *mult;
}


void KnotSet::read(std::istream& is)
{
	int Ord;
	std::cout << "order of Bspline curve";
	is >> Ord;
	int Num;
	std::cout << "number of control points";
	is >> Num;
	Vector<double> Kts(Ord+Num);
	std::cout << "input knots";
	is >> Kts;
	*this = KnotSet(Kts,Ord,Ord+Num);
} 


void KnotSet::writefile(std::ofstream& ofs) 
{
	ofs << "ord, num\n";
	ofs << ord << num;
	ofs << "Knots\n";
	ofs << *kts;
	ofs << "distinct knots\n";
	ofs << *dts;
	ofs << "multiplicities\n";
	ofs << *mult;
}


void KnotSet::readfile(std::ifstream& ifs)
{
	int Ord, Num;
	ifs >> Ord >> Num;
	Vector<double> Kts(Ord+Num);
	ifs >> Kts;
	*this = KnotSet(Kts,Ord,Ord+Num);
} 



void SortCollection::write(std::ostream& os)
{
	os << "Knot Set\n";
	os << kts;
	os << "distinct knots\n";
	os << dts;
	os << "multiplicities\n";
	os << mult;
}


void SortCollection::read(std::istream& is)
{
	
} 


void SortCollection::writefile(std::ofstream& ofs) 
{
	ofs << "Knots\n";
	ofs << kts;
	ofs << "distinct knots\n";
	ofs << dts;
	ofs << "multiplicities\n";
	ofs << mult;
}


void SortCollection::readfile(std::ifstream& ifs)
{
	
} 