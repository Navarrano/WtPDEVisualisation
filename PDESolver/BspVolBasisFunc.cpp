
#include "bspvol.h"
#include "mathfunctions.h"
#include "compbezvol.h"
#include "bspsurf.h"


// CONSTRUCTORS

// default constructor
BspVolBasisFunc::BspVolBasisFunc() : ordu(0), ordv(0), ordw(0), ktsu(), ktsv(), ktsw() {}

// constructor taking orders amd knot sets in u and v
BspVolBasisFunc::BspVolBasisFunc(int Ordu, int Ordv, int Ordw, const Vector<double>& Ktsu, const Vector<double>& Ktsv,
	const Vector<double>& Ktsw) : 
	ordu(Ordu), ordv(Ordv), ordw(Ordw), ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)), ktsw(new Vector<double>(Ktsw)), 
	b(new BspVol<double>(CreateBspVol())) { }


// ACCESS FUNCTIONS

// compute the dimension of the u knot set
int BspVolBasisFunc::ComputeDimU() const
{
	// create KnotSet object 
	KnotSet kset(*ktsu,ordu,ordu+1);
	int num = kset.GetNumDistinct();
	int num_cond=0;
	// find the dimension
	for (int i=1; i<num-1; i++) num_cond = num_cond + ordu - kset.GetMult()[i];
	return ((num-1)*ordu-num_cond);
}

// compute the dimension of the v knot set
int BspVolBasisFunc::ComputeDimV() const
{
	// create KnotSet object
	KnotSet kset(*ktsv,ordv,ordv+1);
	int num = kset.GetNumDistinct();
	int num_cond=0;
	// find the dimension
	for (int i=1; i<num-1; i++) num_cond = num_cond + ordv - kset.GetMult()[i];
	return ((num-1)*ordv-num_cond);
}


// compute the dimension of the v knot set
int BspVolBasisFunc::ComputeDimW() const
{
	// create KnotSet object
	KnotSet kset(*ktsw,ordw,ordw+1);
	int num = kset.GetNumDistinct();
	int num_cond=0;
	// find the dimension
	for (int i=1; i<num-1; i++) num_cond = num_cond + ordw - kset.GetMult()[i];
	return ((num-1)*ordw-num_cond);
}


// evaluate the basis function at u,v
BspVol<double> BspVolBasisFunc::GetBspVol() const
{    
	return *b;
}

// create the BspVol representation of the basis function
BspVol<double> BspVolBasisFunc::CreateBspVol() const
{
	// compute dimension in u and v
	int dimu = ComputeDimU();
	int dimv = ComputeDimV();
	int dimw = ComputeDimW();
	
	Matrix3D<double> cpts(dimu,dimv,dimw,0.0);

	KnotSet ksetu(*ktsu,ordu,ordu+1);
	KnotSet ksetv(*ktsv,ordv,ordv+1);
	KnotSet ksetw(*ktsw,ordw,ordw+1);	
	
	int numu = ksetu.GetNumDistinct();
	int numv = ksetv.GetNumDistinct();
	int numw = ksetw.GetNumDistinct();
	
	// create distinct knot and multiplicity vectors 
	Vector<double> dts1(ksetu.GetDistinctKnots());
	Vector<int> mult1(ksetu.GetMult());
	
	Vector<double> dts2(ksetv.GetDistinctKnots());
	Vector<int> mult2(ksetv.GetMult());

	Vector<double> dts3(ksetw.GetDistinctKnots());
	Vector<int> mult3(ksetw.GetMult());

	// assign multiplicity ordu,ordv at the two ends
	mult1[0]=ordu;
	mult1[numu-1]=ordu;

	mult2[0]=ordv;
	mult2[numv-1]=ordv;

	mult3[0]=ordw;
	mult3[numw-1]=ordw;

	// compute offsets
	int start1=ksetu.GetMult()[0];
	int start2=ksetv.GetMult()[0];
	int start3=ksetw.GetMult()[0];

	KnotSet kset1(dts1,mult1,ordu,numu);
	KnotSet kset2(dts2,mult2,ordv,numv);
	KnotSet kset3(dts3,mult3,ordw,numw);


	// assign 1.0 to appropriate control point
	cpts[ordw-start3][ordu-start1][ordv-start2]=1.0;
	
	//  create and return the BspVol
	return BspVol<double>(cpts,kset1.GetKnots(),kset2.GetKnots(),kset3.GetKnots(),ordu,ordv,ordw,dimu,dimv,dimw);
}



// EVALUATORS

// evaluate the basis function at u,v
double BspVolBasisFunc::Eval(double u, double v, double w) const
{    
	// create the BspVol for the basis function and evaluate it at u,v
	if (u < (*ktsu)[0] || u > (*ktsu)[ordu]) return 0;
	if (v < (*ktsv)[0] || v > (*ktsv)[ordv]) return 0;
	if (w < (*ktsw)[0] || w > (*ktsw)[ordw]) return 0;
	return (*b).Eval1(u,v,w);
}

// evaluate the basis function at u,v
double BspVolBasisFunc::operator()(double u, double v, double w) const
{    
	// create the BspVol for the basis function and evaluate it at u,v
	if (u < (*ktsu)[0] || u > (*ktsu)[ordu]) return 0;
	if (v < (*ktsv)[0] || v > (*ktsv)[ordv]) return 0;
	if (w < (*ktsw)[0] || w > (*ktsw)[ordw]) return 0;
	return (*b).Eval1(u,v,w);
}


// evaluate the basis function at u,v
double BspVolBasisFunc::operator()(int valu, int valv, int valw, double u, double v, double w) const
{    
	return (*b).Derive(valu,valv,valw,u,v,w);	
}


double BspVolBasisFunc::Derive(int m, int n, int p, double u, double v, double w) const
{
	return (*b).Derive(m,n,p,u,v,w);
}

double BspVolBasisFunc::GetLeftLimitU() const
{
	return (*ktsu)[0];
}

double BspVolBasisFunc::GetRightLimitU() const
{
	return (*ktsu)[ordu];
}

double BspVolBasisFunc::GetLeftLimitV() const
{
	return (*ktsv)[0];
}

double BspVolBasisFunc::GetRightLimitV() const
{
	return (*ktsv)[ordv];
}

double BspVolBasisFunc::GetLeftLimitW() const
{
	return (*ktsw)[0];
}

double BspVolBasisFunc::GetRightLimitW() const
{
	return (*ktsw)[ordw];
}


Vector<double> BspVolBasisFunc::GetKnotsU() const
{
	return *ktsu;
}


Vector<double> BspVolBasisFunc::GetKnotsV() const
{
	return *ktsv;
}


Vector<double> BspVolBasisFunc::GetKnotsW() const
{
	return *ktsw;
}



// READ and WRITE

void BspVolBasisFunc::write(std::ostream& os) 
{
	os << "Bspline Volume basis function\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "order in w is " << ordw << "\n";
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are \n";
	os << *ktsv;
	
	os << "\nknots in w are \n";
	os << *ktsw;
	
	os << "\nBspVol representation\n";
	os << *b;
}


void BspVolBasisFunc::read(std::istream& is)
{
	int Ordu, Ordv, Ordw;
	std::cout << "order of Bspline Volume basis in u and v and w";
	is >> Ordu >> Ordv >> Ordw;
	Vector<double> Ktsu(Ordu+1), Ktsv(Ordv+1), Ktsw(Ordw+1);
	std::cout << "knots in u";
	is >> Ktsu;
	std::cout << "knots in v";
	is >> Ktsv;
	std::cout << "knots in w";
	is >> Ktsw;
	
	*this = BspVolBasisFunc(Ordu,Ordv,Ordw,Ktsu,Ktsv,Ktsw);
} 


void BspVolBasisFunc::writefile(std::ofstream& ofs)
{
	ofs << "Bspline Volume basis function\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "order in w is " << ordw << "\n";
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are \n";
	ofs << *ktsv;
	
	ofs << "\nknots in w are \n";
	ofs << *ktsw;
	ofs << "\nBspVol basis reprsentation\n";
	ofs << *b;
}


void BspVolBasisFunc::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw;
	ifs >> Ordu >> Ordv >> Ordw;

	Vector<double> Ktsu(Ordu+1), Ktsv(Ordv+1), Ktsw(Ordw+1);
	ifs >> Ktsu;
	ifs >> Ktsv;
	ifs >> Ktsw;
	
	*this = BspVolBasisFunc(Ordu,Ordv,Ordw,Ktsu,Ktsv,Ktsw);
} 


// CONSTRUCTORS

// default constructor
BspVolBasisFuncSet::BspVolBasisFuncSet() : ordu(0), ordv(0), ordw(0), numu(0), numv(0), numw(0), ktsu(), ktsv(), ktsw(), b() {}


// constructor taking order in u, order in v, knots in u and v
BspVolBasisFuncSet::BspVolBasisFuncSet(int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw, 
const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw) : 
	ordu(Ordu), ordv(Ordv), ordw(Ordw), numu(Numu), numv(Numv), numw(Numw), 
		ktsu(new Vector<double>(Ktsu)), ktsv(new Vector<double>(Ktsv)), ktsw(new Vector<double>(Ktsw)), 
		b(new Matrix3D<BspVolBasisFunc>(Numu-Ordu,Numv-Ordv,Numw-Ordw,BspVolBasisFunc()))
{
	Vector<double> v1(ordu+1), v2(ordv+1), v3(ordw+1);
	for (int k=0; k<numw-ordw; k++) {
		for (int i=0; i<numu-ordu; i++) {
			for (int j=0; j<numv-ordv; j++) {
				for (int l=0; l<=ordu; l++) v1[l] = (*ktsu)[i+l];
				for (int m=0; m<=ordv; m++) v2[m] = (*ktsv)[j+m];
				for (int n=0; n<=ordw; n++) v3[n] = (*ktsw)[k+n];
				(*b)[k][i][j] = BspVolBasisFunc(ordu,ordv,ordw,v1,v2,v3);
			}
		}
	}
}


// EVALUATORS
// evaluate the basis function at u, v. Creates the BezVol
// representation and evaluates it
Matrix3D<double> BspVolBasisFuncSet::Eval(double u, double v, double w) const
{    
	Matrix3D<double> mat(numu-ordu,numv-ordv,numw-ordw);

	for (int k=0; k<numw-ordw; k++)
		for (int i=0; i<numu-ordu; i++) 
			for (int j=0; j<numv-ordv; j++) 
				mat[k][i][j] = (*b)[k][i][j].Eval(u,v,w);

	return mat;
}


// evaluate the basis function at u, v. Creates the BezVol
// representation and evaluates it
Matrix3D<double> BspVolBasisFuncSet::operator()(double u, double v, double w) const
{    
	Matrix3D<double> mat(numu-ordu,numv-ordv,numw-ordw);

	for (int k=0; k<numw-ordw; k++)
		for (int i=0; i<numu-ordu; i++) 
			for (int j=0; j<numv-ordv; j++) 
				mat[k][i][j] = (*b)[k][i][j].Eval(u,v,w);

	return mat;
}

Matrix3D<double> BspVolBasisFuncSet::CreateIntegral(double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<double> res(numu-ordu,numv-ordv,numw-ordw);

	for (int k=0; k<numw-ordw; k++)
		for (int i=0; i<numu-ordu; i++) 
			for (int j=0; j<numv-ordv; j++) {
				res[k][i][j] = (*b)[k][i][j].GetBspVol().Integrate(u1,u2,v1,v2,w1,w2);
			}

	return res;
}


Matrix<double> BspVolBasisFuncSet::CreateMatrixKnotAveragesU() const
{
	return BspCurvBasisFuncSet(*ktsu,ordu,numu).CreateMatrixKnotAverages();
}


Matrix<double> BspVolBasisFuncSet::CreateMatrixKnotAveragesV() const
{
	return BspCurvBasisFuncSet(*ktsv,ordv,numv).CreateMatrixKnotAverages();
}


Matrix<double> BspVolBasisFuncSet::CreateMatrixKnotAveragesW() const
{
	return BspCurvBasisFuncSet(*ktsw,ordw,numw).CreateMatrixKnotAverages();
}


Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisationU(int levu, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(levu,0,0,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(levu);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(0);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(0);
	
	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int l=0; l<derivv.GetNumRows(); l++) {
					for (int m=0; m<derivv.GetNumRows(); m++)
						sum = sum + mat[k][i][l][m][l][j];

					final[k][i][j] = sum;
					sum=0.0;
				}

	return Math::mult5((Math::mult6(Math::transpose(derivu),final)),derivu);
}


Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisationV(int levv, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(0,levv,0,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(0);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(levv);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(0);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int r=0; r<derivw.GetNumRows(); r++) {
					for (int l=0; l<derivv.GetNumRows(); l++)
						for (int m=0; m<derivv.GetNumRows(); m++)
							for (int n=0; n<derivv.GetNumCols(); n++)
								sum = sum + mat[k][i][l][r][j][m]*derivv[l][n]*derivv[m][n];

					final[k][i][j] = sum;
					sum=0.0;
				}

	return final;
}


Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisationV1(int levv, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(0,levv,0,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(0);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(levv);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(0);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<derivw.GetNumRows(); k++) {
					for (int l=0; l<derivv.GetNumRows(); l++)
						for (int m=0; m<derivv.GetNumRows(); m++)
							for (int n=0; n<derivv.GetNumCols(); n++)
								sum = sum + mat[k][i][l][k][j][m]*derivv[l][n]*derivv[m][n];

					final[k][i][j] = sum;
					sum=0.0;
				}

	return final;
}


Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisationW(int levw, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(0,0,levw,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(0);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(0);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(levw);

	Matrix<double> derivwt = Math::transpose(derivw);

	Matrix<double> derivwm = Math::mult1(derivw,derivwt);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++)
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int l=0; l<derivw.GetNumRows(); l++) {
					for (int r=0; r<derivu.GetNumRows(); r++)
						for (int s=0; s<derivw.GetNumRows(); s++)
							sum = sum + mat[l][i][r][s][r][j]*derivwm[l][k];

					final[k][i][j] = sum;
					sum=0.0;
				}

	return final;
}


Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisationUV(int levu, int levv, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(levu,levv,0,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(levu);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(levv);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(0);

	Matrix<double> derivwt = Math::transpose(derivw);

	Matrix<double> derivwm = Math::mult1(derivw,derivwt);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int r=0; r<derivw.GetNumRows(); r++) {
					for (int l=0; l<derivv.GetNumRows(); l++)
						for (int m=0; m<derivv.GetNumRows(); m++)
							for (int n=0; n<derivv.GetNumCols(); n++)
								sum = sum + mat[k][i][l][r][j][m]*derivv[l][n]*derivv[m][n];

					final[k][i][j] = sum;
					sum=0.0;
				}

	return Math::mult5((Math::mult6(Math::transpose(derivu),final)),derivu);
}



Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisationUW(int levu, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(levu,0,levw,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(levu);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(0);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(levw);

	Matrix<double> derivwt = Math::transpose(derivw);

	Matrix<double> derivwm = Math::mult1(derivw,derivwt);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++)
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int l=0; l<derivw.GetNumRows(); l++) {
					for (int r=0; r<derivu.GetNumRows(); r++)
						for (int s=0; s<derivw.GetNumRows(); s++)
							sum = sum + mat[l][i][r][s][r][j]*derivwm[l][k];

					final[k][i][j] = sum;
					sum=0.0;
				}

	return Math::mult5((Math::mult6(Math::transpose(derivu),final)),derivu);
}

Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisationVW(int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(0,levv,levw,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(0);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(levv);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(levw);

	Matrix<double> derivwt = Math::transpose(derivw);

	Matrix<double> derivwm = Math::mult1(derivw,derivwt);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int r=0; r<derivw.GetNumRows(); r++) {
					for (int l=0; l<derivv.GetNumRows(); l++)
						for (int m=0; m<derivv.GetNumRows(); m++)
							for (int n=0; n<derivv.GetNumCols(); n++)
								for (int p=0; p<derivw.GetNumRows(); p++) 
								sum = sum + mat[p][i][l][r][j][m]*derivv[l][n]*derivv[m][n]*derivwm[p][k];

					final[k][i][j] = sum;
					sum=0.0;
				}

	return final;
}


Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisation(int levu, int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(levu,levv,levw,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(levu);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(levv);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(levw);

	Matrix<double> derivwt = Math::transpose(derivw);

	Matrix<double> derivwm = Math::mult1(derivw,derivwt);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int r=0; r<derivw.GetNumRows(); r++) {
					for (int l=0; l<derivv.GetNumRows(); l++)
						for (int m=0; m<derivv.GetNumRows(); m++)
							for (int n=0; n<derivv.GetNumCols(); n++)
								for (int p=0; p<derivw.GetNumRows(); p++) 
									sum = sum + mat[p][i][l][r][j][m]*derivv[l][n]*derivv[m][n]*derivwm[p][k];

					final[k][i][j] = sum;
					sum=0.0;
				}
				std::cerr << final;

	return Math::mult5((Math::mult6(Math::transpose(derivu),final)),derivu);
}


Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisation1(int levu, int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(levu,levv,levw,u1,u2,v1,v2,w1,w2);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(levu);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(levv);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(levw);
	
	Matrix<double> derivwt = Math::transpose(derivw);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),numw-ordw,0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<numw-ordw; k++) {
					for (int l=0; l<derivv.GetNumRows(); l++)
						for (int m=0; m<derivw.GetNumRows(); m++)
							for (int n=0; n<derivv.GetNumRows(); n++)
								for (int p=0; p<derivw.GetNumRows(); p++) 
									for (int q=0; q<derivv.GetNumCols(); q++) 
										sum = sum + mat[m][i][l][p][j][n]*derivv[l][q]*derivv[n][q]*derivw[m][k]*derivw[p][k];

					final[k][i][j] = sum;
					sum=0.0;
				}


	return Math::mult5((Math::mult6(Math::transpose(derivu),final)),derivu);
}

Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisation1(int levu, int levv, int levw) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(levu,levv,levw);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(levu);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(levv);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(levw);
	
	Matrix<double> derivwt = Math::transpose(derivw);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),numw-ordw,0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<numw-ordw; k++) {
					for (int l=0; l<derivv.GetNumRows(); l++)
						for (int m=0; m<derivw.GetNumRows(); m++)
							for (int n=0; n<derivv.GetNumRows(); n++)
								for (int p=0; p<derivw.GetNumRows(); p++) 
									for (int q=0; q<derivv.GetNumCols(); q++) 
										sum = sum + mat[m][i][l][p][j][n]*derivv[l][q]*derivv[n][q]*derivw[m][k]*derivw[p][k];

					final[k][i][j] = sum;
					sum=0.0;
				}


	return Math::mult5((Math::mult6(Math::transpose(derivu),final)),derivu);
}

Matrix3D<double> BspVolBasisFuncSet::CreateMatrixMinimisation(int levu, int levv, int levw) const
{
	Matrix3D<Matrix3D<double> > mat = CreateMatrixIntegral(levu,levv,levw);

	Matrix<double> derivu = KnotSet(*ktsu,ordu,numu).CreateMatrixDeriv(levu);
	Matrix<double> derivv = KnotSet(*ktsv,ordv,numv).CreateMatrixDeriv(levv);
	Matrix<double> derivw = KnotSet(*ktsw,ordw,numw).CreateMatrixDeriv(levw);

	Matrix<double> derivwt = Math::transpose(derivw);

	Matrix<double> derivwm = Math::mult1(derivw,derivwt);

	Matrix3D<double> final(derivu.GetNumRows(),derivu.GetNumRows(),derivw.GetNumRows(),0.0);

	double sum=0.0;

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++) 
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int r=0; r<derivw.GetNumRows(); r++) {
					for (int l=0; l<derivv.GetNumRows(); l++)
						for (int m=0; m<derivv.GetNumRows(); m++)
							for (int n=0; n<derivv.GetNumCols(); n++)
								sum = sum + mat[k][i][l][r][j][m]*derivv[l][n]*derivv[m][n];

					final[k][i][j] = sum;
					sum=0.0;
				}

	for (int i=0; i<derivu.GetNumRows(); i++)
		for (int j=0; j<derivu.GetNumRows(); j++)
			for (int k=0; k<derivw.GetNumRows(); k++) 
				for (int l=0; l<derivw.GetNumRows(); l++) {
					for (int r=0; r<derivu.GetNumRows(); r++)
						for (int s=0; s<derivw.GetNumRows(); s++)
							sum = sum + mat[l][i][r][s][r][j]*derivwm[l][k];

					final[k][i][j] = sum;
					sum=0.0;
				}

	return Math::mult5((Math::mult6(Math::transpose(derivu),final)),derivu);
}

Matrix3D<double> BspVolBasisFuncSet::CreateIntegralNewU(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<double> res(numu-ordu,numv-ordv,numw-ordw,0.0);

	int count=0;
	for (int l=0; l<numw-ordw; l++) 
		for (int j=0; j<numu-ordu; j++) 
			for (int k=0; k<numv-ordv; k++) 
				for (int i=1; i<=s.GetOrdU(); i++) {
					
					BspVol<double> vol = (*b)[l][j][k].GetBspVol().IntegrateU(i);
				
					double y1=vol.GetLeftLimitU();
					double y2=vol.GetRightLimitU();
					
					BspVol<double> s1 = s.DeriveU(i-1).Subdivide(vol.GetLeftLimitU(),vol.GetRightLimitU(),vol.GetLeftLimitV(),vol.GetRightLimitV(),vol.GetLeftLimitW(),vol.GetRightLimitW());
			
					BspSurf<double> c2 = vol.GetIsoparametricU(y2);
					BspSurf<double> c3 = s1.GetIsoparametricU(y2);
					BspSurf<double> c4 = vol.GetIsoparametricU(y1);
					BspSurf<double> c5 = s1.GetIsoparametricU(y1);
			
			
					BspSurf<double> prod1 = c2.Product(c3);

					BspSurf<double> prod2 = c4.Product(c5);

					BspSurf<double> sub = prod1.Subtract(prod2);
						
					res[l][j][k] = res[l][j][k] + pow(-1.0,double(i-1))*sub.Integrate(v1,v2,w1,w2);
	
					count++;
				}
	return res;
}


Matrix3D<double> BspVolBasisFuncSet::CreateIntegralNewV(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<double> res(numu-ordu,numv-ordv,numw-ordw,0.0);

	for (int j=0; j<numu-ordu; j++) 
		for (int k=0; k<numv-ordv; k++) 
			for (int l=0; l<numw-ordw; l++) {
				for (int i=1; i<=s.GetOrdV(); i++) {	
					BspVol<double> vol = (*b)[l][j][k].GetBspVol().IntegrateV(i);
					double y1=vol.GetLeftLimitV();
					double y2=vol.GetRightLimitV();
					BspVol<double> s1 = s.DeriveV(i-1).Subdivide(vol.GetLeftLimitU(),vol.GetRightLimitU(),vol.GetLeftLimitV(),vol.GetRightLimitV(),vol.GetLeftLimitW(),vol.GetRightLimitW());
					BspSurf<double> c2 = vol.GetIsoparametricV(y2);
					BspSurf<double> c3 = s1.GetIsoparametricV(y2);
					BspSurf<double> c4 = vol.GetIsoparametricV(y1);
					BspSurf<double> c5 = s1.GetIsoparametricV(y1);
					BspSurf<double> prod1 = c2.Product(c3);
					BspSurf<double> prod2 = c4.Product(c5);
					BspSurf<double> sub = prod1.Subtract(prod2);
					res[l][j][k] = res[l][j][k] + pow(-1.0,double(i-1))*sub.Integrate(u1,u2,w1,w2);
					
		}
	}


	return res;
}


Matrix3D<double> BspVolBasisFuncSet::CreateIntegralNewW(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const
{
	Matrix3D<double> res(numu-ordu,numv-ordv,numw-ordw,0.0);

	for (int j=0; j<numu-ordu; j++) 
		for (int k=0; k<numv-ordv; k++) 
			for (int l=0; l<numw-ordw; l++) {
				for (int i=1; i<=s.GetOrdW(); i++) {
					BspVol<double> vol = (*b)[l][j][k].GetBspVol().IntegrateW(i);
					double y1=vol.GetLeftLimitW();
					double y2=vol.GetRightLimitW();
					BspVol<double> s1 = s.DeriveW(i-1).Subdivide(vol.GetLeftLimitU(),vol.GetRightLimitU(),vol.GetLeftLimitV(),vol.GetRightLimitV(),vol.GetLeftLimitW(),vol.GetRightLimitW());
					BspSurf<double> c2 = vol.GetIsoparametricW(y2);
					BspSurf<double> c3 = s1.GetIsoparametricW(y2);
					BspSurf<double> c4 = vol.GetIsoparametricW(y1);
					BspSurf<double> c5 = s1.GetIsoparametricW(y1);
					BspSurf<double> prod1 = c2.Product(c3);
					BspSurf<double> prod2 = c4.Product(c5);
					BspSurf<double> sub = prod1.Subtract(prod2);
					res[l][j][k] = res[l][j][k] + pow(-1.0,double(i-1))*sub.Integrate(u1,u2,v1,v2);
		}
	}


	return res;
}


Matrix3D<Matrix3D<double> > BspVolBasisFuncSet::CreateMatrixIntegral(int levu, int levv, int levw, double U1, double U2, double V1, double V2, double W1, double W2) const
{
	KnotSet ksetu = KnotSet(*ktsu,ordu,numu).CreateKnotSetDeriv(levu);
	KnotSet ksetv = KnotSet(*ktsv,ordv,numv).CreateKnotSetDeriv(levv);
	KnotSet ksetw = KnotSet(*ktsw,ordw,numw).CreateKnotSetDeriv(levw);

	Matrix3D<Matrix3D<double> > mat(ksetu.GetNum()-(ordu-levu),ksetv.GetNum()-(ordv-levv),ksetw.GetNum()-(ordw-levw),Matrix3D<double>(ksetu.GetNum()-ordu+levu,ksetv.GetNum()-ordv+levv,ksetw.GetNum()-ordw+levw,0.0));
	
	BspVolBasisFuncSet basis(ordu-levu,ordv-levv,ordw-levw,ksetu.GetNum(),ksetv.GetNum(),ksetw.GetNum(),ksetu.GetKnots(),ksetv.GetKnots(),ksetw.GetKnots());

	// find the first basis function for which the last distinct
	// knot is greater than x1

	int indu1=-1;
	do {
		indu1++;
	} while ((*(basis.b))[0][indu1][0].GetKnotsU()[ordu-levu] <= U1);


	// find the last basis function for which the first distinct
	// knot is less than x2

	int indu2=-1;
	do {
		indu2++;
	} while (indu2 < ksetu.GetNum()-ordu+levu && (*(basis.b))[0][indu2][0].GetKnotsU()[0] < U2);



	// find the first basis function for which the last distinct
	// knot is greater than x1

	int indv1=-1;
	do {
		indv1++;
	} while ((*(basis.b))[0][0][indv1].GetKnotsV()[ordv-levv] <= V1);


	// find the last basis function for which the first distinct
	// knot is less than x2

	int indv2=-1;
	do {
		indv2++;
	} while (indv2 < ksetv.GetNum()-ordv+levv && (*(basis.b))[0][0][indv2].GetKnotsV()[0] < V2);


		// find the first basis function for which the last distinct
	// knot is greater than x1

	int indw1=-1;
	do {
		indw1++;
	} while ((*(basis.b))[indw1][0][0].GetKnotsW()[ordw-levw] <= W1);


	// find the last basis function for which the first distinct
	// knot is less than x2

	int indw2=-1;
	do {
		indw2++;
	} while (indw2 < ksetw.GetNum()-ordw+levw && (*(basis.b))[indw2][0][0].GetKnotsW()[0] < W2);


	for (int i=indu1; i<=indu2-1; i++)
		for (int j=indv1; j<=indv2-1; j++) 
			for (int k=indw1; k<=indw2-1; k++) 
				//	mat[i][j] = Matrix<double>(ksetu.GetNum()-ordu+levu,ksetv.GetNum()-ordv+levv,0.0);
				for (int l=indu1; l<=indu2-1; l++) 
					for (int m=indv1; m<=indv2-1; m++) 
						for (int n=indw1;n<=indw2-1; n++)  {
					// create the two sets representing the two knot sets
		
							Vector<double> temp1((*(basis.b))[k][i][j].GetKnotsU());
							Vector<double> temp2((*(basis.b))[n][l][m].GetKnotsU());
							Vector<double> temp3((*(basis.b))[k][i][j].GetKnotsV());
							Vector<double> temp4((*(basis.b))[n][l][m].GetKnotsV());
							Vector<double> temp5((*(basis.b))[k][i][j].GetKnotsW());
							Vector<double> temp6((*(basis.b))[n][l][m].GetKnotsW());

							std::set<double> s1(temp1.begin(),temp1.end());
							std::set<double> s2(temp2.begin(),temp2.end());
							std::set<double> s3(temp3.begin(),temp3.end());
							std::set<double> s4(temp4.begin(),temp4.end());
							std::set<double> s5(temp5.begin(),temp5.end());
							std::set<double> s6(temp6.begin(),temp6.end());

			
			
							// if there is an intersection
							//if (*(--s2.end()) > *(s1.begin()) && *(--s4.end()) > *(s3.begin())) {
							// form the intersection
							std::set<double> si1, si2, si3;
							std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(si1,si1.begin()));
							std::set_intersection(s3.begin(),s3.end(),s4.begin(),s4.end(),std::inserter(si2,si2.begin()));
							std::set_intersection(s5.begin(),s5.end(),s6.begin(),s6.end(),std::inserter(si3,si3.begin()));			

							if (si1.size() > 1 && si2.size() > 1 && si3.size() > 1) {
								Vector<double> v1(si1.size());
								Vector<double> v2(si2.size());
								Vector<double> v3(si3.size());

								std::set<double>::iterator si = si1.begin();
								std::set<double>::iterator sj = si2.begin();
								std::set<double>::iterator sk = si3.begin();

								// copy the elements into a vector
								for (unsigned int k1=0; k1<si1.size(); k1++) v1[k1] = *si++;
								for (unsigned int k1=0; k1<si2.size(); k1++) v2[k1] = *sj++;
								for (unsigned int k1=0; k1<si3.size(); k1++) v3[k1] = *sk++;

								// create the compbezcurvs
								Matrix3D<BezVol<double> > mat1(si1.size()-1,si2.size()-1,si3.size()-1);
								Matrix3D<BezVol<double> > mat2(si1.size()-1,si2.size()-1,si3.size()-1);

								BspVol<double> b1((*(basis.b))[k][i][j].GetBspVol()), b2((*(basis.b))[n][l][m].GetBspVol());
								// find the segments of intersection
								for (unsigned int k1=0; k1<si1.size()-1; k1++) 
									for (unsigned int l1=0; l1<si2.size()-1; l1++) 
										for (unsigned int m1=0; m1<si3.size()-1; m1++) {
											int segb1 = b1.GetKnotSetU().Find_segment(v1[k1]);
											int segb2 = b1.GetKnotSetV().Find_segment(v2[l1]);
											int segb3 = b1.GetKnotSetW().Find_segment(v3[m1]);
											int segb4 = b2.GetKnotSetU().Find_segment(v1[k1]);
											int segb5 = b2.GetKnotSetV().Find_segment(v2[l1]);
											int segb6 = b2.GetKnotSetW().Find_segment(v3[m1]);
											mat1[m1][k1][l1] = b1.GetVol(segb1,segb2,segb3);
											mat2[m1][k1][l1] = b2.GetVol(segb4,segb5,segb6);
										}

								CompBezVol<double> cb1(mat1,si1.size()-1,si2.size()-1,si3.size()-1,v1,v2,v3);
								CompBezVol<double> cb2(mat2,si1.size()-1,si2.size()-1,si3.size()-1,v1,v2,v3);
								CompBezVol<double> prod = cb1.Product(cb2);
								mat[k][i][j][n][l][m] = prod.ConvertBspVol().Integrate(U1,U2,V1,V2,W1,W2);
							}
						}


	return mat;
}


Matrix3D<Matrix3D<double> > BspVolBasisFuncSet::CreateMatrixIntegral(int levu, int levv, int levw) const
{
	KnotSet ksetu = KnotSet(*ktsu,ordu,numu).CreateKnotSetDeriv(levu);
	KnotSet ksetv = KnotSet(*ktsv,ordv,numv).CreateKnotSetDeriv(levv);
	KnotSet ksetw = KnotSet(*ktsw,ordw,numw).CreateKnotSetDeriv(levw);

	Matrix3D<Matrix3D<double> > mat(ksetu.GetNum()-(ordu-levu),ksetv.GetNum()-(ordv-levv),ksetw.GetNum()-(ordw-levw),Matrix3D<double>(ksetu.GetNum()-ordu+levu,ksetv.GetNum()-ordv+levv,ksetw.GetNum()-ordw+levw,0.0));
	
	BspVolBasisFuncSet basis(ordu-levu,ordv-levv,ordw-levw,ksetu.GetNum(),ksetv.GetNum(),ksetw.GetNum(),ksetu.GetKnots(),ksetv.GetKnots(),ksetw.GetKnots());



	for (int i=0; i<ksetu.GetNum()-(ordu-levu); i++)
		for (int j=0; j<ksetv.GetNum()-(ordv-levv); j++) 
			for (int k=0; k<ksetw.GetNum()-(ordw-levw); k++) 
				for (int l=0; l<ksetu.GetNum()-(ordu-levu); l++) 
					for (int m=0; m<ksetv.GetNum()-(ordv-levv); m++) 
						for (int n=0;n<ksetw.GetNum()-(ordw-levw); n++) {
					// create the two sets representing the two knot sets
		
							Vector<double> temp1((*(basis.b))[k][i][j].GetKnotsU());
							Vector<double> temp2((*(basis.b))[n][l][m].GetKnotsU());
							Vector<double> temp3((*(basis.b))[k][i][j].GetKnotsV());
							Vector<double> temp4((*(basis.b))[n][l][m].GetKnotsV());
							Vector<double> temp5((*(basis.b))[k][i][j].GetKnotsW());
							Vector<double> temp6((*(basis.b))[n][l][m].GetKnotsW());

							std::set<double> s1(temp1.begin(),temp1.end());
							std::set<double> s2(temp2.begin(),temp2.end());
							std::set<double> s3(temp3.begin(),temp3.end());
							std::set<double> s4(temp4.begin(),temp4.end());
							std::set<double> s5(temp5.begin(),temp5.end());
							std::set<double> s6(temp6.begin(),temp6.end());

			
			
							// if there is an intersection
							//if (*(--s2.end()) > *(s1.begin()) && *(--s4.end()) > *(s3.begin())) {
							// form the intersection
							std::set<double> si1, si2, si3;
							std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(si1,si1.begin()));
							std::set_intersection(s3.begin(),s3.end(),s4.begin(),s4.end(),std::inserter(si2,si2.begin()));
							std::set_intersection(s5.begin(),s5.end(),s6.begin(),s6.end(),std::inserter(si3,si3.begin()));			

							if (si1.size() > 1 && si2.size() > 1 && si3.size() > 1) {
								Vector<double> v1(si1.size());
								Vector<double> v2(si2.size());
								Vector<double> v3(si3.size());

								std::set<double>::iterator si = si1.begin();
								std::set<double>::iterator sj = si2.begin();
								std::set<double>::iterator sk = si3.begin();

								// copy the elements into a vector
								for (unsigned int k1=0; k1<si1.size(); k1++) v1[k1] = *si++;
								for (unsigned int k1=0; k1<si2.size(); k1++) v2[k1] = *sj++;
								for (unsigned int k1=0; k1<si3.size(); k1++) v3[k1] = *sk++;

								// create the compbezcurvs
								Matrix3D<BezVol<double> > mat1(si1.size()-1,si2.size()-1,si3.size()-1);
								Matrix3D<BezVol<double> > mat2(si1.size()-1,si2.size()-1,si3.size()-1);

								BspVol<double> b1((*(basis.b))[k][i][j].GetBspVol()), b2((*(basis.b))[n][l][m].GetBspVol());
								// find the segments of intersection
								for (unsigned int k1=0; k1<si1.size()-1; k1++) 
									for (unsigned int l1=0; l1<si2.size()-1; l1++) 
										for (unsigned int m1=0; m1<si3.size()-1; m1++) {
											int segb1 = b1.GetKnotSetU().Find_segment(v1[k1]);
											int segb2 = b1.GetKnotSetV().Find_segment(v2[l1]);
											int segb3 = b1.GetKnotSetW().Find_segment(v3[m1]);
											int segb4 = b2.GetKnotSetU().Find_segment(v1[k1]);
											int segb5 = b2.GetKnotSetV().Find_segment(v2[l1]);
											int segb6 = b2.GetKnotSetW().Find_segment(v3[m1]);
											mat1[m1][k1][l1] = b1.GetVol(segb1,segb2,segb3);
											mat2[m1][k1][l1] = b2.GetVol(segb4,segb5,segb6);
										}

								CompBezVol<double> cb1(mat1,si1.size()-1,si2.size()-1,si3.size()-1,v1,v2,v3);
								CompBezVol<double> cb2(mat2,si1.size()-1,si2.size()-1,si3.size()-1,v1,v2,v3);
								CompBezVol<double> prod = cb1.Product(cb2);
								mat[k][i][j][n][l][m] = prod.ConvertBspVol().Integrate((*ktsu)[ordu-1],(*ktsu)[numu-ordu+1],(*ktsv)[ordv-1],(*ktsv)[numv-ordv+1],(*ktsw)[ordw-1],(*ktsw)[numw-ordw+1]);
							}
						}


	return mat;
}



Matrix<double> BspVolBasisFuncSet::ComputeUBasisMatrix() const
{

	
	Matrix<double> mat(numu-ordu,numu-ordu);
	Vector<double> u(KnotSet((*ktsu),ordu,numu).ComputeKnotSetAver());

	BspCurvBasisFuncSet b((*ktsu),ordu,numu);

	for (int i=0; i<numu-ordu; i++)
			mat[i] = b.ComputeVecBasis(u[i]);

	return mat;
}

Matrix<double> BspVolBasisFuncSet::ComputeVBasisMatrix() const
{
	Matrix<double> mat(numv-ordv,numv-ordv);
	Vector<double> v(KnotSet((*ktsv),ordv,numv).ComputeKnotSetAver());

	BspCurvBasisFuncSet b((*ktsv),ordv,numv);

	for (int i=0; i<numv-ordv; i++)
			mat[i] = b.ComputeVecBasis(v[i]);

	return mat;
}

Matrix<double> BspVolBasisFuncSet::ComputeWBasisMatrix() const
{
	Matrix<double> mat(numw-ordw,numw-ordw);
	Vector<double> v(KnotSet((*ktsw),ordw,numw).ComputeKnotSetAver());

	BspCurvBasisFuncSet b((*ktsw),ordw,numw);

	for (int i=0; i<numw-ordw; i++)
			mat[i] = b.ComputeVecBasis(v[i]);

	return mat;
}



// READ AND WRITE

void BspVolBasisFuncSet::write(std::ostream& os) 
{
	os << "Bspline Volace basis function set\n";
	os << "order in u is " << ordu << "\n";
	os << "order in v is " << ordv << "\n";
	os << "knots in u are\n";
	os << *ktsu;
	os << "\nknots in v are\n";
	os << *ktsv;
	os << "\nknots in w are\n";
	os << *ktsw;
	os << "\nBspVol representation\n";
	for (int i=0; i<numu-ordu; i++)
		for (int j=0; j<numv-ordv; j++) 
			for (int k=0; k<numw-ordw; k++) os << (*b)[k][i][j];
}


void BspVolBasisFuncSet::read(std::istream& is)
{
	int Ordu, Ordv, Ordw, Numu, Numv, Numw;
	std::cout << "order of Bspline Volace  basis function in u and v, w";
	is >> Ordu >> Ordv >> Ordw;
	std::cout << "number of knots in u and v, w\n";
	is >> Numu >> Numv >> Numw;
	Vector<double> Ktsu(Numu);
	Vector<double> Ktsv(Numv);
	Vector<double> Ktsw(Numw);
	is >> Ktsu;
	is >> Ktsv;
	is >> Ktsw;
	*this = BspVolBasisFuncSet(Ordu,Ordv,Ordw,Numu,Numv,Numw,Ktsu,Ktsv,Ktsw);
} 



void BspVolBasisFuncSet::writefile(std::ofstream& ofs)
{
	ofs << "Bspline Volace basis function\n";
	ofs << "order in u is " << ordu << "\n";
	ofs << "order in v is " << ordv << "\n";
	ofs << "order in w is " << ordw << "\n";
	ofs << "knots in u are\n";
	ofs << *ktsu;
	ofs << "\nknots in v are\n";
	ofs << *ktsv;
	ofs << "\nknots in w are\n";
	ofs << *ktsw;
	
	ofs << "\nBspVol representation is\n";
	for (int i=0; i<numu-ordu; i++)
		for (int j=0; j<numv-ordv; j++) 
			for (int k=0; k<numw-ordw; k++) ofs << (*b)[k][i][j];
}


void BspVolBasisFuncSet::readfile(std::ifstream& ifs)
{
	int Ordu, Ordv, Ordw, Numu, Numv, Numw;
	ifs >> Ordu >> Ordv >> Ordw;
	ifs >> Numu >> Numv >> Numw;
	Vector<double> Ktsu(Numu), Ktsv(Numv), Ktsw(Numw);
	ifs >> Ktsu; 
	ifs >> Ktsv;
	ifs >> Ktsw;
	*this = BspVolBasisFuncSet(Ordu,Ordv,Ordw,Numu,Numv,Numw,Ktsu,Ktsv,Ktsw);
} 






