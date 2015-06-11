
#include "pbspvol.h"


BspVolDouble::BspVolDouble() : BspVol<double>() {}
BspVolDouble::BspVolDouble(const BspVol<double>& v) : BspVol<double>(v) {}
BspVolDouble::BspVolDouble(const Matrix3D<double>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw) :
BspVol<double>(Cpts, Ktsu, Ktsv, Ktsw, Ordu, Ordv, Ordw,Numu,Numv,Numw) {}
BspVolDouble::BspVolDouble(const Matrix3D<double>& Cpts, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw) :
BspVol<double>(Cpts,Ordu,Ordv,Ordw,Numu,Numv,Numw) {}





FBspVol::FBspVol() : BspVol<Point1D>() {}
FBspVol::FBspVol(const BspVol<Point1D>& v) : BspVol<Point1D>(v) {}
FBspVol::FBspVol(const Matrix3D<Point1D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw) :
BspVol<Point1D>(Cpts, Ktsu, Ktsv, Ktsw, Ordu, Ordv, Ordw,Numu,Numv,Numw) {}
FBspVol::FBspVol(const Matrix3D<Point1D>& Cpts, int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw) :
BspVol<Point1D>(Cpts,Ordu,Ordv,Ordw,Numu,Numv,Numw) {}



double FBspVol::Poisson3DError(int numx, int numy, int numz)
{
	double sum=0.0;
	double stepx = 1.0/numx;
	double stepy = 1.0/numy;
	double stepz = 1.0/numz;
	double incx = 0;
	double incy = 0;
	double incz = 0;
	double diff = 0;
	double max= 0.0;


	for (int k=0; k<=numz; k++) {
		for (int i=0; i<=numx; i++) {
			for (int j=0; j<=numy; j++) {
				double p = Math::Poisson3D(incx,incy,incz,10);
				double q = (*this)(incx,incy,incz);
				diff =  fabs(q - p);
				if (fabs(p) > 0.00001)
						diff = diff/p;
				if (diff > max) max = diff;
				incy += stepy;
			}
			incy=0.0;
			incx += stepx;
		}
		incx=0.0;
		incy=0.0;
		incz += stepz;
	}
	return max;
}

double FBspVol::Laplace3D1Error(int numx, int numy, int numz)
{
	double sum=0.0;
	double stepx = 1.0/numx;
	double stepy = 1.0/numy;
	double stepz = 1.0/numz;
	double incx = 0;
	double incy = 0;
	double incz = 0;
	double diff = 0;
	double max= 0.0;


	for (int k=0; k<=numz; k++) {
		for (int i=0; i<=numx; i++) {
			for (int j=0; j<=numy; j++) {
				double p = Math::Laplace3D1(incx,incy,incz,50);
				double q = (*this)(incx,incy,incz);
				diff =  fabs(q - p);
				if (fabs(p) > 0.00001)
						diff = diff/p;
				if (diff > max) max = diff;
				incy += stepy;
			}
			incy=0.0;
			incx += stepx;
		}
		incx=0.0;
		incy=0.0;
		incz += stepz;
	}
	return max;
}


double FBspVol::Laplace3D2Error(int numx, int numy, int numz)
{
	double sum=0.0;
	double stepx = 1.0/numx;
	double stepy = 1.0/numy;
	double stepz = 1.0/numz;
	double incx = 0;
	double incy = 0;
	double incz = 0;
	double diff = 0;
	double max= 0.0;


	for (int k=0; k<=numz; k++) {
		for (int i=0; i<=numx; i++) {
			for (int j=0; j<=numy; j++) {
				double p = Math::Laplace3D2(incx,incy,incz,10);
				double q = (*this)(incx,incy,incz);
				diff =  fabs(q - p);
				if (fabs(p) > 0.00001)
						diff = diff/p;
				if (diff > max) max = diff;
				incy += stepy;
			}
			incy=0.0;
			incx += stepx;
		}
		incx=0.0;
		incy=0.0;
		incz += stepz;
	}
	return max;
}


FBspVol FBspVol::ComputeSolution(double L, double W, double H, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, Vector<Matrix3D<double> >& volume, Vector<int>& natgeom,
										Matrix3D<double>& bound) const
{
	BspVolBasisFuncSet b(GetOrdU(),GetOrdV(),GetOrdW(),GetOrdU()+GetNumU(),GetOrdV()+GetNumV(),GetOrdW()+GetNumW(),GetKnotsU(),GetKnotsV(),GetKnotsW());
	BspSurfBasisFuncSet b10(GetOrdU(),GetOrdV(),GetOrdU()+GetNumU(),GetOrdV()+GetNumV(),GetKnotsU(),GetKnotsV());
	BspSurfBasisFuncSet b11(GetOrdU(),GetOrdW(),GetOrdU()+GetNumU(),GetOrdW()+GetNumW(),GetKnotsU(),GetKnotsW());
	BspSurfBasisFuncSet b12(GetOrdV(),GetOrdW(),GetOrdV()+GetNumV(),GetOrdW()+GetNumW(),GetKnotsV(),GetKnotsW());
	
	BspCurvBasisFuncSet b1(GetKnotSetU());
	BspCurvBasisFuncSet b2(GetKnotSetV());
	BspCurvBasisFuncSet b3(GetKnotSetW());
	BspCurvBasisFuncSet b4(GetKnotSetU().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b5(GetKnotSetV().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b6(GetKnotSetW().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b7(GetKnotSetU().CreateKnotSetDeriv(2));
	BspCurvBasisFuncSet b8(GetKnotSetV().CreateKnotSetDeriv(2));
	BspCurvBasisFuncSet b9(GetKnotSetW().CreateKnotSetDeriv(2));
	Matrix<double> mat1 = b1.CreateMatrixMinimisation(1);
	Matrix<double> mat2 = b2.CreateMatrixMinimisation(1);
	Matrix<double> mat3 = b3.CreateMatrixMinimisation(1);
	Matrix<double> mat5 = b1.CreateMatrixMinimisation(0);
	Matrix<double> mat6 = b2.CreateMatrixMinimisation(0);
	Matrix<double> mat7 = b3.CreateMatrixMinimisation(0);

	Matrix<double> vm1, vm2, vm3;
	vm1 = Math::kronecker(mat6,mat1);
	vm2 = Math::kronecker(mat2,mat5);
	vm3 = Math::kronecker(mat6,mat5);
	vm1 = Math::kronecker(mat7,vm1);
	vm2 = Math::kronecker(mat7,vm2);
	vm3 = Math::kronecker(mat3,vm3);
	Matrix<double> matfinal = Math::add(Math::add(vm1,vm2),vm3);

	Matrix3D<double> mat(GetNumU(),GetNumV(),GetNumW(),0.0);
	int j = 0;
	for (int i=0; i<nums[0]; i++) {
		mat = Math::add(mat,Math::mmult(ptforce[j+3],b(ptforce[j],ptforce[j+1],ptforce[j+2])));
		j+=4;
	}

	j = 0;
	int ind = 0;
	
	int laplace = nums[39];
	
	if (!laplace)
	for (int i=0; i<nums[1]; i++) {
	
		mat = Math::mmult(-1.0,Math::add(mat,b.CreateIntegral(disforce[j],disforce[j+1],disforce[j+2],disforce[j+3],disforce[j+4],disforce[j+5])));
		j+=6;
		ind=ind+3;
	}


	// boundary conditions
	// supports
	
	j=0;
	int row_count=0;

	int num_constraint = nums[2];

	// Face 1 (uv at w=0)
	if (nums[3] != 0)
		num_constraint += GetNumU()*GetNumV();

	if (nums[5] != 0)
		if (natgeom[0] == 0)
			num_constraint += (GetNumU()-2)*(GetNumV()-2);

	// Face 2 (uw at v=1)
	if (nums[9] != 0)
		num_constraint += GetNumU()*(GetNumW()-1);

	if (nums[11] != 0)
		if (natgeom[1] == 0)
			num_constraint += (GetNumU()-2)*(GetNumW()-3);

	// Face 3 (uv at w=1)
	if (nums[15] != 0) 
	
		num_constraint += (GetNumU()-2)*(GetNumV()-2);

	if (nums[17] != 0) 
		if (natgeom[2] == 0)
			num_constraint += (GetNumU()-3)*(GetNumV()-4);

	// Face 4 (uw at v=0)
	if (nums[21] != 0)
	
		num_constraint += GetNumU()*(GetNumW()-1);

	if (nums[23] != 0) 
		if (natgeom[3] == 0) 
			num_constraint += (GetNumU()-2)*(GetNumW()-4);

	// Face 5 (vw at u=1)

	if (nums[27] != 0) 
		num_constraint += (GetNumV()-2)*(GetNumW()-1);

	if (nums[29] != 0) 
		if (natgeom[4] == 0)
			num_constraint += (GetNumV()-4)*(GetNumW()-4);


	// Face 6 (vw at u=0)
		if (nums[33] != 0)
			num_constraint += (GetNumV() - 2)*(GetNumW() - 1);
	
	if (nums[35] != 0)
		if (natgeom[5] == 0)
			num_constraint += (GetNumV()-4)*(GetNumW()-4);
	
	Matrix<double> constraints(num_constraint,GetNumU()*GetNumV()*GetNumW()+1,0.0);
	Matrix<double> mu = Math::ComputeIdentityMatrix(GetNumU());
	Matrix<double> mv = Math::ComputeIdentityMatrix(GetNumV());
	Matrix<double> mw = Math::ComputeIdentityMatrix(GetNumW());

	if (nums[3] != 0) {
	// boundary 1
		BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(0),nums[3],nums[4],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());
		for (int i=0; i<GetNumU(); i++) 
			for (int j=0; j<GetNumV(); j++) {
				Vector<double> v11 = b3(0.0);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mv.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v11),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m1(GetNumU(),GetNumV(),0.0);
	if (nums[5] != 0) {
		if (natgeom[0] == 0) {
			BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(1),nums[5],nums[6],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=1; i<GetNumU()-1; i++) 
				for (int j=1; j<GetNumV()-1; j++) {
					Vector<double> v11 = Math::mult4(b6(0.0),GetKnotSetW().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mv.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v11),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b13= PolySurf<double>(bound.GetUV(1),nums[5],nums[6],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).ConvertBspSurf();
			m1 = Math::add(m1,b10.CreateIntegralNewU(b13,b13,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]));
		}
	}
	// boundary 2

	if (nums[9] != 0) {
	// boundary 2
		BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(3),nums[9],nums[10],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());

		for (int i=0; i<GetNumU(); i++) 
			for (int j=1; j<GetNumW(); j++) {
				Vector<double> v11 = b2(GetKnotsV()[GetNumV()]);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m2(GetNumU(),GetNumW(),0.0);
	if (nums[11] != 0) {
		if (natgeom[1] == 0) {
			BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(4),nums[11],nums[12],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
		
			for (int i=1; i<GetNumU()-1; i++) 
				for (int j=2; j<GetNumW()-1; j++) {
					Vector<double> v11 = Math::mult4(b5(GetKnotsV()[GetNumV()]),GetKnotSetV().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b13= PolySurf<double>(bound.GetUV(4),nums[11],nums[12],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]).ConvertBspSurf();
			m2 = Math::add(m2,b11.CreateIntegralNewU(b13,b13,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
		}
	}


	// boundary 3
	if (nums[15] != 0) {
		BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(6),nums[15],nums[16],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());
	
		for (int i=1;/*0;*/ i<GetNumU()-1; i++) 
			for (int j= 1; /*0;*/ j<GetNumV()-1; j++) {
				Vector<double> v11 = b3(GetKnotsW()[GetNumW()]);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mv.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v11),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m3(GetNumU(),GetNumV(),0.0);
	if (nums[17] != 0) {
		if (natgeom[2] == 0) {
			BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(7),nums[17],nums[18],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
		
			for (int i=1; i<GetNumU()-1; i++) 
				for (int j=1; j<GetNumV()-2; j++) {
					Vector<double> v11 = Math::mult4(b6(GetKnotsW()[GetNumW()]),GetKnotSetW().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mv.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v11),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b13= PolySurf<double>(bound.GetUV(7),nums[17],nums[18],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).ConvertBspSurf();
			m3 = Math::add(m3,b10.CreateIntegralNewU(b13,b13,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]));
		}
	}

	// boundary 4

	if (nums[21] != 0) {
	// boundary 4
		BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(9),nums[21],nums[22],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());
	
		for (int i=0; i<GetNumU(); i++) 
			for (int j=1; j<GetNumW();  j++) {
				Vector<double> v11 = b2(0.0);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m4(GetNumU(),GetNumW(),0.0);
	if (nums[23] != 0) {
		if (natgeom[3] == 0) {
			BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(10),nums[23],nums[24],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
		
			for (int i=1; i<GetNumU()-1; i++) 
				for (int j=2; j<GetNumW()-2; j++) {
					Vector<double> v11 = Math::mult4(b5(0.0),GetKnotSetV().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b13= PolySurf<double>(bound.GetUV(10),nums[23],nums[24],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]).ConvertBspSurf();
			m4 = Math::add(m4,b11.CreateIntegralNewU(b13,b13,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
		}
	}

	if (nums[27] != 0) {
	// boundary 5
		BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(12),nums[27],nums[28],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());
	
		
		for (int i=1; i<GetNumV()-1; i++) 
			for (int j=1; j<GetNumW(); j++) {
				Vector<double> v11 = b1(GetKnotsU()[GetNumU()]);
				Vector<double> v21 = mv.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m5(GetNumV(),GetNumW(),0.0);
	if (nums[29] != 0) {
		if (natgeom[4] == 0) {
			BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(13),nums[29],nums[30],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
		
			for (int i=2; i<GetNumV()-2; i++) 
				for (int j=2; j<GetNumW()-2; j++) {
					Vector<double> v11 = Math::mult4(b4(GetKnotsU()[GetNumU()]),GetKnotSetU().CreateMatrixDeriv());
					Vector<double> v21 = mv.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b13= PolySurf<double>(bound.GetUV(13),nums[29],nums[30],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]).ConvertBspSurf();
			m5 = Math::add(m5,b12.CreateIntegralNewU(b13,b13,GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
		}
	}

	if (nums[33] != 0) {

		// boundary 6
		BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(15),nums[33],nums[34],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());
	
		for (int i=1; i<GetNumV()-1; i++) 
			for (int j=1; j<GetNumW(); j++) {
				Vector<double> v11 = b1(0.0);
				Vector<double> v21 = mv.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m6(GetNumV(),GetNumW(),0.0);
	if (nums[35] != 0) {
		if (natgeom[5] == 0) {
			BspSurf<double> b11 = BspSurf<double>(PolySurf<double>(bound.GetUV(16),nums[35],nums[36],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
		
			for (int i=2; i<GetNumV()-2; i++) 
				for (int j=2; j<GetNumW()-2; j++) {
					Vector<double> v11 = Math::mult4(b4(0.0),GetKnotSetU().CreateMatrixDeriv());
					Vector<double> v21 = mv.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b11.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b13= PolySurf<double>(bound.GetUV(16),nums[35],nums[36],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]).ConvertBspSurf();
			m6 = Math::add(m6,b12.CreateIntegralNewU(b13,b13,GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
		}
	}


	j=0;
	for (int i=0; i<nums[2]; i++) {
		double u = supports[j];
		double v = supports[j+1];
		double w = supports[j+2];
		double val = supports[j+3];
		Vector<double> v1 =	b1(u);
		Vector<double> v2 = b2(v);
		Vector<double> v3 = b3(w);
		Matrix<double> mat = Math::kronecker(Matrix<double>(v2),Matrix<double>(v1));
		mat = Math::kronecker(Matrix<double>(v3),mat);
		constraints.InsertRow(mat.GetRow(0),row_count);
		constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = val;
		row_count++;
		j=j+4;
	}


	Matrix3D<double> mat4(GetNumU(),GetNumV(),GetNumW(),0.0);

	// add these explicitly!!!!
	mat4.InsertUV(m1,0);
	mat4.InsertUV(m3,GetNumW()-1);
	mat4.InsertUW(m2,GetNumV()-1);
	mat4.InsertUW(m4,0);
	mat4.InsertVW(m5,GetNumU()-1);
	mat4.InsertVW(m6,0);
		

	mat = Math::subtract(mat,mat4);

	// add them in
	Vector<double> vec1(Math::CreateKroneckerVector(mat));
	vec1 = Math::vmult(-1.0,vec1);
	Vector<int> v2 = Math::GetEliminateIndices(matfinal,constraints);	
	Matrix<double> m10 = Math::EliminateMVariables(matfinal,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(matfinal,constraints,vec1,v2);

    Vector<double> sol = Math::gauss(m10,v1,v1.GetNum());


	Vector<Point1D> nsol(sol.GetNum());
	
	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m30(GetNumU()*GetNumV()*GetNumW(),GetNumU()*GetNumV()*GetNumW(),0.0);
	Vector<double> v3(GetNumU()*GetNumV()*GetNumW());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNumU()*GetNumV()*GetNumW());
		
	for (int i=0; i<GetNumU()*GetNumV()*GetNumW(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<num_constraint; i++)
		for (int j=0; j<GetNumU()*GetNumV()*GetNumW(); j++) 
			m30[i][j] = constraints[i][j];
		

	
	for (int i=num_constraint; i<GetNumU()*GetNumV()*GetNumW(); i++) m30[i][diff[i-num_constraint]] = 1.0;


	// rhs Vector

	for (int i=0; i<num_constraint; i++) v3[i]= constraints[i][GetNumU()*GetNumV()*GetNumW()];//0.0;


	for (int i=num_constraint; i<GetNumU()*GetNumV()*GetNumW(); i++) v3[i] = sol[i-num_constraint];


	// solve

	Vector<double> sol1 = Math::gauss(m30,v3,GetNumU()*GetNumV()*GetNumW());

	
	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);


	return FBspVol(Math::Matrix3DFromVector(nsol1,GetNumU(),GetNumV(),GetNumW()),GetKnotsU(),GetKnotsV(),GetKnotsW(),GetOrdU(),GetOrdV(),GetOrdW(),GetNumU(),GetNumV(),GetNumW());
}


FBspVol FBspVol::ComputeSolution1(double L, double W, double H, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, Vector<Matrix3D<double> >& volume, Vector<int>& natgeom,
										Matrix3D<double>& bound) const
{
	BspVolBasisFuncSet b(GetOrdU(),GetOrdV(),GetOrdW(),GetOrdU()+GetNumU(),GetOrdV()+GetNumV(),GetOrdW()+GetNumW(),GetKnotsU(),GetKnotsV(),GetKnotsW());
	BspSurfBasisFuncSet b101(GetOrdU(),GetOrdV(),GetOrdU()+GetNumU(),GetOrdV()+GetNumV(),GetKnotsU(),GetKnotsV());
	BspSurfBasisFuncSet b111(GetOrdU(),GetOrdW(),GetOrdU()+GetNumU(),GetOrdW()+GetNumW(),GetKnotsU(),GetKnotsW());
	BspSurfBasisFuncSet b121(GetOrdV(),GetOrdW(),GetOrdV()+GetNumV(),GetOrdW()+GetNumW(),GetKnotsV(),GetKnotsW());
	
	BspCurvBasisFuncSet b1(GetKnotSetU());
	BspCurvBasisFuncSet b2(GetKnotSetV());
	BspCurvBasisFuncSet b3(GetKnotSetW());
	BspCurvBasisFuncSet b4(GetKnotSetU().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b5(GetKnotSetV().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b6(GetKnotSetW().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b7(GetKnotSetU().CreateKnotSetDeriv(2));
	BspCurvBasisFuncSet b8(GetKnotSetV().CreateKnotSetDeriv(2));
	BspCurvBasisFuncSet b9(GetKnotSetW().CreateKnotSetDeriv(2));

	Matrix<double> mat1 = b1.CreateMatrixMinimisation(1);
	Matrix<double> mat2 = b2.CreateMatrixMinimisation(1);
	Matrix<double> mat3 = b3.CreateMatrixMinimisation(1);
	Matrix<double> mat5 = b1.CreateMatrixMinimisation(0);
	Matrix<double> mat6 = b2.CreateMatrixMinimisation(0);
	Matrix<double> mat7 = b3.CreateMatrixMinimisation(0);

	Matrix<double> vm1, vm2, vm3;
	vm1 = Math::kronecker(mat6,mat1);
	vm2 = Math::kronecker(mat2,mat5);
	vm3 = Math::kronecker(mat6,mat5);
	vm1 = Math::kronecker(mat7,vm1);
	vm2 = Math::kronecker(mat7,vm2);
	vm3 = Math::kronecker(mat3,vm3);
	Matrix<double> matfinal = Math::add(Math::add(vm1,vm2),vm3);

	Matrix3D<double> mat(GetNumU(),GetNumV(),GetNumW(),0.0);
	int j = 0;
	for (int i=0; i<nums[0]; i++) {
		mat = Math::add(mat,Math::mmult(ptforce[j+3],b(ptforce[j],ptforce[j+1],ptforce[j+2])));
		j+=4;
	}

	j = 0;
	int ind = 0;
	
	int laplace = nums[51];
	if (!laplace)
	for (int i=0; i<nums[1]; i++) {
	
		mat = Math::mmult(-1.0,Math::add(mat,b.CreateIntegral(disforce[j],disforce[j+1],disforce[j+2],disforce[j+3],disforce[j+4],disforce[j+5])));
		j+=6;
		ind=ind+3;
	}

	// boundary conditions
	// supports
	
	j=0;
	int row_count=0;

	int num_constraint = nums[2];
	
	

	int b11 , b12 , b13 , b14 , b15 , b16, b17, b18;
	int b21, b22 , b23 , b24, b25, b26, b27, b28;
	int b31, b32 , b33 , b34, b35, b36, b37, b38;
	int b41, b42 , b43 , b44, b45, b46, b47, b48;
	int b51, b52, b53 , b54, b55, b56, b57, b58;
	int b61, b62 , b63 , b64, b65, b66, b67, b68;

	b11 = b12 = b13 = b14 = b15 = b16= b17= b18=0;
	b21= b22 = b23 = b24= b25= b26= b27= b28=0;
	b31= b32 = b33 = b34= b35= b36= b37= b38=0;
	b41= b42 = b43 = b44= b45= b46= b47= b48=0;
	b51= b52= b53 = b54= b55= b56= b57= b58=0;
	b61= b62 = b63 = b64= b65= b66= b67= b68=0;


	// Face 1 (uv at w=0)
	if (nums[3] != 0) {
		num_constraint += GetNumU()*GetNumV();
		b11 = b12 = b13 = b14 = 0;
	}

	if (nums[5] != 0)
		if (natgeom[0] == 0) {
			num_constraint += (GetNumU())*(GetNumV());
			b15 = b16 = b17 = b18 = 0;
		}

	// Face 2 (uw at v=1)
	if (nums[9] != 0 && nums[5] != 0) {
		num_constraint += GetNumU()*(GetNumW()-2);
		b21 = b22 = 0; b23 = 2; b24 = 0;
	}
	else if (nums[9] != 0 && nums[3] != 0) {
		num_constraint += GetNumU()*(GetNumW()-1);
		b21 = b22 = 0; b23 = 1; b24 = 0;
	}
	else if (nums[9] != 0) {
		num_constraint += GetNumU()*GetNumW();
		b21 = b22 = b23 = b24 = 0;
	}

	if (natgeom[1] == 0) 
		if (nums[11] != 0) {
			if (nums[3] != 0 && nums[5] == 0) {
				num_constraint += GetNumU()*(GetNumW()-1);
				b25 = b26 = 0; b27 = 1; b28 = 0;
			}
			else if (nums[3] != 0 && nums[5] != 0) {
				num_constraint += GetNumU()*(GetNumW()-2);
				b25 = 0; b26 = 0; b27 = 2; b28 = 0;
			}
			else {
				num_constraint += GetNumU()*GetNumW();
				b25 = b26 = b27 = b28 = 0;
			}
		}


	// Face 3 (uv at w=1)
	if (nums[15] != 0 && nums[11] !=0) {
		num_constraint += GetNumU()*(GetNumV()-2);
		b31 = b32 = b33 = 0; b34 = 2;
	}
	else if (nums[15] != 0 && nums[9] != 0) {
		num_constraint += GetNumU()*(GetNumV()-1);
		b31 = b32 = b33 = 0; b34 = 1;
	}
	else if (nums[15] != 0) {
		num_constraint += GetNumU()*GetNumV();
		b31 = b32 = b33 = b34 = 0;
	}

	if (natgeom[2] == 0) 
		if (nums[17] != 0) {
			if (nums[9] != 0 && nums[11] == 0) {
				num_constraint += GetNumU()*(GetNumV()-1);
				b35 = b36 = b37 = 0; b38 = 1;
			}
			else if (nums[9] != 0 && nums[11] != 0) {
				num_constraint += GetNumU()*(GetNumV()-2);
				b35 = b26 = b37 =0; b38 = 1;
			}
			else {
				num_constraint += GetNumU()*GetNumV();
				b35 = b36 = b37 = b38 = 0;
			}
		}
		

	// face 4 (uw at v = 0)
	if (nums[21] != 0 && nums[5] != 0 && nums[17] != 0) {
		num_constraint += GetNumU()*(GetNumW()-4);
		b41 = b42 = 0; b43 = 2; b44 = 2;
	}
	else if (nums[21] != 0 && nums[5] != 0 && nums[15] != 0) {
		num_constraint += GetNumU()*(GetNumW()-3);
		b41 = b42 = 0; b43 = 2; b44 = 1;
	}
	else if (nums[21] != 0 && nums[5] != 0 && nums[15] == 0) {
		num_constraint += GetNumU()*(GetNumW()-2);
		b41 = b42 = 0; b43 = 2; b44 = 0;
	}
	else if (nums[21] != 0 && nums[3] != 0 && nums[17] != 0) {
		num_constraint += GetNumU()*(GetNumW()-3);
		b41 = b42 = 0; b43 = 1; b44 = 2;
	}
	else if (nums[21] != 0 && nums[3] != 0 && nums[15] != 0) {
		num_constraint += GetNumU()*(GetNumW()-2);
		b41 = b42 = 0; b43 = b44 = 1;
	}
	else if (nums[21] != 0 && nums[3] != 0 && nums[15] == 0) {
		num_constraint += GetNumU()*(GetNumW()-1);
		b41 = b42 = 0; b43 = 1; b44 = 0;
	}
	else if (nums[21] != 0 && nums[3] == 0 && nums[15] != 0) {
		num_constraint += GetNumU()*(GetNumW()-1);
		b41 = b42 = b43 = 0; b44 = 1;
	}
	else if (nums[21] != 0 && nums[3] == 0 && nums[17] != 0) {
		num_constraint += GetNumU()*(GetNumW()-2);
		b41 = b42 = b43 = 0; b44 = 2;
	}
	else if (nums[21] != 0 && nums[3] == 0 && nums[15] == 0) {
		num_constraint += GetNumU()*GetNumW();
		b41 = b42 = b43 = b44 = 0;
	}

	if (natgeom[3] == 0) 
		if (nums[23] != 0) {
			if (nums[5] != 0 && nums[17] != 0) {
				num_constraint += GetNumU()*(GetNumW()-4);
				b45 = b46 = 0; b47 = b48 = 2;
			}
			else if (nums[5] != 0 && nums[15] != 0) {
				num_constraint += GetNumU()*(GetNumW()-3);
				b45 = b46 = 0; b47 = 2; b48 = 1;
			}
			else if (nums[5] != 0 && nums[15] == 0) {
				num_constraint += GetNumU()*(GetNumW()-2);
				b45 = b46 = 0; b47 = 2; b48 = 0;
			}
			else if (nums[3] != 0 && nums[17] != 0) {
				num_constraint += GetNumU()*(GetNumU()-3);
				b45 = b46 = 0; b47 = 1; b48 = 2;
			}
			else if (nums[3] != 0 && nums[15] != 0) {
				num_constraint += GetNumU()*(GetNumW()-2);
				b45 = b46 = 0; b47 = b48 = 1;
			}
			else if (nums[3] != 0 && nums[15] == 0) {
				num_constraint += GetNumU()*(GetNumW()-1);
				b45 = b46 = 0; b47 = 1; b48 = 0;
			}
			else if (nums[3] == 0 && nums[17] != 0) {
				num_constraint += GetNumU()*(GetNumW()-2);
				b45 = b46 = b47 = 0; b48 = 2;
			}
			else if (nums[3] == 0 && nums[15] != 0) {
				num_constraint += GetNumU()*(GetNumW()-1);
				b45 = b46 = b47 = 0; b48 = 1;
			}
			else if (nums[3] == 0 && nums[15] == 0) {
				num_constraint += GetNumU()*GetNumW();
				b45 = b46 = b47 = b48 = 0;
			}
		}
		
	// Face 5 (vw at u=1)
	if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-4);
		b51 = b52 = b53 = b54 = 2;
	}

	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-4);
		b51 = 1; b52 = b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b51 = 0; b52 = b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-3);
		b51 = b52 = b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b51 = 1; b52 = 2; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 0; b52 = 2; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b51 = b52 = b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 = 1; b52 = 2; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 0; b52 = 2; b53 = 2; b54 = 0;
	}

	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-4);
		b51 = 2; b52 = 1; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b51 = b52 = 1; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b51 = 0 ; b52 = 1; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b51 = 2; b52 = 1; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = b52 = 1; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 0; b52 = 1; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 0; b52 = 1; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 = 2; b52 = 1; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = b52 = 1; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b51 = 2; b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b51 = 1; b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-4);
		b51 = b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 2; b52 = 0; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 1; b52 = 0; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b51 = b52 = 0; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = 0; b53 = 2; b54 = 0;
	}


	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-3);
		b51 = b52 = 2; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b51 = 1; b52 = 2; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 0; b52 = 2; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b51 = b52 = 2; b53 = 1; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 = 1; b52 = 2; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 0; b52 = 2; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-1);
		b51 = b52 = 2; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b51 = 1; b52 = 2; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 0; b52 = 2; b53 = 1; b54 = 0;
	}

	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b51 = b52 = 2; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 =1; b52 = 2; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 0; b52 = 2; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-1);
		b51 = b52 = 2; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b51 = 1; b52 = 2; b53 = 0; b54 = 1;
	}

	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 0; b52 = 2; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*GetNumW();
		b51 = b52 = 2; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*GetNumW();
		b51 = 1; b52 = 2; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b51 = 0; b52 = 2; b53 = b54 = 0;
	}

	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b51 = 2; b52 = 1; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = b52 = 1; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 0; b52 = 1; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 = 2; b52 = 1; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = b52 = b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 0; b52 = b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b51 = 2; b52 = 1; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = b52 = b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 0; b52 = 1; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 2; b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 1; b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b51 = b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 2; b52 = 0; b53 = 1; b54 = 0;
	}
		
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 1; b52 = 0; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b51 = b52 = 0; b53 = 1; b54 = 0;
	}

	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 =2; b52 = 1; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = b52 = 1; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 0; b52 = 1; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b51 = 2; b52 = 1; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = b52 = 1; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 0; b52 = 1; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*GetNumW();
		b51 = 2; b52 = 1; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b51 = b52 = 1; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b51 = 0; b52 = 1; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 2; b52 = b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 1; b52 = 0; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b51 = b52 = b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b51 = 2; b52 = b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b51 = 1; b52 = b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*GetNumW();
		b51 = b52 = b53 = b54 = 0;
	}

	if (natgeom[4] == 0) 
		if (nums[29] != 0) {
			if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-4);
				b55 = b56 = b57 = b58 = 2;
			}

			else if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-4);
				b55 = 1; b56 = b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b55 = 0; b56 = b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-3);
				b55 = b56 = b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b55 = 1; b56 = 2; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 0; b56 = 2; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b55 = b56 = b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 = 1; b56 = 2; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 0; b56 = 2; b57 = 2; b58 = 0;
			}

			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-4);
				b55 = 2; b56 = 1; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b55 = b56 = 1; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b55 = 0 ; b56 = 1; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b55 = 2; b56 = 1; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = b56 = 1; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 0; b56 = 1; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 0; b56 = 1; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 = 2; b56 = 1; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = b56 = 1; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b55 = 2; b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b55 = 1; b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-4);
				b55 = b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 2; b56 = 0; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 1; b56 = 0; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b55 = b56 = 0; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = 0; b57 = 2; b58 = 0;
			}


			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-3);
				b55 = b56 = 2; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b55 = 1; b56 = 2; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 0; b56 = 2; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b55 = b56 = 2; b57 = 1; b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 = 1; b56 = 2; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 0; b56 = 2; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-1);
				b55 = b56 = 2; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b55 = 1; b56 = 2; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 0; b56 = 2; b57 = 1; b58 = 0;
			}

			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b55 = b56 = 2; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 =1; b56 = 2; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 0; b56 = 2; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-1);
				b55 = b56 = 2; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b55 = 1; b56 = 2; b57 = 0; b58 = 1;
			}

			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 0; b56 = 2; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*GetNumW();
				b55 = b56 = 2; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*GetNumW();
				b55 = 1; b56 = 2; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b55 = 0; b56 = 2; b57 = b58 = 0;
			}
			
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b55 = 2; b56 = 1; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = b56 = 1; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 0; b56 = 1; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 = 2; b56 = 1; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = b56 = b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 0; b56 = b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b55 = 2; b56 = 1; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = b56 = b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 0; b56 = 1; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 2; b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 1; b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b55 = b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 2; b56 = 0; b57 = 1; b58 = 0;
			}
				
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 1; b56 = 0; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b55 = b56 = 0; b57 = 1; b58 = 0;
			}

			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 =2; b56 = 1; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = b56 = 1; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 0; b56 = 1; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b55 = 2; b56 = 1; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = b56 = 1; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 0; b56 = 1; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*GetNumW();
				b55 = 2; b56 = 1; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b55 = b56 = 1; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b55 = 0; b56 = 1; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 2; b56 = b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 1; b56 = 0; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b55 = b56 = b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b55 = 2; b56 = b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b55 = 1; b56 = b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*GetNumW();
				b55 = b56 = b57 = b58 = 0;
			}
		}


			

	// Face 6 (vw at u=0)
	if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-4);
		b61 = b62 = b63 = b64 = 2;
	}

	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-4);
		b61 = 1; b62 = b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b61 = 0; b62 = b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-3);
		b61 = b62 = b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b61 = 1; b62 = 2; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 0; b62 = 2; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b61 = b62 = b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 = 1; b62 = 2; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 0; b62 = 2; b63 = 2; b64 = 0;
	}


	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-4);
		b61 = 2; b62 = 1; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b61 = b62 = 1; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b61 = 0 ; b62 = 1; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b61 = 2; b62 = 1; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = b62 = 1; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 0; b62 = 1; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 0; b62 = 1; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 = 2; b62 = 1; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = b62 = 1; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b61 = 2; b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b61 = 1; b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-4);
		b61 = b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 2; b62 = 0; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 1; b62 = 0; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b61 = b62 = 0; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = 0; b63 = 2; b64 = 0;
	}


	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-3);
		b61 = b62 = 2; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b61 = 1; b62 = 2; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 0; b62 = 2; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b61 = b62 = 2; b63 = 1; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 = 1; b62 = 2; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 0; b62 = 2; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-1);
		b61 = b62 = 2; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b61 = 1; b62 = 2; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 0; b62 = 2; b63 = 1; b64 = 0;
	}
	

	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b61 = b62 = 2; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 =1; b62 = 2; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 0; b62 = 2; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-1);
		b61 = b62 = 2; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b61 = 1; b62 = 2; b63 = 0; b64 = 1;
	}

	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 0; b62 = 2; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*GetNumW();
		b61 = b62 = 2; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*GetNumW();
		b61 = 1; b62 = 2; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b61 = 0; b62 = 2; b63 = b64 = 0;
	}
	
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b61 = 2; b62 = 1; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = b62 = 1; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 0; b62 = 1; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 = 2; b62 = 1; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = b62 = b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 0; b62 = b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b61 = 2; b62 = 1; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = b62 = b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 0; b62 = 1; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 2; b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 1; b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b61 = b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 2; b62 = 0; b63 = 1; b64 = 0;
	}
		
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 1; b62 = 0; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b61 = b62 = 0; b63 = 1; b64 = 0;
	}

	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 =2; b62 = 1; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = b62 = 1; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 0; b62 = 1; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b61 = 2; b62 = 1; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = b62 = 1; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 0; b62 = 1; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*GetNumW();
		b61 = 2; b62 = 1; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b61 = b62 = 1; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b61 = 0; b62 = 1; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 2; b62 = b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 1; b62 = 0; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b61 = b62 = b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b61 = 2; b62 = b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b61 = 1; b62 = b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*GetNumW();
		b61 = b62 = b63 = b64 = 0;
	}

	if (natgeom[5] == 0) 
		if (nums[35] != 0) {
			if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-4);
				b65 = b66 = b67 = b68 = 2;
			}

			else if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-4);
				b65 = 1; b66 = b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b65 = 0; b66 = b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-3);
				b65 = b66 = b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b65 = 1; b66 = 2; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 0; b66 = 2; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b65 = b66 = b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 = 1; b66 = 2; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 0; b66 = 2; b67 = 2; b68 = 0;
			}

			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-4);
				b65 = 2; b66 = 1; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b65 = b66 = 1; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b65 = 0 ; b66 = 1; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b65 = 2; b66 = 1; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = b66 = 1; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 0; b66 = 1; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 0; b66 = 1; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 = 2; b66 = 1; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = b66 = 1; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b65 = 2; b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b65 = 1; b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-4);
				b65 = b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 2; b66 = 0; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 1; b66 = 0; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b65 = b66 = 0; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = 0; b67 = 2; b68 = 0;
			}


			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-3);
				b65 = b66 = 2; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b65 = 1; b66 = 2; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 0; b66 = 2; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b65 = b66 = 2; b67 = 1; b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 = 1; b66 = 2; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 0; b66 = 2; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-1);
				b65 = b66 = 2; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b65 = 1; b66 = 2; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 0; b66 = 2; b67 = 1; b68 = 0;
			}
		
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b65 = b66 = 2; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 =1; b66 = 2; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 0; b66 = 2; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-1);
				b65 = b66 = 2; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b65 = 1; b66 = 2; b67 = 0; b68 = 1;
			}

			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 0; b66 = 2; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*GetNumW();
				b65 = b66 = 2; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*GetNumW();
				b65 = 1; b66 = 2; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b65 = 0; b66 = 2; b67 = b68 = 0;
			}

			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b65 = 2; b66 = 1; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = b66 = 1; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 0; b66 = 1; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 = 2; b66 = 1; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = b66 = b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 0; b66 = b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b65 = 2; b66 = 1; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = b66 = b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 0; b66 = 1; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 2; b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 1; b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b65 = b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 2; b66 = 0; b67 = 1; b68 = 0;
			}
				
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 1; b66 = 0; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b65 = b66 = 0; b67 = 1; b68 = 0;
			}

			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 =2; b66 = 1; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = b66 = 1; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 0; b66 = 1; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b65 = 2; b66 = 1; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = b66 = 1; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 0; b66 = 1; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*GetNumW();
				b65 = 2; b66 = 1; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b65 = b66 = 1; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b65 = 0; b66 = 1; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 2; b66 = b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 1; b66 = 0; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b65 = b66 = b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b65 = 2; b66 = b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b65 = 1; b66 = b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*GetNumW();
				b65 = b66 = b67 = b68 = 0;
			}
		}

	
	Matrix<double> constraints(num_constraint,GetNumU()*GetNumV()*GetNumW()+1,0.0);
	Matrix<double> mu = Math::ComputeIdentityMatrix(GetNumU());
	Matrix<double> mv = Math::ComputeIdentityMatrix(GetNumV());
	Matrix<double> mw = Math::ComputeIdentityMatrix(GetNumW());

	if (nums[3] != 0) {
	// boundary 1
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(0),nums[3],nums[4],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());
		
		for (int i=b11; i<GetNumU()-b12; i++) 
			for (int j=b13; j<GetNumV()-b14; j++) {
				Vector<double> v11 = b3(0.0);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mv.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v11),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m1(GetNumU(),GetNumV(),0.0);
	if (nums[5] != 0) {
		if (natgeom[0] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(1),nums[5],nums[6],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
		
			for (int i=b15; i<GetNumU()-b16; i++) 
				for (int j=b17; j<GetNumV()-b18; j++) {
					Vector<double> v11 = Math::mult4(b6(0.0),GetKnotSetW().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mv.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v11),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		} 
	}
		
	if (natgeom[0] != 0) {
			BspSurf<double> b131= BspSurf<double>(PolySurf<double>(bound.GetUV(18),nums[39],nums[40],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());//ConvertBspSurf();
			m1 = Math::add(m1,b101.CreateIntegralNewU(b131,b131,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]));
	}
	// boundary 2
	
	if (nums[9] != 0) {
	// boundary 2
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(3),nums[9],nums[10],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());
	
		for (int i=b21; i<GetNumU()-b22; i++) 
			for (int j=b23; j<GetNumW()-b24; j++) {
				Vector<double> v11 = b2(GetKnotsV()[GetNumV()]);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m2(GetNumU(),GetNumW(),0.0);
	if (nums[11] != 0) {
		if (natgeom[1] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(4),nums[11],nums[12],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b25; i<GetNumU()-b26; i++) 
				for (int j=b27; j<GetNumW()-b28; j++) {
					Vector<double> v11 = Math::mult4(b5(GetKnotsV()[GetNumV()]),GetKnotSetV().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		}
	}
		
		
	if (natgeom[1] != 0) {	
		BspSurf<double> b131= BspSurf<double>(PolySurf<double>(bound.GetUV(19),nums[41],nums[42],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());//ConvertBspSurf();
		m2 = Math::add(m2,b111.CreateIntegralNewU(b131,b131,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
	}


	// boundary 3
	if (nums[15] != 0) {
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(6),nums[15],nums[16],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());
		
		for (int i=b31; i<GetNumU()-b32; i++) 
			for (int j=b33; j<GetNumV()-b34; j++) {
				Vector<double> v11 = b3(GetKnotsW()[GetNumW()]);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mv.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v11),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m3(GetNumU(),GetNumV(),0.0);
	if (nums[17] != 0) {
		if (natgeom[2] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(7),nums[17],nums[18],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b35; i<GetNumU()-b36; i++) 
				for (int j=b37; j<GetNumV()-b38; j++) {
					Vector<double> v11 = Math::mult4(b6(GetKnotsW()[GetNumW()]),GetKnotSetW().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mv.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v11),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		} 
	}
	

	if (natgeom[2] != 0) {
		BspSurf<double> b131 = BspSurf<double>(PolySurf<double>(bound.GetUV(20),nums[43],nums[44],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());
		m3 = Math::add(m3,b101.CreateIntegralNewU(b131,b131,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]));
	}

	// boundary 4

	if (nums[21] != 0) {
	// boundary 4
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(9),nums[21],nums[22],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());
		for (int i=b41; i<GetNumU()-b42; i++) 
			for (int j=b43; j<GetNumW()-b44; j++) {
				Vector<double> v11 = b2(0.0);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m4(GetNumU(),GetNumW(),0.0);
	if (nums[23] != 0) {
		if (natgeom[3] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(10),nums[23],nums[24],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b45; i<GetNumU()-b46; i++) 
				for (int j=b47; j<GetNumW()-b48; j++) {
					Vector<double> v11 = Math::mult4(b5(0.0),GetKnotSetV().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		}
	}
	
	if (natgeom[3] != 0) {
		BspSurf<double> b131= BspSurf<double>(PolySurf<double>(bound.GetUV(21),nums[45],nums[46],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());//.ConvertBspSurf();
		m4 = Math::add(m4,b111.CreateIntegralNewU(b131,b131,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
	}

	if (nums[27] != 0) {
	// boundary 5
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(12),nums[27],nums[28],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());
		
		for (int i=b51; i<GetNumV()-b52; i++) 
			for (int j=b53; j<GetNumW()-b54; j++) {
				Vector<double> v11 = b1(GetKnotsU()[GetNumU()]);
				Vector<double> v21 = mv.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m5(GetNumV(),GetNumW(),0.0);
	if (nums[29] != 0) {
		if (natgeom[4] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(13),nums[29],nums[30],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b55; i<GetNumV()-b56; i++) 
				for (int j=b57; j<GetNumW()-b58; j++) {
					Vector<double> v11 = Math::mult4(b4(GetKnotsU()[GetNumU()]),GetKnotSetU().CreateMatrixDeriv());
					Vector<double> v21 = mv.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		}
	}
	
	if (natgeom[4] != 0) {
		BspSurf<double> b131= BspSurf<double>(PolySurf<double>(bound.GetUV(22),nums[47],nums[48],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());//.ConvertBspSurf();
		m5 = Math::add(m5,b121.CreateIntegralNewU(b131,b131,GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
	}


	if (nums[33] != 0) {

		// boundary 6
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(15),nums[33],nums[34],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());
	
		for (int i=b61; i<GetNumV()-b62; i++) 
			for (int j=b63; j<GetNumW()-b64; j++) {
				Vector<double> v11 = b1(0.0);
				Vector<double> v21 = mv.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m6(GetNumV(),GetNumW(),0.0);
	if (nums[35] != 0) {
		if (natgeom[5] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(16),nums[35],nums[36],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b65; i<GetNumV()-b66; i++) 
				for (int j=b67; j<GetNumW()-b68; j++) {
					Vector<double> v11 = Math::mult4(b4(0.0),GetKnotSetU().CreateMatrixDeriv());
					Vector<double> v21 = mv.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		}
	}
	
	if (natgeom[5] != 0) {
		BspSurf<double> b131= BspSurf<double>(PolySurf<double>(bound.GetUV(23),nums[49],nums[50],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());//.ConvertBspSurf();
		m6 = Math::add(m6,b121.CreateIntegralNewU(b131,b131,GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
	}

	j=0;
	for (int i=0; i<nums[2]; i++) {
		double u = supports[j];
		double v = supports[j+1];
		double w = supports[j+2];
		double val = supports[j+3];
		Vector<double> v1 =	b1(u);
		Vector<double> v2 = b2(v);
		Vector<double> v3 = b3(w);
		Matrix<double> mat = Math::kronecker(Matrix<double>(v2),Matrix<double>(v1));
		mat = Math::kronecker(Matrix<double>(v3),mat);
		constraints.InsertRow(mat.GetRow(0),row_count);
		constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = val;
		row_count++;
		j=j+4;
	}


	Matrix3D<double> mat4(GetNumU(),GetNumV(),GetNumW(),0.0);


	// add these explicitly!!!!
	mat4.InsertUV(m1,0);

	mat4.InsertUW(m2,GetNumV()-1);
	mat4.InsertUW(m4,0);
	mat4.InsertVW(m5,GetNumU()-1);

	mat4.InsertVW(m6,0);
	mat4.InsertUV(m3,GetNumW()-1);
	mat = Math::subtract(mat,mat4);

	
	// add them in
	Vector<double> vec1(Math::CreateKroneckerVector(mat));
	vec1 = Math::vmult(-1.0,vec1);
	Vector<int> v2 = Math::GetEliminateIndices(matfinal,constraints);
	Matrix<double> m10 = Math::EliminateMVariables(matfinal,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(matfinal,constraints,vec1,v2);
    Vector<double> sol = Math::gauss(m10,v1,v1.GetNum());
	
	Vector<Point1D> nsol(sol.GetNum());
	
	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m30(GetNumU()*GetNumV()*GetNumW(),GetNumU()*GetNumV()*GetNumW(),0.0);
	Vector<double> v3(GetNumU()*GetNumV()*GetNumW());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNumU()*GetNumV()*GetNumW());
		
	for (int i=0; i<GetNumU()*GetNumV()*GetNumW(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<num_constraint; i++)
		for (int j=0; j<GetNumU()*GetNumV()*GetNumW(); j++) 
			m30[i][j] = constraints[i][j];
		

	
	for (int i=num_constraint; i<GetNumU()*GetNumV()*GetNumW(); i++) m30[i][diff[i-num_constraint]] = 1.0;


	// rhs Vector

	for (int i=0; i<num_constraint; i++) v3[i]= constraints[i][GetNumU()*GetNumV()*GetNumW()];//0.0;


	for (int i=num_constraint; i<GetNumU()*GetNumV()*GetNumW(); i++) v3[i] = sol[i-num_constraint];


	// solve
	Vector<double> sol1 = Math::gauss(m30,v3,GetNumU()*GetNumV()*GetNumW());
	
	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);


	return FBspVol(Math::Matrix3DFromVector(nsol1,GetNumU(),GetNumV(),GetNumW()),GetKnotsU(),GetKnotsV(),GetKnotsW(),GetOrdU(),GetOrdV(),GetOrdW(),GetNumU(),GetNumV(),GetNumW());
}


FBspVol FBspVol::ComputeFiniteElementNew(double L, double W, double H, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, BspVol<double>& bvol, Vector<int>& natgeom,
										Matrix3D<double>& bound) const
{
	BspVolBasisFuncSet b(GetOrdU(),GetOrdV(),GetOrdW(),GetOrdU()+GetNumU(),GetOrdV()+GetNumV(),GetOrdW()+GetNumW(),GetKnotsU(),GetKnotsV(),GetKnotsW());
	BspSurfBasisFuncSet b101(GetOrdU(),GetOrdV(),GetOrdU()+GetNumU(),GetOrdV()+GetNumV(),GetKnotsU(),GetKnotsV());
	BspSurfBasisFuncSet b111(GetOrdU(),GetOrdW(),GetOrdU()+GetNumU(),GetOrdW()+GetNumW(),GetKnotsU(),GetKnotsW());
	BspSurfBasisFuncSet b121(GetOrdV(),GetOrdW(),GetOrdV()+GetNumV(),GetOrdW()+GetNumW(),GetKnotsV(),GetKnotsW());
	
	BspCurvBasisFuncSet b1(GetKnotSetU());
	BspCurvBasisFuncSet b2(GetKnotSetV());
	BspCurvBasisFuncSet b3(GetKnotSetW());
	BspCurvBasisFuncSet b4(GetKnotSetU().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b5(GetKnotSetV().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b6(GetKnotSetW().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b7(GetKnotSetU().CreateKnotSetDeriv(2));
	BspCurvBasisFuncSet b8(GetKnotSetV().CreateKnotSetDeriv(2));
	BspCurvBasisFuncSet b9(GetKnotSetW().CreateKnotSetDeriv(2));
	
	Matrix<double> mat1 = b1.CreateMatrixMinimisation(1);
	Matrix<double> mat2 = b2.CreateMatrixMinimisation(1);
	Matrix<double> mat3 = b3.CreateMatrixMinimisation(1);
	Matrix<double> mat5 = b1.CreateMatrixMinimisation(0);
	Matrix<double> mat6 = b2.CreateMatrixMinimisation(0);
	Matrix<double> mat7 = b3.CreateMatrixMinimisation(0);

	Matrix<double> vm1, vm2, vm3;

	vm1 = Math::kronecker(mat6,mat1);
	vm2 = Math::kronecker(mat2,mat5);
	vm3 = Math::kronecker(mat6,mat5);

	vm1 = Math::kronecker(mat7,vm1);
	vm2 = Math::kronecker(mat7,vm2);
	vm3 = Math::kronecker(mat3,vm3);

	Matrix<double> matfinal = Math::add(Math::add(vm1,vm2),vm3);

	Matrix3D<double> mat(GetNumU(),GetNumV(),GetNumW(),0.0);
	int j = 0;
	for (int i=0; i<nums[0]; i++) {
		mat = Math::add(mat,Math::mmult(ptforce[j+3],b(ptforce[j],ptforce[j+1],ptforce[j+2])));
		j+=4;
	}

	j = 0;
	int ind = 0;
	
	for (int i=0; i<nums[1]; i++) {
	
		mat =Math::add(mat,b.CreateIntegralNewU(bvol,bvol,disforce[j],disforce[j+1],disforce[j+2],disforce[j+3],disforce[j+4],disforce[j+5]));
		j+=6;
		ind=ind+3;
	}

	
	// boundary conditions
	// supports
	
	j=0;
	int row_count=0;

	int num_constraint = nums[2];

	int b11, b12 , b13 , b14, b15, b16, b17, b18;
	int b21, b22 , b23 , b24, b25, b26, b27, b28;
	int b31, b32 , b33 , b34, b35, b36, b37, b38;
	int b41, b42 , b43 , b44, b45, b46, b47, b48;
	int b51, b52 , b53 , b54, b55, b56, b57, b58;
	int b61, b62 , b63 , b64, b65, b66, b67, b68;



	// Face 1 (uv at w=0)
	if (nums[3] != 0) {
		num_constraint += GetNumU()*GetNumV();
		b11 = b12 = b13 = b14 = 0;
	}

	if (nums[5] != 0)
		if (natgeom[0] == 0) {
			num_constraint += (GetNumU())*(GetNumV());
			b15 = b16 = b17 = b18 = 0;
		}

	// Face 2 (uw at v=1)
	if (nums[9] != 0 && nums[5] != 0) {
		num_constraint += GetNumU()*(GetNumW()-2);
		b21 = b22 = 0; b23 = 2; b24 = 0;
	}
	else if (nums[9] != 0 && nums[3] != 0) {
		num_constraint += GetNumU()*(GetNumW()-1);
		b21 = b22 = 0; b23 = 1; b24 = 0;
	}
	else if (nums[9] != 0) {
		num_constraint += GetNumU()*GetNumW();
		b21 = b22 = b23 = b24 = 0;
	}

	if (natgeom[1] == 0) 
		if (nums[11] != 0) {
			if (nums[3] != 0 && nums[5] == 0) {
				num_constraint += GetNumU()*(GetNumW()-1);
				b25 = b26 = 0; b27 = 1; b28 = 0;
			}
			else if (nums[3] != 0 && nums[5] != 0) {
				num_constraint += GetNumU()*(GetNumW()-2);
				b25 = 0; b26 = 0; b27 = 2; b28 = 0;
			}
			else {
				num_constraint += GetNumU()*GetNumW();
				b25 = b26 = b27 = b28 = 0;
			}
		}


	// Face 3 (uv at w=1)
	if (nums[15] != 0 && nums[11] !=0) {
		num_constraint += GetNumU()*(GetNumV()-2);
		b31 = b32 = b33 = 0; b34 = 2;
	}
	else if (nums[15] != 0 && nums[9] != 0) {
		num_constraint += GetNumU()*(GetNumV()-1);
		b31 = b32 = b33 = 0; b34 = 1;
	}
	else if (nums[15] != 0) {
		num_constraint += GetNumU()*GetNumV();
		b31 = b32 = b33 = b34 = 0;
	}

	if (natgeom[2] == 0) 
		if (nums[17] != 0) {
			if (nums[9] != 0 && nums[11] == 0) {
				num_constraint += GetNumU()*(GetNumV()-1);
				b35 = b36 = b37 = 0; b38 = 1;
			}
			else if (nums[9] != 0 && nums[11] != 0) {
				num_constraint += GetNumU()*(GetNumV()-2);
				b35 = b26 = b37 =0; b38 = 1;
			}
			else {
				num_constraint += GetNumU()*GetNumV();
				b35 = b36 = b37 = b38 = 0;
			}
		}
		

	// face 4 (uw at v = 0)
	if (nums[21] != 0 && nums[5] != 0 && nums[17] != 0) {
		num_constraint += GetNumU()*(GetNumW()-4);
		b41 = b42 = 0; b43 = 2; b44 = 2;
	}
	else if (nums[21] != 0 && nums[5] != 0 && nums[15] != 0) {
		num_constraint += GetNumU()*(GetNumW()-3);
		b41 = b42 = 0; b43 = 2; b44 = 1;
	}
	else if (nums[21] != 0 && nums[5] != 0 && nums[15] == 0) {
		num_constraint += GetNumU()*(GetNumW()-2);
		b41 = b42 = 0; b43 = 2; b44 = 0;
	}
	else if (nums[21] != 0 && nums[3] != 0 && nums[17] != 0) {
		num_constraint += GetNumU()*(GetNumW()-3);
		b41 = b42 = 0; b43 = 1; b44 = 2;
	}
	else if (nums[21] != 0 && nums[3] != 0 && nums[15] != 0) {
		num_constraint += GetNumU()*(GetNumW()-2);
		b41 = b42 = 0; b43 = b44 = 1;
	}
	else if (nums[21] != 0 && nums[3] != 0 && nums[15] == 0) {
		num_constraint += GetNumU()*(GetNumW()-1);
		b41 = b42 = 0; b43 = 1; b44 = 0;
	}
	else if (nums[21] != 0 && nums[3] == 0 && nums[15] != 0) {
		num_constraint += GetNumU()*(GetNumW()-1);
		b41 = b42 = b43 = 0; b44 = 1;
	}
	else if (nums[21] != 0 && nums[3] == 0 && nums[17] != 0) {
		num_constraint += GetNumU()*(GetNumW()-2);
		b41 = b42 = b43 = 0; b44 = 2;
	}
	else if (nums[21] != 0 && nums[3] == 0 && nums[15] == 0) {
		num_constraint += GetNumU()*GetNumW();
		b41 = b42 = b43 = b44 = 0;
	}

	if (natgeom[3] == 0) 
		if (nums[23] != 0) {
			if (nums[5] != 0 && nums[17] != 0) {
				num_constraint += GetNumU()*(GetNumW()-4);
				b45 = b46 = 0; b47 = b48 = 2;
			}
			else if (nums[5] != 0 && nums[15] != 0) {
				num_constraint += GetNumU()*(GetNumW()-3);
				b45 = b46 = 0; b47 = 2; b48 = 1;
			}
			else if (nums[5] != 0 && nums[15] == 0) {
				num_constraint += GetNumU()*(GetNumW()-2);
				b45 = b46 = 0; b47 = 2; b48 = 0;
			}
			else if (nums[3] != 0 && nums[17] != 0) {
				num_constraint += GetNumU()*(GetNumU()-3);
				b45 = b46 = 0; b47 = 1; b48 = 2;
			}
			else if (nums[3] != 0 && nums[15] != 0) {
				num_constraint += GetNumU()*(GetNumW()-2);
				b45 = b46 = 0; b47 = b48 = 1;
			}
			else if (nums[3] != 0 && nums[15] == 0) {
				num_constraint += GetNumU()*(GetNumW()-1);
				b45 = b46 = 0; b47 = 1; b48 = 0;
			}
			else if (nums[3] == 0 && nums[17] != 0) {
				num_constraint += GetNumU()*(GetNumW()-2);
				b45 = b46 = b47 = 0; b48 = 2;
			}
			else if (nums[3] == 0 && nums[15] != 0) {
				num_constraint += GetNumU()*(GetNumW()-1);
				b45 = b46 = b47 = 0; b48 = 1;
			}
			else if (nums[3] == 0 && nums[15] == 0) {
				num_constraint += GetNumU()*GetNumW();
				b45 = b46 = b47 = b48 = 0;
			}
		}
		
	// Face 5 (vw at u=1)
	if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-4);
		b51 = b52 = b53 = b54 = 2;
	}

	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-4);
		b51 = 1; b52 = b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b51 = 0; b52 = b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-3);
		b51 = b52 = b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b51 = 1; b52 = 2; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 0; b52 = 2; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b51 = b52 = b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 = 1; b52 = 2; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 0; b52 = 2; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b51 = 2; b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b51 = 1; b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-4);
		b51 = b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 =2 ; b52 = 0; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 1; b52 = 0; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b51 = b52 = 0; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = 0; b53 = 2; b54 = 0;
	}


	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-4);
		b51 = 2; b52 = 1; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b51 = b52 = 1; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b51 = 0 ; b52 = 1; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b51 = 2; b52 = 1; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = b52 = 1; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 0; b52 = 1; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 0; b52 = 1; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 = 2; b52 = 1; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = b52 = 1; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b51 = 2; b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b51 = 1; b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-4);
		b51 = b52 = 0; b53 = b54 = 2;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 2; b52 = 0; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 1; b52 = 0; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b51 = b52 = 0; b53 = 2; b54 = 1;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = 2; b54 = 0;
	}
	else if (nums[27] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = 0; b53 = 2; b54 = 0;
	}


	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-3);
		b51 = b52 = 2; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b51 = 1; b52 = 2; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 0; b52 = 2; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b51 = b52 = 2; b53 = 1; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 = 1; b52 = 2; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 0; b52 = 2; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-1);
		b51 = b52 = 2; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b51 = 1; b52 = 2; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 0; b52 = 2; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 2; b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 1; b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b51 = b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 2; b52 = 0; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 1; b52 = 0; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b51 = b52 = 0; b53 = 1; b54 = 0;
	}

	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b51 = b52 = 2; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 =1; b52 = 2; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 0; b52 = 2; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-1);
		b51 = b52 = 2; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b51 = 1; b52 = 2; b53 = 0; b54 = 1;
	}

	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 0; b52 = 2; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*GetNumW();
		b51 = b52 = 2; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*GetNumW();
		b51 = 1; b52 = 2; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b51 = 0; b52 = 2; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = 0; b54 = 2; 
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 2; b52 = 0; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 1; b52 = 0; b53 = 0; b54 = 1; 
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b51 = b52  = b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b51 = 2; b52 = b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b51 = 1; b52 = b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*GetNumW();
		b51 = b52 = b53 = b54 = 0;
	}

	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b51 = 2; b52 = 1; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = b52 = 1; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 0; b52 = 1; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 = 2; b52 = 1; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = b52 = b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 0; b52 = b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b51 = 2; b52 = 1; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = b52 = b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 0; b52 = 1; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b51 = 2; b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b51 = 1; b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b51 = b52 = 0; b53 = 1; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = 0; b53 = b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 2; b52 = 0; b53 = 1; b54 = 0;
	}
		
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 1; b52 = 0; b53 = 1; b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b51 = b52 = 0; b53 = 1; b54 = 0;
	}

	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b51 =2; b52 = 1; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = b52 = 1; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 0; b52 = 1; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b51 = 2; b52 = 1; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = b52 = 1; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 0; b52 = 1; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*GetNumW();
		b51 = 2; b52 = 1; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b51 = b52 = 1; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b51 = 0; b52 = 1; b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b51 = 2; b52 = 0; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b51 = 1; b52 = 0; b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b51 = b52 = b53 = 0; b54 = 2;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b51 = 2; b52 = b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b51 = 1; b52 = 0; b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b51 = b52 = b53 = 0; b54 = 1;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b51 = 2; b52 = b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b51 = 1; b52 = b53 = b54 = 0;
	}
	else if (nums[27] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*GetNumW();
		b51 = b52 = b53 = b54 = 0;
	}

	if (natgeom[4] == 0) 
		if (nums[29] != 0) {
			if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-4);
				b55 = b56 = b57 = b58 = 2;
			}

			else if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-4);
				b55 = 1; b56 = b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b55 = 0; b56 = b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-3);
				b55 = b56 = b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b55 = 1; b56 = 2; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 0; b56 = 2; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b55 = b56 = b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 = 1; b56 = 2; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 0; b56 = 2; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b55 = 2; b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b55 = 1; b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-4);
				b55 = b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 =2 ; b56 = 0; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 1; b56 = 0; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b55 = b56 = 0; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = 0; b57 = 2; b58 = 0;
			}


			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-4);
				b55 = 2; b56 = 1; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b55 = b56 = 1; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b55 = 0 ; b56 = 1; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b55 = 2; b56 = 1; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = b56 = 1; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 0; b56 = 1; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 0; b56 = 1; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 = 2; b56 = 1; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = b56 = 1; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b55 = 2; b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b55 = 1; b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-4);
				b55 = b56 = 0; b57 = b58 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 2; b56 = 0; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 1; b56 = 0; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b55 = b56 = 0; b57 = 2; b58 = 1;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = 2; b58 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = 0; b57 = 2; b58 = 0;
			}


			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-3);
				b55 = b56 = 2; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b55 = 1; b56 = 2; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 0; b56 = 2; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b55 = b56 = 2; b57 = 1; b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 = 1; b56 = 2; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 0; b56 = 2; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-1);
				b55 = b56 = 2; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b55 = 1; b56 = 2; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 0; b56 = 2; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 2; b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 1; b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b55 = b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 2; b56 = 0; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 1; b56 = 0; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b55 = b56 = 0; b57 = 1; b58 = 0;
			}

			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b55 = b56 = 2; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 =1; b56 = 2; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 0; b56 = 2; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-1);
				b55 = b56 = 2; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b55 = 1; b56 = 2; b57 = 0; b58 = 1;
			}

			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 0; b56 = 2; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*GetNumW();
				b55 = b56 = 2; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*GetNumW();
				b55 = 1; b56 = 2; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b55 = 0; b56 = 2; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = 0; b58 = 2; 
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 2; b56 = 0; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 1; b56 = 0; b57 = 0; b58 = 1; 
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b55 = b56  = b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b55 = 2; b56 = b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b55 = 1; b56 = b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*GetNumW();
				b55 = b56 = b57 = b58 = 0;
			}

			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b55 = 2; b56 = 1; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = b56 = 1; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 0; b56 = 1; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 = 2; b56 = 1; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = b56 = b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 0; b56 = b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b55 = 2; b56 = 1; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = b56 = b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 0; b56 = 1; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b55 = 2; b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b55 = 1; b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b55 = b56 = 0; b57 = 1; b58 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = 0; b57 = b58 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 2; b56 = 0; b57 = 1; b58 = 0;
			}
				
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 1; b56 = 0; b57 = 1; b58 = 0;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b55 = b56 = 0; b57 = 1; b58 = 0;
			}

			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b55 =2; b56 = 1; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = b56 = 1; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 0; b56 = 1; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b55 = 2; b56 = 1; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = b56 = 1; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 0; b56 = 1; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*GetNumW();
				b55 = 2; b56 = 1; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b55 = b56 = 1; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b55 = 0; b56 = 1; b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b55 = 2; b56 = 0; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b55 = 1; b56 = 0; b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b55 = b56 = b57 = 0; b58 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b55 = 2; b56 = b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b55 = 1; b56 = 0; b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b55 = b56 = b57 = 0; b58 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b55 = 2; b56 = b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b55 = 1; b56 = b57 = b58 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*GetNumW();
				b55 = b56 = b57 = b58 = 0;
			}
		}


		
			

	// Face 6 (vw at u=0)
	if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-4);
		b61 = b62 = b63 = b64 = 2;
	}

	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-4);
		b61 = 1; b62 = b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b61 = 0; b62 = b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-3);
		b61 = b62 = b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b61 = 1; b62 = 2; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 0; b62 = 2; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b61 = b62 = b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 = 1; b62 = 2; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 0; b62 = 2; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b61 = 2; b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b61 = 1; b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-4);
		b61 = b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 =2 ; b62 = 0; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 1; b62 = 0; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b61 = b62 = 0; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = 0; b63 = 2; b64 = 0;
	}


	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-4);
		b61 = 2; b62 = 1; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b61 = b62 = 1; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b61 = 0 ; b62 = 1; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b61 = 2; b62 = 1; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = b62 = 1; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 0; b62 = 1; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 0; b62 = 1; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 = 2; b62 = 1; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = b62 = 1; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-4);
		b61 = 2; b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-4);
		b61 = 1; b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-4);
		b61 = b62 = 0; b63 = b64 = 2;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 2; b62 = 0; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 1; b62 = 0; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b61 = b62 = 0; b63 = 2; b64 = 1;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = 2; b64 = 0;
	}
	else if (nums[33] != 0 && nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = 0; b63 = 2; b64 = 0;
	}


	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-3);
		b61 = b62 = 2; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b61 = 1; b62 = 2; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 0; b62 = 2; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b61 = b62 = 2; b63 = 1; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 = 1; b62 = 2; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 0; b62 = 2; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-1);
		b61 = b62 = 2; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b61 = 1; b62 = 2; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 0; b62 = 2; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 2; b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 1; b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b61 = b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 2; b62 = 0; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 1; b62 = 0; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b61 = b62 = 0; b63 = 1; b64 = 0;
	}

	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-2);
		b61 = b62 = 2; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 =1; b62 = 2; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 0; b62 = 2; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*(GetNumW()-1);
		b61 = b62 = 2; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b61 = 1; b62 = 2; b63 = 0; b64 = 1;
	}

	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 0; b62 = 2; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-4)*GetNumW();
		b61 = b62 = 2; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-3)*GetNumW();
		b61 = 1; b62 = 2; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b61 = 0; b62 = 2; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = 0; b64 = 2; 
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 2; b62 = 0; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 1; b62 = 0; b63 = 0; b64 = 1; 
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b61 = b62  = b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b61 = 2; b62 = b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b61 = 1; b62 = b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*GetNumW();
		b61 = b62 = b63 = b64 = 0;
	}

	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-3);
		b61 = 2; b62 = 1; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = b62 = 1; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 0; b62 = 1; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 = 2; b62 = 1; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = b62 = b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 0; b62 = b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b61 = 2; b62 = 1; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = b62 = b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 0; b62 = 1; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-3);
		b61 = 2; b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-3);
		b61 = 1; b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-3);
		b61 = b62 = 0; b63 = 1; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = 0; b63 = b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 2; b62 = 0; b63 = 1; b64 = 0;
	}
		
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 1; b62 = 0; b63 = 1; b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b61 = b62 = 0; b63 = 1; b64 = 0;
	}

	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-2);
		b61 =2; b62 = 1; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = b62 = 1; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 0; b62 = 1; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*(GetNumW()-1);
		b61 = 2; b62 = 1; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = b62 = 1; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 0; b62 = 1; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-3)*GetNumW();
		b61 = 2; b62 = 1; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b61 = b62 = 1; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b61 = 0; b62 = 1; b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-2);
		b61 = 2; b62 = 0; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-2);
		b61 = 1; b62 = 0; b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-2);
		b61 = b62 = b63 = 0; b64 = 2;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*(GetNumW()-1);
		b61 = 2; b62 = b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*(GetNumW()-1);
		b61 = 1; b62 = 0; b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
		num_constraint += GetNumV()*(GetNumW()-1);
		b61 = b62 = b63 = 0; b64 = 1;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
		num_constraint += (GetNumV()-2)*GetNumW();
		b61 = 2; b62 = b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
		num_constraint += (GetNumV()-1)*GetNumW();
		b61 = 1; b62 = b63 = b64 = 0;
	}
	else if (nums[33] != 0 && nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
		num_constraint += GetNumV()*GetNumW();
		b61 = b62 = b63 = b64 = 0;
	}


	if (natgeom[5] == 0) 
		if (nums[35] != 0) {
			if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-4);
				b65 = b66 = b67 = b68 = 2;
			}

			else if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-4);
				b65 = 1; b66 = b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b65 = 0; b66 = b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-3);
				b65 = b66 = b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b65 = 1; b66 = 2; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 0; b66 = 2; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b65 = b66 = b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 = 1; b66 = 2; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 0; b66 = 2; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b65 = 2; b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b65 = 1; b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-4);
				b65 = b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 =2 ; b66 = 0; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 1; b66 = 0; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b65 = b66 = 0; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = 0; b67 = 2; b68 = 0;
			}


			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-4);
				b65 = 2; b66 = 1; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b65 = b66 = 1; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b65 = 0 ; b66 = 1; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b65 = 2; b66 = 1; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = b66 = 1; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 0; b66 = 1; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 0; b66 = 1; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 = 2; b66 = 1; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = b66 = 1; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-4);
				b65 = 2; b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-4);
				b65 = 1; b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-4);
				b65 = b66 = 0; b67 = b68 = 2;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 2; b66 = 0; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 1; b66 = 0; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b65 = b66 = 0; b67 = 2; b68 = 1;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = 2; b68 = 0;
			}
			else if (nums[5] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = 0; b67 = 2; b68 = 0;
			}


			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-3);
				b65 = b66 = 2; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b65 = 1; b66 = 2; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 0; b66 = 2; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b65 = b66 = 2; b67 = 1; b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 = 1; b66 = 2; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 0; b66 = 2; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-1);
				b65 = b66 = 2; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b65 = 1; b66 = 2; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 0; b66 = 2; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 2; b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 1; b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b65 = b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 2; b66 = 0; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 1; b66 = 0; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b65 = b66 = 0; b67 = 1; b68 = 0;
			}

			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-2);
				b65 = b66 = 2; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 =1; b66 = 2; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 0; b66 = 2; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*(GetNumW()-1);
				b65 = b66 = 2; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b65 = 1; b66 = 2; b67 = 0; b68 = 1;
			}

			else if (nums[3] == 0 && nums[11] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 0; b66 = 2; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-4)*GetNumW();
				b65 = b66 = 2; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-3)*GetNumW();
				b65 = 1; b66 = 2; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[11] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b65 = 0; b66 = 2; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = 0; b68 = 2; 
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 2; b66 = 0; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 1; b66 = 0; b67 = 0; b68 = 1; 
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b65 = b66  = b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b65 = 2; b66 = b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b65 = 1; b66 = b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[11] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*GetNumW();
				b65 = b66 = b67 = b68 = 0;
			}

			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-3);
				b65 = 2; b66 = 1; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = b66 = 1; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 0; b66 = 1; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 = 2; b66 = 1; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = b66 = b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 0; b66 = b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b65 = 2; b66 = 1; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = b66 = b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 0; b66 = 1; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-3);
				b65 = 2; b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-3);
				b65 = 1; b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-3);
				b65 = b66 = 0; b67 = 1; b68 = 2;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = 0; b67 = b68 = 1;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 2; b66 = 0; b67 = 1; b68 = 0;
			}
				
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 1; b66 = 0; b67 = 1; b68 = 0;
			}
			else if (nums[3] != 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b65 = b66 = 0; b67 = 1; b68 = 0;
			}

			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-2);
				b65 =2; b66 = 1; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = b66 = 1; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 0; b66 = 1; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*(GetNumW()-1);
				b65 = 2; b66 = 1; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = b66 = 1; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 0; b66 = 1; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-3)*GetNumW();
				b65 = 2; b66 = 1; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b65 = b66 = 1; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] != 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b65 = 0; b66 = 1; b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-2);
				b65 = 2; b66 = 0; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-2);
				b65 = 1; b66 = 0; b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[17] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-2);
				b65 = b66 = b67 = 0; b68 = 2;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*(GetNumW()-1);
				b65 = 2; b66 = b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*(GetNumW()-1);
				b65 = 1; b66 = 0; b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] != 0 && nums[21] == 0) {
				num_constraint += GetNumV()*(GetNumW()-1);
				b65 = b66 = b67 = 0; b68 = 1;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[23] != 0) {
				num_constraint += (GetNumV()-2)*GetNumW();
				b65 = 2; b66 = b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] != 0) {
				num_constraint += (GetNumV()-1)*GetNumW();
				b65 = 1; b66 = b67 = b68 = 0;
			}
			else if (nums[3] == 0 && nums[9] == 0 && nums[15] == 0 && nums[21] == 0) {
				num_constraint += GetNumV()*GetNumW();
				b65 = b66 = b67 = b68 = 0;
			}
		}

	
	Matrix<double> constraints(num_constraint,GetNumU()*GetNumV()*GetNumW()+1,0.0);
	Matrix<double> mu = Math::ComputeIdentityMatrix(GetNumU());
	Matrix<double> mv = Math::ComputeIdentityMatrix(GetNumV());
	Matrix<double> mw = Math::ComputeIdentityMatrix(GetNumW());

		if (nums[3] != 0) {
	// boundary 1
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(0),nums[3],nums[4],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());
		for (int i=b11; i<GetNumU()-b12; i++) 
			for (int j=b13; j<GetNumV()-b14; j++) {
				Vector<double> v11 = b3(0.0);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mv.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v11),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m1(GetNumU(),GetNumV(),0.0);
	if (nums[5] != 0) {
		if (natgeom[0] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(1),nums[5],nums[6],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b15; i<GetNumU()-b16; i++) 
				for (int j=b17; j<GetNumV()-b18; j++) {
					Vector<double> v11 = Math::mult4(b6(0.0),GetKnotSetW().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mv.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v11),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b131= PolySurf<double>(bound.GetUV(2),nums[7],nums[8],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).ConvertBspSurf();
			m1 = Math::add(m1,b101.CreateIntegralNewU(b131,b131,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]));
		}
	}
	// boundary 2

	if (nums[9] != 0) {
	// boundary 2
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(3),nums[9],nums[10],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());
	
		for (int i=b21; i<GetNumU()-b22; i++) 
			for (int j=b23; j<GetNumW()-b24; j++) {
				Vector<double> v11 = b2(GetKnotsV()[GetNumV()]);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m2(GetNumU(),GetNumW(),0.0);
	if (nums[11] != 0) {
		if (natgeom[1] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(4),nums[11],nums[12],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b25; i<GetNumU()-b26; i++) 
				for (int j=b27; j<GetNumW()-b28; j++) {
					Vector<double> v11 = Math::mult4(b5(GetKnotsV()[GetNumV()]),GetKnotSetV().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b131= PolySurf<double>(bound.GetUV(5),nums[13],nums[14],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]).ConvertBspSurf();
			m2 = Math::add(m2,b111.CreateIntegralNewU(b131,b131,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
		}
	}


	// boundary 3
	if (nums[15] != 0) {
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(6),nums[15],nums[16],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());
		
		for (int i=b31; i<GetNumU()-b32; i++) 
			for (int j=b33; j<GetNumV()-b34; j++) {
				Vector<double> v11 = b3(GetKnotsW()[GetNumW()]);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mv.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v11),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m3(GetNumU(),GetNumV(),0.0);
	if (nums[17] != 0) {
		if (natgeom[2] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(7),nums[17],nums[18],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b35; i<GetNumU()-b36; i++) 
				for (int j=b37; j<GetNumV()-b38; j++) {
					Vector<double> v11 = Math::mult4(b6(GetKnotsW()[GetNumW()]),GetKnotSetW().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mv.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v31),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v11),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b131 = PolySurf<double>(bound.GetUV(8),nums[19],nums[20],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).ConvertBspSurf();
			m3 = Math::add(m3,b101.CreateIntegralNewU(b131,b131,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]));
		}
	}

	// boundary 4

	if (nums[21] != 0) {
	// boundary 4
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(9),nums[21],nums[22],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());
		for (int i=b41; i<GetNumU()-b42; i++) 
			for (int j=b43; j<GetNumW()-b44; j++) {
				Vector<double> v11 = b2(0.0);
				Vector<double> v21 = mu.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m4(GetNumU(),GetNumW(),0.0);
	if (nums[23] != 0) {
		if (natgeom[3] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(10),nums[23],nums[24],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetU(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b45; i<GetNumU()-b46; i++) 
				for (int j=b47; j<GetNumW()-b48; j++) {
					Vector<double> v11 = Math::mult4(b5(0.0),GetKnotSetV().CreateMatrixDeriv());
					Vector<double> v21 = mu.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b131= PolySurf<double>(bound.GetUV(11),nums[25],nums[26],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]).ConvertBspSurf();
			m4 = Math::add(m4,b111.CreateIntegralNewU(b131,b131,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
		}
	}

	if (nums[27] != 0) {
	// boundary 5
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(12),nums[27],nums[28],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());
		
		for (int i=b51; i<GetNumV()-b52; i++) 
			for (int j=b53; j<GetNumW()-b54; j++) {
				Vector<double> v11 = b1(GetKnotsU()[GetNumU()]);
				Vector<double> v21 = mv.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}

	Matrix<double> m5(GetNumV(),GetNumW(),0.0);
	if (nums[29] != 0) {
		if (natgeom[4] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(13),nums[29],nums[30],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b55; i<GetNumV()-b56; i++) 
				for (int j=b57; j<GetNumW()-b58; j++) {
					Vector<double> v11 = Math::mult4(b4(GetKnotsU()[GetNumU()]),GetKnotSetU().CreateMatrixDeriv());
					Vector<double> v21 = mv.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b131= PolySurf<double>(bound.GetUV(14),nums[31],nums[32],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]).ConvertBspSurf();
			m5 = Math::add(m5,b121.CreateIntegralNewU(b131,b131,GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
		}
	}

	if (nums[33] != 0) {

		// boundary 6
		BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(15),nums[33],nums[34],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());
	
		for (int i=b61; i<GetNumV()-b62; i++) 
			for (int j=b63; j<GetNumW()-b64; j++) {
				Vector<double> v11 = b1(0.0);
				Vector<double> v21 = mv.GetRow(i);
				Vector<double> v31 = mw.GetRow(j);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
				mat = Math::kronecker(Matrix<double>(v31),mat);
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
				row_count++;
			}
	}


	Matrix<double> m6(GetNumV(),GetNumW(),0.0);
	if (nums[35] != 0) {
		if (natgeom[5] == 0) {
			BspSurf<double> b111 = BspSurf<double>(PolySurf<double>(bound.GetUV(16),nums[35],nums[36],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]),GetKnotSetV(),GetKnotSetW());//PolySurf<double>(bound2,nums[4],nums[5],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
			for (int i=b65; i<GetNumV()-b66; i++) 
				for (int j=b67; j<GetNumW()-b68; j++) {
					Vector<double> v11 = Math::mult4(b4(0.0),GetKnotSetU().CreateMatrixDeriv());
					Vector<double> v21 = mv.GetRow(i);
					Vector<double> v31 = mw.GetRow(j);
					Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
					mat = Math::kronecker(Matrix<double>(v31),mat);
					constraints.InsertRow(mat.GetRow(0),row_count);
					constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = b111.GetCPoints()[i][j];
					row_count++;
				}
		} else {
			BspSurf<double> b131= PolySurf<double>(bound.GetUV(17),nums[37],nums[38],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]).ConvertBspSurf();
			m6 = Math::add(m6,b121.CreateIntegralNewU(b131,b131,GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()],GetKnotsW()[GetOrdW()-1],GetKnotsW()[GetNumW()]));
		}
	}

	j=0;
	for (int i=0; i<nums[2]; i++) {
		double u = supports[j];
		double v = supports[j+1];
		double w = supports[j+2];
		double val = supports[j+3];
		Vector<double> v1 =	b1(u);
		Vector<double> v2 = b2(v);
		Vector<double> v3 = b3(w);
		Matrix<double> mat = Math::kronecker(Matrix<double>(v2),Matrix<double>(v1));
		mat = Math::kronecker(Matrix<double>(v3),mat);
		constraints.InsertRow(mat.GetRow(0),row_count);
		constraints[row_count][GetNumU()*GetNumV()*GetNumW()] = val;
		row_count++;
		j=j+4;
	}


	Matrix3D<double> mat4(GetNumU(),GetNumV(),GetNumW(),0.0);


	// add these explicitly!!!!
	mat4.InsertUV(m1,0);
	mat4.InsertUV(m3,GetNumW()-1);
	mat4.InsertUW(m2,GetNumV()-1);
	mat4.InsertUW(m4,0);
	mat4.InsertVW(m5,GetNumU()-1);
	mat4.InsertVW(m6,0);
		
	mat = Math::subtract(mat,mat4);

	// add them in
	Vector<double> vec1(Math::CreateKroneckerVector(mat));
	vec1 = Math::vmult(-1.0,vec1);
	Vector<int> v2 = Math::GetEliminateIndices(matfinal,constraints);	
	Matrix<double> m10 = Math::EliminateMVariables(matfinal,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(matfinal,constraints,vec1,v2);
    Vector<double> sol = Math::gauss(m10,v1,v1.GetNum());
	
	Vector<Point1D> nsol(sol.GetNum());
	
	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m30(GetNumU()*GetNumV()*GetNumW(),GetNumU()*GetNumV()*GetNumW(),0.0);
	Vector<double> v3(GetNumU()*GetNumV()*GetNumW());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNumU()*GetNumV()*GetNumW());
		
	for (int i=0; i<GetNumU()*GetNumV()*GetNumW(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<num_constraint; i++)
		for (int j=0; j<GetNumU()*GetNumV()*GetNumW(); j++) 
			m30[i][j] = constraints[i][j];
		

	
	for (int i=num_constraint; i<GetNumU()*GetNumV()*GetNumW(); i++) m30[i][diff[i-num_constraint]] = 1.0;


	// rhs Vector

	for (int i=0; i<num_constraint; i++) v3[i]= constraints[i][GetNumU()*GetNumV()*GetNumW()];//0.0;


	for (int i=num_constraint; i<GetNumU()*GetNumV()*GetNumW(); i++) v3[i] = sol[i-num_constraint];


	// solve
	Vector<double> sol1 = Math::gauss(m30,v3,GetNumU()*GetNumV()*GetNumW());
	
	
	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);


	return FBspVol(Math::Matrix3DFromVector(nsol1,GetNumU(),GetNumV(),GetNumW()),GetKnotsU(),GetKnotsV(),GetKnotsW(),GetOrdU(),GetOrdV(),GetOrdW(),GetNumU(),GetNumV(),GetNumW());
}


PBspVol3D::PBspVol3D(const FBspVol& v)
{
	Matrix3D<Point4D> m1(v.GetNumU(),v.GetNumV(),v.GetNumW());
	Vector<double> v1(v.GetKnotSetU().ComputeKnotSetAver());
	Vector<double> v2(v.GetKnotSetV().ComputeKnotSetAver());
	Vector<double> v3(v.GetKnotSetV().ComputeKnotSetAver());

	for (int k=0; k<v.GetNumW(); k++) 
		for (int i=0; i<v.GetNumU(); i++) 
			for (int j=0; j<v.GetNumV(); j++) 
				m1[k][i][j] = Point4D(v1[i],v2[j],v3[k],v.GetCPoints()[k][i][j].GetX());

}




