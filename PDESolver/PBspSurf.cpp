#include "pbspsurf.h"
#include "pbspcurv.h"


BspSurfDouble::BspSurfDouble() : BspSurf<double>() {}
BspSurfDouble::BspSurfDouble(const Matrix<double>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv, int Numu, int Numv) :
BspSurf<double>(Cpts,Ktsu,Ktsv,Ordu,Ordv,Numu,Numv) {}
BspSurfDouble::BspSurfDouble(const Matrix<double>& Cpts, int Ordu, int Ordv, int Numu, int Numv) : BspSurf<double>(Cpts,Ordu,Ordv,Numu,Numv) {}
BspSurfDouble::BspSurfDouble(const BspSurf<double>& b) : BspSurf<double>(b.GetCPoints(),
																b.GetKnotsU(),b.GetKnotsV(),b.GetOrdU(),b.GetOrdV(),b.GetNumU(),b.GetNumV())
{
}




FBspSurf::FBspSurf() : BspSurf<Point1D>() {}
FBspSurf::FBspSurf(const BspSurf<Point1D>& b) : BspSurf<Point1D>(b) {}
FBspSurf::FBspSurf(const Matrix<Point1D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv, int Numu, int Numv) :
BspSurf<Point1D>(Cpts,Ktsu,Ktsv,Ordu,Ordv,Numu,Numv) {}
FBspSurf::FBspSurf(const Matrix<Point1D>& Cpts, int Ordu, int Ordv, int Numu, int Numv) : BspSurf<Point1D>(Cpts,Ordu,Ordv,Numu,Numv) {}
FBspSurf::FBspSurf(const BspSurf<double>& b) : BspSurf<Point1D>(Matrix<Point1D>(b.GetNumU(),b.GetNumV()),
																b.GetKnotsU(),b.GetKnotsV(),b.GetOrdU(),b.GetOrdV(),b.GetNumU(),b.GetNumV())
{
	Matrix<Point1D> mat(b.GetNumU(),b.GetNumV());
	for (int i=0; i<b.GetNumU(); i++)
		for (int j=0; j<b.GetNumV(); j++)
			mat[i][j] = Point1D(b.GetCPoints()[i][j]);

	*this = BspSurf<Point1D>(mat,b.GetKnotsU(),b.GetKnotsV(),b.GetOrdU(),b.GetOrdV(),b.GetNumU(),b.GetNumV());	
}



FBspSurf FBspSurf::ComputePDESolution(double L, double W, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, Vector<Matrix<double> >& surface, Vector<int>& natgeom, Matrix<double>& bound) const
{
	BspSurfBasisFuncSet b(GetOrdU(),GetOrdV(),GetOrdU()+GetNumU(),GetOrdV()+GetNumV(),GetKnotsU(),GetKnotsV());
	BspCurvBasisFuncSet b1(GetKnotSetU());
	BspCurvBasisFuncSet b2(GetKnotSetV());
	BspCurvBasisFuncSet b3(GetKnotSetU().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b4(GetKnotSetV().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b5(GetKnotSetU().CreateKnotSetDeriv(2));
	BspCurvBasisFuncSet b6(GetKnotSetV().CreateKnotSetDeriv(2));

	Matrix<double> mat1 = b1.CreateMatrixMinimisation(1);
	Matrix<double> mat2 = b2.CreateMatrixMinimisation(1);
	Matrix<double> mat5 = b1.CreateMatrixMinimisation(0);
	Matrix<double> mat6 = b2.CreateMatrixMinimisation(0);

	mat1 = Math::kronecker(mat6,mat1);
	mat2 = Math::kronecker(mat2,mat5);

	Matrix<double> mat3 = Math::add(mat1,mat2);

	Matrix<double> mat(GetNumU(),GetNumV(),0.0);

	int j = 0;
	for (int i=0; i<nums[0]; i++) {
		mat = Math::add(mat,Math::mmult(ptforce[j+2],b(ptforce[j],ptforce[j+1])));
		j+=3;
	}

	j = 0;
	int ind=0;
	for (int i=0; i<nums[1]; i++) {
		BspSurf<double> bsurf = BspSurf<double>(PolySurf<double>(surface[i],orders[ind],orders[ind+1],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetU(),GetKnotSetV());
	
		mat =Math::add(mat,b.CreateIntegralNewU(bsurf,bsurf,disforce[j],disforce[j+1],disforce[j+2],disforce[j+3]));
		ind = ind+2;
		j+=4;
	}

	
	// boundary conditions
	// supports
	j=0;
	
	int row_count=0;

	int num_constraint = nums[2];
	if (nums[3] != 0) num_constraint += GetNumV();
	if (natgeom[0] == 0) 
		if (nums[4] != 0) num_constraint += GetNumV()-2; 
	if (nums[5] != 0) num_constraint += GetNumV();

	if (nums[6] != 0) num_constraint += GetNumU()-1;
	if (natgeom[1] == 0) 
		if (nums[7] != 0) num_constraint += GetNumU()-3;
	if (nums[8] != 0) num_constraint += GetNumU()-1;

	if (nums[9] != 0) num_constraint += GetNumV()-1;
	if (natgeom[2] == 0) 
		if (nums[10] != 0) num_constraint += GetNumV()-3;
	if (nums[11] != 0) num_constraint += GetNumV()-1;

	if (nums[12] != 0) num_constraint += GetNumU()-2;
	if (natgeom[3] == 0) 
		if (nums[13] != 0) num_constraint += GetNumU()-4;
	if (nums[14] != 0) num_constraint += GetNumU()-2;

	Matrix<double> constraints(num_constraint,GetNumU()*GetNumV()+1,0.0);

	Matrix<double> mi1 = Math::ComputeIdentityMatrix(GetNumV());
	Matrix<double> mi2 = Math::ComputeIdentityMatrix(GetNumU());

	if (nums[3] != 0) {
		BspCurv<double> b11 = BspCurv<double>(PolyCurv<double>(bound.GetRow(0),nums[3],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetV());
		for (int i=0; i<GetNumV(); i++) {
			Vector<double> v11 = b1(GetKnotsU()[GetOrdU() - 1]);  // change
			Vector<double> v21 = mi1.GetRow(i);
			Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
			constraints.InsertRow(mat.GetRow(0),row_count);
			constraints[row_count][GetNumU()*GetNumV()] = b11.GetCPoints()[i];
			row_count++;
		}
	}


	Vector<double> vecb1(GetNumV(),0.0);
	if (nums[4] != 0) {
		if (natgeom[0] == 0) {
			BspCurv<double> b11 = PolyCurv<double>(bound.GetRow(1),nums[4],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
		
			for (int i=1; i<mi1.GetNumRows()-1; i++) {
				Vector<double> v11 = Math::mult4(b3(GetKnotsU()[GetOrdU() - 1]), GetKnotSetU().CreateMatrixDeriv()); // change
				Vector<double> v21 = mi1.GetRow(i);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()] = b11.GetCPoints()[i];//v1[i];
				row_count++;
			}
		} else {
			BspCurv<double> b13= PolyCurv<double>(bound.GetRow(1),nums[4],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).ConvertBspCurv();
			vecb1 = Math::vadd(vecb1,b2.CreateIntegralNew(b13,GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]));
		}
	}

	// boundary 2
	
	if (nums[6] != 0) {
		BspCurv<double> b21 = BspCurv<double>(PolyCurv<double>(bound.GetRow(3),nums[6],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()]),GetKnotSetU());
		for (int i=1; i<GetNumU(); i++) {
			Vector<double> v11 = b2(GetKnotsV()[GetOrdV() - 1]); //change
			Vector<double> v21 = mi2.GetRow(i);
			Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
			constraints.InsertRow(mat.GetRow(0),row_count);
			constraints[row_count][GetNumU()*GetNumV()] = b21.GetCPoints()[i];
			row_count++;
		}
	}


	Vector<double> vecb2(GetNumU(),0.0);
	if (nums[7] != 0) {
		if (natgeom[1] == 0) {
			BspCurv<double> b11 = PolyCurv<double>(bound.GetRow(4),nums[7],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()]).Elevate(GetOrdU()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetU());
		
			for (int i=2; i<mi2.GetNumRows()-1; i++) {
				Vector<double> v11 = Math::mult4(b4(GetKnotsV()[GetOrdV() - 1]), GetKnotSetV().CreateMatrixDeriv());  // change
				Vector<double> v21 = mi2.GetRow(i);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()] = b11.GetCPoints()[i];//v1[i];
				row_count++;
			}
		} else {
			BspCurv<double> b12 = PolyCurv<double>(bound.GetRow(4),nums[7],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()]).ConvertBspCurv();
			vecb2 = Math::vadd(vecb2,b1.CreateIntegralNew(b12,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()]));
		}
	}

	// boundary 3

	if (nums[9] != 0) {
		BspCurv<double> b31 = BspCurv<double>(PolyCurv<double>(bound.GetRow(6),nums[9],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]),GetKnotSetV());
		for (int i=1; i<GetNumV(); i++) {
			Vector<double> v11 = b1(GetKnotsU()[GetNumU()]);
			Vector<double> v21 = mi1.GetRow(i);
			Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
			constraints.InsertRow(mat.GetRow(0),row_count);
			constraints[row_count][GetNumU()*GetNumV()] = b31.GetCPoints()[i];
			row_count++;
		}
	}

	Vector<double> vecb3(GetNumV(),0.0);
	if (nums[10] != 0) {
		if (natgeom[2] == 0) {
			BspCurv<double> b11 = PolyCurv<double>(bound.GetRow(7),nums[10],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).Elevate(GetOrdV()-1).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetV());
		
			for (int i=2; i<mi1.GetNumRows()-1; i++) {
				Vector<double> v11 = Math::mult4(b3(GetKnotsU()[GetNumU()]), GetKnotSetU().CreateMatrixDeriv());
				Vector<double> v21 = mi1.GetRow(i);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v21),Matrix<double>(v11));
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()] = b11.GetCPoints()[i];//v1[i];
				row_count++;
			}
		} else {
			BspCurv<double> b13= PolyCurv<double>(bound.GetRow(7),nums[10],GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]).ConvertBspCurv();
			vecb3 = Math::vadd(vecb3,b2.CreateIntegralNew(b13,GetKnotsV()[GetOrdV()-1],GetKnotsV()[GetNumV()]));
		}
	}


	// boundary 4

	
	if (nums[12] != 0) {
		BspCurv<double> b41 = BspCurv<double>(PolyCurv<double>(bound.GetRow(9),nums[12],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()]),GetKnotSetU());
		for (int i=1; i<GetNumU()-1; i++) {
			Vector<double> v11 = b2(GetKnotsV()[GetNumV()]);
			Vector<double> v21 = mi2.GetRow(i);
			Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
			constraints.InsertRow(mat.GetRow(0),row_count);
			constraints[row_count][GetNumU()*GetNumV()] = b41.GetCPoints()[i];
			row_count++;
		}
	}

	Vector<double> vecb4(GetNumU(),0.0);
	if (nums[13] != 0) {
		if (natgeom[3] == 0) {
			BspCurv<double> b11 = PolyCurv<double>(bound.GetRow(10),nums[13],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()]).Elevate(GetOrdU()-1-nums[13]).ConvertBspCurv().MakeKnotSetCompatable(GetKnotSetU());
		
			for (int i=2; i<mi2.GetNumRows()-2; i++) {
				Vector<double> v11 = Math::mult4(b4(GetKnotsV()[GetNumV()]), GetKnotSetV().CreateMatrixDeriv());
				Vector<double> v21 = mi2.GetRow(i);
				Matrix<double> mat = Math::kronecker(Matrix<double>(v11),Matrix<double>(v21));
				constraints.InsertRow(mat.GetRow(0),row_count);
				constraints[row_count][GetNumU()*GetNumV()] = b11.GetCPoints()[i];//v1[i];
				row_count++;
			}
		} else {
			BspCurv<double> b14 = PolyCurv<double>(bound.GetRow(10),nums[13],GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()]).ConvertBspCurv();
			vecb4 = Math::vadd(vecb4,b1.CreateIntegralNew(b14,GetKnotsU()[GetOrdU()-1],GetKnotsU()[GetNumU()]));
		}
	}


	Matrix<double> mat4(GetNumU(),GetNumV(),0.0);

	mat4.InsertCol(vecb2,0);
	mat4.InsertCol(vecb4,GetNumV()-1);

	mat = Math::subtract(mat,mat4);

	for (int i=0; i<GetNumV(); i++) mat[0][i] = mat[0][i]-vecb1[i];
	for (int i=0; i<GetNumV(); i++) mat[GetNumU()-1][i] = mat[GetNumU()-1][i]-vecb3[i];
	Vector<double> vec1(mat.CreateKroneckerVector());
	vec1 = Math::vmult(-1.0,vec1);

	j=0;
	for (int i=0; i<nums[2]; i++) {
		double u = supports[j];
		double v = supports[j+1];
		double val = supports[j+2];
		Vector<double> v1 =	b1(u);
		Vector<double> v2 = b2(v);
		Matrix<double> mat = Math::kronecker(Matrix<double>(v2),Matrix<double>(v1));
		constraints.InsertRow(mat.GetRow(0),row_count);
		constraints[row_count][GetNumU()*GetNumV()] = val;
		row_count++;
		j=j+3;
	}

	Vector<int> v2 = Math::GetEliminateIndices(mat3,constraints);	

	Matrix<double> m10 = Math::EliminateMVariables(mat3,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(mat3,constraints,vec1,v2);

	Vector<double> sol = Math::gauss(m10,v1,v1.GetNum());

	Vector<Point1D> nsol(sol.GetNum());

	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m30(GetNumU()*GetNumV(),GetNumU()*GetNumV(),0.0);
	Vector<double> v3(GetNumU()*GetNumV());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNumU()*GetNumV());
		
	for (int i=0; i<GetNumU()*GetNumV(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<num_constraint; i++)
		for (int j=0; j<GetNumU()*GetNumV(); j++) 
			m30[i][j] = constraints[i][j];
		
	
	for (int i=num_constraint; i<GetNumU()*GetNumV(); i++) m30[i][diff[i-num_constraint]] = 1.0;


	// rhs Vector

	for (int i=0; i<num_constraint; i++) v3[i]= constraints[i][GetNumU()*GetNumV()];//0.0;


	for (int i=num_constraint; i<GetNumU()*GetNumV(); i++) v3[i] = sol[i-num_constraint];


	// solve
	Vector<double> sol1 = Math::gauss(m30,v3,GetNumU()*GetNumV());

	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);


	return FBspSurf(Math::MatrixFromVector(nsol1,GetNumU(),GetNumV()),GetKnotsU(),GetKnotsV(),GetOrdU(),GetOrdV(),GetNumU(),GetNumV());
}



PBspSurf3D::PBspSurf3D()  {}
PBspSurf3D::PBspSurf3D(const BspSurf<Point3D>& b) : BspSurf<Point3D>(b)
{
}


PBspSurf3D::PBspSurf3D(const BspSurf<Point4D>& b) : BspSurf<Point3D>(b.GetCPoints(), b.GetKnotsU(), b.GetKnotsV(), b.GetOrdU(), b.GetOrdV(), b.GetNumU(), b.GetNumV())
{
}

PBspSurf3D::PBspSurf3D(const Matrix<Point3D>& Cpts, const Vector<double>& Ktsu, const Vector<double>& Ktsv, int Ordu, int Ordv, int Numu, int Numv) :
BspSurf<Point3D>(Cpts, Ktsu, Ktsv, Ordu, Ordv, Numu, Numv) {}
PBspSurf3D::PBspSurf3D(const Matrix<Point3D>& Cpts, int Ordu, int Ordv, int Numu, int Numv) : BspSurf<Point3D>(Cpts, Ordu, Ordv, Numu, Numv) {}


PBspSurf3D::PBspSurf3D(const FBspSurf& s1, const FBspSurf& s2, const FBspSurf& s3)
{
	Matrix<Point3D> m1(s1.GetNumU(), s1.GetNumV());

	for (int i = 0; i<s1.GetNumU(); i++)
	for (int j = 0; j<s1.GetNumV(); j++)
		m1[i][j] = Point3D(s1.GetCPoints()[i][j].GetX(), s2.GetCPoints()[i][j].GetX(), s3.GetCPoints()[i][j].GetX());

	*this = PBspSurf3D(m1, s1.GetKnotsU(), s1.GetKnotsV(), s1.GetOrdU(), s1.GetOrdV(), s1.GetNumU(), s1.GetNumV());
}

PBspSurf3D::PBspSurf3D(const FBspSurf& s) : BspSurf<Point3D>()
{
	Matrix<Point3D> m1(s.GetNumU(), s.GetNumV());
	Vector<double> v1(s.GetKnotSetU().ComputeKnotSetAver());
	Vector<double> v2(s.GetKnotSetV().ComputeKnotSetAver());


	for (int i = 0; i<s.GetNumU(); i++)
	for (int j = 0; j<s.GetNumV(); j++)
		m1[i][j] = Point3D(v1[i], v2[j], s.GetCPoints()[i][j].GetX());

	*this = PBspSurf3D(m1, s.GetKnotsU(), s.GetKnotsV(), s.GetOrdU(), s.GetOrdV(), s.GetNumU(), s.GetNumV());
}




