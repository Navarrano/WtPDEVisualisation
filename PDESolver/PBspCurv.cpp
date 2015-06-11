
#include "pbspcurv.h"
#include "curves.h"


BspCurvDouble::BspCurvDouble() { }
BspCurvDouble::BspCurvDouble(const BspCurv<double>& b) : BspCurv<double>(b) { }
BspCurvDouble::BspCurvDouble(const Vector<double>& Cpts, const Vector<double>& Kts, int Ord, int Num) :
BspCurv<double>(Cpts, Kts, Ord, Num) {}

BspCurvDouble::BspCurvDouble(const Vector<double>& Cpts, int Ord, int Num) : BspCurv<double>(Cpts, Ord, Num) {}


FBspCurv::FBspCurv() { }
FBspCurv::FBspCurv(const BspCurv<Point1D>& b) : BspCurv<Point1D>(b) { }
FBspCurv::FBspCurv(const Vector<Point1D>& Cpts, const Vector<double>& Kts, int Ord, int Num) :
BspCurv<Point1D>(Cpts, Kts, Ord, Num) {}

FBspCurv::FBspCurv(const Vector<Point1D>& Cpts, int Ord, int Num) : BspCurv<Point1D>(Cpts, Ord, Num) {}
FBspCurv::FBspCurv(const BspCurv<double>& b) : BspCurv<Point1D>(Vector<Point1D>(b.GetNum()),
																b.GetKnots(),b.GetOrd(),b.GetNum())
{
	Vector<Point1D> vec(b.GetNum());
	for (int i=0; i<b.GetNum(); i++)
		vec[i] = Point1D(b.GetCPoints()[i]);

	*this = BspCurv<Point1D>(vec,b.GetKnots(),b.GetOrd(),GetNum());	
}



FBspCurv FBspCurv::ComputeDiffEquation1(double L, double C, Vector<int>& nums, Vector<double>& supports, Vector<double>& clamped, Vector<double>& ptforce, Vector<double>& disforce) const
{
	BspCurvBasisFuncSet b(GetKnots(),GetOrd(),GetOrd()+GetNum());
	
	
	Matrix<double> mat = Math::mmult(C,b.CreateMatrixMinimisation(1,0,L));
	mat = Math::add(mat,Math::mmult(-C,b.CreateMatrixMinimisation(0,0,L)));

	
	Vector<double> vec1(GetNum(),0.0);
	int j = 0;
	for (int i=0; i<nums[1]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(ptforce[j+1],b(ptforce[j])));
		j++;
	}

	j = 0;
	Linear l(0.0,1.0,1.0,0.0);

	for (int i=0; i<nums[2]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(disforce[j+2],b.CreateIntegralProduct(l,disforce[j],disforce[j+1])));
		vec1 = Math::vsubtract(vec1,Math::vmult(disforce[j+2],b.CreateIntegralIntegral(disforce[j],disforce[j+1])));
		vec1 = Math::vmult(1.0,vec1);
		j+=3;
	}


	// boundary conditions
	// supports

	Matrix<double> constraints(nums[3]+2*nums[0],GetNum()+1);
	for (int i=0; i<nums[3]; i++) {
		Vector<double> v = GetKnotSet().CreateVectorInterp(supports[i]);
		constraints.InsertRow(v,i);
	}


	// clamped
	j = 0;
	std::cerr << GetKnotSet();
	for (int i=0; i<nums[0]; i++) {
		std::cerr << clamped[i] << "\n";
		Vector<double> v1 = GetKnotSet().CreateVectorInterp(clamped[i]);
		Vector<double> v2 = GetKnotSet().CreateVectorInterpDeriv(2,clamped[i]);

		constraints.InsertRow(v1,nums[3]+j);
		constraints.InsertRow(v2,nums[3]+j+1);
		j+=2;
	}


	Vector<double> rhs(nums[3]+2*nums[0],0.0);

	constraints.InsertCol(rhs,GetNum());


	Vector<int> v2 = Math::GetEliminateIndices(mat,constraints);	
	Matrix<double> m1 = Math::EliminateMVariables(mat,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(mat,constraints,vec1,v2);


	Vector<double> sol = Math::gauss(m1,v1,v1.GetNum());
	Vector<Point1D> nsol (sol.GetNum());

	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m3(GetNum(),GetNum(),0.0);
	Vector<double> v3(GetNum());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNum());
		
	for (int i=0; i<GetNum(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<nums[3]+2*nums[0]; i++)
		for (int j=0; j<GetNum(); j++)
			m3[i][j] = constraints[i][j];

	for (int i=nums[3]+2*nums[0]; i<GetNum(); i++) m3[i][diff[i-nums[3]-2*nums[0]]] = 1.0;

	// rhs Vector

	for (int i=0; i<nums[3]+2*nums[0]; i++) v3[i]=0.0;
	for (int i=nums[3]+2*nums[0]; i<GetNum(); i++) v3[i] = sol[i-nums[3]-2*nums[0]];


	// solve, found remaining variables
	Vector<double> sol1 = Math::gauss(m3,v3,GetNum());

	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);

	return FBspCurv(nsol1,GetKnots(),GetOrd(),GetNum());
}

FBspCurv FBspCurv::ComputeDiffEquation2(double L, double C, Vector<int>& nums, Vector<double>& supports, Vector<double>& clamped, Vector<double>& ptforce, Vector<double>& disforce) const
{
	BspCurvBasisFuncSet b(GetKnots(),GetOrd(),GetOrd()+GetNum());
	Matrix<double> mat = Math::mmult(C,b.CreateMatrixMinimisation(1,0,L));
	Vector<double> vec1(GetNum(),0.0);
	int j = 0;
	for (int i=0; i<nums[1]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(ptforce[j+1],b(ptforce[j])));
		j++;
	}

	j = 0;
	Quadratic l1(0.0,4.0,0.25,0.0,1.0);
	Linear l2(0.0,4.0,0.5,0.0);

	for (int i=0; i<nums[2]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(disforce[j+2],b.CreateIntegralProduct(l1,disforce[j],disforce[j+1])));
		vec1 = Math::vsubtract(vec1,Math::vmult(disforce[j+2],b.CreateIntegralIntegralProduct(l2,disforce[j],disforce[j+1])));
		vec1 = Math::vadd(vec1,Math::vmult(0.5*disforce[j+2],b.CreateIntegralIntegralIntegral(disforce[j],disforce[j+1])));
		vec1 = Math::vmult(-1.0,vec1);
		j+=3;
	}


	// boundary conditions
	// supports

	Matrix<double> constraints(nums[3]+2*nums[0],GetNum()+1);
	for (int i=0; i<nums[3]; i++) {
		Vector<double> v = GetKnotSet().CreateVectorInterp(supports[i]);
		constraints.InsertRow(v,i);
	}


	// clamped
	j = 0;
	for (int i=0; i<nums[0]; i++) {
		Vector<double> v1 = GetKnotSet().CreateVectorInterp(clamped[i]);
		Vector<double> v2 = GetKnotSet().CreateVectorInterpDeriv(2,clamped[i]);

		constraints.InsertRow(v1,nums[3]+j);
		constraints.InsertRow(v2,nums[3]+j+1);
		j+=2;
	}


	Vector<double> rhs(nums[3]+2*nums[0],0.0);

	constraints.InsertCol(rhs,GetNum());

	Vector<int> v2 = Math::GetEliminateIndices(mat,constraints);	
	Matrix<double> m1 = Math::EliminateMVariables(mat,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(mat,constraints,vec1,v2);


	Vector<double> sol = Math::gauss(m1,v1,v1.GetNum());
	Vector<Point1D> nsol (sol.GetNum());

	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m3(GetNum(),GetNum(),0.0);
	Vector<double> v3(GetNum());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNum());
		
	for (int i=0; i<GetNum(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<nums[3]+2*nums[0]; i++)
		for (int j=0; j<GetNum(); j++)
			m3[i][j] = constraints[i][j];

	for (int i=nums[3]+2*nums[0]; i<GetNum(); i++) m3[i][diff[i-nums[3]-2*nums[0]]] = 1.0;

	// rhs Vector

	for (int i=0; i<nums[3]+2*nums[0]; i++) v3[i]=0.0;
	for (int i=nums[3]+2*nums[0]; i<GetNum(); i++) v3[i] = sol[i-nums[3]-2*nums[0]];


	// solve
	Vector<double> sol1 = Math::gauss(m3,v3,GetNum());

	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);

	return FBspCurv(nsol1,GetKnots(),GetOrd(),GetNum());
}

FBspCurv FBspCurv::ComputeDiffEquation3(double L, double C, Vector<int>& nums, Vector<double>& supports, Vector<double>& clamped, Vector<double>& ptforce, Vector<double>& disforce) const
{
	BspCurvBasisFuncSet b(GetKnots(),GetOrd(),GetOrd()+GetNum());
	Matrix<double> mat = Math::mmult(C,b.CreateMatrixMinimisation(2,0,L));
	Vector<double> vec1(GetNum(),0.0);
	int j = 0;
	for (int i=0; i<nums[1]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(ptforce[j+1],b(ptforce[j])));
		j++;
	}

	j = 0;
	

	// boundary conditions
	// supports

	j=0;
	Vector<double> v5(nums[3]+2*nums[0],0.0);
	Matrix<double> constraints(nums[3]+2*nums[0],GetNum()+1);
	for (int i=0; i<nums[3]; i++) {
		Vector<double> v = GetKnotSet().CreateVectorInterp(supports[j]);
		v5[i] = supports[j+1];
		constraints.InsertRow(v,i);
		j+=2;
	}

	
	// clamped
	j = 0;
	std::cerr << GetKnotSet();
	for (int i=0; i<nums[0]; i++) {
		std::cerr << clamped[i] << "\n";
		Vector<double> v1 = GetKnotSet().CreateVectorInterp(clamped[i]);
		Vector<double> v2 = GetKnotSet().CreateVectorInterpDeriv(2,clamped[i]);

		constraints.InsertRow(v1,nums[3]+j);
		constraints.InsertRow(v2,nums[3]+j+1);
		j+=2;
	}


	Vector<double> rhs(nums[3]+2*nums[0],0.0);

	constraints.InsertCol(v5,GetNum());


	Vector<int> v2 = Math::GetEliminateIndices(mat,constraints);	
	Matrix<double> m1 = Math::EliminateMVariables(mat,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(mat,constraints,vec1,v2);


	Vector<double> sol = Math::gauss(m1,v1,v1.GetNum());
	Vector<Point1D> nsol (sol.GetNum());

	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m3(GetNum(),GetNum(),0.0);
	Vector<double> v3(GetNum());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNum());
		
	for (int i=0; i<GetNum(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<nums[3]+2*nums[0]; i++)
		for (int j=0; j<GetNum(); j++)
			m3[i][j] = constraints[i][j];

	for (int i=nums[3]+2*nums[0]; i<GetNum(); i++) m3[i][diff[i-nums[3]-2*nums[0]]] = 1.0;

	// rhs Vector

	for (int i=0; i<nums[3]+2*nums[0]; i++) v3[i]=v5[i];//0.0;
	for (int i=nums[3]+2*nums[0]; i<GetNum(); i++) v3[i] = sol[i-nums[3]-2*nums[0]];


	// solve
	
	
	Vector<double> sol1 = Math::gauss(m3,v3,GetNum());

	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);

	return FBspCurv(nsol1,GetKnots(),GetOrd(),GetNum());
}
	
FBspCurv FBspCurv::ComputeFiniteElementNew(double L, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, Vector<Vector<double> >& curve, Vector<double>& bound) const
{
	BspCurvBasisFuncSet b1(GetKnots(), GetOrd(),GetOrd()+GetNum());
	BspCurvBasisFuncSet b2(GetKnotSet().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b3(GetKnotSet().CreateKnotSetDeriv(2));
	Matrix<double> mat1 = Math::mmult(L,b1.CreateMatrixMinimisation(2));
	Matrix<double> mat2 = b1.CreateMatrixMinimisation(1);
	Matrix<double> mat3 = b1.CreateMatrixMinimisation(0);
	

	Vector<double> vec1(GetNum(),0.0);

	int j = 0;
	for (int i=0; i<nums[0]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(ptforce[j+1],b1(ptforce[j])));
		j+=2;
	}

	j = 0;
	for (int i=0; i<nums[1]; i++) {
		BspCurv<double> bcurv = BspCurv<double>(PolyCurv<double>(curve[i],orders[i],GetKnots()[GetOrd()-1],GetKnots()[GetNum()]),GetKnotSet());
		vec1 =Math::vadd(vec1,b1.CreateIntegralNew(bcurv,disforce[j],disforce[j+1]));
		j+=2;
	}

	
	// boundary conditions
	// supports
	j=0;
	
	int row_count=0;

	int num_constraint = nums[2];
	if (nums[3] != 0) num_constraint += 1;
	if (nums[4] != 0) num_constraint += 1; 
	if (nums[5] != 0) num_constraint += 1;

	if (nums[6] != 0) num_constraint += 1;
	if (nums[7] != 0) num_constraint += 1;
	if (nums[8] != 0) num_constraint += 1;

	Matrix<double> constraints(num_constraint,GetNum()+1,0.0);

	Matrix<double> mi = Math::ComputeIdentityMatrix(GetNum());

	if (nums[3] != 0) {
		Vector<double> v1 = b1(0.0);
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[0];
		row_count++;
	}

	if (nums[4] != 0) {
		Vector<double> v1 = Math::mult4(b2(0.0),GetKnotSet().CreateMatrixDeriv());
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[1];
		row_count++;
	}
	
	if (nums[5] != 0) {
		Vector<double> v1 = Math::mult4(b3(0.0),GetKnotSet().CreateMatrixDeriv(2));
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[2];
		row_count++;
	}
	

	if (nums[6] != 0) {
		Vector<double> v1 = b1(GetKnots()[GetNum()]);
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[3];
		row_count++;
	}

	if (nums[7] != 0) {
		Vector<double> v1 = Math::mult4(b2(GetKnots()[GetNum()]),GetKnotSet().CreateMatrixDeriv());
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[4];
		row_count++;
	}
	
	if (nums[8] != 0) {
		Vector<double> v1 = Math::mult4(b3(GetKnots()[GetNum()]),GetKnotSet().CreateMatrixDeriv(2));
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[5];
		row_count++;
	}
	

	j=0;
	for (int i=0; i<nums[2]; i++) {
		double u = supports[j];
		double val = supports[j+1];
		Vector<double> v1 =	b1(u);
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = val;
		row_count++;
		j=j+2;
	}

	vec1 = Math::vmult(-1.0,vec1);


	Vector<int> v2 = Math::GetEliminateIndices(mat1,constraints);	
	Matrix<double> m10 = Math::EliminateMVariables(mat1,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(mat1,constraints,vec1,v2);

	Vector<double> sol = Math::gauss(m10,v1,v1.GetNum());
	Vector<Point1D> nsol(sol.GetNum());

	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m30(GetNum(),GetNum(),0.0);
	Vector<double> v3(GetNum());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNum());
		
	for (int i=0; i<GetNum(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<num_constraint; i++)
		for (int j=0; j<GetNum(); j++) 
			m30[i][j] = constraints[i][j];
		
	
	for (int i=num_constraint; i<GetNum(); i++) m30[i][diff[i-num_constraint]] = 1.0;


	// rhs Vector

	for (int i=0; i<num_constraint; i++) v3[i]= constraints[i][GetNum()];//0.0;

	for (int i=num_constraint; i<GetNum(); i++) v3[i] = sol[i-num_constraint];

	// solve
	Vector<double> sol1 = Math::gauss(m30,v3,GetNum());

	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);

	return FBspCurv(nsol1,GetKnots(),GetOrd(),GetNum());
}

FBspCurv FBspCurv::ComputeFiniteElementNew(double L, Vector<int>& nums, Vector<double>& supports, Vector<double>& ptforce, Vector<double>& disforce, Vector<int>& orders, BspCurv<double>& curve, Vector<double>& bound) const
{
	BspCurvBasisFuncSet b1(GetKnots(), GetOrd(),GetOrd()+GetNum());
	BspCurvBasisFuncSet b2(GetKnotSet().CreateKnotSetDeriv(1));
	BspCurvBasisFuncSet b3(GetKnotSet().CreateKnotSetDeriv(2));
	Matrix<double> mat1 = Math::mmult(L,b1.CreateMatrixMinimisation(2));
	Matrix<double> mat2 = b1.CreateMatrixMinimisation(1);
	Matrix<double> mat3 = b1.CreateMatrixMinimisation(0);
	

	Vector<double> vec1(GetNum(),0.0);

	int j = 0;
	for (int i=0; i<nums[0]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(ptforce[j+1],b1(ptforce[j])));
		j+=2;
	}

	j = 0;
	for (int i=0; i<nums[1]; i++) {
		vec1 =Math::vadd(vec1,b1.CreateIntegralNew(curve,disforce[j],disforce[j+1]));
		j+=2;
	}

	
	// boundary conditions
	// supports
	j=0;
	
	int row_count=0;

	int num_constraint = nums[2];
	if (nums[3] != 0) num_constraint += 1;
	if (nums[4] != 0) num_constraint += 1; 
	if (nums[5] != 0) num_constraint += 1;

	if (nums[6] != 0) num_constraint += 1;
	if (nums[7] != 0) num_constraint += 1;
	if (nums[8] != 0) num_constraint += 1;

	Matrix<double> constraints(num_constraint,GetNum()+1,0.0);

	Matrix<double> mi = Math::ComputeIdentityMatrix(GetNum());

	if (nums[3] != 0) {
		Vector<double> v1 = b1(0.0);
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[0];
		row_count++;
	}

	if (nums[4] != 0) {
		Vector<double> v1 = Math::mult4(b2(0.0),GetKnotSet().CreateMatrixDeriv());
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[1];
		row_count++;
	}
	
	if (nums[5] != 0) {
		Vector<double> v1 = Math::mult4(b3(0.0),GetKnotSet().CreateMatrixDeriv(2));
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[2];
		row_count++;
	}
	

	if (nums[6] != 0) {
		Vector<double> v1 = b1(GetKnots()[GetNum()]);
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[3];
		row_count++;
	}

	if (nums[7] != 0) {
		Vector<double> v1 = Math::mult4(b2(GetKnots()[GetNum()]),GetKnotSet().CreateMatrixDeriv());
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[4];
		row_count++;
	}
	
	if (nums[8] != 0) {
		Vector<double> v1 = Math::mult4(b3(GetKnots()[GetNum()]),GetKnotSet().CreateMatrixDeriv(2));
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = bound[5];
		row_count++;
	}
	

	j=0;
	for (int i=0; i<nums[2]; i++) {
		double u = supports[j];
		double val = supports[j+1];
		Vector<double> v1 =	b1(u);
		constraints.InsertRow(v1,row_count);
		constraints[row_count][GetNum()] = val;
		row_count++;
		j=j+2;
	}

	vec1 = Math::vmult(-1.0,vec1);


	Vector<int> v2 = Math::GetEliminateIndices(mat1,constraints);	
	Matrix<double> m10 = Math::EliminateMVariables(mat1,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(mat1,constraints,vec1,v2);

	Vector<double> sol = Math::gauss(m10,v1,v1.GetNum());
	Vector<Point1D> nsol(sol.GetNum());

	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m30(GetNum(),GetNum(),0.0);
	Vector<double> v3(GetNum());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNum());
		
	for (int i=0; i<GetNum(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<num_constraint; i++)
		for (int j=0; j<GetNum(); j++) 
			m30[i][j] = constraints[i][j];
		
	
	for (int i=num_constraint; i<GetNum(); i++) m30[i][diff[i-num_constraint]] = 1.0;


	// rhs Vector

	for (int i=0; i<num_constraint; i++) v3[i]= constraints[i][GetNum()];//0.0;

	for (int i=num_constraint; i<GetNum(); i++) v3[i] = sol[i-num_constraint];

	// solve
	Vector<double> sol1 = Math::gauss(m30,v3,GetNum());

	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);

	return FBspCurv(nsol1,GetKnots(),GetOrd(),GetNum());
}




FBspCurv FBspCurv::ComputeFiniteElementNew(double L, double C, Vector<int>& nums, Vector<double>& supports, Vector<double>& clamped, Vector<double>& ptforce, Vector<double>& disforce) const
{
	BspCurvBasisFuncSet b(GetKnots(),GetOrd(),GetOrd()+GetNum());

	Matrix<double> mat = Math::mmult(C,b.CreateMatrixMinimisation(2,0,L));
	std::cerr << mat;
	Vector<double> vec1(GetNum(),0.0);

	int j = 0;
	for (int i=0; i<nums[1]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(ptforce[j+1],b(ptforce[j])));
		j+=2;
	}

	j = 0;
	for (int i=0; i<nums[2]; i++) {
		vec1 = Math::vadd(vec1,Math::vmult(disforce[j+2],b.CreateIntegral(disforce[j],disforce[j+1])));
		j+=3;
	}


	// boundary conditions
	// supports

	j=0;
	Matrix<double> constraints(2*nums[3]+2*nums[0],GetNum()+1);
	for (int i=0; i<nums[3]; i++) {
		Vector<double> v1 = GetKnotSet().CreateVectorInterp(supports[i]);
		Vector<double> v2 = GetKnotSet().CreateVectorInterpDeriv(2,supports[i]);
		constraints.InsertRow(v1,j);
		constraints.InsertRow(v2,j+1);
		j+=2;
	}


	// clamped
	j = 0;
	for (int i=0; i<nums[0]; i++) {
		std::cerr << clamped[i] << "\n";
		Vector<double> v1 = GetKnotSet().CreateVectorInterp(clamped[i]);
		Vector<double> v2 = GetKnotSet().CreateVectorInterpDeriv(1,clamped[i]);

		constraints.InsertRow(v1,2*nums[3]+j);
		constraints.InsertRow(v2,2*nums[3]+j+1);
		j+=2;
	}


	Vector<double> rhs(2*nums[3]+2*nums[0],0.0);

	constraints.InsertCol(rhs,GetNum());


	Vector<int> v2 = Math::GetEliminateIndices(mat,constraints);	
	Matrix<double> m1 = Math::EliminateMVariables(mat,constraints,v2);
	Vector<double> v1 = Math::EliminateVVariables(mat,constraints,vec1,v2);


	Vector<double> sol = Math::gauss(m1,v1,v1.GetNum());
	Vector<Point1D> nsol(sol.GetNum());

	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);


	// find the values of the other variables

	Matrix<double> m3(GetNum(),GetNum(),0.0);
	Vector<double> v3(GetNum());

	std::set<int> set1(v2.begin(),v2.end());
	Vector<int> v4(GetNum());
		
	for (int i=0; i<GetNum(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());


	for (int i=0; i<2*nums[3]+2*nums[0]; i++)
		for (int j=0; j<GetNum(); j++)
			m3[i][j] = constraints[i][j];

	for (int i=2*nums[3]+2*nums[0]; i<GetNum(); i++) m3[i][diff[i-2*nums[3]-2*nums[0]]] = 1.0;

	// rhs Vector

	for (int i=0; i<2*nums[3]+2*nums[0]; i++) v3[i]=0.0;
	for (int i=2*nums[3]+2*nums[0]; i<GetNum(); i++) v3[i] = sol[i-2*nums[3]-2*nums[0]];


	// solve
	Vector<double> sol1 = Math::gauss(m3,v3,GetNum());


	// write away
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);

	return FBspCurv(nsol1,GetKnots(),GetOrd(),GetNum());
}


FBspCurv FBspCurv::ComputeFiniteElement(double C, double L, double Q, double QL, double F, double FL) const
{
	BspCurvBasisFuncSet b(GetKnots(),GetOrd(),GetOrd()+GetNum());

	Matrix<double> mat = b.CreateMatrixMinimisation(2,0,L);

	Vector<double> vec1 = b.CreateIntegral(QL,L);

	Vector<double> vec2 = b(FL);

	mat = Math::mmult(0.5*C,mat);
	vec1 = Math::vmult(Q,vec1);
	vec2 = Math::vmult(F,vec2);

	// boundary conditions

	Vector<double> b1 = GetKnotSet().CreateVectorInterp(0.0);
	Vector<double> b2 = GetKnotSet().CreateVectorInterp(2.0*L/3.0);
	Vector<double> b3 = GetKnotSet().CreateVectorInterp(L);
	Vector<double> b4 = GetKnotSet().CreateVectorInterpDeriv(1,0.0);

	Matrix<double> m(4,GetNum()+1);

	Vector<double> v(4,0.0);

	m.InsertRow(b1,0);
	m.InsertRow(b4,1);
	m.InsertRow(b2,2);
	m.InsertRow(b3,3);
	m.InsertCol(v,GetNum());


	Vector<int> v2 = Math::GetEliminateIndices(mat,m);
	
	Matrix<double> m1 = Math::EliminateMVariables(mat,m,v2);
	
	Vector<double> v1 = Math::EliminateVVariables(mat,m,Math::vadd(vec1,vec2),v2);

	Vector<double> sol = Math::gauss(m1,v1,v1.GetNum());
	Vector<Point1D> nsol (sol.GetNum());

	for (int i=0; i<sol.GetNum(); i++) nsol[i] = Point1D(sol[i]);

	Matrix<double> m3(GetNum(),GetNum(),0.0);

	Vector<double> v3(GetNum());

	std::set<int> set1(v2.begin(),v2.end());
		
	Vector<int> v4(GetNum());
		
	for (int i=0; i<GetNum(); i++) v4[i] = i;

	std::set<int> set2(v4.begin(),v4.end());
	std::set<int> set3;

	std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
	Vector<int> diff(set3.size());
	std::copy(set3.begin(),set3.end(),diff.begin());

	m3.InsertRow(b1,0);
	m3.InsertRow(b4,1);
	m3.InsertRow(b2,2);
	m3.InsertRow(b3,3);

	for (int i=4; i<GetNum(); i++) m3[i][diff[i-4]] = 1.0;


	for (int i=0; i<4; i++) v3[i]=0.0;
	for (int i=4; i<GetNum(); i++) v3[i] = sol[i-4];

	std::cerr << m3;
	std::cerr << v3;


	Vector<double> sol1 = Math::gauss(m3,v3,GetNum());

	std::cerr << sol1;
	Vector<Point1D> nsol1 (sol1.GetNum());

	for (int i=0; i<sol1.GetNum(); i++) nsol1[i] = Point1D(sol1[i]);

	return FBspCurv(nsol1,GetKnots(),GetOrd(),GetNum());
}





PBspCurv2D::PBspCurv2D() { }
PBspCurv2D::PBspCurv2D(const BspCurv<Point2D>& b) : BspCurv<Point2D>(b) { }
PBspCurv2D::PBspCurv2D(const Vector<Point2D>& Cpts, const Vector<double>& Kts, int Ord, int Num) :
BspCurv<Point2D>(Cpts, Kts, Ord, Num) {}



PBspCurv2D::PBspCurv2D(const Vector<Point2D>& Cpts, int Ord, int Num) : BspCurv<Point2D>(Cpts, Ord, Num) {}
PBspCurv2D::PBspCurv2D(const FBspCurv& b)
{
	Vector<Point2D> v1(b.GetNum());
	Vector<double> v2(b.GetKnotSet().ComputeKnotSetAver());

	for (int i=0; i<b.GetNum(); i++) 
		v1[i] = Point2D(v2[i],b.GetCPoints()[i].GetX());

	*this = PBspCurv2D(v1,b.GetKnots(),b.GetOrd(),b.GetNum());
}




PBspCurv3D::PBspCurv3D() { }
PBspCurv3D::PBspCurv3D(const BspCurv<Point3D>& b) : BspCurv<Point3D>(b) { }
PBspCurv3D::PBspCurv3D(const Vector<Point3D>& Cpts, const Vector<double>& Kts, int Ord, int Num) :
BspCurv<Point3D>(Cpts, Kts, Ord, Num) {}


PBspCurv3D::PBspCurv3D(const Vector<Point3D>& Cpts, int Ord, int Num) : BspCurv<Point3D>(Cpts, Ord, Num) {}
PBspCurv3D::PBspCurv3D(const FBspCurv& b)
{
	Vector<Point3D> v1(b.GetNum());
	Vector<double> v2(b.GetKnotSet().ComputeKnotSetAver());

	for (int i=0; i<b.GetNum(); i++) 
		v1[i] = Point3D(v2[i],b.GetCPoints()[i].GetX(),0.0);

	*this = PBspCurv3D(v1,b.GetKnots(),b.GetOrd(),b.GetNum());
}

PBspCurv3D::PBspCurv3D(const PBspCurv2D& b)
{
	Vector<Point3D> v(b.GetNum());

	for (int i=0; i<b.GetNum(); i++) 
		v[i] = Point3D(b.GetCPoints()[i]);

	*this = PBspCurv3D(v,b.GetKnots(),b.GetOrd(),b.GetNum());
}



