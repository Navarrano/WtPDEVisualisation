
#ifndef BSPVOLBASISFUNC
#define BSPVOLBASISFUNC

/* template<class T>
class BspVol;*/

#include "object.h"
#include "BspVol.h"
//#include "KnotSet.h"
#include "Matrix3d.h"
//#include "Vector.h"

//#include "headers.h"


class BspVolBasisFunc : public Vol<double> {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int ordw;	// order in w
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<BspVol<double> > b;	// BspVol representation of basis function

	// private functions
	BspVol<double> CreateBspVol() const;	// creates the BspVol
	virtual ObjectID Identity() const { return std::string("BspVolBasisFunc"); }
public:
	// constructors
	BspVolBasisFunc();
	BspVolBasisFunc(int Ordu, int Ordv, int Ordw, const Vector<double>& Ktsu, 
		const Vector<double>& Ktsv, const Vector<double>& Ktsw);
	BspVolBasisFunc(BuildFlag) { }
	static Ptr<FileObject> Build();

	// access functions
	int ComputeDimU() const;
	int ComputeDimV() const;
	int ComputeDimW() const;
	BspVol<double> GetBspVol() const;		// get BspVol representation

	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	Vector<double> GetKnotsW() const;

	// evaluators
	double Eval(double u, double v, double w) const;	// evaluate basis function
	virtual double operator()(double u, double v, double w) const;
	virtual double operator()(int, int, int, double u, double v, double w) const;
	
	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void * Read(bifstream &is);
	virtual void Write(bofstream &os);
	
	virtual double Derive(int, int , int, double, double, double) const;
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
	virtual double GetLeftLimitW() const;
	virtual double GetRightLimitW() const;
};


class BspVolBasisFuncSet : public TextObject, public FTextObject, public FileObject {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int ordw;
	int numu;
	int numv;
	int numw;
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<Matrix3D<BspVolBasisFunc> > b;
	virtual ObjectID Identity() const { return std::string("BspVolBasisFuncSet"); }
public:
	// constructors
	BspVolBasisFuncSet();
	BspVolBasisFuncSet(int Ordu, int Ordv, int Ordw, int Numu, int Numv, int Numw, const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw);
	BspVolBasisFuncSet(BuildFlag) { }
	static Ptr<FileObject> Build();

	// evaluators
	Matrix3D<double> Eval(double u, double v, double w) const;	// evaluate basis function
	Matrix3D<double> operator()(double u, double v, double w) const;
	Matrix3D<double> operator()(int, int, int, double u, double v, double w) const;
	Matrix3D<double> Derive(int, int, int, double u, double v, double w) const;
	
	Matrix<double> ComputeUBasisMatrix() const;
	Matrix<double> ComputeVBasisMatrix() const;
	Matrix<double> ComputeWBasisMatrix() const;
	Matrix3D<double> CreateIntegral(double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateIntegralNewU(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateIntegralNewV(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateIntegralNewW(const BspVol<double>& s, const BspVol<double>& norm, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<Matrix3D<double> > CreateMatrixIntegral(int levu, int levv, int levw) const;
	Matrix3D<Matrix3D<double> > CreateMatrixIntegral(int levu, int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisation1(int levu, int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisation1(int levu, int levv, int levw) const;
	Matrix3D<double> CreateMatrixMinimisation(int levu, int levv, int levw) const;
	Matrix3D<double> CreateMatrixMinimisation(int levu, int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationUV(int levu, int levv, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationUW(int levu, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationVW(int levv, int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationU(int levu, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationV(int levv, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationV1(int levv, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix3D<double> CreateMatrixMinimisationW(int levw, double u1, double u2, double v1, double v2, double w1, double w2) const;
	Matrix<double> CreateMatrixKnotAveragesU() const;
	Matrix<double> CreateMatrixKnotAveragesV() const;
	Matrix<double> CreateMatrixKnotAveragesW() const;
	Vector<double> GetKnotsU() const;
	Vector<double> GetKnotsV() const;
	Vector<double> GetKnotsW() const;

	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void * Read(bifstream &is);
	virtual void Write(bofstream &os);
};


#endif