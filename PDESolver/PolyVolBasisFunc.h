
#ifndef POLYVOLBASISFUNC
#define POLYVOLBASISFUNC

#include "PolyVol.h"
//#include "matrix3d.h"*/

//#include "headers.h"


class PolyVolBasisFunc : public Vol<double> {

	// data
	int ordu;	 // order of basis function
	int ordv;	// order in v
	int ordw;
	int indexu;
	int indexv;
	int indexw;
	double leftlimitu;  // ranges in u
	double rightlimitu;
	double leftlimitv; // rnages in v
	double rightlimitv;
	double leftlimitw;
	double rightlimitw;
	Ptr<PolyVol<double> > b;	// PolyVol representation of basis function

	// private functions
	PolyVol<double> CreatePolyVol() const;	// creates the PolyVol
	virtual ObjectID Identity() const { return std::string("PolyVolBasisFunc"); }
public:

	// constructors
	PolyVolBasisFunc();
	PolyVolBasisFunc(int Ordu, int Ordv, int Ordw, int Indexu, int Indexv, int Indexw, double LeftLimitU=0.0, double RightLimitU=1.0, 
		double LeftLimitV=0.0, double RightLimitV=1.0, double LeftLimitW=0.0, double RightLimitW=1.0);
	PolyVolBasisFunc(BuildFlag) { }
	static Ptr<FileObject> Build();
	

	// evaluators
	double Eval(double u, double v, double w) const;
	virtual double operator()(double u, double v, double w) const;
	virtual double operator()(int, int, int, double u, double v, double w) const;
	// access functions
	PolyVol<double> GetPolyVol() const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
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


class PolyVolBasisFuncSet : public TextObject, public FTextObject, public FileObject {
	// data
	int ordu;	 // order of basis function
	int ordv;
	int ordw;
	double leftlimitu;
	double rightlimitu;
	double leftlimitv;
	double rightlimitv;
	double leftlimitw;
	double rightlimitw;

	Ptr<Matrix3D<PolyVolBasisFunc> > b;	// PolyVol representation of basis function
	virtual ObjectID Identity() const { return std::string("PolyVolBasisFuncSet"); }
public:
	// constructors
	PolyVolBasisFuncSet();
	PolyVolBasisFuncSet(int Ordu, int Ordv, int Ordw, double LU, double RU, double LV, double RV, double LW, double RW);
	PolyVolBasisFuncSet(BuildFlag) { }
	static Ptr<FileObject> Build();

	// evaluators
	Matrix3D<double> Eval(double u, double v, double w) const;
	Matrix3D<double> operator()(double u, double v, double w) const;

	// read and write
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void * Read(bifstream &is);
	virtual void Write(bofstream &os);
};

#endif