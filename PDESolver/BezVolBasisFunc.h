
#ifndef BEZVOLBASISFUNC
#define BEZVOLBASISFUNC


//#include "headers.h"

/* template<class T>
class BezVol;*/


#include "object.h"
#include "BezVol.h"
//#include "KnotSet.h"
#include "Matrix3D.h"

class BezVolBasisFunc : public Vol<double> {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int ordw;	// order in w
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<BezVol<double> > b;	// BezVol representation of basis function
	
	// private functions
	BezVol<double> CreateBezVol() const;	// creates the BezVol
	virtual ObjectID Identity() const { return std::string("BezVolBasisFunc"); }
public:

	// constructors
	BezVolBasisFunc();
	BezVolBasisFunc(int Ordu, int Ordv, int Ordw, const Vector<double>& Ktsu, const Vector<double>& Ktsv,
			const Vector<double>& Ktsw);
	BezVolBasisFunc(BuildFlag) { }
	static Ptr<FileObject> Build();

	// evaluators
	double Eval(double u, double v, double w) const;	// evaluate basis function
	virtual double operator()(double u, double v, double w) const;
	virtual double operator()(int, int, int, double u, double v, double w) const;

	// access functions
	BezVol<double> GetBezVol() const;		// get BezVol representation
	int ComputeDimU() const;
	int ComputeDimV() const;
	int ComputeDimW() const;

	
	virtual double Derive(int, int , int, double, double, double) const;
	virtual double GetLeftLimitU() const;
	virtual double GetRightLimitU() const;
	virtual double GetLeftLimitV() const;
	virtual double GetRightLimitV() const;
	virtual double GetLeftLimitW() const;
	virtual double GetRightLimitW() const;
	
	// read and write

	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void * Read(bifstream &is);
	virtual void Write(bofstream &os);
};


class BezVolBasisFuncSet : public TextObject, public FTextObject, public FileObject {

	// data
	int ordu;	 // order of basis function in u
	int ordv;	// order of basis function in v
	int ordw;
	Ptr<Vector<double> > ktsu;		// knots in u
	Ptr<Vector<double> > ktsv;		// knots in v
	Ptr<Vector<double> > ktsw;
	Ptr<Matrix3D<BezVolBasisFunc> > b;
	virtual ObjectID Identity() const { return std::string("BezVolBasisFuncSet"); }
public:

	// constructors
	BezVolBasisFuncSet();
	BezVolBasisFuncSet(int Ordu, int Ordv, int Ordw, 
		const Vector<double>& Ktsu, const Vector<double>& Ktsv, const Vector<double>& Ktsw);
	BezVolBasisFuncSet(BuildFlag) { }
	static Ptr<FileObject> Build();

	// evaluators
	Matrix3D<double> Eval(double u, double v, double w) const;	// evaluate basis function
	Matrix3D<double> operator()(double u, double v, double w) const;
	Matrix3D<double> operator()(int, int, int, double u, double v, double w) const;
	Matrix3D<double> Derive(int, int, int, double u, double v, double w) const;
	
	// read and write
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs);
	virtual void * Read(bifstream &is);
	virtual void Write(bofstream &os);
};

#endif