#ifndef CURVES
#define CURVES

#include "object.h"

class SinePi : public Curve<double> {

	double x1, x2;
	virtual ObjectID Identity() const { return std::string("class SinePi"); }
public:
	SinePi(double X1, double X2);
	virtual double operator()(double) const;
	virtual double operator()(int, double) const;
	virtual double Derive(int, double) const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs) ;
	virtual void write(std::ostream& os) ;
	virtual void read(std::istream &is);
};




class Sine : public Curve<double> {

	double x1, x2;
	virtual ObjectID Identity() const { return std::string("class Sine"); }
public:
	Sine(double X1, double X2);
	virtual	double operator()(double) const;
	virtual double operator()(int, double) const;
	virtual double Derive(int, double) const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs) ;
	virtual void write(std::ostream& os) ;
	virtual void read(std::istream &is);
};



// cosine function class

class Cosine : public Curve<double> {

	double x1, x2;
	virtual ObjectID Identity() const { return std::string("class Cosine"); }
public:
	Cosine(double X1, double X2);
	virtual double operator()(double) const;
	virtual double operator()(int, double) const;
	virtual double Derive(int, double) const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs) ;
	virtual void write(std::ostream& os) ;
	virtual void read(std::istream &is);
};




// quadratic function class


class Quadratic : public Curve<double> {

	double x1, x2;
	double a, b, c;
	virtual ObjectID Identity() const { return std::string("class Quadratic"); }
public:
	Quadratic(double X1, double X2, double A, double B, double C);
	virtual double operator()(double) const;
	virtual double operator()(int, double) const;
	virtual double Derive(int, double) const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs) ;
	virtual void write(std::ostream& os) ;
	virtual void read(std::istream &is);
};


class Linear : public Curve<double> {

	double x1, x2;
	double a, b;
	virtual ObjectID Identity() const { return std::string("class Linear"); }
public:
	Linear(double X1, double X2, double A, double B);
	virtual double operator()(double) const;
	virtual double operator()(int, double) const;
	virtual double Derive(int, double) const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs) ;
	virtual void write(std::ostream& os) ;
	virtual void read(std::istream &is);
};


class Circle : public Curve<Point3D> {

	double x1, x2;
	double r;
	virtual ObjectID Identity() const { return std::string("class Circle"); }
public:
	Circle(double X1, double X2, double R);
	virtual Point3D operator()(double x) const;
	virtual Point3D operator()(int, double) const;
	virtual Point3D Derive(int, double) const;
	virtual double GetLeftLimit() const;
	virtual double GetRightLimit() const;
	virtual double GetRadius() const;
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs)  ;
	virtual void write(std::ostream& os);
	virtual void read(std::istream &is);
};



#endif