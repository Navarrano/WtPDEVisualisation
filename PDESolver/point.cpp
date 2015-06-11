
#include "point.h"

Point1D::Point1D() : x(0.0) {}
Point1D::Point1D(double X) : x(X) {}
Point1D::Point1D(const Point2D& pt) : x(pt.GetX()) {}
Point1D::Point1D(const Point3D& pt) : x(pt.GetX()) {}
Point1D::Point1D(const Vector<double>& v) : x(v[0]) {}
Point1D::Point1D(const FVector<double>& v) : x(v[v.GetMin()]) {}
inline double Point1D::GetX() const { return x; }

Point1D& Point1D::operator=(const Point1D& p)
{
	x = p.x;
	return *this;
}

Point1D& Point1D::operator=(double d)
{
	x = d;
	return *this;
}


Point1D operator*(const Point1D& p1, const Point1D& p2) {
	return Point1D(p1.GetX()*p2.GetX());
}

Point1D operator*(double d, const Point1D& p) {
	return Point1D(d*p.GetX());
}

Point1D operator*(const Point1D& p, double d) {
	return Point1D(d*p.GetX());
}


Point1D operator+(double d, const Point1D& p) {
	return Point1D(d+p.GetX());
}

Point1D operator+(const Point1D& p, double d) {
	return Point1D(d+p.GetX());
}

Point1D operator+(const Point1D& p1, const Point1D& p2) {
	return Point1D(p1.GetX()+p2.GetX());
}

Point1D operator-(double d, const Point1D& p) {
	return Point1D(d-p.GetX());
}

Point1D operator-(const Point1D& p, double d) {
	return Point1D(p.GetX()-d);
}

Point1D operator-(const Point1D& p1, const Point1D&p2) {
	return Point1D(p1.GetX()-p2.GetX());
}

Point1D operator/(const Point1D& p, double d) {
	return Point1D(p.GetX()/d);
}

bool operator>(const Point1D& p, double d) {
	if (fabs(p.x) > d) return true;
	else return false;
}

bool operator>(const Point1D& p1, const Point1D& p2) {
	if (p1.x > p2.x) return true;
	else return false;
}

Point1D::operator double() const
{
	return x;
}

Point1D::operator FVector<double>() const
{
	FVector<double> vec(1,3,0.0);
	vec[1] = x;
	vec[2] = 0.0;
	vec[3] = 0.0;
	return vec;
}

void Point1D::write(std::ostream& os)
{
	os << GetX();
}



void Point1D::read(std::istream& is)
{
	std::cout << "input x\n";
	is >> x;
}

void Point1D::writefile(std::ofstream& ofs)
{
	ofs << GetX();
}

void Point1D::readfile(std::ifstream& ifs)
{
	ifs >> x;
}



Point2D::Point2D() : x(0.0), y(0.0) {}
Point2D::Point2D(double d) : x(d), y(d) {}
Point2D::Point2D(double X, double Y) : x(X), y(Y) {}
Point2D::Point2D(const Point1D& pt) : x(pt.GetX()), y(0.0) {}
Point2D::Point2D(const Point3D& pt) : x(pt.GetX()), y(pt.GetY()) {}
Point2D::Point2D(const Vector<double>& v) : x(v[0]), y(v[1]) {}
Point2D::Point2D(const FVector<double>& v) : x(v[v.GetMin()]), y(v[v.GetMin()+1]) {}
inline double Point2D::GetX() const { return x; }
inline double Point2D::GetY() const { return y; }

Point2D& Point2D::operator=(const Point2D& p)
{
	x = p.x;
	y = p.y;
	return *this;
}

Point2D& Point2D::operator=(const Point3D& p)
{
	x = p.GetX();
	y = p.GetY();
	return *this;
}


Point2D& Point2D::operator=(double d)
{
	x = d;
	y = d;
	return *this;
}


Point2D operator*(const Point2D& p1, const Point2D& p2) {
	return Point2D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetY());
}

Point2D operator*(const Point2D& p1, const Point1D& p2) {
	return Point2D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetX());
}

Point2D operator*(double d, const Point2D& p) {
	return Point2D(d*p.GetX(),d*p.GetY());
}

Point2D operator*(const Point2D& p, double d) {
	return Point2D(d*p.GetX(),d*p.GetY());
}


Point2D operator+(double d, const Point2D& p) {
	return Point2D(d+p.GetX(),d+p.GetY());
}

Point2D operator+(const Point2D& p, double d) {
	return Point2D(d+p.GetX(),d+p.GetY());
}

Point2D operator+(const Point2D& p1, const Point2D& p2) {
	return Point2D(p1.GetX()+p2.GetX(),p1.GetY()+p2.GetY());
}

Point2D operator-(double d, const Point2D& p) {
	return Point2D(d-p.GetX(),d-p.GetY());
}

Point2D operator-(const Point2D& p, double d) {
	return Point2D(p.GetX()-d,p.GetY()-d);
}

Point2D operator-(const Point2D& p1, const Point2D&p2) {
	return Point2D(p1.GetX()-p2.GetX(),p1.GetY()-p2.GetY());
}

Point2D operator/(const Point2D& p, double d) {
	return Point2D(p.GetX()/d,p.GetY()/d);
}

bool operator>(const Point2D& p, double d) {
	if (fabs(p.x) > d || fabs(p.y) > d) return true;
	else return false;
}

bool operator>(const Point2D& p1, const Point2D& p2) {
	if (sqrt(p1.x*p1.x+p1.y*p1.y+p1) > sqrt(p2.x*p2.x+p2.y*p2.y)) return true;
	else return false;
}



Point2D::operator double() const
{
	return sqrt(x*x+y*y);
}

Point2D::operator FVector<double>() const
{
	FVector<double> vec(1,3,0.0);
	vec[1] = x;
	vec[2] = y;
	vec[3] = 0.0;
	return vec;
}


void Point2D::write(std::ostream& os)
{
	os << GetX() << GetY();
}



void Point2D::read(std::istream& is)
{
	std::cout << "input x, y\n";
	is >> x >> y;
}

void Point2D::writefile(std::ofstream& ofs)
{
	ofs << GetX() << " " << GetY() << "\n";
}

void Point2D::readfile(std::ifstream& ifs)
{
	ifs >> x >> y;
}




Point3D::Point3D() : x(0.0), y(0.0), z(0.0) {}
Point3D::Point3D(double d) : x(d), y(d), z(d) {}
Point3D::Point3D(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
Point3D::Point3D(const FVector<double>& v) : x(v[v.GetMin()]), y(v[v.GetMin()+1]), z(v[v.GetMin()+2]) {}
Point3D::Point3D(const Vector<double>& v) : x(v[0]), y(v[1]), z(v[2]) {}
Point3D::Point3D(const Point1D& pt) : x(pt.GetX()), y(0.0), z(0.0) {}
Point3D::Point3D(const Point2D& pt) : x(pt.GetX()), y(pt.GetY()), z(0.0) {}
Point3D::Point3D(const Point4D& pt) : x(pt.GetX()/pt.GetW()), y(pt.GetY()/pt.GetW()), z(pt.GetZ()/pt.GetW()) {}
double Point3D::GetX() const { return x; }
double Point3D::GetY() const { return y; }
double Point3D::GetZ() const { return z; }

Point3D& Point3D::operator=(const Point3D& p)
{
	x = p.x;
	y = p.y;
	z = p.z;
	return *this;
}

Point3D& Point3D::operator=(const FVector<double>& v)
{
	x=v[v.GetMin()];
	y=v[v.GetMin()+1];
	z=v[v.GetMin()+2];
	return *this;
}

Point3D& Point3D::operator=(const Point4D& p)
{
	x = p.GetX()/p.GetW();
	y = p.GetY()/p.GetW();
	z = p.GetZ()/p.GetW();
	return *this;
}


Point3D& Point3D::operator=(double d)
{
	x = d;
	y = d;
	z = d;
	return *this;
}


Point3D operator*(const Point3D& p1, const Point3D& p2) {
	return Point3D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetY(),p1.GetZ()*p2.GetZ());
}

Point3D operator*(const Point3D& p1, const Point2D& p2) {
	return Point3D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetY(),0.0);
}


Point3D operator*(const Point3D& p1, const Point1D& p2) {
	return Point3D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetX(),p1.GetZ()*p2.GetX());
}

Point3D operator*(double d, const Point3D& p) {
	return Point3D(d*p.GetX(),d*p.GetY(),d*p.GetZ());
}

Point3D operator*(const Point3D& p, double d) {
	return Point3D(d*p.GetX(),d*p.GetY(),d*p.GetZ());
}


Point3D operator+(double d, const Point3D& p) {
	return Point3D(d+p.GetX(),d+p.GetY(),d+p.GetZ());
}

Point3D operator+(const Point3D& p, double d) {
	return Point3D(d+p.GetX(),d+p.GetY(),d+p.GetZ());
}

Point3D operator+(const Point3D& p1, const Point3D& p2) {
	return Point3D(p1.GetX()+p2.GetX(),p1.GetY()+p2.GetY(),p1.GetZ()+p2.GetZ());
}

Point3D operator+=(const Point3D& p1, const Point3D& p2) {
	return Point3D(p1.GetX()+p2.GetX(),p1.GetY()+p2.GetY(),p1.GetZ()+p2.GetZ());
}


Point3D operator-(double d, const Point3D& p) {
	return Point3D(d-p.GetX(),d-p.GetY(),d-p.GetZ());
}

Point3D Point3D::operator-() {
	return Point3D(-x,-y,-z);
}

Point3D operator-(const Point3D& p, double d) {
	return Point3D(p.GetX()-d,p.GetY()-d,p.GetZ()-d);
}

Point3D operator-(const Point3D& p1, const Point3D&p2) {
	return Point3D(p1.GetX()-p2.GetX(),p1.GetY()-p2.GetY(),p1.GetZ()-p2.GetZ());
}

Point3D operator-=(const Point3D& p1, const Point3D&p2) {
	return Point3D(p1.GetX()-p2.GetX(),p1.GetY()-p2.GetY(),p1.GetZ()-p2.GetZ());
}

Point3D operator/(const Point3D& p, double d) {
	return Point3D(p.GetX()/d,p.GetY()/d,p.GetZ()/d);
}


bool operator>(const Point3D& p, double d) {
	if (fabs(p.x) > d || fabs(p.y) > d || fabs(p.z) > d) return true;
	else return false;
}

bool operator>(const Point3D& p1, const Point3D& p2) {
	if (sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z) > sqrt(p2.x*p2.x+p2.y*p2.y+p2.z*p2.z)) return true;
	else return false;
}

double Point3D::Product(const Point3D& p) const
{
	return (x*p.GetX()+y*p.GetY()+z*p.GetZ());
}

Point3D::operator double() const
{
	return sqrt(x*x+y*y+z*z);
}

double Point3D::Length() const
{
	return sqrt(x*x+y*y+z*z);
}

Point3D::operator FVector<double>() const
{
	FVector<double> vec(1,3,0.0);
	vec[1] = x;
	vec[2] = y;
	vec[3] = z;
	return vec;
}

void Point3D::write(std::ostream& os)
{
	os << GetX() << " " << GetY() << " " << GetZ() << "\n";
}



void Point3D::read(std::istream& is)
{
	std::cout << "input x, y, z\n";
	is >> x >> y >> z;
}

void Point3D::writefile(std::ofstream& ofs)
{
	ofs << GetX() << " " << GetY() << " " << GetZ() << "\n"; 
}

void Point3D::readfile(std::ifstream& ifs)
{
	ifs >> x >> y >> z;
}



Point4D::Point4D() : x(0.0), y(0.0), z(0.0), w(1.0) {}
Point4D::Point4D(double d) : x(d), y(d), z(d), w(d) {}
Point4D::Point4D(double X, double Y, double Z, double W) : x(X), y(Y), z(Z), w(W) {}
Point4D::Point4D(const FVector<double>& v) : x(v[v.GetMin()]), y(v[v.GetMin()+1]), z(v[v.GetMin()+2]), w(v[v.GetMin()+3]) {}
Point4D::Point4D(const Vector<double>& v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}
Point4D::Point4D(const Point1D& pt) : x(pt.GetX()), y(0.0), z(0.0), w(1.0) {}
Point4D::Point4D(const Point2D& pt) : x(pt.GetX()), y(pt.GetY()), z(0.0), w(1.0) {}
Point4D::Point4D(const Point3D& pt) : x(pt.GetX()), y(pt.GetY()), z(pt.GetZ()), w(1.0) {}
double Point4D::GetX() const { return x; }
double Point4D::GetY() const { return y; }
double Point4D::GetZ() const { return z; }
double Point4D::GetW() const { return w; }

Point4D& Point4D::operator=(const Point4D& p)
{
	x = p.x;
	y = p.y;
	z = p.z;
	w = p.w;
	return *this;
}

Point4D& Point4D::operator=(const Point3D& p)
{
	x = p.GetX();
	y = p.GetY();
	z = p.GetZ();
	w = 1.0;
	return *this;
}

Point4D& Point4D::operator=(double d)
{
	x = d;
	y = d;
	z = d;
	w = d;
	return *this;
}


Point4D operator*(const Point4D& p1, const Point4D& p2) {
	return Point4D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetY(),p1.GetZ()*p2.GetZ(),p1.GetW()*p2.GetW());
}

Point4D operator*(const Point4D& p1, const Point3D& p2) {
	return Point4D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetY(),p1.GetZ()*p2.GetZ(),0.0);
}


Point4D operator*(const Point4D& p1, const Point2D& p2) {
	return Point4D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetY(),0.0,0.0);
}


Point4D operator*(const Point4D& p1, const Point1D& p2) {
	return Point4D(p1.GetX()*p2.GetX(),p1.GetY()*p2.GetX(),p1.GetZ()*p2.GetX(),p1.GetW()*p2.GetX());
}

Point4D operator*(double d, const Point4D& p) {
	return Point4D(d*p.GetX(),d*p.GetY(),d*p.GetZ(),d*p.GetW());
}

Point4D operator*(const Point4D& p, double d) {
	return Point4D(d*p.GetX(),d*p.GetY(),d*p.GetZ(),d*p.GetW());
}


Point4D operator+(double d, const Point4D& p) {
	return Point4D(d+p.GetX(),d+p.GetY(),d+p.GetZ(),d+p.GetW());
}

Point4D operator+(const Point4D& p, double d) {
	return Point4D(d+p.GetX(),d+p.GetY(),d+p.GetZ(),d+p.GetW());
}

Point4D operator+(const Point4D& p1, const Point4D& p2) {
	return Point4D(p1.GetX()+p2.GetX(),p1.GetY()+p2.GetY(),p1.GetZ()+p2.GetZ(),p1.GetW()+p2.GetW());
}

Point4D operator+=(const Point4D& p1, const Point4D& p2) {
	return Point4D(p1.GetX()+p2.GetX(),p1.GetY()+p2.GetY(),p1.GetZ()+p2.GetZ(),p1.GetW()+p2.GetW());
}


Point4D operator-(double d, const Point4D& p) {
	return Point4D(d-p.GetX(),d-p.GetY(),d-p.GetZ(),d-p.GetW());
}

Point4D operator-(const Point4D& p, double d) {
	return Point4D(p.GetX()-d,p.GetY()-d,p.GetZ()-d,p.GetW()-d);
}

Point4D operator-(const Point4D& p1, const Point4D&p2) {
	return Point4D(p1.GetX()-p2.GetX(),p1.GetY()-p2.GetY(),p1.GetZ()-p2.GetZ(),p1.GetW()-p2.GetW());
}

Point4D operator-=(const Point4D& p1, const Point4D&p2) {
	return Point4D(p1.GetX()-p2.GetX(),p1.GetY()-p2.GetY(),p1.GetZ()-p2.GetZ(),p1.GetW()-p2.GetW());
}

Point4D operator/(const Point4D& p, double d) {
	return Point4D(p.GetX()/d,p.GetY()/d,p.GetZ()/d,p.GetW()/d);
}


bool operator>(const Point4D& p, double d) {
	if (fabs(p.x) > d || fabs(p.y) > d || fabs(p.z) > d || fabs(p.w) > d) return true;
	else return false;
}

bool operator>(const Point4D& p1, const Point4D& p2) {
	if (sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z+p1.w*p1.w) > sqrt(p2.x*p2.x+p2.y*p2.y+p2.z*p2.z+p2.w*p2.w)) return true;
	else return false;
}

double Point4D::Product(const Point4D& p) const
{
	return (x*p.GetX()+y*p.GetY()+z*p.GetZ()+w*p.GetW());
}

Point4D::operator double() const
{
	return sqrt(x*x+y*y+z*z+w*w);
}

Point4D::operator FVector<double>() const
{
	FVector<double> vec(0,3,0.0);
	vec[0] = x;
	vec[1] = y;
	vec[2] = z;
	vec[3] = w;
	return vec;
}

void Point4D::write(std::ostream& os)
{
	os << GetX() << " " << GetY() << " " << GetZ() << " " << GetW() << "\n";
}



void Point4D::read(std::istream& is)
{
	std::cout << "input x, y, z\n";
	is >> x >> y >> z >> w;
}

void Point4D::writefile(std::ofstream& ofs)
{
	ofs << GetX() << GetY() << GetZ() << GetW(); 
}

void Point4D::readfile(std::ifstream& ifs)
{
	ifs >> x >> y >> z >> w;
}

