
#include "curves.h"
#include "point.h"


void SinePi::read(std::istream& is)
{
	std::cout << "input the left and right limits for the sine curve\n";
	is >> x1 >> x2;
}


void SinePi::readfile(std::ifstream& ifs)
{
	ifs >> x1 >> x2;
}


void SinePi::write(std::ostream& os) 
{
	os << "the curve is the sine curve\n";
	os << "the left and right limits are\n";
	os << x1 << " " << x2 << "\n";
}


void SinePi::writefile(std::ofstream& ofs) 
{
	ofs << "the curve is the sine curve\n";
	ofs << "the left and right limits are\n";
	ofs << x1 << " " << x2 << "\n";
}


double SinePi::GetLeftLimit() const
{
	return x1;
}


double SinePi::GetRightLimit() const
{
	return x2;
}


SinePi::SinePi(double X1, double X2) : x1(X1), x2(X2)
{
}


double SinePi::operator()(double x) const
{
	return sin(Pi*x);
}


double SinePi::operator()(int val, double x) const
{    
	return Derive(val,x);
}


double SinePi::Derive(int level, double x) const
{
	if (level == 0) return sin(Pi*x);
	else if (level == 1) return Pi*cos(Pi*x);
	else if (level == 2) return -Pi*Pi*sin(Pi*x);
	else return -Pi*Pi*Pi*cos(Pi*x);
}


void Sine::read(std::istream& is)
{
	std::cout << "input the left and right limits for the sine curve\n";
	is >> x1 >> x2;
}

void Sine::readfile(std::ifstream& ifs)
{
	ifs >> x1 >> x2;
}


void Sine::write(std::ostream& os) 
{
	os << "the curve is the sine curve\n";
	os << "the left and right limits are\n";
	os << x1 << " " << x2 << "\n";
}


void Sine::writefile(std::ofstream& ofs) 
{
	ofs << "the curve is the sine curve\n";
	ofs << "the left and right limits are\n";
	ofs << x1 << " " << x2 << "\n";
}


double Sine::GetLeftLimit() const
{
	return x1;
}


double Sine::GetRightLimit() const
{
	return x2;
}


Sine::Sine(double X1, double X2) : x1(X1), x2(X2)
{
}

double Sine::operator()(double x) const
{
	return sin(x);
}

double Sine::operator()(int val, double x) const
{    
	return Derive(val,x);
}

double Sine::Derive(int level, double x) const
{
	if (level == 0) return sin(x);
	else if (level == 1) return cos(x);
	else if (level == 2) return -sin(x);
	else return -cos(x);
}


double Cosine::GetLeftLimit() const
{
	return x1;
}


double Cosine::GetRightLimit() const
{
	return x2;
}


Cosine::Cosine(double X1, double X2) : x1(X1), x2(X2)
{
}


double Cosine::operator()(double x) const
{
	return cos(x);
}

double Cosine::operator()(int val, double x) const
{    
	return Derive(val,x);
}

double Cosine::Derive(int level, double x) const
{
	if (level == 0) return cos(x);
	else if (level == 1) return -sin(x);
	else if (level == 2) return -cos(x);
	else return sin(x);
}



void Cosine::read(std::istream& is)
{
	std::cout << "input the left and right limits for the sine curve\n";
	is >> x1 >> x2;
}


void Cosine::readfile(std::ifstream& ifs)
{
	ifs >> x1 >> x2;
}


void Cosine::write(std::ostream& os) 
{
	os << "the curve is the cosine curve\n";
	os << "the left and right limits are\n";
	os << x1 << " " << x2 << "\n";
}


void Cosine::writefile(std::ofstream& ofs)
{
	ofs << "the curve is the cosine curve\n";
	ofs << "the left and right limits are\n";
	ofs << x1 << " " << x2 << "\n";
}



double Circle::GetLeftLimit() const
{
	return x1;
}


double Circle::GetRightLimit() const
{
	return x2;
}


double Circle::GetRadius() const
{
	return r;
}


Circle::Circle(double X1, double X2, double R) : x1(X1), x2(X2), r(R)
{
}

Point3D Circle::operator()(int val, double x) const
{    
	return Derive(val,x);
}

Point3D Circle::operator()(double x) const
{
	return Point3D(r*cos(x),r*sin(x),0.0);
}



Point3D Circle::Derive(int level, double x) const
{
	if (level == 0) return Point3D(r*cos(x),r*sin(x),0.0);
	else if (level == 1) return Point3D(-r*sin(x),r*cos(x),0.0);
	else return Point3D(-r*cos(x),-r*sin(x),0.0);
}



void Circle::read(std::istream& is)
{
	std::cout<< "input the left and right limits for the curve\n";
	is >> x1 >> x2;
	std::cout << " input the radius ";
	is >> r;
}


void Circle::readfile(std::ifstream& ifs)
{
	ifs >> x1 >> x2;
	ifs >> r;
}


void Circle::write(std::ostream& os) 
{
	os << "left and right limits for the curve\n";
	os << x1 << x2;
	os << "radius\n";
	os << r;
}


void Circle::writefile(std::ofstream& ofs) 
{
	ofs << "left and right limits for the curve\n";
	ofs << x1 << x2;
	ofs << "radius\n";
	ofs << r;
}




double Quadratic::GetLeftLimit() const
{
	return x1;
}


double Quadratic::GetRightLimit() const
{
	return x2;
}


Quadratic::Quadratic(double X1, double X2, double A, double B, double C) : x1(X1), x2(X2), a(A), b(B), c(C)
{
}


double Quadratic::operator()(double x) const
{
	return a*x*x+b*x+c;
}

double Quadratic::operator()(int val, double x) const
{    
	return Derive(val,x);
}

double Quadratic::Derive(int level, double x) const
{
	if (level == 0) return a*x*x+b*x+c;
	else if (level == 1) return 2*a*x+b;
	else return 2*a;
}




void Quadratic::read(std::istream& is)
{
	std::cout << "input the left and right limits for the curve\n";
	is >> x1 >> x2;
	std::cout << "input the coeffs\n";
	is >> a >> b >> c;
}


void Quadratic::readfile(std::ifstream& ifs)
{
	ifs >> x1 >> x2;
	ifs >> a >> b >> c;
}


void Quadratic::write(std::ostream& os) 
{
	os << "left and right limits for the curve\n";
	os << x1 << x2;
	os << "coeffs \n";
	os << a << " " << b << " " << c;
}


void Quadratic::writefile(std::ofstream& ofs) 
{
	ofs << "left and right limits for the curve\n";
	ofs << x1 << x2;
	ofs << "coeffs \n";
	ofs << a << " " << b  << " " << c;
}




double Linear::GetLeftLimit() const
{
	return x1;
}


double Linear::GetRightLimit() const
{
	return x2;
}


Linear::Linear(double X1, double X2, double A, double B) : x1(X1), x2(X2), a(A), b(B)
{
}


double Linear::operator()(double x) const
{
	return a*x+b;
}

double Linear::operator()(int val, double x) const
{    
	return Derive(val,x);
}

double Linear::Derive(int level, double x) const
{
	if (level == 0) return a*x+b;
	else if (level == 1) return a;
	else return 0.0;
}




void Linear::read(std::istream& is)
{
	std::cout << "input the left and right limits for the curve\n";
	is >> x1 >> x2;
	std::cout << "input the coeffs\n";
	is >> a >> b;
}


void Linear::readfile(std::ifstream& ifs)
{
	ifs >> x1 >> x2;
	ifs >> a >> b;
}


void Linear::write(std::ostream& os) 
{
	os << "left and right limits for the curve\n";
	os << x1 << x2;
	os << "coeffs \n";
	os << a << " " << b  << "\n";
}


void Linear::writefile(std::ofstream& ofs) 
{
	ofs << "left and right limits for the curve\n";
	ofs << x1 << x2;
	ofs << "coeffs \n";
	ofs << a << " " << b << "\n";
}

