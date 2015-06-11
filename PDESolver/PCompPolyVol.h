
#ifndef PCOMPPOLYVOL
#define PCOMPPOLYVOL


#include "comppolyvol.h"


class CompPolyVolDouble : public CompPolyVol<double>
{
	virtual ObjectID Identity() const { return std::string("class CompPolyVolDouble"); }
public:
	CompPolyVolDouble();
	CompPolyVolDouble(const CompPolyVol<double>&);
};

class FCompPolyVol : public CompPolyVol<Point1D>
{
	virtual ObjectID Identity() const { return std::string("class FCompPolyVol"); }
public:
	FCompPolyVol();
	FCompPolyVol(const CompPolyVol<Point1D>&);
};


class PCompPolyVol3D : public CompPolyVol<Point3D>
{
	virtual ObjectID Identity() const { return std::string("class PCompPolyVol3D"); }
public:
	PCompPolyVol3D();
	PCompPolyVol3D(const CompPolyVol<Point3D>&);
};

#endif

