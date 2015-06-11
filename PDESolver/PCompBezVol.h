
#ifndef PCOMPBEZVOL
#define PCOMPBEZVOL

#include "compbezvol.h"


class CompBezVolDouble : public CompBezVol<double>
{
	virtual ObjectID Identity() const { return std::string("class CompBezVolDouble"); }
public:
	CompBezVolDouble();
	CompBezVolDouble(const CompBezVol<double>&);
};

class FCompBezVol : public CompBezVol<Point1D>
{
	virtual ObjectID Identity() const { return std::string("class FCompBezVol"); }
public:
	FCompBezVol();
	FCompBezVol(const CompBezVol<Point1D>&);

};

class PCompBezVol3D : public CompBezVol<Point3D>
{
	virtual ObjectID Identity() const { return std::string("class PCompBezVol3D"); }
public:
	PCompBezVol3D();
	PCompBezVol3D(const CompBezVol<Point3D>&);

};

#endif