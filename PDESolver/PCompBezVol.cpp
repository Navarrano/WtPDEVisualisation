
#include "pcompbezvol.h"


CompBezVolDouble::CompBezVolDouble() : CompBezVol<double>() {}
CompBezVolDouble::CompBezVolDouble(const CompBezVol<double>& v) : CompBezVol<double>(v) {}

FCompBezVol::FCompBezVol() : CompBezVol<Point1D>() {}
FCompBezVol::FCompBezVol(const CompBezVol<Point1D>& v) : CompBezVol<Point1D>(v) {}




PCompBezVol3D::PCompBezVol3D() : CompBezVol<Point3D>() {}
PCompBezVol3D::PCompBezVol3D(const CompBezVol<Point3D>& v) : CompBezVol<Point3D>(v) {}




