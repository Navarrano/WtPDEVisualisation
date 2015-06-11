
#include "pcomppolyvol.h"

CompPolyVolDouble::CompPolyVolDouble() : CompPolyVol<double>() {}

CompPolyVolDouble::CompPolyVolDouble(const CompPolyVol<double>& v) : CompPolyVol<double>(v) {}



FCompPolyVol::FCompPolyVol() : CompPolyVol<Point1D>() {}

FCompPolyVol::FCompPolyVol(const CompPolyVol<Point1D>& v) : CompPolyVol<Point1D>(v) {}




PCompPolyVol3D::PCompPolyVol3D() : CompPolyVol<Point3D>() {}


PCompPolyVol3D::PCompPolyVol3D(const CompPolyVol<Point3D>& v) : CompPolyVol<Point3D>(v) {}

