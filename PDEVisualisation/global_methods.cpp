#include "global_methods.h"

#include <iomanip>
#include <sstream>

std::string toString(const double val, short precision)
{
	std::stringstream str;
	str << std::fixed << std::setprecision(precision) << val;
	return str.str();
}