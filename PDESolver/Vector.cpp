#include "vector.h"

std::ofstream & operator<< (std::ofstream &ofs, FTextObject& r)
{
	r.writefile(ofs);
	return ofs;
}


std::ofstream & operator<< (std::ofstream &ofs, PFTextObject p)
{
	p->writefile(ofs);
	return ofs;
}

std::ifstream & operator>> (std::ifstream &ifs, RFTextObject r)
{
	r.readfile(ifs);
	return ifs;
}

std::ifstream & operator>> (std::ifstream& ifs, PFTextObject p)
{
	p->readfile(ifs);
	return ifs;
}

std::ostream & operator<< (std::ostream &os, TextObject& r)
{
	r.write(os);
	return os;
}


std::ostream & operator<< (std::ostream &os, PTextObject p)
{
	p->write(os);
	return os;
}

std::istream & operator>> (std::istream &is, RTextObject r)
{
	r.read(is);
	return is;
}

std::istream & operator>> (std::istream &is, PTextObject p)
{
	p->read(is);
	return is;
}
