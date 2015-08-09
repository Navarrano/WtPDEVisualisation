#ifndef POLYNOMIAL_PARSER_H
#define POLYNOMIAL_PARSER_H

#include <string>
#include <iostream>
#include <algorithm>
#include <array>
#include <iterator>
#include <sstream>
#include <vector>
#include <map>

#include "infix_iterator.h"

class SingleEquation;

struct PolynomialData
{
	PolynomialData() {};
	PolynomialData(short order, std::vector<double> coeffs) : coeffs_(coeffs), order_(order) {}
	std::vector<double> coeffs_;
	short order_;
};

class PolynomialParser
{
public:
	PolynomialParser(std::string equation = "");
	PolynomialParser(std::string equation, const char invalidVariable, size_t variableNo = 2);
	~PolynomialParser();
	void validate();
	PolynomialData parse();

	friend std::ostream& operator<<(std::ostream& o, const PolynomialParser& p);
private:
	void removeSpaces(std::string& str);
	void removeMultiplicationSign(std::string& str);
	std::string convertMinusSign(std::string& str);
	size_t getSingleEquationIndex(SingleEquation& se);

	void printAllowedChars(std::ostream& o = std::cout);
	bool isCharAllowed(const char& c) const;
	bool isVariableChar(const char& c) const;

	std::string eq_;
	static const char allowedChars_[];
	char invalidVariable_;
	const size_t variablesNo_;
	size_t eqOrder_;
	size_t coeffArraySize_;
};

class SingleEquation
{
public:
	SingleEquation(std::string str);
	~SingleEquation();
	short getOrder();
	double getCoeff() { return coeff_; };
	short getPower(const char c) { return powers_[c]; };
	size_t variablesNumber();
	
	friend std::ostream& operator<<(std::ostream& o, const SingleEquation& v);
private:
	bool isVariableChar(const char& c) const;
	

	double coeff_;
	std::map<char, short> powers_;
};

#endif