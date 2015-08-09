#include "PolynomialParser.h"

#include <cctype>

const char PolynomialParser::allowedChars_[] = { '+', '-', '*', '.', '^', 'u', 'v', 'w' };

PolynomialParser::PolynomialParser(std::string equation) : PolynomialParser(equation, 'f', 3)
{
}

PolynomialParser::PolynomialParser(std::string equation, const char invalidVariable, size_t variableNo) : eq_(equation), invalidVariable_(invalidVariable), variablesNo_(variableNo), eqOrder_(1), coeffArraySize_(0)
{
	removeSpaces(eq_);
	removeMultiplicationSign(eq_);
	eq_ = convertMinusSign(eq_);
}

void PolynomialParser::removeSpaces(std::string& str)
{
	str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
}

void PolynomialParser::removeMultiplicationSign(std::string& str)
{
	str.erase(std::remove(str.begin(), str.end(), '*'), str.end());
}

std::string PolynomialParser::convertMinusSign(std::string& str)
{
	std::stringstream ss;
	for (int i = 0; i < eq_.size(); i++)
	{
		const char c = eq_[i];
		if (c == '-' && i != 0 && eq_[i] != '^')
		{
			ss << '+';
		}
		ss << c;
	}
	return ss.str();
}

void PolynomialParser::validate()
{
	for (int i = 0; i < eq_.size(); i++)
	{
		const char c = eq_[i];
		
		if (!isCharAllowed(c) || c == invalidVariable_)
		{	
			std::stringstream ss;
			ss << "Invalid character " << c << " at column " << i << ". Allowed characters are: [";
			printAllowedChars(ss);
			ss << "]\n";
			throw std::invalid_argument(ss.str());
		}
	}
	for (int i = 0; i < eq_.size(); i++)
	{
		const char c = eq_[i];
		std::stringstream ss;
		if (c == '^')
		{
			if (i == 0 || !isVariableChar(eq_[i - 1]))
			{
				ss << "^ at column " << i << " must be preceded by variable character [u,v,w]";
				throw std::invalid_argument(ss.str());
			}
			if (i == eq_.size() || !isdigit(eq_[i + 1]))
			{
				ss << "^ at column " << i << " must be followed by digit";
				throw std::invalid_argument(ss.str());
			}

		}
	}
}

PolynomialData PolynomialParser::parse()
{
	std::map<char, bool> vars = { { 'u', false }, { 'v', false }, { 'w', false } };
	size_t prevPos = 0;
	size_t pos;
	std::vector<SingleEquation> equations;
	do
	{
		pos = eq_.find('+', prevPos + 1);
		SingleEquation v(eq_.substr(prevPos, pos - prevPos));
		eqOrder_ = v.getOrder() > eqOrder_ ? v.getOrder() : eqOrder_;
		equations.push_back(v);
		prevPos = pos;	
	} while (pos != std::string::npos);
	//std::copy(equations.begin(), equations.end(), std::ostream_iterator<SingleEquation>(std::cout, "\n"));

	std::vector<double> result(std::pow(eqOrder_, variablesNo_));
	std::fill(result.begin(), result.end(), 0.0);
	for (int i = 0; i < equations.size(); i++)
	{
		size_t idx = getSingleEquationIndex(equations[i]);
		result[idx] = equations[i].getCoeff();
	}

	return PolynomialData(eqOrder_, result);
}

size_t PolynomialParser::getSingleEquationIndex(SingleEquation& se)
{
	short first = 0;
	short second = 0;
	short third = 0;
	if (variablesNo_ == 2)
	{
		if (invalidVariable_ == 'u')
		{
			first = se.getPower('v');
			second = se.getPower('w');
		}
		else if (invalidVariable_ == 'v')
		{
			first = se.getPower('u');
			second = se.getPower('w');
		}
		else if (invalidVariable_ == 'w')
		{
			first = se.getPower('u');
			second = se.getPower('v');
		}
		return second * eqOrder_ + first;
	}
	else
	{
		first = se.getPower('u');
		second = se.getPower('v');
		third = se.getPower('w');

		return (second * eqOrder_ + first) + (third * std::pow(eqOrder_, 2));
	}
}

void PolynomialParser::printAllowedChars(std::ostream& o)
{
	std::vector<char> chars(std::begin(allowedChars_), std::end(allowedChars_));
	chars.erase(std::remove(chars.begin(), chars.end(), invalidVariable_), chars.end());
	std::copy(std::begin(chars), std::end(chars), infix_ostream_iterator<char>(o, ","));
}

bool PolynomialParser::isCharAllowed(const char& c) const
{
	const char* ch = std::find(std::begin(allowedChars_), std::end(allowedChars_), c);
	if (ch != std::end(allowedChars_))
		return true;
	return isdigit(c);
}

bool PolynomialParser::isVariableChar(const char& c) const
{
	if (c == 'u' || c == 'v' || c == 'w')
		return true;
	return false;
}


PolynomialParser::~PolynomialParser()
{
}

std::ostream& operator<<(std::ostream& o, const PolynomialParser& p)
{
	o << p.eq_ << std::endl;
	o << "Number of variables: " << p.variablesNo_ << std::endl;
	o << "Order: " << p.eqOrder_ << std::endl;
	o << "Array size: " << p.coeffArraySize_ << std::endl;
	return o;
}


//////// SINGLE EQUATION

SingleEquation::SingleEquation(std::string str) : coeff_(0)
{
	powers_['u'] = 0;
	powers_['v'] = 0;
	powers_['w'] = 0;

	short signCoeff = 1;
	bool isPowerNumber = false;
	char var;
	int start = 0;
	bool inNumber = false;
	bool inVariable = false;

	std::stringstream ss;
	for (int i = 0; i < str.size(); i++)
	{
		const char c = str[i];
		if (std::isdigit(c) || c == '.'){
			inNumber = true;
			ss << c;
		}
		else
		{	
			if (inNumber)
			{	
				double number = std::stod(ss.str());
				ss.str("");
				if (!isPowerNumber)
					coeff_ += signCoeff * number;
				else
					powers_[var] += signCoeff * number;
				signCoeff = 1;
				isPowerNumber = false;
				inNumber = false;
			}
			if (c == '-')
			{
				signCoeff = -1;
			}
			if (isVariableChar(c) && str[i+1] != '^')
			{
				powers_[c]++;
			}
			if (c == '^')
			{
				isPowerNumber = true;
				var = str[i - 1];
			}
		}
	}
	if (ss.str().size() > 0)
	{
		double number = std::stod(ss.str());
		ss.clear();
		if (!isPowerNumber)
			coeff_ = signCoeff * number;
		else
			powers_[var] = signCoeff * number;
	}
	if (coeff_ == 0 && (powers_['u'] != 0 || powers_['v'] != 0 || powers_['w'] != 0))
		coeff_ = 1;
}

bool SingleEquation::isVariableChar(const char& c) const
{
	if (c == 'u' || c == 'v' || c == 'w')
		return true;
	return false;
}

short SingleEquation::getOrder()
{
	short max = 0;
	for (std::map<char, short>::iterator it = powers_.begin(); it != powers_.end(); ++it)
	{
		if (it->second > max)
		{
			max = it->second;
		}
	}
	return max + 1;
}

size_t SingleEquation::variablesNumber()
{
	size_t n = 0;
	for (std::map<char, short>::iterator it = powers_.begin(); it != powers_.end(); ++it)
	{
		if (it->second > 0)
		{
			n++;
		}
	}
	return n;
}

std::ostream& operator<<(std::ostream& o, const SingleEquation& v)
{
	std::map<char, short> p = v.powers_;

	o << "Coeff: " << v.coeff_ << std::endl;
	o << "Power of u: " << p['u'] << std::endl;
	o << "Power of v: " << p['v'] << std::endl;
	o << "Power of w: " << p['w'] << std::endl;

	return o;
}

SingleEquation::~SingleEquation()
{

}