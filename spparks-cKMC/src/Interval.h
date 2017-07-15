// Interval.h: interface for the Interval class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GENERALIZED_INTERVAL_H_)
#define _GENERALIZED_INTERVAL_H_


#include "IVException.h"

#include <iostream>
using namespace std;

//enum BOOLEAN { FALSE, TRUE };

class Interval  
{
private:
	double inf;
//	double nominal;
	double sup;

public:
	static const double POS_INF;
	static const double NEG_INF;

public:
	//Arithmetic operations
	Interval operator+(const Interval& right);
	Interval operator-(const Interval& right);
	Interval operator*(const Interval& right);
	Interval operator/(const Interval& right);
	Interval& operator=(const Interval& right);
	Interval power(double exp);

	Interval dual();
	Interval prop();
	Interval impr();

	//mathmatical operations/functions
	friend Interval operator+(double left, const Interval& right);
	friend Interval operator-(double left, const Interval& right);
	friend Interval operator*(double left, const Interval& right);
	friend Interval operator/(double left, const Interval& right);
	friend Interval operator-(const Interval& right);
	friend Interval sin(const Interval& p);
	friend Interval cos(const Interval& p);

	//Set operations
	Interval union_with(const Interval& right);
	Interval intersect_with(const Interval& right);

	//relations
	int operator==(const Interval &right);
	int operator>=(const Interval &right);
	int operator>(const Interval &right);
	int operator<=(const Interval &right);
	int operator<(const Interval &right);

	int EQ(const Interval &right);
	int GT_EQ(const Interval &right);
	int GT(const Interval &right);
	int LT_EQ(const Interval &right);
	int LT(const Interval &right);
	int S_GT_EQ(const Interval &right);
	int S_GT(const Interval &right);
	int S_LT_EQ(const Interval &right);
	int S_LT(const Interval &right);
	int INCLUDE(const Interval &right);

	//Modifiers
//	void setNominal(double nominal);
	void setInf(double L);
	void setSup(double R);
	void setLBound(double l);
	void setUBound(double u);

	//Access methods
	double Inf();
	double Sup();
	double LBound();
	double UBound();
//	double Nominal();
	double Mid();
	double Width();
	int			isProperInterval();
	int			isImproperInterval();
	int			isPointInterval();
	Interval	SET();

	//copy constructor
	Interval(const Interval &src);

	//constructors and destructors
	Interval(double value);
	Interval(double l, double u);
	Interval(double l, double n, double u);
	Interval();
	virtual ~Interval();

	friend ostream &operator<<(ostream &stream, const Interval &itv);
};

#ifndef _CONSTANT_PI_
#define _CONSTANT_PI_
const double PI = 3.1415926535897932384626;
#endif


#endif // !defined(_GENERALIZED_INTERVAL_H_)
