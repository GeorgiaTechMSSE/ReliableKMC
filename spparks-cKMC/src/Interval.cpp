// Interval.cpp: implementation of the Interval class.
//
//////////////////////////////////////////////////////////////////////

#include "Interval.h"

#include <math.h>

const double Interval::POS_INF = 1e200;
const double Interval::NEG_INF = -1e200;


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Interval::Interval()
{
	inf = 0.0;
//	nominal = 0.0;
	sup = 0.0;
}

Interval::Interval(double l, double r)
{
//	if(l>u)	throw IVException("IVException: Construct empty interval.");
	inf = l;
//	nominal = (l+r)/2;
	sup = r;
}
/*
Interval::Interval(double l, double n, double r)
{
	if(n<l && n<r || n>l && n>r)
	{
		cerr<<"["<<l<<", "<<n<<", "<<r<<"]"<<endl;
		throw IVException("IVException: Construct invalid interval - nominal value.");
	}
	inf = l;
	nominal = n;
	sup = r;
}
*/
Interval::Interval(double value)
{
	inf = value;
//	nominal = value;
	sup = value;
}

Interval::~Interval()
{

}

/////////////////////
// copy constructor
/////////////////////
Interval::Interval(const Interval &src)
{
	inf = src.inf;
//	nominal = src.nominal;
	sup = src.sup;
}

//////////////////////////////////////////////////////////////////////
// Access functions
//////////////////////////////////////////////////////////////////////

double Interval::Sup()
{
	return sup;
}

double Interval::Inf()
{
	return inf;
}
/*
double Interval::Nominal()
{
	return nominal;
}
*/
double Interval::Mid()
{
	return (sup+inf)/2.0;
}

double Interval::Width()
{
	return fabs(sup-inf);
}

int Interval::isProperInterval()
{
	return inf<=sup ? 1 : 0;
}

int Interval::isImproperInterval()
{
	return inf>=sup ? 1 : 0;
}

int Interval::isPointInterval()
{
	return inf==sup ? 1 : 0;
}

Interval Interval::SET()
{
//	return Interval(this->LBound(), this->Nominal(), this->UBound());
	return Interval(this->LBound(), this->UBound());
}

double Interval::UBound()
{
	return inf >= sup ? inf : sup;
}

double Interval::LBound()
{
	return inf <= sup ? inf : sup;
}

//////////////////////////////////////////////////////////////////////
// Modifiers
//////////////////////////////////////////////////////////////////////

void Interval::setSup(double r)
{
	sup = r;
}

void Interval::setInf(double l)
{
	inf = l;
}

void Interval::setUBound(double r)
{
	sup = r;
}

void Interval::setLBound(double l)
{
	inf = l;
}
/*
void Interval::setNominal(double nominal)
{
	this->nominal = nominal;
}
*/

//////////////////////////////////////////////////////////////////////
// I/O overload
//////////////////////////////////////////////////////////////////////

ostream &operator <<(ostream &stream, const Interval &itv)
{
//	stream << "[" << itv.inf <<", "<< itv.nominal <<", "<< itv.sup << "]";
	stream << "[" << itv.inf <<", "<< itv.sup << "]";
	return stream;
}


//////////////////////////////////////////////////////////////////////
// Operator overloads
//////////////////////////////////////////////////////////////////////

Interval Interval::operator +(const Interval &right)
{
	Interval result = right;
//	if(result.isEmpty())	
//		throw IVException("IVException: add an empty interval.");
	result.setInf(inf+right.inf);
//	result.setNominal(nominal+right.nominal);
	result.setSup(sup+right.sup);
	return result;
}

Interval Interval::operator -(const Interval &right)
{
	Interval result = right;
//	if(result.isEmpty())	
//		throw IVException("IVException: subtract an empty interval.");
	result.setInf(inf-right.sup);
//	result.setNominal(nominal-right.nominal);
	result.setSup(sup-right.inf);
	return result;
}

Interval Interval::operator *(const Interval &right)
{
//	Interval result = right;
//	if(result.isEmpty())	
//		throw IVException("IVException: multiply an empty interval.");
	double a1b1, a1b2, a2b1, a2b2;
	a1b1 = inf*right.inf;
	a1b2 = inf*right.sup;
	a2b1 = sup*right.inf;
	a2b2 = sup*right.sup;
/*	int min = 0, max = 0;
	for(int i=1; i<4; i++)
	{
		if(m[i]<m[min]) min=i;
		if(m[i]>m[max]) max=i;
	}
	result.setL(m[min]);
	result.setNominal(nominal*right.nominal);
	result.setR(m[max]);
	return result;
*/
	if ( inf>=0 && sup>=0 && right.inf>=0 && right.sup>=0 )
		return Interval(a1b1, a2b2);
	else if ( inf>=0 && sup>=0 && right.inf>=0 && right.sup<0 )
		return Interval(a1b1, a1b2);
	else if ( inf>=0 && sup>=0 && right.inf<0 && right.sup>=0 )
		return Interval(a2b1, a2b2);
	else if ( inf>=0 && sup>=0 && right.inf<0 && right.sup<0 )
		return Interval(a2b1, a1b2);
	else if ( inf>=0 && sup<0 && right.inf>=0 && right.sup>=0 )
		return Interval(a1b1, a2b1);
	else if ( inf>=0 && sup<0 && right.inf>=0 && right.sup<0 )
		return Interval((a2b2>=a1b1?a2b2:a1b1), (a2b1<=a1b2?a2b1:a1b2));
	else if ( inf>=0 && sup<0 && right.inf<0 && right.sup>=0 )
		return Interval(0.0,0.0);
	else if ( inf>=0 && sup<0 && right.inf<0 && right.sup<0 )
		return Interval(a2b2, a1b2);
	else if ( inf<0 && sup>=0 && right.inf>=0 && right.sup>=0 )
		return Interval(a1b2, a2b2);
	else if ( inf<0 && sup>=0 && right.inf>=0 && right.sup<0 )
		return Interval(0.0, 0.0);
	else if ( inf<0 && sup>=0 && right.inf<0 && right.sup>=0 )
		return Interval((a1b2<=a2b1?a1b2:a2b1), (a1b1>=a2b2?a1b1:a2b2));
	else if ( inf<0 && sup>=0 && right.inf<0 && right.sup<0 )
		return Interval(a2b1, a1b1);
	else if ( inf<0 && sup<0 && right.inf>=0 && right.sup>=0 )
		return Interval(a1b2, a2b1);
	else if ( inf<0 && sup<0 && right.inf>=0 && right.sup<0 )
		return Interval(a2b2, a2b1);
	else if ( inf<0 && sup<0 && right.inf<0 && right.sup>=0 )
		return Interval(a1b2, a1b1);
	else
		return Interval(a2b2, a1b1);
	
}

Interval Interval::operator/(const Interval& right)
{
	Interval a = right;
	double a1b1, a1b2, a2b1, a2b2;

//	if(a.isEmpty())	
//		throw IVException("IVException: divide an empty interval.");
	// dividend does not include zero
/*	if(!a.INCLUDE(0))
	{
		Interval div(1.0/right.R, 1.0/right.nominal, 1.0/right.inf);
		return this->operator *(div);
	}
	else
	{
		// dividend includes zero
		if(right.inf==0 && right.R==0)
			return Interval(Interval::NEG_INF, 
							0,
							Interval::POS_INF);
		if(this->sup<=0 && right.sup==0)
			return Interval(this->sup/right.inf, 
							this->sup/right.inf,
							Interval::POS_INF);
		if(this->sup<=0 && right.sup>0 && right.inf<0)
			return Interval(Interval::NEG_INF, 
							0,
							Interval::POS_INF);
		if(this->sup<=0 && right.inf==0)
			return Interval(Interval::NEG_INF,
							this->sup/right.sup, 
							this->sup/right.sup);
		if(this->inf<0 && this->sup>0)
			return Interval(Interval::NEG_INF, 
							0,
							Interval::POS_INF);
		if(this->inf>=0 && right.sup==0)
			return Interval(Interval::NEG_INF, 
							this->inf/right.inf,
							this->inf/right.inf);
		if(this->inf>=0 && right.sup>0 && right.inf<0)
			return Interval(Interval::NEG_INF, 
							0,
							Interval::POS_INF);
		if(this->inf>=0 && right.inf==0)
			return Interval(this->inf/right.sup, 
							this->inf/right.sup,
							Interval::POS_INF);
		throw IVException("IVException: Divid by zero");
	}
*/
	if(a.INCLUDE(0))
		throw IVException("IVException: Divid by zero");
	else
	{
		a1b1 = inf/right.inf;
		a1b2 = inf/right.sup;
		a2b1 = sup/right.inf;
		a2b2 = sup/right.sup;

	if ( inf>=0 && sup>=0 && right.inf>0 && right.sup>0 )
		return Interval(a1b2, a2b1);
	else if ( inf>=0 && sup>=0 && right.inf<0 && right.sup<0 )
		return Interval(a2b2, a1b1);
	else if ( inf>=0 && sup<0 && right.inf>0 && right.sup>0 )
		return Interval(a1b2, a2b2);
	else if ( inf>=0 && sup<0 && right.inf<0 && right.sup<0 )
		return Interval(a2b1, a1b1);
	else if ( inf<0 && sup>=0 && right.inf>0 && right.sup>0 )
		return Interval(a1b1, a2b1);
	else if ( inf<0 && sup>=0 && right.inf<0 && right.sup<0 )
		return Interval(a2b2, a1b2);
	else if ( inf<0 && sup<0 && right.inf>0 && right.sup>0 )
		return Interval(a1b1, a2b2);
	else if ( inf<0 && sup<0 && right.inf<0 && right.sup<0 )
		return Interval(a2b1, a1b2);
	else
		throw IVException("IVException: Invalid division");
	}


}


Interval& Interval::operator=(const Interval& right)
{
	if(this == &right) return *this;
	this->inf = right.inf;
//	this->nominal = right.nominal;
	this->sup = right.sup;
	return *this;
}


Interval Interval::power(double exp)
{
	double l, u;
	double a, b;

	if(exp>=0)
	{
/*		
		if(	sup<=0 && l<0 && u<=0			||
			inf<0 && sup>0 && l<0			||
			inf>=0						
			) return Interval(l, n, u);
		else if( sup<=0 && l>0 && u>=0
			) return Interval(u, n, l);
		else if( inf<0 && sup>0 && l>=u && u>0
			) return Interval(0, n, l);
		else if( inf<0 && sup>0 && u>l && l>0
			) return Interval(0, n, u);
		else return Interval(l<u?l:u, n, l>u?l:u);
*/
		if(exp>1)
		{
			int i = exp/2;  
			if(i == exp/2.0) //if exp is even
			{
				if(inf>=0 && sup>=0)
				{
					l = pow(inf,exp);
					u = pow(sup,exp);
				}
				if(inf<0 && sup<0)
				{
					l = pow(sup,exp);
					u = pow(inf,exp);
				}
				if(inf<0 && sup>=0) 
				{
					a = pow(fabs(inf),exp);
					b = pow(fabs(sup),exp);
					l = 0.0;
					u = a>b?a:b;
				}
				if(inf>=0 && sup<0)
				{
					a = pow(fabs(inf),exp);
					b = pow(fabs(sup),exp);
					l = a>b?a:b;
					u = 0.0;
				}
			}
			else //if exp is odd
			{
				l = pow(inf,exp);
				u = pow(sup,exp);
			}
			return Interval(l,u);
		}
		else
		{
			if(inf>=0 && sup>=0)
			{
				l = pow(inf,exp);
				u = pow(sup,exp);
				return Interval(l,u);
			}
			else
				return Interval(1,1);
		}


	}
	
	return 1/this->power(-exp);
}

Interval Interval::dual()
{
//	return Interval(sup, nominal, inf);
	return Interval(sup, inf);
}

Interval Interval::prop()
{
	if (inf <= sup)
		return Interval(inf, sup);
	else
		return Interval(sup, inf);
}

Interval Interval::impr()
{
	if (inf >= sup)
		return Interval(inf, sup);
	else
		return Interval(sup, inf);
}


//////////////////////////////////////////////////////////////////////
//Friend mathmatical operations/functions
//////////////////////////////////////////////////////////////////////
Interval operator+(double left, const Interval& right)
{
	Interval l(left);
	return l+right;
}

Interval operator-(double left, const Interval& right)
{
	Interval l(left);
	return l-right;
}

Interval operator*(double left, const Interval& right)
{
	Interval l(left);
	return l*right;
}

Interval operator/(double left, const Interval& right)
{
	Interval l(left);
	return l/right;
}

Interval operator-(const Interval& right)
{
//	return Interval(-right.sup, -right.nominal, -right.inf);
	return Interval(-right.sup, -right.inf);
}

Interval sin(const Interval& p)
{
	double v_2pi = 2*PI;
//	double v_n = sin(p.nominal);
	if(fabs(p.sup-p.inf)>= v_2pi)
//		return Interval(-1,v_n,1);
		return Interval(-1,1);
	else
	{
		double v_pi2 = PI/2;
		double v_3pi2 = 3*PI/2;
		double v_u = sin(p.sup);
		double v_l = sin(p.inf);
		int addi = 0;
		if(p.inf>0) addi = 0;
		if(p.sup<0) addi = -1;
		double u = p.sup-v_2pi*((int)(p.sup/v_2pi)+addi);
		double l = p.inf-v_2pi*((int)(p.inf/v_2pi)+addi);
		if(u>v_3pi2) u -= v_2pi;
		if(l>v_3pi2) l -= v_2pi;
		if(l<=u  &&  l<v_pi2  &&  u<v_pi2)
//			return Interval(v_l, v_n, v_u);
			return Interval(v_l, v_u);
		if(l<=u  &&  l<v_pi2  &&  u>=v_pi2)
			if(v_l <= v_u)
//				return Interval(v_l, v_n, 1);
				return Interval(v_l, 1);
			else 
//				return Interval(v_u, v_n, 1);
				return Interval(v_u, 1);
		if(l<=u  &&  l>=v_pi2)
//			return Interval(v_u, v_n, v_l);
			return Interval(v_u, v_l);
		if(l>u  &&  l<v_pi2  &&  u<v_pi2)
//			return Interval(-1, v_n, 1);
			return Interval(v_l, v_u);
		if(l>u  &&  l>=v_pi2  &&  u<v_pi2)
			if(v_l <= v_u)
//				return Interval(-1, v_n, v_l);
				return Interval(-1, v_l);
			else 
//				return Interval(-1, v_n, v_u);
				return Interval(-1, v_u);
		if(l>u  &&  u>=v_pi2)
//			return Interval(-1, v_n, 1);
			return Interval(v_u, v_l);
		throw new IVException("IVException: unknown at sin(Interval)");
	}
}
/*
Interval sin(const Interval& p)
{
	double v_2pi = 2*PI;
//	double v_n = sin(p.nominal);
	if(fabs(p.sup-p.inf)>= v_2pi)
//		return Interval(-1,v_n,1);
		return Interval(-1,1);
	else
	{
		double v_pi2 = PI/2;
		double v_3pi2 = 3*PI/2;
		double v_u = sin(p.sup);
		double v_l = sin(p.inf);
		int addi = 0;
		if(p.inf>0) addi = 0;
		if(p.sup<0) addi = -1;
		double u = p.sup-v_2pi*((int)(p.sup/v_2pi)+addi);
		double l = p.inf-v_2pi*((int)(p.inf/v_2pi)+addi);
		if(u>v_3pi2) u -= v_2pi;
		if(l>v_3pi2) l -= v_2pi;
		if(l<=u  &&  l<v_pi2  &&  u<v_pi2)
//			return Interval(v_l, v_n, v_u);
			return Interval(v_l, v_u);
		if(l<=u  &&  l<v_pi2  &&  u>=v_pi2)
			if(v_l <= v_u)
//				return Interval(v_l, v_n, 1);
				return Interval(v_l, 1);
			else 
//				return Interval(v_u, v_n, 1);
				return Interval(v_u, 1);
		if(l<=u  &&  l>=v_pi2)
//			return Interval(v_u, v_n, v_l);
			return Interval(v_u, v_l);
		if(l>u  &&  l<v_pi2  &&  u<v_pi2)
//			return Interval(-1, v_n, 1);
			return Interval(-1, 1);
		if(l>u  &&  l>=v_pi2  &&  u<v_pi2)
			if(v_l <= v_u)
//				return Interval(-1, v_n, v_l);
				return Interval(-1, v_l);
			else 
//				return Interval(-1, v_n, v_u);
				return Interval(-1, v_u);
		if(l>u  &&  u>=v_pi2)
//			return Interval(-1, v_n, 1);
			return Interval(-1, 1);
		throw new IVException("IVException: unknown at sin(Interval)");
	}
}
*/

Interval cos(const Interval& p)
{
	double v_2pi = 2*PI;
//	double v_n = cos(p.nominal);
	if(fabs(p.sup-p.inf)>= v_2pi)
//		return Interval(-1,v_n,1);
		return Interval(-1,1);
	else
	{
		double v_u = cos(p.sup);
		double v_l = cos(p.inf);
		int addi = 0;
		if(p.inf>0) addi = 0;
		if(p.sup<0) addi = -1;
		double u = p.sup-v_2pi*((int)(p.sup/v_2pi)+addi);
		double l = p.inf-v_2pi*((int)(p.inf/v_2pi)+addi);
		if(l<=u  &&  l<PI  &&  u<PI)
//			return Interval(v_u, v_n, v_l);
			return Interval(v_u, v_l);
		if(l<=u  &&  l<PI  &&  u>=PI)
			if(v_l <= v_u)
//				return Interval(-1, v_n, v_u);
				return Interval(-1, v_u);
			else
//				return Interval(-1, v_n, v_l);
				return Interval(-1, v_l);
		if(l<=u  &&  l>=PI)
//			return Interval(v_l, v_n, v_u);
			return Interval(v_l, v_u);
		if(l>u  &&  l<PI  &&  u<PI)
//			return Interval(-1, v_n, 1);
			return Interval(v_u, v_l);
		if(l>u  &&  l>=PI  &&  u<PI)
			if(v_l <= v_u)
//				return Interval(v_l, v_n, 1);
				return Interval(v_l, 1);
			else 
//				return Interval(v_u, v_n, 1);
				return Interval(v_u, 1);
		if(l>u  &&  u>=PI)
//			return Interval(-1, v_n, 1);
			return Interval(v_l, v_u);
		throw new IVException("IVException: unknown at cos(Interval)");
	}
}
/*
Interval cos(const Interval& p)
{
	double v_2pi = 2*PI;
//	double v_n = cos(p.nominal);
	if(fabs(p.sup-p.inf)>= v_2pi)
//		return Interval(-1,v_n,1);
		return Interval(-1,1);
	else
	{
		double v_u = cos(p.sup);
		double v_l = cos(p.inf);
		int addi = 0;
		if(p.inf>0) addi = 0;
		if(p.sup<0) addi = -1;
		double u = p.sup-v_2pi*((int)(p.sup/v_2pi)+addi);
		double l = p.inf-v_2pi*((int)(p.inf/v_2pi)+addi);
		if(l<=u  &&  l<PI  &&  u<PI)
//			return Interval(v_u, v_n, v_l);
			return Interval(v_u, v_l);
		if(l<=u  &&  l<PI  &&  u>=PI)
			if(v_l <= v_u)
//				return Interval(-1, v_n, v_u);
				return Interval(-1, v_u);
			else
//				return Interval(-1, v_n, v_l);
				return Interval(-1, v_l);
		if(l<=u  &&  l>=PI)
//			return Interval(v_l, v_n, v_u);
			return Interval(v_l, v_u);
		if(l>u  &&  l<PI  &&  u<PI)
//			return Interval(-1, v_n, 1);
			return Interval(-1, 1);
		if(l>u  &&  l>=PI  &&  u<PI)
			if(v_l <= v_u)
//				return Interval(v_l, v_n, 1);
				return Interval(v_l, 1);
			else 
//				return Interval(v_u, v_n, 1);
				return Interval(v_u, 1);
		if(l>u  &&  u>=PI)
//			return Interval(-1, v_n, 1);
			return Interval(-1, 1);
		throw new IVException("IVException: unknown at cos(Interval)");
	}
}
*/
//////////////////////////////////////////////////////////////////////
//Set operations
//////////////////////////////////////////////////////////////////////
Interval Interval::union_with(const Interval& right)
{
	Interval result = right;
/*	if(result.isEmpty())
	{
		result.setL(this->L);
		result.setNominal(this->nominal);
		result.setR(this->R);
	}
	else
*/	{
		if(this->inf < right.inf)
			result.setInf(this->inf);
		else
			result.setInf(right.inf);
		if(this->sup > right.sup)
			result.setSup(this->sup);
		else
			result.setSup(right.sup);
//		result.setNominal((this->nominal+right.nominal)/2);
	}
	return result;
}

Interval Interval::intersect_with(const Interval& right)
{
	Interval result = right;
/*	if(result.isEmpty())
	{
		result.setL(right.L);
		result.setNominal(right.nominal);
		result.setR(right.R);
	}
	else
*/	{
		if(this->inf > right.inf)
			result.setInf(this->inf);
		else
			result.setInf(right.inf);
		if(this->sup < right.sup)
			result.setSup(this->sup);
		else
			result.setSup(right.sup);
//		result.setNominal(result.Mid());
	}
	return result;
}


//////////////////////////////////////////////////////////////////////
// Relations
//////////////////////////////////////////////////////////////////////
int Interval::operator ==(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if( inf==right.inf && 
//		nominal==right.nominal && 
		sup==right.sup)
		return 1;
	else
		return 0;
}

int Interval::operator >=(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup>=right.sup &&
		inf>=right.inf)
		return 1;
	else
		return 0;
}

int Interval::operator >(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup>right.sup &&
		inf>right.inf)
		return 1;
	else
		return 0;
}

int Interval::operator <=(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<=right.sup &&
		inf<=right.inf)
		return 1;
	else
		return 0;
}

int Interval::operator <(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<right.sup &&
		inf<right.inf)
		return 1;
	else
		return 0;
}


//////////////////////////////////////////////////////////////////////
// Comparison and relations
//////////////////////////////////////////////////////////////////////

int Interval::EQ(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if( inf==right.inf && 
//		nominal==right.nominal && 
		sup==right.sup)
		return 1;
	else
		return 0;
}

int Interval::GT_EQ(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup>=right.sup &&
		inf>=right.inf)
		return 1;
	else
		return 0;
}

int Interval::GT(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup>right.sup &&
		inf>right.inf)
		return 1;
	else
		return 0;
}

int Interval::LT_EQ(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<=right.sup &&
		inf<=right.inf)
		return 1;
	else
		return 0;
}

int Interval::LT(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<right.sup &&
		inf<right.inf)
		return 1;
	else
		return 0;
}


int Interval::S_GT_EQ(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(inf>=right.sup)
		return 1;
	else
		return 0;
}

int Interval::S_GT(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(inf>right.sup)
		return 1;
	else
		return 0;
}

int Interval::S_LT_EQ(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<=right.inf)
		return 1;
	else
		return 0;
}

int Interval::S_LT(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/	if(sup<right.inf)
		return 1;
	else
		return 0;
}

int Interval::INCLUDE(const Interval &right)
{
/*	Interval a = right;
	if(a.isEmpty())	
		throw IVException("IVException: compare to an empty interval.");
*/
//	if( sup>=right.sup && inf<=right.inf  || sup<=right.sup && inf>=right.inf )
	if( sup>=right.sup && inf<=right.inf )
		return 1;
	else
		return 0;
}
