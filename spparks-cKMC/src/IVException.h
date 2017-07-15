// IVException.h: interface for the Exception class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_EXCEPTION_H__F386C9E0_B462_11D6_B9B6_006008BFAB46__INCLUDED_)
#define AFX_EXCEPTION_H__F386C9E0_B462_11D6_B9B6_006008BFAB46__INCLUDED_

#include <string.h>

class IVException  
{
	char message[64];
public:
	char* getMessage();

	IVException(const char* msg);
	virtual ~IVException();

};

#endif // !defined(AFX_EXCEPTION_H__F386C9E0_B462_11D6_B9B6_006008BFAB46__INCLUDED_)
