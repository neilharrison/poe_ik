// Elapsed.h: interface for the CElapsed class.
//
// CElapsed class wraps up windows queryperformancecounter functionalities
// Provides microseconds timing accuracy.
//
//////////////////////////////////////////////////////////////////////

#if !defined(DEF_ELAPSED)
#define DEF_ELAPSED

#pragma once

//#if _MSC_VER > 1000
//#pragma once
//#endif // _MSC_VER > 1000
#include <chrono>

class CElapsed
{
public:
	CElapsed();
	virtual ~CElapsed();

	bool Begin();    // start timing
	double End();
	bool Available();
	long long int GetFreq();
	std::chrono::time_point<std::chrono::high_resolution_clock> m_iBeginTime;
	double m_dElapsed;
private :
    int m_iInitialized;
    long long int m_iFrequency;
    

};

#endif
