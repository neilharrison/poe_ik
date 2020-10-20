// Elapsed.cpp: implementation of the CElapsed class.
//
//////////////////////////////////////////////////////////////////////
//#include <windows.h>
#include "Elapsed.h"
#include <inttypes.h>
#include <iostream>
//#ifdef _DEBUG
//#undef THIS_FILE
//static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
//#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CElapsed::CElapsed()
{
	// get the frequency of the counter
//	m_iInitialized = QueryPerformanceFrequency( (LARGE_INTEGER *)&m_iFrequency );
//std::cout<<"Init elapsed"<<std::endl;
m_iInitialized = 1;
}

CElapsed::~CElapsed()
{

}

bool CElapsed::Begin()    // start timing
{
	//if (!m_iInitialized)
	//	return 0;   // error - couldn't get frequency

	// get the starting counter value
	m_iBeginTime = std::chrono::high_resolution_clock::now();
	//std::cout<<"Hello"<<std::endl;
	return 1 ;
//	return QueryPerformanceCounter( (LARGE_INTEGER *)&m_iBeginTime );
}

double CElapsed::End()    // stop timing and get elapsed time in seconds
{
	//if (!m_iInitialized )
		//return 0.0; // error - couldn't get frequency

	// get the ending counter value
	//QueryPerformanceCounter( (LARGE_INTEGER *)&endtime );
	std::chrono::time_point<std::chrono::high_resolution_clock> endtime = std::chrono::high_resolution_clock::now();
	// determine the elapsed counts
	std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(endtime - m_iBeginTime);
	//std::cout<<elapsed.count()<<std::endl;
	m_dElapsed = elapsed.count();
	// convert counts to time in seconds and return it
	return m_dElapsed;
}
bool CElapsed::Available()  // returns true if the perf counter is available
{ 
	return m_iInitialized; 
}

long long int CElapsed::GetFreq() // return perf counter frequency as large int
{
	return m_iFrequency; 
}