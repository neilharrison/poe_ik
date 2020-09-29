#if !defined __PHV_PKM__
#define __PHV_PKM__

#pragma once

#include "SerialChain.h"


class CPkm
{
public:
	CPkm(void);
	~CPkm(void);

//variables
	char Name[16];				// PKM name
	int numLegs;				// number of indipendent serial chains
	Matrix4d base;					// PKM base matrix;
	Matrix4d tool;					// PKM tool matrix;
	Eigen::SparseMatrix<double> pathMatrix;		// PKM path matrix;

// Serial chain definitions
	std::vector<CSerialChain*> Leg;
	std::vector<CSerialChain*>::iterator Leg_it;
};
#endif // __PHV_PKM__
