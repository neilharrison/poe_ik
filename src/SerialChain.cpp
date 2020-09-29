#include "SerialChain.h"
#include <iostream>
CSerialChain::CSerialChain(t_kinParam* lp, VectorXd q):legParam(lp), jointAngles(q)
{
	dof = jointAngles.size();
	if(dof != legParam->extDenHart.rows())
		std::cout << "Number of joints mismatch" << std::endl;

}

CSerialChain::~CSerialChain(void)
{
	// delete[] link;
}
