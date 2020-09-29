#if !defined __PHV_SERIALCHAIN__
#define __PHV_SERIALCHAIN__

#pragma once
#include "Common.h"
#include "Link.h"
	
class CSerialChain
{
public:
	CSerialChain(t_kinParam*, VectorXd);
	virtual ~CSerialChain(void);

	char Name[16];	// chain name

//structure variables
	int dof;				// chain Degrees Of Freedom
	t_kinParam* legParam;   // leg kinematic parameters
	//Matd extDenHart;		// extended Denavit-Hartenberg parameters (8 elements)
	//Mat4d baseFrame;		// chain base reference frame
	//Mat4d toolFrame;		// chain tool reference frame
	
// state variables
	VectorXd jointAngles;		// current joint angles
	Matrix4d Ton;              // base-tool transformation;
	//Clink* link;	// chain link

};
#endif // __PHV_SERIALCHAIN__