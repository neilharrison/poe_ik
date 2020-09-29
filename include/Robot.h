#if !defined __PHV_ROBOT__
#define __PHV_ROBOT__

#pragma once

#include "Common.h"
#include "Link.h"

class CRobot
{
public:
	CRobot(void);
	~CRobot(void);

//variables
	int dof;		// robot Degrees Of Freedom
	char Name[16];	// robot name
	double* q;		// current joint angles
	Clink* link;	// robot link
	Matrix4d base;		// robot base matrix;
	Matrix4d tool;		// robot tool matrix;
};
#endif // __PHV_ROBOT__
