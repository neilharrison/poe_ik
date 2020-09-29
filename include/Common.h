#if !defined __COMMON__
#define __COMMON__
#pragma once
#ifdef NDEBUG
	#undef _DEBUG
#endif
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
//#include "vl/VLfd.h"
#include <Eigen/Eigen>

#include <string>
#include <stdlib.h>

#include <stdio.h>
//#include <tchar.h>


using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;

//typedef Mat4d test4;


struct t_kinParam{
	MatrixXd extDenHart;		// extended Denavit-Hartenberg parameters (8 elements)
	Matrix4d baseFrame;		// chain base reference frame
	Matrix4d toolFrame;		// chain tool reference frame
};


#endif // __COMMON__