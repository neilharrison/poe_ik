// VerneKinByPoE.cpp : Defines the entry point for the console application.
//
#define POE_BLOCK
#define _DEBUG

#include <iostream>
#include <fstream>
//#include <windows.h>
#include "VerneKinByPoE.h"

// SVD stuff
//#include "./SVD/ap.h"
//#include "./SVD/reflections.h"
//#include "./SVD/bidiagonal.h"
//#include "./SVD/qr.h"
//#include "./SVD/lq.h"
//#include "./SVD/blas.h"
//#include "./SVD/rotations.h"
//#include "./SVD/bdsvd.h"
#include "./SVD/svd.h"

#include <Eigen/Eigen>
#include <algorithm>



#include <sstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;



// #define NUMB_OF_LEGS 4
// #define NUMB_OF_JOINTS 23
#define NUMB_OF_LEGS 1
#define NUMB_OF_JOINTS 6

// Denavit-Hartenberg parameters are stored in (n, 7) matrices,
// where rows are structured in this manner:
//
// [ alpha, A, theta, D, sigma, mdh, offset]
//
// alpha: link twist angle
// A: link length
// theta: initial link rotation angle (variable if joint is revolute)
// D: initial link offset distance (variable if joint is prismatic)
// sigma: 0 if joint is revolute, 1 if is prismatic
// mdh: 0 if DH notation is standard, 1 if is modified (Craig notation)
// offset: joint coordinate offset

// #define PATH_MATRIX    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
// 					   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
// 					   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
// 					   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
#define PATH_MATRIX 	1.0,1.0,1.0,1.0,1.0,1.0
					

#define DH_PARAM_11    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,        \
					  -M_PI_2, 0.0, 0.0, 0.08, 0.0, 1.0, 0.0,    \
					  -M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, -M_PI_2, \
					   M_PI_2, 0.0, 0.0, 0.85, 0.0, 1.0, 0.0,    \
					   M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, -M_PI_2, \
					   M_PI_2, 0.0, 0.0, 0.13, 0.0, 1.0, 0.0

#define DH_PARAM_12   0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,        \
					 -M_PI_2, 0.0, 0.0, -0.08, 0.0, 1.0, 0.0,   \
					 -M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, -M_PI_2, \
					  M_PI_2, 0.0, 0.0, 0.85, 0.0, 1.0, 0.0,    \
					  M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, -M_PI_2, \
					  M_PI_2, 0.0, 0.0, -0.13, 0.0, 1.0, 0.0

#define DH_PARAM_II   0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,       \
					  0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,       \
					 -M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, M_PI_2, \
					  M_PI_2, 0.0, 0.0, 0.95, 0.0, 1.0, M_PI,  \
					  M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, M_PI_2, \
					  M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0

#define DH_PARAM_III   0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,       \
					   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,       \
					  -M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, M_PI_2, \
					   M_PI_2, 0.0, 0.0, 0.95, 0.0, 1.0, M_PI,  \
					   M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, M_PI_2, \
					   M_PI_2, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0

#define DH_PARAM_UR5   M_PI_2, 0.0,   0.0, 0.089159, 0.0, 0.0, 0.0,       \
					   0.0,  -0.425,  0.0, 0.0,      0.0, 0.0, 0.0,       \
					   0.0, -0.39225, 0.0, 0.0,      0.0, 0.0, 0.0, \
					   M_PI_2, 0.0,   0.0, 0.10915,  0.0, 0.0, 0.0,  \
					  -M_PI_2, 0.0,   0.0, 0.09465,  0.0, 0.0, 0.0, \
					   0.0,    0.0,   0.0, 0.0823,   0.0, 0.0, 0.0


//--------------------------------------------------------------
// Baase and (initial) tool transformation matrices for each leg
//--------------------------------------------------------------

#define BASE_11     1.0, 0.0, 0.0, 0.55, \
					0.0, 1.0, 0.0, 0.0,   \
					0.0, 0.0, 1.0, 1.2859, \
					0.0, 0.0, 0.0, 1.0

#define TOOL_11		1.0, 0.0, 0.0, -0.23,  \
					0.0, 0.0, 1.0, 0.0,  \
					0.0, -1.0, 0.0, 0.0,  \
					0.0, 0.0, 0.0, 1.0

#define BASE_12		1.0, 0.0, 0.0, 0.55, \
					0.0, 1.0, 0.0, 0.0,   \
					0.0, 0.0, 1.0, 1.2859, \
					0.0, 0.0, 0.0, 1.0

#define TOOL_12		1.0, 0.0, 0.0, -0.23,  \
					0.0, 0.0, 1.0, 0.0,  \
					0.0, -1.0, 0.0, 0.0,  \
					0.0, 0.0, 0.0, 1.0

#define BASE_II		1.0, 0.0, 0.0, -0.285,  \
					0.0, 1.0, 0.0, -0.565,  \
					0.0, 0.0, 1.0, 1.3392,   \
					0.0, 0.0, 0.0, 1.0

#define TOOL_II		1.0, 0.0, 0.0, 0.045, \
					0.0, 1.0, 0.0, 0.19,  \
					0.0, 0.0, 1.0, 0.0,  \
					0.0, 0.0, 0.0, 1.0

#define BASE_III	1.0, 0.0, 0.0, -0.285,  \
					0.0, 1.0, 0.0, 0.565,  \
					0.0, 0.0, 1.0, 1.3392,   \
					0.0, 0.0, 0.0, 1.0

#define TOOL_III	1.0, 0.0, 0.0, 0.045, \
					0.0, 1.0, 0.0, -0.19,  \
					0.0, 0.0, 1.0, 0.0,  \
					0.0, 0.0, 0.0, 1.0

#define BASE_UR5	1.0, 0.0, 0.0, 0.0,  \
					0.0, 1.0, 0.0, 0.0,  \
					0.0, 0.0, 1.0, 0.0,   \
					0.0, 0.0, 0.0, 1.0

#define TOOL_UR5	1.0, 0.0, 0.0, 0.0, \
					0.0, 1.0, 0.0, 0.0,  \
					0.0, 0.0, 1.0, 0.0,  \
					0.0, 0.0, 0.0, 1.0

// Initial joint values for each leg
#define LEG_11_INIT_CONF 0.0, -1.1841, 0.0589, 0.0, 0.0589, -1.1841
#define LEG_12_INIT_CONF 0.0, -1.1841, -0.0589, 0.0, -0.0589, -1.1841
#define LEG_II_INIT_CONF 0.0, 1.0015, 1.0830, 0.0, -1.0830, -1.0015
#define LEG_III_INIT_CONF 0.0, -1.0015, 1.0830, 0.0, -1.0830, 1.0015

#define LEG_UR5_INIT_CONF 0.0,-M_PI_2, M_PI_2, 0.0, M_PI_2, 0.0

// starting trajectory at point 777
//#define LEG_11_INIT_CONF -0.4365, -0.5948, -0.5613, 0.3664, -0.2922, -0.7492
//#define LEG_12_INIT_CONF -0.4365, -0.4603, -0.6719, 0.4218, -0.4363, -0.6862
//#define LEG_II_INIT_CONF 0.1836, 4.3627, 1.4600, 2.6811, -1.2458, -1.6898
//#define LEG_III_INIT_CONF -0.5689, -1.6116, 0.3762, 0.0243, -0.8027, 1.6254

#define TRAJECTORY_FILE_NAME "VERNE_test_th_1mm.txt"
#define MAX_ITERATION_NUMBER 200
#define ITERATION_TOL 1.0e-4

int POINT_CNT  = 0;  // total number of trajectory points

MatrixXd test_trajectory(MatrixXd::Zero(3500,4)); // TCP points

// path matrix
MatrixXd path_matrix((MatrixXd(NUMB_OF_LEGS, NUMB_OF_JOINTS)<<PATH_MATRIX).finished());	
//<< NUMB_OF_LEGS, NUMB_OF_JOINTS
// joint angles
MatrixXd joint_angles_11(MatrixXd::Zero(3500,6));
MatrixXd joint_angles_12(MatrixXd::Zero(3500,6));
MatrixXd joint_angles_II(MatrixXd::Zero(3500,6));
MatrixXd joint_angles_III(MatrixXd::Zero(3500,6));


// Create matrices for DH parameters for each leg
MatrixXd DH_11((MatrixXd(6,7)<<DH_PARAM_11).finished());
MatrixXd DH_12((MatrixXd(6,7)<<DH_PARAM_12).finished());
MatrixXd DH_II((MatrixXd(6,7)<<DH_PARAM_II).finished());
MatrixXd DH_III((MatrixXd(6,7)<<DH_PARAM_III).finished());

MatrixXd DH_UR5((MatrixXd(6,7)<<DH_PARAM_UR5).finished());



// Base and Tool transformation for each leg
Matrix4d Base_11((Matrix4d()<<BASE_11).finished());
Matrix4d Tool_11((Matrix4d()<<TOOL_11).finished());
Matrix4d Base_12((Matrix4d()<<BASE_12).finished());
Matrix4d Tool_12((Matrix4d()<<TOOL_12).finished());
Matrix4d Base_II((Matrix4d()<<BASE_II).finished());
Matrix4d Tool_II((Matrix4d()<<TOOL_II).finished());
Matrix4d Base_III((Matrix4d()<<BASE_III).finished());
Matrix4d Tool_III((Matrix4d()<<TOOL_III).finished());

Matrix4d Base_UR5((Matrix4d()<<BASE_UR5).finished());
Matrix4d Tool_UR5((Matrix4d()<<TOOL_UR5).finished());


// Initial configuration
VectorXd current_joint_angles_11((VectorXd(6)<< LEG_11_INIT_CONF).finished());
VectorXd current_joint_angles_12((VectorXd(6)<< LEG_12_INIT_CONF).finished());
VectorXd current_joint_angles_II((VectorXd(6)<< LEG_II_INIT_CONF).finished());
VectorXd current_joint_angles_III((VectorXd(6)<< LEG_III_INIT_CONF).finished());

VectorXd current_joint_angles_UR5((VectorXd(6)<< LEG_UR5_INIT_CONF).finished());



int ToString(char *dest, const char *fmt, double x, int count )
{
	int ret = snprintf(dest, count, fmt, x );
	return ret;
}

void saveTestData (const char *fileName, MatrixXd &data)
{
	char dest[15];
	char *d = &dest[0];

	
	std::string line;
	std::ofstream data_file;

	data_file.open(fileName);
	if(data_file.fail())            //is it ok?
   { std::cout<<"Input file did not open please check it\n";}

	int numCols = data.cols();

	/* Open file in text mode: */
 	if(data_file.is_open())  {
		for(int i = 0; i < POINT_CNT; i++) {
			for(int j = 0; j < numCols; j++)  {
				//std::cout<< (data(i,j)) <<std::endl;
				data_file  << (data(i,j));
				data_file << "  ";

			}
			data_file<< "\n";

		}
	}
	else
		std::cout << "Test data file not opened\n" << std::endl;
	data_file.close();
}

void Initialize()
{
	// store initial joint angles
	joint_angles_11.row(0) = current_joint_angles_11;
	joint_angles_12.row(0) = current_joint_angles_12;
	joint_angles_II.row(0) = current_joint_angles_II;
	joint_angles_III.row(0) = current_joint_angles_III;

	/*
	std::string line;
	std::ifstream ftraj;

	ftraj.open(TRAJECTORY_FILE_NAME);

	if(ftraj.fail())            //is it ok?
   		{ std::cout<<"Input file did not open please check it\n";}

	int column_cnt = 0;
	std::string::size_type sz;
	std::string delimiter = "  ";
	std::string token;

	// Cycle until end of file reached:
	while ( getline (ftraj,line) ) {
		column_cnt=0;
		size_t pos = 0;
		while ((pos = line.find(delimiter)) != std::string::npos) {
    		token = line.substr(0, pos);
    		//std::cout << std::stod(token, &sz) << std::endl;
			test_trajectory(POINT_CNT,column_cnt) = std::stod(line, &sz);
			column_cnt++;
    		line.erase(0, pos + delimiter.length());
		}
		//std::cout<<test_trajectory.row(POINT_CNT)<<std::endl;
		POINT_CNT++;
	}
	ftraj.close();

	// reverse the trajectory if points are too few
	if(POINT_CNT < 40) {
		for(int i = 0; i < POINT_CNT; i++)
			for(int j = 0; j < 4; j++)
				test_trajectory(POINT_CNT + i,j) = test_trajectory(POINT_CNT - i - 1,j);
		POINT_CNT += POINT_CNT;
	}
	*/
}

MatrixXd PreCalcDHparams(MatrixXd &dh_params)
{
	MatrixXd dh_precalc(MatrixXd::Zero(6,8));
	for(int i = 0; i < 6; i++)  {
		dh_precalc(i,0) = cos(dh_params(i,0));  // calculate cos(alpha)
		dh_precalc(i,1) = sin(dh_params(i,0));  // calculate sin(alpha)
		for(int j = 2; j < 8; j++)  
			dh_precalc(i,j) = dh_params(i,j-1);
	}
	return dh_precalc;
}



void SVDecomposition(MatrixXd &A, MatrixXd &U, MatrixXd&V, VectorXd &s)
{
	int m = U.rows();
	int n = U.cols();
    ap::real_2d_array a;
    ap::real_2d_array u;
    ap::real_2d_array vt;
    ap::real_1d_array w;
    
	a.setbounds(1, m, 1, n);
	for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)
			a(i+1,j+1) = A(i,j);

	if( !svddecomposition(a, m, n, 2, 2, 2, w, u, vt) )
    {
        printf("SVD decomposition failed!\n");
		exit(1);
    }

	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			U(i,j) = u(i+1,j+1);
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			V(i,j) = vt(j+1,i+1);
	for(int i = 0; i < n; i++)
		s(i) = w(i+1);

}

// Calculation of pseudo-inverse using SVD decomposition
MatrixXd pinv(MatrixXd M)
{
	double tol;
	if(M.cols() > M.rows())
		M.transposeInPlace();
	int m = M.rows();
	int n = M.cols();
	MatrixXd U(MatrixXd::Zero(m,n));
	MatrixXd V(MatrixXd::Zero(n,n));
	VectorXd s(VectorXd::Zero(n));
	MatrixXd X;

	cout<<M<<endl;
	SVDecomposition(M, U, V, s);
	cout<<"U"<<U<<endl;
	cout<<"V"<<V<<endl;
	cout<<"s"<<s<<endl;
	// Find the subvector of s with elements > tol
	// values of s are in decreasing order
	tol = m * s(0) * EPS;
	int r = 0;
	while(s(r) > tol && r < n-1)
		r++;
	// if r != 0 invert all element of s,
	// create a diagonal matrix and multiply by V and U
	if(r)  {
		MatrixXd S(MatrixXd::Zero(r, r));
		for(int i = 0; i < r; i++)
			S(i,i) = 1.0/s(i);
		X = V.block(0,0,n,r) * (S*((U.block(0,0,m,r)).transpose()));
	}
	else  {
		X.Zero(n, m);
	}

	return X;
}

// Calculation of pseudo-inverse using SVD factorization
// Neil didnt bother updating

// Matd pinv2(Matd &M)
// {
// 	double tol;
// 	if(M.Cols() > M.Rows())
// 		M = trans(M);
// 	int m = M.Rows();
// 	int n = M.Cols();
// 	Matd U(m, n, vl_0);
// 	Matd V(n ,n, vl_0);
// 	Vecd s(n, vl_0);
// 	Matd X;
	
// 	/// THIS OPTIONAL ALGORITHM IS SLOWER
// 	SVDFactorization(M, U, V, s);
// 	// find max value of s; this a type of norm(M)
// 	double max_s = s[0];
// 	for(int i = 0; i < n; i++)
// 		max_s = max(max_s, s[i]);
// 	tol = m * max_s * EPS;
// 	Matd S(n, n, vl_0);
// 	for(int i = 0; i < n; i++)
// 		if(s[i] > tol)
// 			S[i][i] = 1.0/s[i];
// 		else {
// 			cout << "SMALL s" << endl;
// 			X.SetSize(n, m);
// 			X = vl_0;
// 		}
// 	X = V * (S * trans(U));

// 	return X;
// }


Matrix4d fkine(MatrixXd &dh_param, Matrix4d &base_frame, Matrix4d &tool_frame, VectorXd &current_joint_angles)
{
	Matrix4d base_tool_transform(Matrix4d::Zero());
	Matrix4d temp_matrix(Matrix4d::Zero());
	temp_matrix = base_frame;
	
	for(int i = 0; i < 6; i++)  { //6 should be variable? (number of links in leg?)
		//base_tool_transform = temp_matrix * linkTransformMod(dh_param.row(i), current_joint_angles(i));
		base_tool_transform = temp_matrix * linkTransform(dh_param.row(i), current_joint_angles(i));
		temp_matrix = base_tool_transform;
	}
	base_tool_transform = temp_matrix * tool_frame;
	//cout<<base_tool_transform<<endl;
	return base_tool_transform;
}

MatrixXd BstarS(MatrixXd &dh, Matrix4d &baseFrame, VectorXd &q, int &leg_index)
{
	int dof = dh.rows();
	
	// Calculates B*S and replace as "B" matrix 
	// Vectors are built up as row vectors
	// B need to be transposed before being used
	MatrixXd B(MatrixXd::Zero(VERNE.pathMatrix.cols(), dof));
	VectorXd s(VectorXd::Zero(6)); // Twist local coordinates
	Matrix4d Toi = baseFrame;
	
	int i = 0;
	//SVIterd j_it;
	//for (j_it.Begin(VERNE.pathMatrix.row(leg_index); !j_it.AtEnd(); j_it.Inc()) {
	Eigen::SparseVector<double> pm_row= VERNE.pathMatrix.row(leg_index);
	
	for (Eigen::SparseVector<double>::InnerIterator it(pm_row);it;++it){
		//cout<<it.index()<<endl;
		//cout<<dq_all(it.index())<<endl;		
		//Toi = Toi * linkTransformMod(dh.row(i), q(i));
		
		s(2) = 1.0 - dh(i,5);
		s(5) = dh(i,5);
		if (i ==0) {
			B.row(it.index()) = s;
		}
		else {
			Toi = Toi * linkTransform(dh.row(i-1), q(i-1));
			B.row(it.index()) = Adjnt(Toi)*s;
		}
		//cout<<it.index()<<endl;
		//cout<<B.row(it.index())<<endl;
		i++;
	}

	//cout << B.transpose() << endl;
	return B.transpose(); //*S;
}

VectorXd poseDiff(Eigen::Ref<Eigen::Matrix3d> Rnew, Eigen::Ref<Eigen::Matrix3d> R, Eigen::Ref<Eigen::Vector3d> xyz_new, Eigen::Ref<Eigen::Vector3d> xyz)
{
	Matrix3d Rtrans = R.transpose();
	//cout<<R<<endl;
	//cout<<xyz<<endl;
	// Calculate inv(Ton)*Ton_new
	Matrix3d Rtot = Rtrans * Rnew;
	Vector3d xyz_tot = xyz_new - xyz;
	xyz_tot = Rtrans * xyz_tot;

	Matrix3d omega(Matrix3d::Zero());
	Vector3d omg;
	double theta;
	double s_theta;
	/////////// Calculate log(inv(Ton)*Ton_new)
	double tr_Rtot = Rtot.trace();
	if (tr_Rtot >= 3){
		// its already zero  /// omega = Matrix3d::Zero();
		//VectorXd D = (Vectorxd(6)<<0,0,0,xyz_tot).finished();
		return (VectorXd(6)<<0,0,0,xyz_tot).finished();
	}
	else if(tr_Rtot <= -1)  {
        if ((1+Rtot(2,2))!=0) {

			omg = (1 / sqrt(2 * (1 + Rtot(2, 2)))) * (Vector3d()<<Rtot(0, 2), Rtot(1, 2), 1 + Rtot(2, 2)).finished()*M_PI;
		}
		else if ((1+Rtot(1,1))!=0) {
			omg = (1 / sqrt(2 * (1 + Rtot(1, 1)))) * (Vector3d()<<Rtot(0, 1), 1 + Rtot(1, 1), Rtot(2, 1)).finished()*M_PI;
		}
		else {
			omg = (1 / sqrt(2 * (1 + Rtot(0, 0)))) * (Vector3d()<<1 + Rtot(0, 0), Rtot(1, 0), Rtot(2, 0)).finished()*M_PI;
		}
	
		omega = (Matrix3d()<<0,-omg(2),omg(1),omg(2),0,-omg(0),-omg(1),omg(0),0).finished();
		theta = M_PI;
		s_theta = sin(theta);
	//TO do - check what this does to bottom of A fraction ->0 -> inf??
	}
	else{
		theta = acos((tr_Rtot-1.0)/2.0);
		s_theta = sin(theta);
        omega = (theta/2.0/s_theta)*(Rtot-Rtot.transpose());
	}
	
    
	Matrix3d A(Matrix3d::Identity());
    // Calculate the adjoint rotational component of log(inv(Ton)*Ton_new)
    // omega = log(Rtot) = (fi/2sin(fi))(Rtot-Rtot')
	//if(abs(fi) > EPS) {
	
        // Calculate the adjoint translational component A*d of log(inv(Ton)*Ton_new)
    A = A - (0.5 * omega) + (2.0*s_theta-theta*(1+cos(theta)))/(2.0*theta*theta*s_theta) * (omega * omega);
	//}

    Vector3d v = A * xyz_tot;
    
    // Adjoint representation of log(inv(Ton)*Ton_new)
    Matrix4d Adj_log(Matrix4d::Zero());
	for(int i = 0; i < 3; i++)  {
		for(int j = 0; j < 3; j++)
			Adj_log(i,j) = omega(i,j);
			Adj_log(i,3) = v(i);
	}

	// Rearrange Adj_log in a 6x1 vector (pose difference vector)
	VectorXd D = quickTr2Diff(Adj_log);
	//cout<<"D\n"<<D<<endl;
	return D;

}

// Calculation of Body Manipulator Jacobian
MatrixXd bm_jacobian(MatrixXd &dh, Matrix4d &baseFrame, VectorXd &q, Eigen::Ref<Eigen::Matrix3d> &Ron, Eigen::Ref<Eigen::Vector3d> &xyz)
{
	Matrix3d Ron_trans = Ron.transpose();
	int dof = dh.rows();
	MatrixXd J(MatrixXd::Zero(dof, 6));  // this is the transpose of J, convenient for calcs

	// Build Delta
	Matrix3d Dlt((Matrix3d()<< 0.0, -xyz(2), xyz(1),
				   xyz(2), 0.0, -xyz(0),
				  -xyz(1), xyz(0), 0.0).finished());

	Dlt = -Ron_trans * Dlt;
			   
	MatrixXd inv_Adj_Ton((MatrixXd(6,6)<< Ron_trans(0,0), Ron_trans(0,1), Ron_trans(0,2), 0.0, 0.0, 0.0,
						   Ron_trans(1,0), Ron_trans(1,1), Ron_trans(1,2), 0.0, 0.0, 0.0,
						   Ron_trans(2,0), Ron_trans(2,1), Ron_trans(2,2), 0.0, 0.0, 0.0,
						   Dlt(0,0), Dlt(0,1), Dlt(0,2), Ron_trans(0,0), Ron_trans(0,1), Ron_trans(0,2),
						   Dlt(1,0), Dlt(1,1), Dlt(1,2), Ron_trans(1,0), Ron_trans(1,1), Ron_trans(1,2),
						   Dlt(2,0), Dlt(2,1), Dlt(2,2), Ron_trans(2,0), Ron_trans(2,1), Ron_trans(2,2)).finished());
	
	// Build the Adjoint matrix for each link
	Matrix4d Toi = baseFrame ;//new Matrix4d;
	Eigen::Ref<Matrix3d> Roi = Toi.block(0,0,3,3);

	VectorXd s(VectorXd::Zero(6));
	Matrix3d Dlt_oi(Matrix3d::Zero());
	MatrixXd Adj_Toi(MatrixXd::Zero(6,6));

	for(int i = 0; i < dof; i++)  {
		Toi =  Toi * linkTransformMod(dh.row(i), q(i));
		//cout << Toi << endl;
		s(2) = 1.0 - dh(i,5);
		s(5) = dh(i,5);
		// Build Delta matrix for link i
		Dlt_oi(0,0) = 0.0; Dlt_oi(0,1) = -Toi(2,3); Dlt_oi(0,2) = Toi(1,3);
		Dlt_oi(1,0) = Toi(2,3); Dlt_oi(1,1) = 0.0; Dlt_oi(1,2) = -Toi(0,3),
		Dlt_oi(2,0) = -Toi(1,3); Dlt_oi(2,1) = Toi(0,3); Dlt_oi(2,2) = 0.0;
		Dlt_oi = Dlt_oi * Roi;
		//cout << Dlt_oi << endl;
		//cout << "\n" << endl;
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++)  {
				Adj_Toi(j,k) = Roi(j,k);
				Adj_Toi(j+3,k) = Dlt_oi(j,k);
				Adj_Toi(j+3,k+3) = Roi(j,k);
			}
		J.row(i) = inv_Adj_Ton * (Adj_Toi * s);
	}

	return J.transpose();
}

VectorXd ikinVERNE_POE(MatrixXd &dh, Matrix4d &baseFrame, Eigen::Ref<Eigen::Matrix3d> &Ron_new, Eigen::Ref<Eigen::Matrix3d> &Ron, Eigen::Ref<Eigen::Vector3d> &xyz_new, Eigen::Ref<Eigen::Vector3d> &xyz, VectorXd &q)
{
	Matrix3d Ron_trans = Ron.transpose();
	
	// Calculate inv(Ton)*Ton_new
	Matrix3d Rtot = Ron_trans * Ron_new;
	Vector3d xyz_tot = xyz_new - xyz;
	xyz_tot = Ron_trans * xyz_tot;

	/////////// Calculate log(inv(Ton)*Ton_new)
	double tr_Rtot = Rtot.trace();
	if(tr_Rtot == -1)  { //modern robotics pg 85 - singularity like euler
        std::cout << "Trace of Rtot = -1" << std::endl;
        exit(2);
	}
	double fi = acos((tr_Rtot-1.0)/2.0);
    Matrix3d omega(Matrix3d::Zero());
	Matrix3d A(Matrix3d::Identity());
    // Calculate the adjoint rotational component of log(inv(Ton)*Ton_new)
    // omega = log(Rtot) = (fi/2sin(fi))(Rtot-Rtot')
	if(abs(fi) > 1e-17) {
		double s_fi = sin(fi);
		//MR - 1/2sin(theta)*(R-R^T) - why extra fi?
        omega = (fi/2.0/s_fi)*(Rtot-Rtot.transpose());
        // Calculate the adjoint translational component A*d of log(inv(Ton)*Ton_new)
        A = A - (0.5 * omega) + (2.0*s_fi-fi*(1+cos(fi)))/(2.0*fi*fi*s_fi) * (omega * omega);
	}
    //else
        //omega = zeros(3);
        //A= eye(3);
        //cout << "warning: NULL fi" << endl;
    
    Vector3d v = A * xyz_tot;
    //std::cout<<"v"<<v<<std::endl;
    // Adjoint representation of log(inv(Ton)*Ton_new)
    Matrix4d Adj_log(Matrix4d::Zero());
	for(int i = 0; i < 3; i++)  {
		for(int j = 0; j < 3; j++)
			Adj_log(i,j) = omega(i,j);
		Adj_log(i,3) = v(i);
	}
	//std::cout<<Adj_log<<std::endl;
	// Rearrange Adj_log in a 6x1 vector (pose difference vector)
	// body twist Vd? 

	//Neil trying to recreate///

	//log(inv(Ton)*Ton_new)
	// Matrix4d Ton; //This is dumb, remove it
	// for(int i = 0; i < 3; i++)  {
	// 	for(int j = 0; j < 3; j++)
	// 		Ton(i,j) = Rtot(i,j);
	// 	Ton(i,3) = xyz_tot(i);
	// }
	
	//std::cout<<MatrixLog6(Ton)<<std::endl;
	//Gives same as original! good?



	VectorXd D = quickTr2Diff(Adj_log);
	//std::cout<<D<<std::endl;
	// Calculate the body manipulator Jacobian
	MatrixXd J = bm_jacobian(dh, baseFrame, q, Ron, xyz);

	MatrixXd J_star = pinv(J);  // J pseudoinverse

	//Eigen pinv apparently unstable? better to use svd solve as below
	MatrixXd Eigenpinv = J.completeOrthogonalDecomposition().pseudoInverse();
#ifdef _DEBUG
	//std::cout << "D\n" << D << std::endl;
	//std::cout << "J\n" << J << std::endl;

	// Solve
    // Matd J_star = pinv(bm_jacobian(dh, baseFrame, q, Ron, xyz));  // J pseudoinverse
	//std::cout << "J_star\n" << J_star << std::endl;
	//std::cout<<"EigenPinv\n"<<Eigenpinv<<endl;
#else
	
	//Matd J = bm_jacobian(dh, baseFrame, q, Ron, xyz);
    //cout << "J\n" << J << endl;
	//Matd J_star = pinv(J);
	//Vecd dq = J_star * D;
	//exit(0);

#endif
	//std::cout<<q<<std::endl;
	//cout<<"D\n"<<D<<endl;
	
    // Update position
	//VectorXd dq = pinv(bm_jacobian(dh, baseFrame, q, Ron, xyz)) * D;
	//VectorXd dq = Eigenpinv * D;
	VectorXd dq = J.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(D);
	//cout<<"dq\n"<<dq<<endl;

	q += dq;
	
 	return D;
}

void iksolve(Matrix4d &Ton_new, int &n_it) { 

	MatrixXd dh;
	Matrix4d baseFrame;
	Matrix4d toolFrame;

	Matrix4d Ton(Matrix4d::Identity());
	// extract rotation matrices
	Eigen::Ref<Matrix3d> Ron_new = Ton_new.block(0,0, 3, 3);   
	Eigen::Ref<Matrix3d> Ron = Ton.block(0,0, 3, 3);
	
	Eigen::Ref<Vector3d> xyz = Ton.col(3).segment(0,3);
	Eigen::Ref<Vector3d> xyz_new = Ton_new.col(3).segment(0,3);

	// body manipulator jacobian 
	// Should really define size here
	MatrixXd J;

	// Adjoint matrix of inv(Ton)
	MatrixXd A(MatrixXd::Zero(6*NUMB_OF_LEGS, 6*NUMB_OF_LEGS));
	// Path adjoint matrix (pre-multiplied by S)
	MatrixXd B(MatrixXd::Zero(6*VERNE.pathMatrix.rows(), VERNE.pathMatrix.cols()));
	// Pose difference vector
	VectorXd D(VectorXd::Zero(6*NUMB_OF_LEGS));
	// local temp variables
	VectorXd current_joint_angles;
	VectorXd dq_all(VectorXd::Zero(NUMB_OF_JOINTS));
	VectorXd dq_all_test(VectorXd::Zero(NUMB_OF_JOINTS));

	int leg_index,leg_dof;
	double error;

	error = 1.0;
	 n_it = 0;
	
	
	while(error > ITERATION_TOL)  {
		leg_index = 0;
		// prev_leg_dofs = 0;
		A = MatrixXd::Zero(6*NUMB_OF_LEGS, 6*NUMB_OF_LEGS);
		B = MatrixXd::Zero(6*VERNE.pathMatrix.rows(), VERNE.pathMatrix.cols());
		D = VectorXd::Zero(6*NUMB_OF_LEGS);
		for((VERNE.Leg_it) = VERNE.Leg.begin(); (VERNE.Leg_it) <  VERNE.Leg.end(); ++(VERNE.Leg_it))  {
			// acquire leg data
			Ton = (*(VERNE.Leg_it))->Ton;
			changeSubMatd(B, BstarS((*(VERNE.Leg_it))->legParam->extDenHart, (*(VERNE.Leg_it))->legParam->baseFrame, (*(VERNE.Leg_it))->jointAngles, leg_index),
							leg_index*6, 0, 6, VERNE.pathMatrix.cols());

			changeSubMatd(A, AdjInvH(Ton), leg_index*6, leg_index*6, 6, 6);
							
			changeSubVecd(D, poseDiff(Ron_new, Ron, xyz_new, xyz), leg_index*6, 6);
			
			leg_index++;
		}
		// calculate all joints increments
		
		J = A * B;
		//cout<<"J\n"<<J<<endl;
		//cout<<"D"<<D<<endl;

		//dq_all_test = pinv(J) * D;
		dq_all = J.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(D);
		
		//cout<<"dq\n"<<dq_all<<endl;
		
		leg_index = 0;
		
		// update joint angles for each leg
		error = 0.0;
		for((VERNE.Leg_it) = VERNE.Leg.begin(); (VERNE.Leg_it) <  VERNE.Leg.end(); ++(VERNE.Leg_it))  {
			leg_dof = (*(VERNE.Leg_it))->dof;
			
			// update joint angles
			
			Eigen::SparseVector<double> pm_row= VERNE.pathMatrix.row(leg_index);
			int j=0;
			for (Eigen::SparseVector<double>::InnerIterator it(pm_row);it;++it){
				//cout<<"joint angles"<<(*(VERNE.Leg_it))->jointAngles(j)<<endl;
				(*(VERNE.Leg_it))->jointAngles(j) += dq_all(it.index());
				//cout<<"joint angles"<<(*(VERNE.Leg_it))->jointAngles(j)<<endl;
				j++;
			}

			// update kinematics
			
			Ton = fkine((*(VERNE.Leg_it))->legParam->extDenHart, (*(VERNE.Leg_it))->legParam->baseFrame, (*(VERNE.Leg_it))->legParam->toolFrame, (*(VERNE.Leg_it))->jointAngles);
			(*(VERNE.Leg_it))->Ton = Ton;
			//cout<<Ton<<endl;
			changeSubVecd(D, poseDiff(Ron_new, Ron, xyz_new, xyz), leg_index*6, 6);

			// evaluate max tcp calculation error across all legs
			//tr2rpy(rpy, Ton);
			//tcp(0) = xyz(0); tcp(1) = xyz(1); tcp(2) = xyz(2);
			//tcp(3) = rpy(0); tcp(4) = rpy(1); tcp(5) = rpy(2);
			//double nm = (xyz_new - xyz).norm();
			//max_tcp_error = max(max_tcp_error, len(xyz_new - xyz));
			//double nm = len(tcp_new - tcp);
			//error = max(error, len(xyz_new - xyz));
			//max_tcp_error = max(max_tcp_error, len(tcp_new - tcp));
			leg_index++;
		}
		/*error = dq_all[0];
		for(int i = 1; i < D.Elts(); i++)
			error = max(error, dq_all[i]);*/

		if (n_it>MAX_ITERATION_NUMBER)
			{
				cout<<"Maybe stuck"<<endl;
				break;
			}
		error = D.norm();
		cout<<"error\n"<<error<<endl;
		//error = len(dq_all);
		//error = max_tcp_error;
		n_it++;
	}
}


#ifdef POE_BLOCK

int main(int argc, char **argv)
{
	
	//Setting up legs - add to pkm(VERNE) class
	t_kinParam leg_11, leg_12, leg_II, leg_III;

	// leg_11.extDenHart = PreCalcDHparams(DH_11);
	// leg_11.baseFrame = Base_11;
	// leg_11.toolFrame = Tool_11;
	// //PKM.push_back(new CSerialChain(&leg_11, current_joint_angles_11));
	// VERNE.Leg.push_back(new CSerialChain(&leg_11, current_joint_angles_11));

	// leg_12.extDenHart = PreCalcDHparams(DH_12);
	// leg_12.baseFrame = Base_12;
	// leg_12.toolFrame = Tool_12;
	// //PKM.push_back(new CSerialChain(&leg_12, current_joint_angles_12));
	// VERNE.Leg.push_back(new CSerialChain(&leg_12, current_joint_angles_12));

	// leg_II.extDenHart = PreCalcDHparams(DH_II);
	// leg_II.baseFrame = Base_II;
	// leg_II.toolFrame = Tool_II;
	// //PKM.push_back(new CSerialChain(&leg_II, current_joint_angles_II));
	// VERNE.Leg.push_back(new CSerialChain(&leg_II, current_joint_angles_II));
	// leg_III.extDenHart = PreCalcDHparams(DH_III);
	// leg_III.baseFrame = Base_III;
	// leg_III.toolFrame = Tool_III;
	// //PKM.push_back(new CSerialChain(&leg_III, current_joint_angles_III));
	// VERNE.Leg.push_back(new CSerialChain(&leg_III, current_joint_angles_III));

	t_kinParam leg_UR5;

	leg_UR5.extDenHart = PreCalcDHparams(DH_UR5);
	leg_UR5.baseFrame = Base_UR5;
	leg_UR5.toolFrame = Tool_UR5;
	// //PKM.push_back(new CSerialChain(&leg_11, current_joint_angles_11));
	VERNE.Leg.push_back(new CSerialChain(&leg_UR5, current_joint_angles_UR5));


	//path matrix to sparse matrix inside pkm class
	VERNE.pathMatrix.resize(path_matrix.rows(), path_matrix.cols());
	// // Set path matrix as VERNE sparse matrix
	 for(int i = 0; i < path_matrix.rows(); i++){
	 	for(int j = 0; j < path_matrix.cols(); j++){
		 	if(path_matrix(i,j)!=0) {
	 		VERNE.pathMatrix.insert(i,j) = path_matrix(i,j);
			 }
		 }
	 }


	int max_it = MAX_ITERATION_NUMBER;
	Initialize();

	Matrix4d Ton_new(Matrix4d::Identity());
	Matrix4d Ton(Matrix4d::Identity());
	// extract rotation matrices
	Eigen::Ref<Matrix3d> Ron_new = Ton_new.block(0,0, 3, 3);   
	Eigen::Ref<Matrix3d> Ron = Ton.block(0,0, 3, 3);
	// extract translational vectors
	Eigen::Ref<Vector3d> xyz = Ton.col(3).segment(0,3);
	Eigen::Ref<Vector3d> xyz_new = Ton_new.col(3).segment(0,3);
	

	Vector3d rpy_new(Vector3d::Zero());
	Vector3d rpy(Vector3d::Zero());
	VectorXd tcp_new(VectorXd::Zero(6));
	VectorXd tcp(VectorXd::Zero(6));
	
	
	// data storage vectors and matrices
	// [calc error, number of iterations, elapsed time]
	MatrixXd test_data_leg11(MatrixXd::Zero(POINT_CNT, 3));
	MatrixXd test_data_leg12(MatrixXd::Zero(POINT_CNT, 3));
	MatrixXd test_data_legII(MatrixXd::Zero(POINT_CNT, 3));
	MatrixXd test_data_legIII(MatrixXd::Zero(POINT_CNT, 3));
	MatrixXd test_data_global(MatrixXd::Zero(POINT_CNT, 3));
	MatrixXd joint_values(MatrixXd::Zero(POINT_CNT, NUMB_OF_JOINTS));


	int prev_leg_dofs, leg_dof;
	int leg_index, n_it;
	double error, max_tcp_error;

	int point_cnt = 0;

	// Initialise forward kinematics
	for(VERNE.Leg_it = VERNE.Leg.begin(); VERNE.Leg_it <  VERNE.Leg.end(); ++(VERNE.Leg_it))  {
		// acquire leg data
		
		Ton = fkine((*(VERNE.Leg_it))->legParam->extDenHart, (*(VERNE.Leg_it))->legParam->baseFrame, (*(VERNE.Leg_it))->legParam->toolFrame, (*(VERNE.Leg_it))->jointAngles);
		(*(VERNE.Leg_it))->Ton = Ton;
	}
	homeTon = Ton;
	cout<<"Home\n"<<homeTon<<endl;
	ros::init(argc, argv, "poe_ik_node");
	//ros::NodeHandle n;

	//ros::Publisher joint_angles_pub = n.advertise<std_msgs::Float64MultiArray>("/arm/joint_group_position_controller/command", 1000);
	//ros::Rate loop_rate(10);

	//ros::Subscriber sub = n.subscribe("/ik_goal", 1000, ikCallback);
	SubscribeAndPublish subpub;
	ros::spin();


	// }
	// main loop
	for(int i = 0; i < POINT_CNT; i++)  {
		point_cnt++;
		// cout << "move count: " << point_cnt << "\n" << endl;
		// set new tcp point
		rpy_new(2) = test_trajectory(i,3);
		// First update Ton_new with new rpy angles
		rpy2tr(Ton_new, rpy_new);  
		// then update Ton_new with new xyz values
		tcp_new(0) = xyz_new(0) = test_trajectory(i,0);
		tcp_new(1) = xyz_new(1) = test_trajectory(i,1);
		tcp_new(2) = xyz_new(2) = test_trajectory(i,2);
		tcp_new(3) = rpy_new(0); tcp_new(4) = rpy_new(1); tcp_new(5) = rpy_new(2);

		m_el.Begin(); // start calc of elapsed time here

		//Function solves ik to move to Ton_new and updates joint angles in class
		Ton_new = (Matrix4d()<<  0.0000  , -1.0000  , -0.0000   , 0.0946,
   -1.0000   ,      0  , -0.0000 ,  -0.1091,
    0.0000   , 0.0000 ,  -1.0000  ,  0.0396,
         0     ,    0      ,   0 ,   1.0000).finished();
		cout<<Ton_new<<endl;

		Vector4d quat = (Vector4d()<<0,0,0,1).finished();
		Vector3d pose = (Vector3d()<<0.0856,-0.1091,0.0396).finished();

		Ton_new = transformFromPose(quat, pose);

		iksolve(Ton_new,n_it);

		m_el.End(); // end of elapsed time
		cout<<n_it<<endl;
		// evaluate max tcp calculation error across all legs
		max_tcp_error = 0.0;
		leg_index = 0;
		for((VERNE.Leg_it) = VERNE.Leg.begin(); (VERNE.Leg_it) <  VERNE.Leg.end(); ++(VERNE.Leg_it))  {
			Ton = (*(VERNE.Leg_it))->Ton;
			// tr2rpy(rpy, Ton);
			/*tcp[0] = xyz[0]; tcp[1] = xyz[1]; tcp[2] = xyz[2];
			tcp[5] = rpy[2];*/
			max_tcp_error = std::max(max_tcp_error, (xyz_new - xyz).norm());
			// store joint angles
			int j = 0;
			//for(int j = 0; j < leg_dof ; j++){
			Eigen::SparseVector<double> pm_row= VERNE.pathMatrix.row(leg_index);
			for (Eigen::SparseVector<double>::InnerIterator it(pm_row);it;++it){
				//	cout<<it.index()<<endl;
					//cout<<dq_all(it.index())<<endl;
				if((*(VERNE.Leg_it))->legParam->extDenHart(j,5))
					//joint_values(i,j_it.Index()) = 1000.0 * (*(VERNE.Leg_it))->jointAngles(j);
					joint_values(i,it.index()) = 1000.0 * (*(VERNE.Leg_it))->jointAngles(j);

				else
					//joint_values(i,j_it.Index()) = RAD_TO_DEG * (*(VERNE.Leg_it))->jointAngles(j);
					joint_values(i,it.index()) = RAD_TO_DEG * (*(VERNE.Leg_it))->jointAngles(j);
				j++;
			}
			leg_index++;
		}
		
		// store global test data
		test_data_global(i,0) = max_tcp_error;
		test_data_global(i,1) = n_it;
		test_data_global(i,2) = m_el.m_dElapsed;
		//cout<<joint_values.row(i)<<endl;
	}
	saveTestData("testData.txt", test_data_global);
	saveTestData("jointValues.txt", joint_values);
	
	//out = dest + ' ';
	//cout << out << dest << endl;
	//num_char = ToString(dest, "%e", test_data_global[0][0], 25);
	//cout << num_char << endl;
	//cout << dest << endl;
	//num_char = ToString(dest, "%e", test_data_global[0][1], 25);
	//cout << num_char << endl;
	//cout << dest << endl;
	//num_char = ToString(d, "%e", test_data_global[0][2], 15);
	//cout << num_char << endl;
	//cout << dest << endl;
	
}

#else

int main()
{
	std::cout << "Please change program priority to Realtime from Task Manager" << std::endl;
	//system("pause");

	t_kinParam leg_11, leg_12, leg_II, leg_III;

	leg_11.extDenHart = PreCalcDHparams(DH_11);
	leg_11.baseFrame = Base_11;
	leg_11.toolFrame = Tool_11;
	VERNE.Leg.push_back(new CSerialChain(&leg_11, current_joint_angles_11));

	leg_12.extDenHart = PreCalcDHparams(DH_12);
	leg_12.baseFrame = Base_12;
	leg_12.toolFrame = Tool_12;
	VERNE.Leg.push_back(new CSerialChain(&leg_12, current_joint_angles_12));

	leg_II.extDenHart = PreCalcDHparams(DH_II);
	leg_II.baseFrame = Base_II;
	leg_II.toolFrame = Tool_II;
	VERNE.Leg.push_back(new CSerialChain(&leg_II, current_joint_angles_II));

	leg_III.extDenHart = PreCalcDHparams(DH_III);
	leg_III.baseFrame = Base_III;
	leg_III.toolFrame = Tool_III;
	VERNE.Leg.push_back(new CSerialChain(&leg_III, current_joint_angles_III));


	VERNE.pathMatrix.resize(path_matrix.rows(), path_matrix.cols());
// // Set path matrix as VERNE sparse matrix
	for(int i = 0; i < path_matrix.rows(); i++){
	for(int j = 0; j < path_matrix.cols(); j++){
		if(path_matrix(i,j)!=0) {
		VERNE.pathMatrix.insert(i,j) = path_matrix(i,j);
			}
		}
	}

	int max_it = MAX_ITERATION_NUMBER;
	Initialize();

	Matrix4d Ton_new(Matrix4d::Identity());
	Matrix4d Ton(Matrix4d::Identity());
	// extract rotation matrices
	Eigen::Ref<Matrix3d> Ron_new = Ton_new.block(0,0, 3, 3);   
	Eigen::Ref<Matrix3d> Ron = Ton.block(0,0, 3, 3);

	// extract translational vectors
	Eigen::Ref<Vector3d> xyz = Ton.col(3).segment(0,3);
	Eigen::Ref<Vector3d> xyz_new = Ton_new.col(3).segment(0,3);
	

	Vector3d rpy_new(Vector3d::Zero());
	Vector3d rpy(Vector3d::Zero());
	VectorXd tcp_new(VectorXd::Zero(6));
	VectorXd tcp(VectorXd::Zero(6));
	
	// local temp variables
	VectorXd current_joint_angles(VectorXd::Zero(6));
	MatrixXd dh;
	Matrix4d baseFrame;
	Matrix4d toolFrame;
	VectorXd dq;
	VectorXd D;


	// data storage vectors and matrices
	// [calc error, number of iterations, elapsed time]
	MatrixXd test_data_leg11(MatrixXd::Zero(POINT_CNT,3));
	MatrixXd test_data_leg12(MatrixXd::Zero(POINT_CNT,3));
	MatrixXd test_data_legII(MatrixXd::Zero(POINT_CNT,3));
	MatrixXd test_data_legIII(MatrixXd::Zero(POINT_CNT,3));
	MatrixXd test_data_global(MatrixXd::Zero(POINT_CNT,3));
	MatrixXd joint_values(MatrixXd::Zero(POINT_CNT, NUMB_OF_JOINTS));


	int leg_index, n_it;
	double error, tcp_error;

	int point_cnt = 0;
	// main loop
	for(int i = 0; i < POINT_CNT; i++)  {
		point_cnt++;
		//cout << "move count: " << point_cnt << "\n" << endl;
		// set new tcp point
		rpy_new(2) = test_trajectory(i,3);
		//std::cout<<test_trajectory(i,0)<<std::endl;
		rpy2tr(Ton_new, rpy_new);  // update Ton_new with new rpy angles
		// update Ton_new with new xyz values
		//Ron_new = Ton_new.block(0,0,3,3);
		//xyz_new = Ton_new.block(0,3,3,1);
		//std::cout<<Ton_new<<std::endl;
		tcp_new(0) = xyz_new(0) = test_trajectory(i,0);
		tcp_new(1) = xyz_new(1) = test_trajectory(i,1);
		tcp_new(2) = xyz_new(2) = test_trajectory(i,2);
		tcp_new(3) = rpy_new(0); tcp_new(4) = rpy_new(1); tcp_new(5) = rpy_new(2);
		//std::cout<<"Ton_new"<<Ton_new<<std::endl;
		leg_index = 0;
		
		m_el.Begin(); // start calc of elapsed time here
		for((VERNE.Leg_it) = VERNE.Leg.begin(); (VERNE.Leg_it) <  VERNE.Leg.end(); ++(VERNE.Leg_it))  {
			current_joint_angles = (*(VERNE.Leg_it))->jointAngles;
			//std::cout<<current_joint_angles<<std::endl;
			dh = (*(VERNE.Leg_it))->legParam->extDenHart;
			baseFrame = (*(VERNE.Leg_it))->legParam->baseFrame;
			toolFrame = (*(VERNE.Leg_it))->legParam->toolFrame;
			error = 1.0;
			tcp_error = 1.0;
			n_it = 0;
			// Initial guess
			Ton = fkine(dh, baseFrame, toolFrame, current_joint_angles);
			//Ron = Ton.block(0,0,3,3);
			//xyz = Ton.block(0,3,3,1);
#ifdef _DEBUG
			//if(point_cnt > 300) {
				//std::cout << "LEG " << leg_index << "\n" << std::endl;
				//std::cout << "Initial joint angles: " << current_joint_angles << std::endl;
				//std::cout << "Initial guess\n" << Ton << std::endl;
				//std::cout << "New tcp matrix:\n" << Ton_new << std::endl;
			//}
#endif
			while(error > ITERATION_TOL)  {
				D = ikinVERNE_POE(dh, baseFrame, Ron_new, Ron,
							  xyz_new, xyz, current_joint_angles);
				Ton = fkine(dh, baseFrame, toolFrame, current_joint_angles);
				//std::cout<<"Ton"<<Ton<<std::endl;
				//std::cout<<Ron<<std::endl;
				//std::cout<<current_joint_angles<<std::endl;

				//Ron = Ton.block(0,0,3,3);
				//xyz = Ton.block(0,3,3,1);
				tr2rpy(rpy, Ton);
				tcp(0) = xyz(0); tcp(1) = xyz(1); tcp(2) = xyz(2);
				tcp(5) = rpy(2);
				tcp_error = (xyz_new-xyz).norm();

				error = (D).norm();
				//std::cout<<"error"<<error<<std::endl;

				//error = tcp_error;
				std::cout<<n_it<<std::endl;
				if (n_it>MAX_ITERATION_NUMBER)
				{
					cout<<"Maybe stuck"<<endl;
					break;
				}
				n_it++;
#ifdef _DEBUG
				/*cout << "Ton:" << endl;
				cout << Ton << endl;*/
#endif
			}
			
			// update leg joint angles
			(*(VERNE.Leg_it))->jointAngles = current_joint_angles;
			cout<<m_el.End()<<endl;
			// store test data for each leg
			switch(leg_index)  {
				case 0:
					test_data_leg11(i,0) = tcp_error;
					test_data_leg11(i,1) = n_it;
					test_data_leg11(i,2) = m_el.m_dElapsed;
					//std::cout<<test_data_leg11(i,0)<<std::endl;
					break;
				case 1:
					test_data_leg12(i,0) = tcp_error;
					test_data_leg12(i,1) = n_it;
					test_data_leg12(i,2) = m_el.m_dElapsed;
					break;
				case 2:
					test_data_legII(i,0) = tcp_error;
					test_data_legII(i,1) = n_it;
					test_data_legII(i,2) = m_el.m_dElapsed;
					break;
				case 3:
					test_data_legIII(i,0) = tcp_error;
					test_data_legIII(i,1) = n_it;
					test_data_legIII(i,2) = m_el.m_dElapsed;
					break;
				default:
					std::cout << "leg_index error\n" << std::endl;
					exit(1);
			}

			int j = 0;
			//for(int j = 0; j < leg_dof ; j++){
			Eigen::SparseVector<double> pm_row= VERNE.pathMatrix.row(leg_index);
			for (Eigen::SparseVector<double>::InnerIterator it(pm_row);it;++it){
				//	cout<<it.index()<<endl;
					//cout<<dq_all(it.index())<<endl;
				if((*(VERNE.Leg_it))->legParam->extDenHart(j,5))
					//joint_values(i,j_it.Index()) = 1000.0 * (*(VERNE.Leg_it))->jointAngles(j);
					joint_values(i,it.index()) = 1000.0 * (*(VERNE.Leg_it))->jointAngles(j);

				else
					//joint_values(i,j_it.Index()) = RAD_TO_DEG * (*(VERNE.Leg_it))->jointAngles(j);
					joint_values(i,it.index()) = RAD_TO_DEG * (*(VERNE.Leg_it))->jointAngles(j);
				//(*(VERNE.Leg_it))->jointAngles(j) += dq_all(it.index());
				//cout<<(*(VERNE.Leg_it))->jointAngles(j)<<endl;
				j++;
			}
			leg_index++;

#ifdef _DEBUG
			if(point_cnt > 300) {
				//std::cout << "AFTER" << std::endl;
				//std::cout << "iteration error: " << error << "\n" << std::endl;
				//std::cout << "Ton:" << std::endl;
				//std::cout << Ton << "\n" << std::endl;
				// std::cout << "joint angles\n" << current_joint_angles << "\n" << std::endl;
				//std::cout << "tcp: " << tcp << "\n" << std::endl;
				// std::cout << "number of iterations: " << n_it << "\n" << std::endl;
				// std::cout << "\n" <<std::endl;
			}
#endif
		}
		 // end of elapsed time 
		// store global test data
		double max_error = test_data_leg11(i,0);
		max_error = std::max(max_error, test_data_leg12(i,0));
		max_error = std::max(max_error, test_data_legII(i,0));
		test_data_global(i,0) = std::max(max_error, test_data_legIII(i,0));
		test_data_global(i,1) = test_data_leg11(i,1) +
								 test_data_leg12(i,1) +
								 test_data_legII(i,1) +
								 test_data_legIII(i,1);
		test_data_global(i,2) = m_el.m_dElapsed;
	// store joint angles
		


	}
	//std::cout<<test_data_leg11<<std::endl;
	saveTestData("testDataLeg11.txt", test_data_leg11);
	//saveTestData("testDataLeg12.txt", test_data_leg12);
	//saveTestData("testDataLegII.txt", test_data_legII);
	//saveTestData("testDataLegIII.txt", test_data_legIII);
	saveTestData("testData.txt", test_data_global);
	saveTestData("jointValues.txt", joint_values);
	
	//out = dest + ' ';
	//cout << out << dest << endl;
	//num_char = ToString(dest, "%e", test_data_global[0][0], 25);
	//cout << num_char << endl;
	//cout << dest << endl;
	//num_char = ToString(dest, "%e", test_data_global[0][1], 25);
	//cout << num_char << endl;
	//cout << dest << endl;
	//num_char = ToString(d, "%e", test_data_global[0][2], 15);
	//cout << num_char << endl;
	//cout << dest << endl;
	
}

#endif // POE_BLOCK
