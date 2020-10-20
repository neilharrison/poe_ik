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

	Matrix4d fkine()
	{
		Matrix4d base_tool_transform(Matrix4d::Zero());
		Matrix4d temp_matrix(Matrix4d::Zero());
		temp_matrix = legParam->baseFrame;
		
		for(int i = 0; i < legParam->extDenHart.rows(); i++)  { 
			//base_tool_transform = temp_matrix * linkTransformMod(dh_param.row(i), current_joint_angles(i));
			base_tool_transform = temp_matrix * linkTransform(legParam->extDenHart.row(i), jointAngles(i));
			temp_matrix = base_tool_transform;
		}
		//std::cout<<legParam->toolFrame<<std::endl;
		//std::cout<<base_tool_transform<<std::endl;
		base_tool_transform = temp_matrix * legParam->toolFrame;
		//std::cout<<base_tool_transform<<std::endl;

		//cout<<base_tool_transform<<endl;
		//Ton = base_tool_transform;
		return base_tool_transform;
	}

	MatrixXd BstarS(Matrix4d &parentFrame)
	{
		//int dof = dh.rows();
		
		// Calculates B*S and replace as "B" matrix 
		// Vectors are built up as row vectors
		// B need to be transposed before being used
		MatrixXd B(MatrixXd::Zero(legParam->extDenHart.rows(), 6));
		VectorXd s(VectorXd::Zero(6)); // Twist local coordinates
		Matrix4d Toi = parentFrame * legParam->baseFrame;
		
		int i = 0;
		//SVIterd j_it;
		//for (j_it.Begin(VERNE.pathMatrix.row(leg_index); !j_it.AtEnd(); j_it.Inc()) {
		//Eigen::SparseVector<double> pm_row = VERNE.pathMatrix.row(leg_index);
		
		//for (Eigen::SparseVector<double>::InnerIterator it(pm_row);it;++it){
		for (int it = 0; it<legParam->extDenHart.rows();it++) {
			//cout<<it.index()<<endl;
			//cout<<dq_all(it.index())<<endl;		
			//Toi = Toi * linkTransformMod(dh.row(i), q(i));
			
			s(2) = 1.0 - legParam->extDenHart(i,5);
			s(5) = legParam->extDenHart(i,5);
			if (i ==0) {
				B.row(it) = Adjnt(parentFrame)*s;
			}
			else {
				Toi = Toi * linkTransform(legParam->extDenHart.row(i-1), jointAngles(i-1));
				B.row(it) = Adjnt(Toi)*s;
			}
			//cout<<it.index()<<endl;
			//cout<<B.row(it.index())<<endl;
			i++;
		}

		//cout << B.transpose() << endl;
		return B.transpose(); //*S;
	}

//structure variables
	int dof;				// chain Degrees Of Freedom
	t_kinParam* legParam;   // leg kinematic parameters
	//Matd extDenHart;		// extended Denavit-Hartenberg parameters (8 elements)
	//Mat4d baseFrame;		// chain base reference frame
	//Mat4d toolFrame;		// chain tool reference frame
	std::vector<std::string> parent;
	std::vector<std::string> children;
	
// state variables
	VectorXd jointAngles;		// current joint angles
	Matrix4d Ton;              // base-tool transformation;
	//Clink* link;	// chain link

};
#endif // __PHV_SERIALCHAIN__