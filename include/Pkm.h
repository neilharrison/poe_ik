#if !defined __PHV_PKM__
#define __PHV_PKM__

#pragma once

#include "SerialChain.h"


class CPkm
{
public:
	CPkm(void);
	~CPkm(void);

	void fkinAll() 
	{
		Leg_Map["Base"]->Ton = Leg_Map["Base"]->fkine();
		fkinRecur("Base");
	}
	void fkinRecur(std::string p)
	{
		for (std::string &c:Leg_Map[p]->children){
			Leg_Map[c]->Ton = Leg_Map[p]->Ton*Leg_Map[c]->fkine();
			std::cout<<c<<std::endl;
			std::cout<<Leg_Map[c]->Ton<<std::endl;
			fkinRecur(c);
		}
	}


	MatrixXd Jacobian(Matrix4d &Ton_new) 
	{
		//Also calculates pose difference (D) from Ton_new as convienient to do here
		int numLegs = pathMatrix.rows();
		A = MatrixXd::Zero(6*numLegs, 6*numLegs);
		B = MatrixXd::Zero(6*numLegs,pathMatrix.cols());
		D = VectorXd::Zero(6*numLegs);

		int leg_index = 0;
		int pm_index = 0;
		
		ABDrecur(End_effectors[0],leg_index,pm_index,Ton_new);

		return A*B;
	}

	void ABDrecur(std::string p, int &legindex,int &pm_index, Matrix4d &Ton_new) 
	{
		for (std::string & p:Leg_Map[p]->parent)
			{
				if (p != "Base")
				{	
					if (pm_index==0){ //At end-effector leg - now calc A matrix and D Vector
						A.block(legindex*6,legindex*6,6,6) = AdjInvH(Leg_Map[p]->Ton);
						D.segment(legindex*6,6) = poseDiff(Ton_new,Leg_Map[p]->Ton);
					}
					//This bit calcs B*S Matrix
					// 100% sure there are efficiencies to be made here - lots of exta variable/iterations
					MatrixXd temp_B = Leg_Map[p]->BstarS(Leg_Map[Leg_Map[p]->parent[0]]->Ton);
					Eigen::SparseVector<double> pm_row = pathMatrix.row(legindex);
					
					int i = 0;
					int j = 0;

					for (Eigen::SparseVector<double>::ReverseInnerIterator it(pm_row);it;--it){
						if (j>=pm_index){ //skip through to where we left off
							if (i <temp_B.cols()){
								B.block(legindex*6,it.index(),6,1) = temp_B.col(temp_B.cols()-i-1);
								i++;}
							else{break;}
							pm_index++;
						}
						j++;						
					}
					ABDrecur(p,legindex,pm_index,Ton_new);
				}
				else 
				{	
					pm_index = 0;
					legindex++;
				}
			}
	}

	VectorXd recalcD(Matrix4d &Ton_new)
	{
		D = VectorXd::Zero(6*pathMatrix.rows());
		int legindex = 0;
		for (std::string & p:Leg_Map[End_effectors[0]]->parent)
		{
			D.segment(legindex*6,6) = poseDiff(Ton_new,Leg_Map[p]->Ton);
			legindex++;
		}
		return D;
	}

	void updateJointAngles(VectorXd &dq)
	{
		int leg_index = 0;
		int pm_index = 0;
		jointAngleRecur(End_effectors[0],leg_index,pm_index,dq);
	}

	void jointAngleRecur(std::string p, int &legindex,int &pm_index,VectorXd &dq)
	{
		for (std::string & p:Leg_Map[p]->parent)
			{
				if (p != "Base")
				{	
					Eigen::SparseVector<double> pm_row = pathMatrix.row(legindex);
					int i = 0;
					int j = 0;

					for (Eigen::SparseVector<double>::ReverseInnerIterator it(pm_row);it;--it){
						if (j>=pm_index){ //skip through to where we left off
							if (i <Leg_Map[p]->jointAngles.rows()){
								//B.block(legindex*6,it.index(),6,1) = temp_B.col(temp_B.cols()-i-1);
								//std::cout<<p<<Leg_Map[p]->jointAngles.rows()-i-1<<"-"<< dq(it.index())<<std::endl;
								Leg_Map[p]->jointAngles(Leg_Map[p]->jointAngles.rows()-i-1) += dq(it.index());
								i++;}
							else{break;}
							pm_index++;
						}
						j++;						
					}
					jointAngleRecur(p,legindex,pm_index,dq);
				}
				else 
				{	
					pm_index = 0;
					legindex++;
				}
			}
	}

	void iksolve(Matrix4d &Ton_new) {

		MatrixXd J = MatrixXd::Zero(pathMatrix.rows()*6,pathMatrix.cols());
		VectorXd dq_all(VectorXd::Zero(pathMatrix.cols()));
		int n_it = 0;
		double error = 1.0;
		fkinAll();
		//recalcD(Ton_new);
		//error = D.norm();

		while(error > 1E-3)  {
			J = Jacobian(Ton_new);
			//std::cout<<J<<std::endl;
			//std::cout<<"J\n"<<J<<std::endl;
			//std::cout<<"D\n"<<D<<std::endl;
			dq_all = J.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(D);
			//Put dq_all into legs
			//std::cout<<"dq\n"<<dq_all<<std::endl;
			//std::cout<<"pre\n"<<Leg_Map["UR5_1"]->jointAngles<<std::endl;
			updateJointAngles(dq_all);
			//std::cout<<"post\n"<<Leg_Map["UR5_1"]->jointAngles<<std::endl;
			fkinAll();
			recalcD(Ton_new);
			if (n_it>5000)
			{
				std::cout<<"Maybe stuck"<<std::endl;
				break;
			}
			//std::cout<<"D\n"<<D<<std::endl;
			error = D.norm();
			std::cout<<n_it<<std::endl;
			n_it++;
		}
	}


//variables
	char Name[16];				// PKM name
	int numLegs;				// number of independent serial chains
	Matrix4d base;					// PKM base matrix;
	Matrix4d tool;					// PKM tool matrix;
	Eigen::SparseMatrix<double> pathMatrix;		// PKM path matrix;

	MatrixXd B;
	MatrixXd A;
	VectorXd D;

// Serial chain definitions
	std::vector<CSerialChain*> Leg;
	std::vector<CSerialChain*>::iterator Leg_it;
	std::map<std::string,CSerialChain*> Leg_Map;
	std::vector<std::string> End_effectors;
};
#endif // __PHV_PKM__
