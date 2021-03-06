#if !defined __PHV_VERNEKINBYPOE__
#define __PHV_VERNEKINBYPOE__

#pragma inline_depth (254)

#define EPS   2.2204e-016//1e-17

#define RAD_TO_DEG 57.2957795130823209
//#include "Elapsed.h"
//#include "SerialChain.h"
#include "Pkm.h"
#include "Elapsed.h"

#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "std_msgs/Float64MultiArray.h"
#include "std_msgs/Float64.h"
#include "geometry_msgs/PoseStamped.h"

using std::cout;
using std::endl;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;
using Eigen::Vector4d;

std::vector<MatrixXd> DH_Parameters;
std::vector<MatrixXd>::iterator DH_Parameters_it;

CElapsed m_el;

std::vector<CSerialChain*> PKM ;
std::vector<CSerialChain*>::iterator PKM_it;

CPkm VERNE;

void saveTestData (char *, MatrixXd &d);

//Matd addDiagBlock(Matd &, Matd &);

MatrixXd BstarS(MatrixXd &dh, Matrix4d &baseFrame, VectorXd &q, int &leg_index);
MatrixXd bm_jacobian(MatrixXd &dh, Matrix4d &baseFrame, VectorXd &q, Matrix3d &Ron, Vector3d &xyz);
VectorXd ikinVERNE_POE(MatrixXd &dh, Matrix4d &baseFrame, Matrix3d &Ron_new, Matrix3d &Ron, Vector3d &xyz_new, Vector3d &xyz, VectorXd &q);
MatrixXd PreCalcDHparams(MatrixXd &dh_params);
//void rpy2tr(Mat4d &, Vec3d &);
//void tr2rpy(Vec3d &,Mat4d &);
void SVDecomposition(MatrixXd &A, MatrixXd &U, MatrixXd&V, VectorXd &s);
Matrix4d fkine(MatrixXd &dh_param, Matrix4d &base_frame, Matrix4d &tool_frame, VectorXd &current_joint_angles);
void Initialize(void);

MatrixXd pinv(MatrixXd M);
// Matd pinv2(Matd &);


Matrix4d transformFromPose(Vector4d &q, Vector3d &p);

void iksolve(Matrix4d &Ton_new, int &n_it);
Matrix4d homeTon;

//Ros subscriber/publisher 
//This does everything! subscribes to ik goals, calls ik solve, and publishes joint angles
class SubscribeAndPublish
{
public:
  	SubscribeAndPublish()
  	{
    	//Topic you want to publish
    	pub_UR1 = n_.advertise<std_msgs::Float64MultiArray>("/UR1/joint_group_position_controller/command", 1000);
		pub_UR2 = n_.advertise<std_msgs::Float64MultiArray>("/UR2/joint_group_position_controller/command", 1000);

		pub_joint_A = n_.advertise<std_msgs::Float64>("/wrist/joint1_position_controller/command", 1000);
		pub_joint_B = n_.advertise<std_msgs::Float64>("/wrist/joint2_position_controller/command", 1000);


    	//Topic you want to subscribe
    	sub_ = n_.subscribe("/ik_goal", 1000, &SubscribeAndPublish::ikCallback, this);
  	}

 
  	void ikCallback(const geometry_msgs::PoseStamped::ConstPtr& msg) {

		Vector3d pos = (Vector3d()<<msg->pose.position.x,msg->pose.position.y,msg->pose.position.z).finished();
		Vector4d quat = (Vector4d()<<msg->pose.orientation.x,msg->pose.orientation.y,msg->pose.orientation.z,msg->pose.orientation.w).finished();

		//Making sure pos is within reach
		cout<<pos<<endl;
		//Matrix4d temptest =(Matrix4d()<<1,0,0,0.1,0,1,0,0,0,0,1,0,0,0,0,1).finished();
		Matrix4d Ton_new = homeTon * transformFromPose(quat,pos) ;
		//Matrix4d Ton_new = homeTon * temptest;
		/*
		/*
		Eigen::Ref<Vector3d> pos_new = Ton_new.col(3).segment(0,3);
		//0.0767 from wrist geometry - sqrt(2*L^2 -2*L^2*cos(angle))
		pos_new = 0.076693586680288*(pos_new/pos_new.norm());
		cout<<pos_new<<endl;
		
		cout<<Ton_new<<endl;
		*/

		int n_it=0;
		m_el.Begin();
		//iksolve(Ton_new,n_it);
		VERNE.iksolve(Ton_new);
		cout<<m_el.End()<<endl;

		std_msgs::Float64MultiArray joints_UR1;
		std_msgs::Float64MultiArray joints_UR2;

		std_msgs::Float64 joint_a;
		std_msgs::Float64 joint_b;
	
	
		joints_UR1.data.clear();
		joints_UR2.data.clear();

		for (int i = 0;i<VERNE.Leg_Map["UR5_1"]->jointAngles.rows();i++){
			joints_UR1.data.push_back(VERNE.Leg_Map["UR5_1"]->jointAngles(i));
			std::cout<<"UR1"<<VERNE.Leg_Map["UR5_1"]->jointAngles(i)<<std::endl;
		}
		for (int i = 0;i<VERNE.Leg_Map["UR5_2"]->jointAngles.rows();i++){
			joints_UR2.data.push_back(VERNE.Leg_Map["UR5_2"]->jointAngles(i));
			std::cout<<"UR2"<<VERNE.Leg_Map["UR5_2"]->jointAngles(i)<<std::endl;

		}


		// joint_a.data = ja.data[0];
		// joint_b.data = ja.data[4];
		// cout<<joint_a.data<<endl;
		// cout<<joint_b.data<<endl;
		pub_UR1.publish(joints_UR1);
		pub_UR2.publish(joints_UR2);

		//pub_joint_A.publish(joint_a);
		//pub_joint_B.publish(joint_b);
		//ros::spinOnce();
		//cout<<n_it<<endl;
  	}


private:
  	ros::NodeHandle n_; 
  	ros::Publisher pub_UR1;
	ros::Publisher pub_UR2;

	ros::Publisher pub_joint_A;
	ros::Publisher pub_joint_B;
  	ros::Subscriber sub_;

};//End of class SubscribeAndPublish


inline VectorXd quickTr2Diff(Matrix4d &m)
{
	VectorXd d(VectorXd::Zero(6));
	d(0) = m(2,1);
	d(1) = m(0,2);
	d(2) = m(1,0);
	d(3) = m(0,3);
	d(4) = m(1,3);
	d(5) = m(2,3);
	return d;
}

inline MatrixXd Adjnt(Matrix4d &T)
{
	Matrix3d R((Matrix3d()<< T(0,0), T(0,1), T(0,2),
			T(1,0), T(1,1), T(1,2),
			T(2,0), T(2,1), T(2,2)).finished());

	Matrix3d Dlt((Matrix3d()<< 0.0, -T(2,3), T(1,3),
			  T(2,3), 0.0, -T(0,3),
			 -T(1,3), T(0,3), 0.0).finished());

	MatrixXd Adj_T(MatrixXd::Zero(6, 6));

	Dlt = Dlt * R;
	for(int j = 0; j < 3; j++)
		for(int k = 0; k < 3; k++)  {
			Adj_T(j,k) = R(j,k);
			Adj_T(j+3,k) = Dlt(j,k);
			Adj_T(j+3,k+3) = R(j,k);
		}
	
	return Adj_T;
}

inline MatrixXd AdjInvH(Matrix4d &T)
{
	// Matrix3d R_trans(T[0][0], T[1][0], T[2][0],
	// 				T[0][1], T[1][1], T[2][1],
	// 				T[0][2], T[1][2], T[2][2]);

	Matrix3d R_trans = T.block(0,0,3,3).transpose();

	Matrix3d Dlt((Matrix3d()<< 0.0, -T(2,3), T(1,3),
			 				  T(2,3), 0.0, -T(0,3),
							 -T(1,3), T(0,3), 0.0).finished());



	MatrixXd Adj_inv_T(MatrixXd::Zero(6, 6));

	Dlt = -R_trans * Dlt;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)  {
			Adj_inv_T(i,j) = R_trans(i,j);
			Adj_inv_T(i+3,j) = Dlt(i,j);
			Adj_inv_T(i+3,j+3) = R_trans(i,j);
		}
	
	return Adj_inv_T;
}
inline void changeSubMatd(MatrixXd &Y, MatrixXd X, int row, int col, int elts_row, int elts_col)
{
	if(X.rows() != elts_row) 
		printf("changeSubMatd:: rows mismatch");
	if(X.cols() != elts_col)
		printf("changeSubMatd:: columns mismatch");
	Y.block(row,col,elts_row,elts_col) = X;
	//SubMatd SubOfY = sub(Y, row, col, elts_row, elts_col);
	//SubOfY = X;
}

inline void changeSubVecd(VectorXd &Y, VectorXd X, int row, int elts_row)
{
	if(X.size() != elts_row) 
		printf("changeSubMatd:: rows mismatch");
	Y.segment(row,elts_row) = X;
	//SubVecd SubOfY = sub(Y, row, elts_row);
	//SubOfY = X;
}
	
// inline void addRowBlock(MatrixXd &Y, MatrixXd &X)
// {
// 	MatrixXd Ytemp;
// 	int my = Y.rows(); int mx = X.rows();
// 	int ny = Y.cols(); int nx = X.cols();

// 	if(my && ny)  {
// 		Ytemp = Y;
// 		if (mx != my)
// 			printf("addRowBlock::number of row mismatch\n");
// 	}

// 	Y.SetSize(mx, (ny+nx));

// 	for(int i = 0; i < mx; i++) {
// 		if(ny)
// 			for(int j = 0; j < ny; j++)
// 				Y[i][j] = Ytemp[i][j];
// 		for(int j = 0; j < nx; j++)
// 			Y[i][j+ny] = X[i][j];
// 	}
// }

// inline void addColBlock(Matd &Y, Matd &X)
// {
// 	Matd Ytemp;
// 	int my = Y.Rows(); int mx = X.Rows();
// 	int ny = Y.Cols(); int nx = X.Cols();

// 	if(my && ny)  {
// 		Ytemp = Y;
// 		if (nx != ny)
// 			printf("addColBlock::number of columns mismatch\n");
// 	}

// 	Y.SetSize((mx+my), nx);

// 	for(int j = 0; j < nx; j++) {
// 		if(my)
// 			for(int i = 0; i < my; i++)
// 				Y[i][j] = Ytemp[i][j];
// 		for(int i = 0; i < mx; i++)
// 			Y[i+my][j] = X[i][j];
// 	}
// }
// inline void addVecBlock(Vecd &Y, Vecd &X)
// {
// 	Vecd Ytemp;
// 	int nx = X.Elts();
// 	int ny = Y.Elts();
// 	if(ny)
// 		Ytemp = Y;

// 	Y.SetSize(ny + nx);

// 	if(ny)
// 		for(int i = 0; i < ny; i++)
// 				Y[i] = Ytemp[i];
// 	for(int i = 0; i < nx; i++)
// 		Y[i+ny] = X[i];

// }

// inline void addDiagBlock(Matd &Y, Matd &X)
// {
// 	Matd Ytemp;
	
// 	int my = Y.Rows(); int mx = X.Rows();
// 	int ny = Y.Cols(); int nx = X.Cols();
	
// 	if(my && ny)
// 		Ytemp = Y;
	
// 	Y.SetSize((my + mx), (ny + nx));
// 	Y = vl_0;
// 	if(my && ny)
// 		for(int i = 0; i < my; i++)
// 			for(int j = 0; j < ny; j++)
// 				Y[i][j] = Ytemp[i][j];
// 	for(int i = 0; i < mx; i++)
// 		for(int j = 0; j < nx; j++)
// 			Y[my+i][ny+j] = X[i][j];
// }

inline void rpy2tr(Matrix4d &T, Vector3d &rpy)
{
	Matrix3d m;
	m = Eigen::AngleAxisd(rpy(0), Vector3d::UnitZ())
    * (Eigen::AngleAxisd(rpy(1), Vector3d::UnitY())
    * Eigen::AngleAxisd(rpy(2), Vector3d::UnitX()));

	T.block(0,0,3,3) = m;
}

inline void tr2rpy(Vector3d &rpy, Matrix4d &T)
{
	if( abs(T(0,0)) < EPS && abs(T(1,0)) < EPS)  {
		// singularity
		rpy(0) = 0;
		rpy(1) = atan2(-T(2,0), T(0,0));
		rpy(2) = atan2(-T(1,2), T(1,1));

	}
	else  {
		rpy(0) = atan2(T(1,0), T(0,0));
		double sp = sin(rpy(0));
		double cp = cos(rpy(0));
		rpy(1) = atan2(-T(2,0), (cp * T(0,0) + sp * T(1,0)));
		rpy(2) = atan2((sp * T(0,2) - cp * T(1,2)), (cp * T(1,1) - sp * T(0,1)));
	}
}

inline Matrix4d linkTransformMod(VectorXd dh_link_param, double current_joint_angle)
{
	// Link to link transformation matrix 
	// based on Craig notation
	
	// remember: dh_param has one more element at the beggining (because of cos(alpha) and sin(alpha))
	double DSa, DCa;
	double theta, Ct, St;
	Matrix4d T(Matrix4d::Identity());

	double cja = current_joint_angle + dh_link_param(7);
	if(dh_link_param(5)) {
		// Prismatic joint
		theta = dh_link_param(3);
		DSa = (dh_link_param(4) + cja) * dh_link_param(1);
		DCa = (dh_link_param(4) + cja) * dh_link_param(0);
	}
	else  {
		// Revolute joint
		theta = cja;
		DSa = dh_link_param(4) * dh_link_param(1);
		DCa = dh_link_param(4) * dh_link_param(0);
	}
	Ct = cos(theta);
	St = sin(theta);
	T(0,0) = Ct; T(0,1) = -St; T(0,2) = 0.0; T(0,3) = dh_link_param(2);
	T(1,0) = St*dh_link_param(0); T(1,1) = Ct*dh_link_param(0); T(1,2) = -dh_link_param(1); T(1,3) = -DSa;
	T(2,0) = St*dh_link_param(1); T(2,1) = Ct*dh_link_param(1); T(2,2) = dh_link_param(0); T(2,3) = DCa;

	return T;
}

inline Matrix4d linkTransform(VectorXd dh_link_param, double current_joint_angle)
{
	// Link to link transformation matrix 
	// based on standard notation
	// remember: dh_param has one more element at the beggining (because of cos(alpha) and sin(alpha))
	//std::cout<<current_joint_angle<<std::endl;
	double theta, Ct, St, D;
	Matrix4d T(Matrix4d::Identity());

	double cja = current_joint_angle + dh_link_param(7);
	if(dh_link_param(5)) {
		// Prismatic joint
		theta = dh_link_param(3);
		D = dh_link_param(4) + cja;
	}
	else  {
		// Revolute joint
		theta = cja;
		D = dh_link_param(4);
	}
	Ct = cos(theta);
	St = sin(theta);
	T(0,0) = Ct; T(0,1) = -St * dh_link_param(0); T(0,2) = St * dh_link_param(1); T(0,3) = dh_link_param(2) * Ct;
	T(1,0) = St; T(1,1) = Ct * dh_link_param(0); T(1,2) = -Ct * dh_link_param(1); T(1,3) = dh_link_param(2) * St;
	T(2,0) = 0.0; T(2,1) = dh_link_param(1); T(2,2) = dh_link_param(0); T(2,3) = D;
	
	return T;
}


inline Matrix4d transformFromPose(Vector4d &q, Vector3d &p){

	Matrix4d T = (Matrix4d()<< 1.0f - 2.0f*q(1)*q(1) - 2.0f*q(2)*q(2), 2.0f*q(0)*q(1) - 2.0f*q(2)*q(3), 2.0f*q(0)*q(2) + 2.0f*q(1)*q(3), p(0),
    				2.0f*q(0)*q(1) + 2.0f*q(2)*q(3), 1.0f - 2.0f*q(0)*q(0) - 2.0f*q(2)*q(2), 2.0f*q(1)*q(2) - 2.0f*q(0)*q(3), p(1),
    				2.0f*q(0)*q(2) - 2.0f*q(1)*q(3), 2.0f*q(1)*q(2) + 2.0f*q(0)*q(3), 1.0f - 2.0f*q(0)*q(0) - 2.0f*q(1)*q(1), p(2),
    				0.0f, 0.0f, 0.0f, 1.0f).finished();
	return T;
}


//From Modern Robotics textbook  - (git)

	bool NearZero(const double val) {
		return (std::abs(val) < .000001);
	}
	std::vector<Eigen::MatrixXd> TransToRp(const Eigen::MatrixXd& T) {
		std::vector<Eigen::MatrixXd> Rp_ret;
		Eigen::Matrix3d R_ret;
		// Get top left 3x3 corner
		R_ret = T.block<3, 3>(0, 0);

		Eigen::Vector3d p_ret(T(0, 3), T(1, 3), T(2, 3));

		Rp_ret.push_back(R_ret);
		Rp_ret.push_back(p_ret);

		return Rp_ret;
	}

	Eigen::Matrix3d VecToso3(const Eigen::Vector3d& omg) {
		Eigen::Matrix3d m_ret;
		m_ret << 0, -omg(2), omg(1),
			omg(2), 0, -omg(0),
			-omg(1), omg(0), 0;
		return m_ret;
	}


	Eigen::Matrix3d MatrixLog3(const Eigen::Matrix3d& R) {
		double acosinput = (R.trace() - 1) / 2.0;
		Eigen::MatrixXd m_ret = Eigen::MatrixXd::Zero(3, 3);
		if (acosinput >= 1)
			return m_ret;
		else if (acosinput <= -1) {
			Eigen::Vector3d omg;
			if (!NearZero(1 + R(2, 2)))
				omg = (1.0 / std::sqrt(2 * (1 + R(2, 2))))*Eigen::Vector3d(R(0, 2), R(1, 2), 1 + R(2, 2));
			else if (!NearZero(1 + R(1, 1)))
				omg = (1.0 / std::sqrt(2 * (1 + R(1, 1))))*Eigen::Vector3d(R(0, 1), 1 + R(1, 1), R(2, 1));
			else
				omg = (1.0 / std::sqrt(2 * (1 + R(0, 0))))*Eigen::Vector3d(1 + R(0, 0), R(1, 0), R(2, 0));
			m_ret = VecToso3(M_PI * omg);
			return m_ret;
		}
		else {
			double theta = std::acos(acosinput);
			m_ret = theta / 2.0 / sin(theta)*(R - R.transpose());
			return m_ret;
		}
	}


	Eigen::MatrixXd MatrixLog6(const Eigen::MatrixXd& T) {
		Eigen::MatrixXd m_ret(4, 4);
		auto rp = TransToRp(T);
		Eigen::Matrix3d omgmat = MatrixLog3(rp.at(0));
		Eigen::Matrix3d zeros3d = Eigen::Matrix3d::Zero(3, 3);
		if (NearZero(omgmat.norm())) {
			m_ret << zeros3d, rp.at(1),
				0, 0, 0, 0;
		}
		else {
			double theta = std::acos((rp.at(0).trace() - 1) / 2.0);
			Eigen::Matrix3d logExpand1 = Eigen::MatrixXd::Identity(3, 3) - omgmat / 2.0;
			Eigen::Matrix3d logExpand2 = (1.0 / theta - 1.0 / std::tan(theta / 2.0) / 2)*omgmat*omgmat / theta;
			Eigen::Matrix3d logExpand = logExpand1 + logExpand2;
			m_ret << omgmat, logExpand*rp.at(1),
				0, 0, 0, 0;
		}
		return m_ret;
	}





#endif // __PHV_VERNEKINBYPOE__