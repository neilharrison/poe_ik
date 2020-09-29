#if !defined __PHV_LINK__
#define __PHV_LINK__

#pragma once

class Clink
{
public:
	Clink(void);
	~Clink(void);

// DH parameters
	double dh_alpha;
	double dh_a;
	double dh_theta;
	double dh_d;

	double q; // current joint angle
	unsigned int dh_sigma; // 0 for revolute joints, 1 for prismatic
	unsigned int mdh; // DH convention: 0 for standard, 1 for modified (Craig convention)
};

#endif // __PHV_LINK__
