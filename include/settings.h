/* �������ñ��� */

#include <string>
#include "CArray.h"


#pragma region ��������

// ά������
extern int D;
extern int L;

// DVD���
extern int N;   // ��������
extern double* x_n;   // ��������x
extern double* y_n;   // ��������y
extern double* z_n;   // ��������z
// DVM���
extern double dvm_bg;   // dvm
extern double dvm_st;
extern double* dvm;
extern int M;
extern double* uk;
extern double* vk;
extern double* wk;
extern double* U2k;
// BGK���
extern double flag_collide;
extern double kappa;
// ���������
extern double dx;
extern double dy;
extern double CFLMax;

//2D���������ʼ����
/*  |-------|(2*x_num-1,2*y_num-1)
    |0,1|1,1|
    |-------|
    |0,0|1,0|
    |-------|(2*x_num-1,0)  */
extern double r_init[2][2];
extern double u_init[2][2];
extern double v_init[2][2];
extern double p_init[2][2];
extern double E_init[2][2];

// ������ʱ���
extern double t_current;
extern int step_count;    // ��-1��ʼ

// ���㾫��
extern double REDEQUIL_ERROR;

// ��ָ�ʽ
extern std::string SCHEME;
// ��ɢ��ѡȡ����
extern std::string dvm_strategy;

// �����õ�����
extern int save_max;            // ��󴢴����
extern double save_interval;    // ����ʱ����������ʱ��Ϊ������-1��*������
extern int Save_Count;

// �����ô�����
extern CArray MF;       //ǰһʱ������(x��y������dvm)
extern CArray MF_c;     //��������(x��y������dvm)
extern CArray Meq;      //ƽ��̬
extern CArray macro;    //�����(x��y����rho,u,v,e��)

// ����DVM
void SetDVM();

// ����u,v,w,U^2
void SetUVW();

// ���÷���
void SetDirections();

void Set3DNormalVec(double* v1, double* v2, double* v3, double x, double y, double z);
#pragma endregion