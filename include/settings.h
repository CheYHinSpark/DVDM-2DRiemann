/* 声明设置变量 */

#include <string>
#include "CArray.h"


#pragma region 各种设置

// 维度设置
extern int D;
extern int L;

// DVD相关
extern int N;   // 方向数量
extern double* x_n;   // 各个方向x
extern double* y_n;   // 各个方向y
extern double* z_n;   // 各个方向z
// DVM相关
extern double dvm_bg;   // dvm
extern double dvm_st;
extern double* dvm;
extern int M;
extern double* uk;
extern double* vk;
extern double* wk;
extern double* U2k;
// BGK相关
extern double flag_collide;
extern double kappa;
// 网格差分相关
extern double dx;
extern double dy;
extern double CFLMax;

//2D黎曼问题初始设置
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

// 计算总时间等
extern double t_current;
extern int step_count;    // 从-1开始

// 计算精度
extern double REDEQUIL_ERROR;

// 差分格式
extern std::string SCHEME;
// 离散点选取方法
extern std::string dvm_strategy;

// 保存用的数组
extern int save_max;            // 最大储存次数
extern double save_interval;    // 储存时间间隔，最终时间为（次数-1）*储存间隔
extern int Save_Count;

// 计算用大数组
extern CArray MF;       //前一时刻数据(x，y，方向，dvm)
extern CArray MF_c;     //复制数据(x，y，方向，dvm)
extern CArray Meq;      //平衡态
extern CArray macro;    //宏观量(x，y，（rho,u,v,e）)

// 设置DVM
void SetDVM();

// 设置u,v,w,U^2
void SetUVW();

// 设置方向
void SetDirections();

void Set3DNormalVec(double* v1, double* v2, double* v3, double x, double y, double z);
#pragma endregion