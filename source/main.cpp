// twoD_Riemann.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
#include <iostream>
#include <string>
#include <stdlib.h>

#include "../include/twoDRP.h"

const string VERSION_NAME = "22.03.14";

using namespace std;
int main()
{
    // TEST
    //dvm_bg = -144;
    //dvm_st = 1;
    //dvm_strategy = "sqrt-central";
    //dvm_strategy = "common-difference";
    //SetDirections();
    //SetDVM();
    //SetUVW();
    //double* rhoeq = new double[N] {};
    //double* ueq = new double[N]{};
    //double sig2eq;
    //double* dvddvm1 = new double[N * M]{};
    //double* dvddvm2 = new double[N * M]{};
    //double r = 1;
    //double u = 0.4;
    //double v = 0.1;
    //double e = 1;
    //double rt = (2 * e - u * u - v * v) / 3;
    //double a[4] = { -1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, -1 / rt };  // 分别为a, b1, b2, c
    //RedEquil_DVD(r, u, v, e, a, rhoeq, ueq, &sig2eq);
    //double aa[4] = { -1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, -1 / rt };  // 分别为a, b1, b2, c
    //RedEquil_DVM(r, u, v, e, aa, dvddvm1);
    //for (int n = 0;n < N;++n)
    //{
    //    RedEquil_DVDDVM(rhoeq[n], ueq[n], 0.5 * (ueq[n] * ueq[n] + sig2eq), &dvddvm2[n * M]);
    //}
    //double diff = 0;
    //for (int k = M * N - 1;k >= 0;--k)
    //{
    //    diff += abs(dvddvm1[k] - dvddvm2[k]);
    //}
    //cout << diff << endl;

    // 测试内容结束

    MF.~CArray();       //前一时刻数据(x，y，方向，dvm)
    MF_c.~CArray();     //复制数据(x，y，方向，dvm)
    Meq.~CArray();      //平衡态
    macro.~CArray();    //宏观量(x，y，（rho,u,v,e）)

    cout << "2D黎曼问题正式开始！！！当前版本号" << VERSION_NAME << endl;
    cout << endl;

#pragma region 输入配置
    string init_setting;
    cout << "输入初始值设置文件名，或乱输入使用默认设置：" << endl;
    cin >> init_setting;
    string* nlist = new string[20]{};
    string* vlist = new string[20]{};
    if (CSVRead(init_setting, 20, nlist, vlist))
    {
        for (int i = 0;i < 20;++i)
        {
            string n = nlist[i];
            string v_str = vlist[i];
            double v_dbl = atof(v_str.c_str());
            if (n == "r00") { r_init[0][0] = v_dbl; }
            if (n == "r10") { r_init[1][0] = v_dbl; }
            if (n == "r01") { r_init[0][1] = v_dbl; }
            if (n == "r11") { r_init[1][1] = v_dbl; }
            if (n == "u00") { u_init[0][0] = v_dbl; }
            if (n == "u10") { u_init[1][0] = v_dbl; }
            if (n == "u01") { u_init[0][1] = v_dbl; }
            if (n == "u11") { u_init[1][1] = v_dbl; }
            if (n == "v00") { v_init[0][0] = v_dbl; }
            if (n == "v10") { v_init[1][0] = v_dbl; }
            if (n == "v01") { v_init[0][1] = v_dbl; }
            if (n == "v11") { v_init[1][1] = v_dbl; }
            if (n == "p00") { p_init[0][0] = v_dbl; }
            if (n == "p10") { p_init[1][0] = v_dbl; }
            if (n == "p01") { p_init[0][1] = v_dbl; }
            if (n == "p11") { p_init[1][1] = v_dbl; }
        }
    }
    else
    {
        cout << "将使用默认值" << endl;
        init_setting = "C0";
    }
    delete[] nlist;
    delete[] vlist;
    cin.clear();
    cin.sync();

    string compute_setting;
    cout << "输入计算配置文件名，或乱输入使用默认配置：" << endl;
    cin >> compute_setting;
    nlist = new string[20]{};
    vlist = new string[20]{};
    if (CSVRead(compute_setting, 20, nlist, vlist))
    {
        for (int i = 0;i < 20;++i)
        {
            string n = nlist[i];
            string v_str = vlist[i];
            double v_dbl = atof(v_str.c_str());
            if (n == "Nlen") { N = (int)v_dbl; }
            if (n == "dx") { dx = v_dbl; }
            if (n == "dy") { dy = v_dbl; }
            if (n == "flag_collide") { flag_collide = v_dbl; }
            if (n == "kappa") { kappa = v_dbl; }
            if (n == "save_interval") { save_interval = v_dbl; }
            if (n == "save_max") { save_max = (int)v_dbl; }
            if (n == "CFLMax") { CFLMax = v_dbl; }
            if (n == "dvm_bg") { dvm_bg = v_dbl; }
            if (n == "dvm_st") { dvm_st = v_dbl; }
            if (n == "REDEQUIL_ERROR") { REDEQUIL_ERROR = v_dbl; }
            if (n == "SCHEME") { SCHEME = v_str; }
            if (n == "dvm_strategy") { dvm_strategy = v_str; }
        }
    }
    else
    {
        cout << "将使用默认值" << endl;
        compute_setting = "default";
    }
    delete[] nlist;
    delete[] vlist;
#pragma endregion

    string mode;
    cout << "输入计算模式，1为dvm，2为dvddvm，3为dvm_gh，4为dvddvm_gh，其他为矩方法" << endl;
    cin >> mode;
    
    if (mode == "1") { twoDRP_DVM(init_setting, compute_setting); }
    else if (mode == "2") { twoDRP_DVDDVM(init_setting, compute_setting); }
    else if (mode == "3") { twoDRP_DVM_gh(init_setting, compute_setting); }
    else if (mode == "4") { twoDRP_DVDDVM_gh(init_setting, compute_setting); }
    else { twoDRP_EQMOM(init_setting, compute_setting); }

    cout << "程序运行完毕，按任意键继续……";

    return 0;
}

