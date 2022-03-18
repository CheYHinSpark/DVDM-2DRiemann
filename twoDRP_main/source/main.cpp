// twoD_Riemann.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
#include <iostream>
#include <string>
#include <cstring>
#include <stdlib.h>

#include "../include/twoDRP.h"

const string VERSION_NAME = "22.03.19";

using namespace std;
int main(int argc, char* argv[])
{
    /* TEST
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

     测试内容结束*/

    MF.~CArray();       //前一时刻数据(x，y，方向，dvm)
    MF_c.~CArray();     //复制数据(x，y，方向，dvm)
    MF_m.~CArray();     //DUGKS专用的中间态数据(x，y，方向，dvm)
    Meq.~CArray();      //平衡态
    macro.~CArray();    //宏观量(x，y，（rho,u,v,e）)

    cout << "2D黎曼问题正式开始！！！当前版本号" << VERSION_NAME << endl;
    cout << endl;

    string mode;
    string compute_setting;
    string init_setting;

    if (argc == 4)
    {
        for (int i = 1;i < argc;i++)
        {
            char* pchar = argv[i];
            switch (pchar[0])
            {
            case '-':
            {
                switch (pchar[1])
                {
                case 'h':
                { // 帮助
                    break;
                }
                }
                break;
            }
            case 'i':
            {
                init_setting = &pchar[2];
                init_setting = ReadSettings(init_setting, "C0");
                break;
            }
            case 'c':
            {
                compute_setting = &pchar[2];
                compute_setting = ReadSettings(compute_setting, "default");
                break;
            }
            case 'm':
            {
                mode = &pchar[2];
                break;
            }
            default: { break; }
            }
        }
    }
    else
    {
        // 需要输入
        cout << "输入初始值设置文件名，或乱输入使用默认设置：" << endl;
        cin >> init_setting;
        init_setting = ReadSettings(init_setting, "C0");

        cout << "输入计算配置文件名，或乱输入使用默认配置：" << endl;
        cin >> compute_setting;
        compute_setting = ReadSettings(compute_setting, "default");

        cout << "输入计算模式，1为dvm，2为dvddvm，3为dvm_gh，4为dvddvm_gh，其他为矩方法" << endl;
        cin >> mode;
    }

    if (mode == "1") { twoDRP_DVM(init_setting, compute_setting); }
    else if (mode == "2") { twoDRP_DVDDVM(init_setting, compute_setting); }
    else if (mode == "3") { twoDRP_DVDDVM_gh(init_setting, compute_setting, 3); }
    else if (mode == "4") { twoDRP_DVDDVM_gh(init_setting, compute_setting, 4); }
    else { twoDRP_EQMOM(init_setting, compute_setting); }

    cout << "程序运行完毕，按任意键继续……";

    return 0;
}
