#include "../include/twoDRP.h"

void twoDRP_EQMOM(string init_setting, string compute_setting)
{

#pragma region 初始化计算参数，创建需要的数组

    // 设置各个方向，要避免平行于壁面的方向
    SetDirections();

    // 以当前时间确定一个文件夹名
    string foldername = "EQMOM" + to_string(N) + "d"
        + init_setting + compute_setting + CurrentLocalTime();
    string command;
    command = "mkdir " + foldername;
    system(command.c_str());


    // 创建各种数组
    double* uk = new double[M * N];
    double* vk = new double[M * N];
    double* wk = new double[M * N];
    for (int n = 0;n < N;++n)
    {
        for (int m = 0;m < M;++m)
        {
            uk[n * M + m] = x_n[n] * dvm[m];
            vk[n * M + m] = y_n[n] * dvm[m];
            wk[n * M + m] = z_n[n] * dvm[m];
        }
    }

    // 半隐式格式需要的系数
    double tau_pre = kappa;

    // 差分网格尺寸
    int x_num = int(0.5 / dx);    //  x轴格点数的一半
    int y_num = int(0.5 / dy);    //  y轴格点数的一半

    // 创建各种数组

    MF.reConstruct(N, 2 * x_num, 2 * y_num, 5);        //前一时刻数据(x，y，方向，dvm)
    // 实际上程序只需要存一组MF数据就够了，并且因为使用EQMOM，5阶矩完全没有用到，故略去
    Meq.reConstruct(N, 2 * x_num, 2 * y_num, 5);       //平衡态
    macro.reConstruct(2 * x_num, 2 * y_num, 4);
    CArray para_f(N, 2 * x_num, 2 * y_num, 5);      //双节点EQMOM的系数
    CArray alpha(2 * x_num, 2 * y_num, 4);          //平衡态参数
    //CArray alpha(2 * x_num, 2 * y_num, 5);          //平衡态参数

    // F
    CArray flux_p(N, 2 * x_num, 2 * y_num, 5);      // 正半轴积分，1到5阶
    CArray flux_n(N, 2 * x_num, 2 * y_num, 5);      // 负半轴积分，1到5阶
    CArray Flux_x(N, 2 * x_num + 1, 2 * y_num, 5);  // x方向多一格，各个格子对左右的影响，1到5阶
    CArray Flux_y(N, 2 * x_num, 2 * y_num + 1, 5);  // y方向多一格，各个格子对上下的影响，1到5阶

    // 保存用的
    double* save_t = new double[save_max] {};
    // 注意这里就是要把x和y反过来
    CArray save_r(2 * y_num, 2 * x_num);      // rho
    CArray save_u(2 * y_num, 2 * x_num);      // u
    CArray save_v(2 * y_num, 2 * x_num);      // v
    CArray save_e(2 * y_num, 2 * x_num);      // e
#pragma endregion

    cout << "计算方法EQMOM，方向数N=" << N << "，碰撞参数flag_collide=" << flag_collide << "，kappa=" << kappa <<
        "，网格数" << 2 * x_num << "*" << 2 * y_num << "，CFLMax=" << CFLMax << endl;
    cout << "将计算至t=" << (save_max - 1) * save_interval << endl;
    cout << "即将开始计算" << endl;

#pragma region 处理初始值
    // 初始能量
    for (int i = 0;i < 2;++i)
    {
        for (int j = 0;j < 2;++j)
        {
            E_init[i][j] = 0.5 * (pow(u_init[i][j], 2) + pow(v_init[i][j], 2) + 3 * p_init[i][j] / r_init[i][j]);
        }
    }
    double* req_i = new double[N] {};
    double* ueq_i = new double[N] {};
    double sig2eq_i = 0;
    // I,J大区坐标
    for (int I = 0;I < 2;++I)
    {
        for (int J = 0;J < 2;++J)
        {
            double u = u_init[I][J];
            double v = v_init[I][J];
            double rt = p_init[I][J] / r_init[I][J];
            double a[4] = { -1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, -1 / rt };  // 分别为a, b1, b2, c
            //double a[5] = { - 1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, 0, -1 / rt };  // 分别为a, b1, b2, c

            RedEquil_DVD(r_init[I][J], u_init[I][J], v_init[I][J], E_init[I][J],
                a, req_i, ueq_i, &sig2eq_i);
            for (int i = 0;i < x_num;++i)
            {
                for (int j = 0;j < y_num;++j)
                {
                    for (int k = 0;k < 4;++k)
                        //for (int k = 0;k < 5;++k)
                    {
                        alpha(i + I * x_num, j + J * y_num, k) = a[k];
                    }

                    for (int k = 0;k < N;++k)
                    {
                        // 设置初始MF
                        double* m = &MF(k, i + I * x_num, j + J * y_num, 0);
                        m[0] = req_i[k];
                        m[1] = req_i[k] * ueq_i[k];
                        m[2] = ueq_i[k] * m[1] + sig2eq_i * m[0];
                        m[3] = ueq_i[k] * m[2] + 2 * sig2eq_i * m[1];
                        m[4] = ueq_i[k] * m[3] + 3 * sig2eq_i * m[2];
                        // 初始双节点EQMOM
                        invEQ2(m[0], m[1], m[2], m[3], m[4], &para_f(k, i + I * x_num, j + J * y_num, 0));
                    }
                    // 初始宏观量
                    macro(i + I * x_num, j + J * y_num, 0) = r_init[I][J];
                    macro(i + I * x_num, j + J * y_num, 1) = u_init[I][J];
                    macro(i + I * x_num, j + J * y_num, 2) = v_init[I][J];
                    macro(i + I * x_num, j + J * y_num, 3) = E_init[I][J];

                    save_r(j + J * y_num, i + I * x_num) = r_init[I][J];
                    save_u(j + J * y_num, i + I * x_num) = u_init[I][J];
                    save_v(j + J * y_num, i + I * x_num) = v_init[I][J];
                    save_e(j + J * y_num, i + I * x_num) = E_init[I][J];
                }
            }
        }
    }
    delete[] req_i;
    delete[] ueq_i;

    // 第一次储存
    CSVWriteMatrix(foldername + "\\" + "r0.csv", &save_r(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "u0.csv", &save_u(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "v0.csv", &save_v(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "e0.csv", &save_e(0, 0), 2 * y_num, 2 * x_num);

#pragma endregion

    cout << "初始值处理完毕，主循环开始" << endl;

#pragma region 主循环
    double clock_1, clock_2, clock_3, clock_4, clock_5, clock_begin, clock_end;
    clock_begin = clock();
    while (Save_Count < save_max)
    {
        step_count++;
        clock_1 = clock();

#pragma region 第1步，确定时间步长
        // 模仿之前的程序
        double Umax = 0.;
#pragma omp parallel for
        for (int k = 0;k < N;++k)
        {
            double Uijmax;
            for (int i = 0;i < 2 * x_num;++i)
            {
                for (int j = 0;j < 2 * y_num;++j)
                {
                    Uijmax = max(abs(para_f(k, i, j, 1)) + 1.8 * sqrt(2 * para_f(k, i, j, 4)),
                        abs(para_f(k, i, j, 3)) + 1.8 * sqrt(2 * para_f(k, i, j, 4)));
                    Umax = max(Umax, Uijmax);
                }
            }
        }
        // CFL确定时间步长，参考Fox08年文献
        double dt;
        if (flag_collide == 1)
        {
            if (kappa == 0)
            {
                dt = CFLMax * dx / max(Umax, 0.1);
            }
            else
            {
                dt = min(CFLMax * dx / max(Umax, 0.1), kappa / 10);
            }
        }
        else
        {
            dt = CFLMax * dx / max(Umax, 0.1);
        }
        t_current += dt;    // 当前时刻
        // cout << dt << endl;
#pragma endregion

        clock_2 = clock();

#pragma region 第2步，根据宏观量计算出各方向的平衡态
        // 如果是无碰撞情形可以跳过这一步
        if (flag_collide == 1)
        {
#pragma omp parallel for
            for (int i = 0;i < 2 * x_num;++i)
            {
                double* req = new double[N] {};
                double* ueq = new double[N] {};
                double sig2eq = 0;
                for (int j = 0;j < 2 * y_num;++j)
                {
                    RedEquil_DVD(macro(i, j, 0), macro(i, j, 1), macro(i, j, 2), macro(i, j, 3),
                        &alpha(i, j, 0), req, ueq, &sig2eq);
                    for (int k = 0;k < N;++k)
                    {
                        double* m = &Meq(k, i, j, 0);
                        m[0] = req[k];
                        m[1] = req[k] * ueq[k];
                        m[2] = ueq[k] * m[1] + sig2eq * m[0];
                        m[3] = ueq[k] * m[2] + 2 * sig2eq * m[1];
                        m[4] = ueq[k] * m[3] + 3 * sig2eq * m[2];
                    }
                }
                delete[] req;
                delete[] ueq;
            }
        }
#pragma endregion

        clock_3 = clock();

#pragma region 第3步，计算Flux
#pragma omp parallel for
        for (int k = 0;k < N;++k)
        {
            double u_f_n[6][2]{};   // 计算Flux时用到的临时变量
            double u_f_p[6][2]{};
            double w1, u1, w2, u2, sig2;
            bool is_cosk_p = x_n[k] > 0;
            bool is_sink_p = y_n[k] > 0;
            // 先计算每个网格对两边的影响
            for (int i = 0;i < 2 * x_num;++i)
            {
                for (int j = 0;j < 2 * y_num;++j)
                {
                    w1 = para_f(k, i, j, 0);
                    u1 = para_f(k, i, j, 1);
                    w2 = para_f(k, i, j, 2);
                    u2 = para_f(k, i, j, 3);
                    sig2 = para_f(k, i, j, 4);

                    // <v ^ 0>
                    u_f_p[0][0] = 0.5 * erfc(-u1 / sqrt(2 * sig2));
                    u_f_n[0][0] = 0.5 * erfc(u1 / sqrt(2 * sig2));
                    u_f_p[0][1] = 0.5 * erfc(-u2 / sqrt(2 * sig2));
                    u_f_n[0][1] = 0.5 * erfc(u2 / sqrt(2 * sig2));
                    // <v ^ 1>
                    u_f_p[1][0] = u1 * u_f_p[0][0] + sqrt(sig2 / 2 / M_PI) * exp(-u1 * u1 / 2 / sig2);
                    u_f_n[1][0] = u1 * u_f_n[0][0] - sqrt(sig2 / 2 / M_PI) * exp(-u1 * u1 / 2 / sig2);
                    u_f_p[1][1] = u2 * u_f_p[0][1] + sqrt(sig2 / 2 / M_PI) * exp(-u2 * u2 / 2 / sig2);
                    u_f_n[1][1] = u2 * u_f_n[0][1] - sqrt(sig2 / 2 / M_PI) * exp(-u2 * u2 / 2 / sig2);
                    flux_p(k, i, j, 0) = w1 * u_f_p[1][0] + w2 * u_f_p[1][1];
                    flux_n(k, i, j, 0) = w1 * u_f_n[1][0] + w2 * u_f_n[1][1];
                    // <v ^ m>
                    for (int m = 2;m < 6;++m)
                    {
                        u_f_p[m][0] = u1 * u_f_p[m - 1][0] + (m - 1) * sig2 * u_f_p[m - 2][0];
                        u_f_n[m][0] = u1 * u_f_n[m - 1][0] + (m - 1) * sig2 * u_f_n[m - 2][0];
                        u_f_p[m][1] = u2 * u_f_p[m - 1][1] + (m - 1) * sig2 * u_f_p[m - 2][1];
                        u_f_n[m][1] = u2 * u_f_n[m - 1][1] + (m - 1) * sig2 * u_f_n[m - 2][1];
                        flux_p(k, i, j, m - 1) = w1 * u_f_p[m][0] + w2 * u_f_p[m][1];
                        flux_n(k, i, j, m - 1) = w1 * u_f_n[m][0] + w2 * u_f_n[m][1];
                    }
                }
            }
            // 这里上下的循环必须分开
            // 计算Flux
            for (int m = 0;m < 5;++m)
            {
                for (int j = 0;j < 2 * y_num;++j)
                {
                    //计算Flux_x中间(2*x_num-1)*2*y_num个的
                    for (int i = 1;i < 2 * x_num;++i)
                    {
                        if (is_cosk_p)
                        {
                            Flux_x(k, i, j, m) = flux_p(k, i - 1, j, m) + flux_n(k, i, j, m);
                        }
                        else
                        {
                            Flux_x(k, i, j, m) = flux_n(k, i - 1, j, m) + flux_p(k, i, j, m);
                        }
                    }
                    // 边上两格，Neumann边界条件，直接取为边上格子的正负相加
                    Flux_x(k, 0, j, m) = flux_p(k, 0, j, m) + flux_n(k, 0, j, m);
                    Flux_x(k, 2 * x_num, j, m) = flux_p(k, 2 * x_num - 1, j, m) + flux_n(k, 2 * x_num - 1, j, m);
                }
                for (int i = 0;i < 2 * x_num;++i)
                {
                    // 计算Flux_y中间2*x_num*(2*y_num-1)个的
                    for (int j = 1;j < 2 * y_num;++j)
                    {
                        if (is_sink_p)
                        {
                            Flux_y(k, i, j, m) = flux_p(k, i, j - 1, m) + flux_n(k, i, j, m);
                        }
                        else
                        {
                            Flux_y(k, i, j, m) = flux_n(k, i, j - 1, m) + flux_p(k, i, j, m);
                        }
                    }
                    // 边上两格，Neumann边界条件，直接取为边上格子的正负相加
                    Flux_y(k, i, 0, m) = flux_p(k, i, 0, m) + flux_n(k, i, 0, m);
                    Flux_y(k, i, 2 * y_num, m) = flux_p(k, i, 2 * y_num - 1, m) + flux_n(k, i, 2 * y_num - 1, m);
                }
            }
        }
#pragma endregion

        clock_4 = clock();

#pragma region 第4步，求解差分方程，更新双节点EQMOM参数，处理宏观量等
#pragma omp parallel for
        for (int i = 0;i < 2 * x_num;++i)
        {
            for (int j = 0;j < 2 * y_num;++j)
            {
                // 计算本地温度
                double rt_local = macro(i, j, 3) - 0.5 * (pow(macro(i, j, 1), 2) + pow(macro(i, j, 2), 2));
                // 此时本轮计算宏观量已经没有用了
                macro(i, j, 0) = 0;
                macro(i, j, 1) = 0;
                macro(i, j, 2) = 0;
                macro(i, j, 3) = 0;
                for (int k = 0;k < N;++k)
                {
                    // 得到下一时刻各阶矩
                    if (flag_collide == 1 && kappa == 0)
                    {
                        // 用平衡态计算
                        for (int m = 0;m < 5;++m)
                        {
                            MF(k, i, j, m) = Meq(k, i, j, m)
                                - (x_n[k] * dt / dx) * (Flux_x(k, i + 1, j, m) - Flux_x(k, i, j, m))
                                - (y_n[k] * dt / dy) * (Flux_y(k, i, j + 1, m) - Flux_y(k, i, j, m));
                        }
                    }
                    else
                    {
                        if (flag_collide == 1)
                        {
                            // 补上碰撞项
                            // 计算出本地的tau
                            double tau = tau_pre / sqrt(rt_local);
                            for (int m = 0;m < 5;++m)
                            {
                                // 使用指数型半隐式格式，注意0到4阶都要加
                                MF(k, i, j, m) = exp(-dt / tau) * (
                                    MF(k, i, j, m)
                                    - (x_n[k] * dt / dx) * (Flux_x(k, i + 1, j, m) - Flux_x(k, i, j, m))
                                    - (y_n[k] * dt / dy) * (Flux_y(k, i, j + 1, m) - Flux_y(k, i, j, m))
                                    ) + (1 - exp(-dt / tau)) * Meq(k, i, j, m);
                            }
                        }
                        else
                        {
                            // 无碰撞
                            for (int m = 0;m < 5;++m)
                            {
                                MF(k, i, j, m) = MF(k, i, j, m)
                                    - (x_n[k] * dt / dx) * (Flux_x(k, i + 1, j, m) - Flux_x(k, i, j, m))
                                    - (y_n[k] * dt / dy) * (Flux_y(k, i, j + 1, m) - Flux_y(k, i, j, m));
                            }
                        }
                    }
                    // 0到4阶已经算完
                    // 还要计算出5阶矩，顺便得到新的双节点EQMOM参数
                    invEQ2(MF(k, i, j, 0), MF(k, i, j, 1), MF(k, i, j, 2),
                        MF(k, i, j, 3), MF(k, i, j, 4), &para_f(k, i, j, 0));

#pragma region 计算下一时刻的宏观量
                    // 计算下一时刻的宏观量
                    macro(i, j, 0) += MF(k, i, j, 0);                   // rho
                    macro(i, j, 1) += MF(k, i, j, 1) * x_n[k];   // rho u
                    macro(i, j, 2) += MF(k, i, j, 1) * y_n[k];   // rho v
                    macro(i, j, 3) += MF(k, i, j, 2) * 0.5;             // rho E
                }
                macro(i, j, 1) /= macro(i, j, 0);           // u
                macro(i, j, 2) /= macro(i, j, 0);           // v
                macro(i, j, 3) /= macro(i, j, 0);           // E
#pragma endregion

            }
        }
#pragma endregion

        clock_5 = clock();
        // 一次循环结束，输出当前时间
        cout << "第" << step_count + 1 << "步，" << "t=" << t_current
            << "，计算耗时" << (clock_5 - clock_1) / CLOCKS_PER_SEC << "秒"
            /* << "，算步长" << (clock_2 - clock_1) / CLOCKS_PER_SEC << "秒"
             << "，算平衡态" << (clock_3 - clock_2) / CLOCKS_PER_SEC << "秒"
             << "，算Flux" << (clock_4 - clock_3) / CLOCKS_PER_SEC << "秒"
             << "，求差分方程等" << (clock_5 - clock_4) / CLOCKS_PER_SEC << "秒"*/
            << endl;

#pragma region 看看是不是需要储存
        if (t_current >= Save_Count * save_interval)
        {
            save_t[Save_Count] = t_current;
#pragma omp parallel for
            for (int i = 0;i < 2 * x_num;++i)
            {
                for (int j = 0;j < 2 * y_num;++j)
                {
                    save_r(j, i) = macro(i, j, 0);
                    save_u(j, i) = macro(i, j, 1);
                    save_v(j, i) = macro(i, j, 2);
                    save_e(j, i) = macro(i, j, 3);
                }
            }

            CSVWriteMatrix(foldername + "\\" + "r" + to_string(Save_Count) + ".csv", &save_r(0, 0), 2 * y_num, 2 * x_num);
            CSVWriteMatrix(foldername + "\\" + "u" + to_string(Save_Count) + ".csv", &save_u(0, 0), 2 * y_num, 2 * x_num);
            CSVWriteMatrix(foldername + "\\" + "v" + to_string(Save_Count) + ".csv", &save_v(0, 0), 2 * y_num, 2 * x_num);
            CSVWriteMatrix(foldername + "\\" + "e" + to_string(Save_Count) + ".csv", &save_e(0, 0), 2 * y_num, 2 * x_num);

            clock_end = clock();
            cout << "第" << Save_Count << "次储存，当前t=" << t_current
                << "，计算已用时" << (clock_end - clock_begin) / CLOCKS_PER_SEC << "秒"
                << endl;
            Save_Count++;
        }
#pragma endregion

    }
    clock_end = clock();
#pragma endregion

    cout << "主循环结束，总用时" << (clock_end - clock_begin) / CLOCKS_PER_SEC << "秒" << endl;

#pragma region 储存数据

    CSVWriteArray(foldername + "\\" + "t.csv", &save_t[1], save_max - 1);
    // 最后储存一些基本信息
    ofstream outFile;
    outFile.open(foldername + "\\" + "setting.csv", ios::out); // 打开模式可省略
    outFile << "N," << N << endl;
    outFile << "dx," << dx << endl;
    outFile << "dy," << dy << endl;
    outFile << "flag_collide," << flag_collide << endl;
    outFile << "kappa," << kappa << endl;
    outFile.close();

#pragma endregion

    cout << "数据储存完毕" << endl;
}