#include "../include/twoDRP.h"

void twoDRP_EQMOM(string init_setting, string compute_setting)
{

#pragma region ��ʼ�����������������Ҫ������

    // ���ø�������Ҫ����ƽ���ڱ���ķ���
    SetDirections();

    // �Ե�ǰʱ��ȷ��һ���ļ�����
    string foldername = "EQMOM" + to_string(N) + "d"
        + init_setting + compute_setting + CurrentLocalTime();
    string command;
    command = "mkdir " + foldername;
    system(command.c_str());


    // ������������
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

    // ����ʽ��ʽ��Ҫ��ϵ��
    double tau_pre = kappa;

    // �������ߴ�
    int x_num = int(0.5 / dx);    //  x��������һ��
    int y_num = int(0.5 / dy);    //  y��������һ��

    // ������������

    MF.reConstruct(N, 2 * x_num, 2 * y_num, 5);        //ǰһʱ������(x��y������dvm)
    // ʵ���ϳ���ֻ��Ҫ��һ��MF���ݾ͹��ˣ�������Ϊʹ��EQMOM��5�׾���ȫû���õ�������ȥ
    Meq.reConstruct(N, 2 * x_num, 2 * y_num, 5);       //ƽ��̬
    macro.reConstruct(2 * x_num, 2 * y_num, 4);
    CArray para_f(N, 2 * x_num, 2 * y_num, 5);      //˫�ڵ�EQMOM��ϵ��
    CArray alpha(2 * x_num, 2 * y_num, 4);          //ƽ��̬����
    //CArray alpha(2 * x_num, 2 * y_num, 5);          //ƽ��̬����

    // F
    CArray flux_p(N, 2 * x_num, 2 * y_num, 5);      // ��������֣�1��5��
    CArray flux_n(N, 2 * x_num, 2 * y_num, 5);      // ��������֣�1��5��
    CArray Flux_x(N, 2 * x_num + 1, 2 * y_num, 5);  // x�����һ�񣬸������Ӷ����ҵ�Ӱ�죬1��5��
    CArray Flux_y(N, 2 * x_num, 2 * y_num + 1, 5);  // y�����һ�񣬸������Ӷ����µ�Ӱ�죬1��5��

    // �����õ�
    double* save_t = new double[save_max] {};
    // ע���������Ҫ��x��y������
    CArray save_r(2 * y_num, 2 * x_num);      // rho
    CArray save_u(2 * y_num, 2 * x_num);      // u
    CArray save_v(2 * y_num, 2 * x_num);      // v
    CArray save_e(2 * y_num, 2 * x_num);      // e
#pragma endregion

    cout << "���㷽��EQMOM��������N=" << N << "����ײ����flag_collide=" << flag_collide << "��kappa=" << kappa <<
        "��������" << 2 * x_num << "*" << 2 * y_num << "��CFLMax=" << CFLMax << endl;
    cout << "��������t=" << (save_max - 1) * save_interval << endl;
    cout << "������ʼ����" << endl;

#pragma region �����ʼֵ
    // ��ʼ����
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
    // I,J��������
    for (int I = 0;I < 2;++I)
    {
        for (int J = 0;J < 2;++J)
        {
            double u = u_init[I][J];
            double v = v_init[I][J];
            double rt = p_init[I][J] / r_init[I][J];
            double a[4] = { -1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, -1 / rt };  // �ֱ�Ϊa, b1, b2, c
            //double a[5] = { - 1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, 0, -1 / rt };  // �ֱ�Ϊa, b1, b2, c

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
                        // ���ó�ʼMF
                        double* m = &MF(k, i + I * x_num, j + J * y_num, 0);
                        m[0] = req_i[k];
                        m[1] = req_i[k] * ueq_i[k];
                        m[2] = ueq_i[k] * m[1] + sig2eq_i * m[0];
                        m[3] = ueq_i[k] * m[2] + 2 * sig2eq_i * m[1];
                        m[4] = ueq_i[k] * m[3] + 3 * sig2eq_i * m[2];
                        // ��ʼ˫�ڵ�EQMOM
                        invEQ2(m[0], m[1], m[2], m[3], m[4], &para_f(k, i + I * x_num, j + J * y_num, 0));
                    }
                    // ��ʼ�����
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

    // ��һ�δ���
    CSVWriteMatrix(foldername + "\\" + "r0.csv", &save_r(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "u0.csv", &save_u(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "v0.csv", &save_v(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "e0.csv", &save_e(0, 0), 2 * y_num, 2 * x_num);

#pragma endregion

    cout << "��ʼֵ������ϣ���ѭ����ʼ" << endl;

#pragma region ��ѭ��
    double clock_1, clock_2, clock_3, clock_4, clock_5, clock_begin, clock_end;
    clock_begin = clock();
    while (Save_Count < save_max)
    {
        step_count++;
        clock_1 = clock();

#pragma region ��1����ȷ��ʱ�䲽��
        // ģ��֮ǰ�ĳ���
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
        // CFLȷ��ʱ�䲽�����ο�Fox08������
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
        t_current += dt;    // ��ǰʱ��
        // cout << dt << endl;
#pragma endregion

        clock_2 = clock();

#pragma region ��2�������ݺ����������������ƽ��̬
        // ���������ײ���ο���������һ��
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

#pragma region ��3��������Flux
#pragma omp parallel for
        for (int k = 0;k < N;++k)
        {
            double u_f_n[6][2]{};   // ����Fluxʱ�õ�����ʱ����
            double u_f_p[6][2]{};
            double w1, u1, w2, u2, sig2;
            bool is_cosk_p = x_n[k] > 0;
            bool is_sink_p = y_n[k] > 0;
            // �ȼ���ÿ����������ߵ�Ӱ��
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
            // �������µ�ѭ������ֿ�
            // ����Flux
            for (int m = 0;m < 5;++m)
            {
                for (int j = 0;j < 2 * y_num;++j)
                {
                    //����Flux_x�м�(2*x_num-1)*2*y_num����
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
                    // ��������Neumann�߽�������ֱ��ȡΪ���ϸ��ӵ��������
                    Flux_x(k, 0, j, m) = flux_p(k, 0, j, m) + flux_n(k, 0, j, m);
                    Flux_x(k, 2 * x_num, j, m) = flux_p(k, 2 * x_num - 1, j, m) + flux_n(k, 2 * x_num - 1, j, m);
                }
                for (int i = 0;i < 2 * x_num;++i)
                {
                    // ����Flux_y�м�2*x_num*(2*y_num-1)����
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
                    // ��������Neumann�߽�������ֱ��ȡΪ���ϸ��ӵ��������
                    Flux_y(k, i, 0, m) = flux_p(k, i, 0, m) + flux_n(k, i, 0, m);
                    Flux_y(k, i, 2 * y_num, m) = flux_p(k, i, 2 * y_num - 1, m) + flux_n(k, i, 2 * y_num - 1, m);
                }
            }
        }
#pragma endregion

        clock_4 = clock();

#pragma region ��4��������ַ��̣�����˫�ڵ�EQMOM����������������
#pragma omp parallel for
        for (int i = 0;i < 2 * x_num;++i)
        {
            for (int j = 0;j < 2 * y_num;++j)
            {
                // ���㱾���¶�
                double rt_local = macro(i, j, 3) - 0.5 * (pow(macro(i, j, 1), 2) + pow(macro(i, j, 2), 2));
                // ��ʱ���ּ��������Ѿ�û������
                macro(i, j, 0) = 0;
                macro(i, j, 1) = 0;
                macro(i, j, 2) = 0;
                macro(i, j, 3) = 0;
                for (int k = 0;k < N;++k)
                {
                    // �õ���һʱ�̸��׾�
                    if (flag_collide == 1 && kappa == 0)
                    {
                        // ��ƽ��̬����
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
                            // ������ײ��
                            // ��������ص�tau
                            double tau = tau_pre / sqrt(rt_local);
                            for (int m = 0;m < 5;++m)
                            {
                                // ʹ��ָ���Ͱ���ʽ��ʽ��ע��0��4�׶�Ҫ��
                                MF(k, i, j, m) = exp(-dt / tau) * (
                                    MF(k, i, j, m)
                                    - (x_n[k] * dt / dx) * (Flux_x(k, i + 1, j, m) - Flux_x(k, i, j, m))
                                    - (y_n[k] * dt / dy) * (Flux_y(k, i, j + 1, m) - Flux_y(k, i, j, m))
                                    ) + (1 - exp(-dt / tau)) * Meq(k, i, j, m);
                            }
                        }
                        else
                        {
                            // ����ײ
                            for (int m = 0;m < 5;++m)
                            {
                                MF(k, i, j, m) = MF(k, i, j, m)
                                    - (x_n[k] * dt / dx) * (Flux_x(k, i + 1, j, m) - Flux_x(k, i, j, m))
                                    - (y_n[k] * dt / dy) * (Flux_y(k, i, j + 1, m) - Flux_y(k, i, j, m));
                            }
                        }
                    }
                    // 0��4���Ѿ�����
                    // ��Ҫ�����5�׾أ�˳��õ��µ�˫�ڵ�EQMOM����
                    invEQ2(MF(k, i, j, 0), MF(k, i, j, 1), MF(k, i, j, 2),
                        MF(k, i, j, 3), MF(k, i, j, 4), &para_f(k, i, j, 0));

#pragma region ������һʱ�̵ĺ����
                    // ������һʱ�̵ĺ����
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
        // һ��ѭ�������������ǰʱ��
        cout << "��" << step_count + 1 << "����" << "t=" << t_current
            << "�������ʱ" << (clock_5 - clock_1) / CLOCKS_PER_SEC << "��"
            /* << "���㲽��" << (clock_2 - clock_1) / CLOCKS_PER_SEC << "��"
             << "����ƽ��̬" << (clock_3 - clock_2) / CLOCKS_PER_SEC << "��"
             << "����Flux" << (clock_4 - clock_3) / CLOCKS_PER_SEC << "��"
             << "�����ַ��̵�" << (clock_5 - clock_4) / CLOCKS_PER_SEC << "��"*/
            << endl;

#pragma region �����ǲ�����Ҫ����
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
            cout << "��" << Save_Count << "�δ��棬��ǰt=" << t_current
                << "����������ʱ" << (clock_end - clock_begin) / CLOCKS_PER_SEC << "��"
                << endl;
            Save_Count++;
        }
#pragma endregion

    }
    clock_end = clock();
#pragma endregion

    cout << "��ѭ������������ʱ" << (clock_end - clock_begin) / CLOCKS_PER_SEC << "��" << endl;

#pragma region ��������

    CSVWriteArray(foldername + "\\" + "t.csv", &save_t[1], save_max - 1);
    // ��󴢴�һЩ������Ϣ
    ofstream outFile;
    outFile.open(foldername + "\\" + "setting.csv", ios::out); // ��ģʽ��ʡ��
    outFile << "N," << N << endl;
    outFile << "dx," << dx << endl;
    outFile << "dy," << dy << endl;
    outFile << "flag_collide," << flag_collide << endl;
    outFile << "kappa," << kappa << endl;
    outFile.close();

#pragma endregion

    cout << "���ݴ������" << endl;
}