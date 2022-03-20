#include "../include/twoDRP.h"

void twoDRP_DVDDVM(string init_setting, string compute_setting, int mode)
{

#pragma region ��ʼ�����������������Ҫ������

    // ���ø�������Ҫ����ƽ���ڱ���ķ���
    SetDirections();

    // �������dvm
    SetDVM();

    // �Ե�ǰʱ��ȷ��һ���ļ�����
    string foldername = "";
    if (mode == 1)
    { 
        foldername += "DVM";
    }
    else if (mode == 2)
    {
        foldername += "DVDDVM";
    }
    foldername += to_string(N) + "x" + to_string(M)
        + init_setting + compute_setting + CurrentLocalTime() + SCHEME;
    string command;
    command = "mkdir " + foldername;
    system(command.c_str());

    // �������ߴ�
    int x_num = int(0.5 / dx);    //  x��������һ��
    int y_num = int(0.5 / dy);    //  y��������һ��
    int grid_num = 4 * x_num * y_num;

    // ������������
    SetUVW();

    MF.reConstruct(2 * x_num, 2 * y_num, N, M);     //ǰһʱ������(x��y������dvm)
    MF_c.reConstruct(2 * x_num, 2 * y_num, N, M);   //��������(x��y������dvm)
    Meq.reConstruct(2 * x_num, 2 * y_num, N, M);    //ƽ��̬
    CArray alpha(2 * x_num, 2 * y_num, 4);          //ƽ��̬����
    macro.reConstruct(2 * x_num, 2 * y_num, 4);     //�����
    if (flag_collide == 1 && kappa != 0 && SCHEME == "DUGKS")
    {
        MF_m.reConstruct(2 * x_num + 1, 2 * y_num + 1, N, M);
    }

    // �����õ�
    double* save_t = new double[save_max] {};
    // ע���������Ҫ��x��y������
    CArray save_r(2 * y_num, 2 * x_num);      // rho
    CArray save_u(2 * y_num, 2 * x_num);      // u
    CArray save_v(2 * y_num, 2 * x_num);      // v
    CArray save_e(2 * y_num, 2 * x_num);      // e
#pragma endregion

    if (mode == 1)
    {
        cout << "���㷽��DVM��";
    }
    else if (mode == 2)
    {
        cout << "���㷽��DVDDVM��";
    }
    cout << "������N=" << N << "����ɢ�ٶȸ���" << M << "����ײ����flag_collide=" << flag_collide <<
        "��kappa=" << kappa << "��������" << 2 * x_num << "*" << 2 * y_num << "��CFLMax=" << CFLMax << endl;
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
    double* MFi = new double[M * N]{};
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

            if (mode == 1)
            {
                RedEquil_DVM(r_init[I][J], u_init[I][J], v_init[I][J], E_init[I][J],
                    a, MFi);
            }
            else if (mode == 2)
            {
                RedEquil_DVD(r_init[I][J], u_init[I][J], v_init[I][J], E_init[I][J],
                    a, req_i, ueq_i, &sig2eq_i);
                for (int n = 0;n < N;++n)
                {
                    RedEquil_DVDDVM(req_i[n], ueq_i[n], 0.5 * (ueq_i[n] * ueq_i[n] + sig2eq_i),
                        &MFi[n * M]);
                }
            }
            for (int i = 0;i < x_num;++i)
            {
                for (int j = 0;j < y_num;++j)
                {
                    for (int k = 0;k < 4;++k)
                        //for (int k = 0;k < 5;++k)
                    {
                        alpha(i + I * x_num, j + J * y_num, k) = a[k];
                    }
                    for (int k = M * N - 1;k >= 0;--k)
                    {
                        MF(i + I * x_num, j + J * y_num, 0, k) = MFi[k];
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
    delete[] MFi;

    // ��һ�δ���
    CSVWriteMatrix(foldername + "\\" + "r0.csv", &save_r(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "u0.csv", &save_u(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "v0.csv", &save_v(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "e0.csv", &save_e(0, 0), 2 * y_num, 2 * x_num);

#pragma endregion

    //exit(0);
    cout << "��ʼֵ������ϣ���ѭ����ʼ" << endl;

#pragma region ��ѭ��
    double clock_1, clock_2, clock_3, clock_4, clock_5, clock_begin, clock_end;
    clock_begin = clock();
    while (Save_Count < save_max)
    {
        step_count++;
        clock_1 = clock();

#pragma region ��1����ȷ��ʱ�䲽��
        // Guo, ÿ�����Զ���
        double dt = 0.5* min(dx, dy) / abs(dvm_bg);
        t_current += dt;    // ��ǰʱ��
        // cout << dt << endl;
#pragma endregion

        clock_2 = clock();

#pragma region ��2�������ݺ����������������ƽ��̬
        // ���������ײ���ο���������һ��
        if (flag_collide == 1)
        {
            int finish_count = 0;
            if (mode == 1)
            {
#pragma omp parallel for
                for (int i = grid_num - 1; i >= 0;i--)
                {
                    RedEquil_DVM(macro(0, i, 0), macro(0, i, 1), macro(0, i, 2), macro(0, i, 3),
                        &alpha(0, i, 0), &Meq(0, i, 0, 0));
                    finish_count++;
                    printf("����ƽ��̬������� %d / %d\r", finish_count, grid_num);
                }
            }
            else if (mode == 2)
            {
#pragma omp parallel for
                for (int i = grid_num - 1; i >= 0;i--)
                {
                    double* req = new double[N] {};
                    double* ueq = new double[N] {};
                    double sig2eq = 0;
                    RedEquil_DVD(macro(0, i, 0), macro(0, i, 1), macro(0, i, 2), macro(0, i, 3),
                        &alpha(0, i, 0), req, ueq, &sig2eq);
                    for (int n = 0;n < N;++n)
                    {
                        RedEquil_DVDDVM(req[n], ueq[n], 0.5 * (ueq[n] * ueq[n] + sig2eq), &Meq(0, i, n, 0));
                    }
                    delete[] req;
                    delete[] ueq;
                    finish_count++;
                    printf("����ƽ��̬������� %d / %d\r", finish_count, grid_num);
                }
            }
        }
#pragma endregion

        clock_3 = clock();

#pragma region ��3��������

        swapPtr(&MF.ptr, &MF_c.ptr);

#pragma endregion

        clock_4 = clock();

#pragma region ��4��������ַ��̣�����˫�ڵ�EQMOM����������������

        // ӭ���ʽ
        if (SCHEME == "wind-up") { D_wind_up(2 * x_num, 2 * y_num, dt); }
        // DUGKS��ʽ
        if (SCHEME == "DUGKS") { DUGKS(2 * x_num, 2 * y_num, dt); }
        // TODO ����DUGKS��ʽ

#pragma endregion

        clock_5 = clock();
        // һ��ѭ�������������ǰʱ��
        cout << "��" << step_count + 1 << "����" << "t=" << t_current
            << "�������ʱ" << (clock_5 - clock_1) / CLOCKS_PER_SEC << "��"
            /* << "���㲽��" << (clock_2 - clock_1) / CLOCKS_PER_SEC << "��"
             << "����ƽ��̬" << (clock_3 - clock_2) / CLOCKS_PER_SEC << "��"
             << "����Flux" << (clock_4 - clock_3) / CLOCKS_PER_SEC << "��"
             << "�����ַ��̵�" << (clock_5 - clock_4) / CLOCKS_PER_SEC << "��" */
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
    CSVWriteArray(foldername + "\\" + "dvm.csv", &dvm[0], M);
    // ��󴢴�һЩ������Ϣ
    ofstream outFile;
    outFile.open(foldername + "\\" + "setting.csv", ios::out); // ��ģʽ��ʡ��
    outFile << "N," << N << endl;
    outFile << "dvm_strategy," << dvm_strategy << endl;
    outFile << "dvm_bg," << dvm_bg << endl;
    outFile << "dvm_st," << dvm_st << endl;
    outFile << "dx," << dx << endl;
    outFile << "dy," << dy << endl;
    outFile << "flag_collide," << flag_collide << endl;
    outFile << "kappa," << kappa << endl;
    outFile.close();

#pragma endregion

    cout << "���ݴ������" << endl;
}
