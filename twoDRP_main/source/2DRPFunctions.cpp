#include <iostream>
#include "../include/2DRPFunctions.h"


// 交换指针
void swapPtr(double** a, double** b)
{
	double* c = *a;
	*a = *b;
	*b = c;
}

// 向量点乘
double VdotV(double* a, double* b, int n) 
{
	double result = 0;
	for (int i = 0;i < n;++i)
	{
		result += a[i] * b[i];
	}
	return result;
}

// 方阵乘向量
void MtimesV(double* A, double* b, int n, double* result) 
{
	memset(result, 0, sizeof(double) * n);
	for (int i = 0;i < n;++i)
	{
		for (int j = 0;j < n;++j)
		{
			result[i] += A[i * n + j] * b[j];
		}
	}
}

// 数的符号
int sign(double n) 
{
	if (n > 0) { return 1; }
	else if (n < 0) { return -1; }
	return 0;
}

bool CSVRead(string filename, const int Max, string* name_list, string* value_list)
{
	filename += ".txt";
	ifstream csv(filename);
	if (csv.good())
	{
		string temp, var_name, var_value;
		int i = 0;
		while (getline(csv, temp))
		{
			string field;
			istringstream iss(temp);
			getline(iss, field, ',');
			name_list[i] = field;
			getline(iss, field);
			value_list[i] = field;
			i++;
			if (i >= Max) { break; }
		}
		csv.close();
		return true;
	}
	else
	{
		cout << "打开 " << filename << "失败" << endl;
		return false;
	}
}

string CurrentLocalTime()
{
	time_t tt = time(0);
	struct tm* lcl = new tm;
	localtime_s(lcl, &tt);
	return to_string(((lcl->tm_mday) * 100 + lcl->tm_hour) * 100 + lcl->tm_min);
}

void CSVWriteArray(string filename, double* varray, int X)
{
	ofstream outFile;
	outFile.open(filename, ios::out); // 打开模式可省略
	for (int x = 0;x < X;++x)
	{
		outFile << varray[x] << endl;
	}
	outFile.close();
}

void CSVWriteMatrix(string filename, double* varray, int Y, int X)
{
	ofstream outFile;
	outFile.open(filename, ios::out); // 打开模式可省略
	for (int y = 0;y < Y;++y)
	{
		outFile << varray[y * X];
		for (int x = 1;x < X;++x)
		{
			outFile << ',' << varray[y * X + x];
		}
		outFile << endl;
	}
	outFile.close();
}

string ReadSettings(string csv_name, string default_name)
{
	string* nlist = new string[40]{};
	string* vlist = new string[40]{};
	if (CSVRead(csv_name, 40, nlist, vlist))
	{
		for (int i = 0;i < 40;++i)
		{
			string n = nlist[i];
			string v_str = vlist[i];
			double v_dbl = atof(v_str.c_str());
			// 初始设置
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
			// 计算参数
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
		csv_name = default_name;
	}
	delete[] nlist;
	delete[] vlist;
	return csv_name;
}