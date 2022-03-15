#ifndef FUNC_H
#define FUNC_H

#include <string>
#include <cstring>

#include <fstream>
#include <sstream>

#include <time.h>

#include "settings.h"

using namespace std;

// 交换指针
void swapPtr(double** a, double** b);

// 向量点乘
double VdotV(double* a, double* b, int n);

// 方阵乘向量
void MtimesV(double* A, double* b, int n, double* result);

// 数的符号
int sign(double n);

bool CSVRead(string filename, const int Max, string* name_list, string* value_list);

string CurrentLocalTime();

void CSVWriteArray(string filename, double* varray, int X);

void CSVWriteMatrix(string filename, double* varray, int Y, int X);

string ReadSettings(string csv_name, string default_name);
#endif