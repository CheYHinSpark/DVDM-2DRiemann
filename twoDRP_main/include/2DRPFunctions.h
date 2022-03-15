#ifndef FUNC_H
#define FUNC_H

#include <string>
#include <cstring>

#include <fstream>
#include <sstream>

#include <time.h>

#include "settings.h"

using namespace std;

// ����ָ��
void swapPtr(double** a, double** b);

// �������
double VdotV(double* a, double* b, int n);

// ���������
void MtimesV(double* A, double* b, int n, double* result);

// ���ķ���
int sign(double n);

bool CSVRead(string filename, const int Max, string* name_list, string* value_list);

string CurrentLocalTime();

void CSVWriteArray(string filename, double* varray, int X);

void CSVWriteMatrix(string filename, double* varray, int Y, int X);

string ReadSettings(string csv_name, string default_name);
#endif