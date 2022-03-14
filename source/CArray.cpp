#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <ctime>

#include "../include/CArray.h"
using namespace std;

CArray::CArray(int d, int size_new[]) : dim(d)
{
	size_num = 1;
	size = new int[d] {};
	for (int i = 0; i < d; ++i)
	{
		size[i] = size_new[i];
		size_num *= size_new[i];
	}
	ptr = new double[size_num] {};
	CheckMem();
}

CArray::CArray(CArray& marray)
{
	dim = marray.dim;
	size_num = marray.size_num;
	size = new int[dim] {};
	ptr = new double[size_num] {};
	CheckMem();
	int i;
	for (i = 0; i < dim; ++i)	size[i] = marray.size[i];
	for (i = 0; i < size_num; ++i)	ptr[i] = marray.ptr[i];
}

CArray::CArray(int d0, int d1)
{
	reConstruct(d0, d1);
}

CArray::CArray(int d0, int d1, int d2)
{
	reConstruct(d0, d1, d2);
}

CArray::CArray(int d0, int d1, int d2, int d3)
{
	reConstruct(d0, d1, d2, d3);
}

CArray::CArray(int d0, int d1, int d2, int d3, int d4)
{
	reConstruct(d0, d1, d2, d3, d4);
}

CArray::CArray(int d0, int d1, int d2, int d3, int d4, int d5)
{
	reConstruct(d0, d1, d2, d3, d4, d5);
}

void CArray::reConstruct(int d0, int d1)
{
	size_num = d0 * d1;
	dim = 2;
	size = new int[2]{};
	size[0] = d0;
	size[1] = d1;
	ptr = new double[size_num] {};
	CheckMem();
}

void CArray::reConstruct(int d0, int d1, int d2)
{
	size_num = d0 * d1 * d2;
	dim = 3;
	size = new int[3]{};
	size[0] = d0;
	size[1] = d1;
	size[2] = d2;
	ptr = new double[size_num] {};
	CheckMem();
}

void CArray::reConstruct(int d0, int d1, int d2, int d3)
{
	size_num = d0 * d1 * d2 * d3;
	dim = 4;
	size = new int[4]{};
	size[0] = d0;
	size[1] = d1;
	size[2] = d2;
	size[3] = d3;
	ptr = new double[size_num] {};
	CheckMem();
}

void CArray::reConstruct(int d0, int d1, int d2, int d3, int d4)
{
	size_num = d0 * d1 * d2 * d3 * d4;
	dim = 5;
	size = new int[5]{};
	size[0] = d0;
	size[1] = d1;
	size[2] = d2;
	size[3] = d3;
	size[4] = d4;
	ptr = new double[size_num] {};
	CheckMem();
}

void CArray::reConstruct(int d0, int d1, int d2, int d3, int d4, int d5)
{
	size_num = d0 * d1 * d2 * d3 * d4 * d5;
	dim = 6;
	size = new int[6]{};
	size[0] = d0;
	size[1] = d1;
	size[2] = d2;
	size[3] = d3;
	size[4] = d4;
	size[5] = d5;
	ptr = new double[size_num] {};
	CheckMem();
}

void CArray::CheckDim(int n)
{
	if (dim < n)
	{
		cout << "wrong dimension for CArray" << endl;
		exit(1);
	}
}

void CArray::CheckMem()
{
	if (!ptr)
	{
		cerr << "Error: out of memory!" << endl;
		exit(1);
	}
}