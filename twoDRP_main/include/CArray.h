#ifndef CARRAY_H
#define CARRAY_H

using namespace std;

class CArray			// Definition of a multidimensional array. The dimension has to be 2,3,4,5,6.
{
public:
	int dim;			// dim:dimension, and it has to be 2,3,4,5,6.
	int* size;			// size array
	int size_num;		// total size of the array
	double* ptr;		// Pointing to a dy_namic allocated array

	CArray(int d, int size_new[]);				// Constructor
	CArray(CArray& marray);						// Copy constructor
	CArray(int d0, int d1);
	CArray(int d0, int d1, int d2);
	CArray(int d0, int d1, int d2, int d3);
	CArray(int d0, int d1, int d2, int d3, int d4);
	CArray(int d0, int d1, int d2, int d3, int d4, int d5);
	// 需要先保证已经清除了内存
	void reConstruct(int d0, int d1);
	void reConstruct(int d0, int d1, int d2);
	void reConstruct(int d0, int d1, int d2, int d3);
	void reConstruct(int d0, int d1, int d2, int d3, int d4);
	void reConstruct(int d0, int d1, int d2, int d3, int d4, int d5);
	~CArray() { delete[] size; delete[] ptr; }	// Destructor

	double& operator()(int n0, int n1)									// Any index starts from 0.
	{
		CheckDim(2);
		return ptr[n0 * size[1] + n1];
	}
	
	double& operator()(int n0, int n1, int n2)							// Any index starts from 0.
	{
		CheckDim(3);
		return ptr[(n0 * size[1] + n1) * size[2] + n2];
	}
	
	double& operator()(int n0, int n1, int n2, int n3)					// Any index starts from 0.
	{
		CheckDim(4);
		return ptr[((n0 * size[1] + n1) * size[2] + n2) * size[3] + n3];
	}
	
	double& operator()(int n0, int n1, int n2, int n3, int n4)			// Any index starts from 0.
	{
		CheckDim(5);
		return ptr[(((n0 * size[1] + n1) * size[2] + n2) * size[3] + n3) * size[4] + n4];
	}
	
	double& operator()(int n0, int n1, int n2, int n3, int n4, int n5)	// Any index starts from 0.
	{
		CheckDim(6);
		return ptr[((((n0 * size[1] + n1) * size[2] + n2) * size[3] + n3) * size[4] + n4) * size[5] + n5];
	}

private:
	void CheckDim(int n);
	void CheckMem();
};

#endif