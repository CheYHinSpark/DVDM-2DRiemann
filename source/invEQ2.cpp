#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <ctime>

#include "../include/invEQ2.h"

using namespace std;

void invEQ2(double m0, double m1, double m2, double m3, double m4, double mpara[])
{
	// 这两个参数要专门设置！！
	double dulim = 1.;
	double dumax = 8;//20.;

	int i;
	double U, Th, q, yta;
	double c1, c2, c3, ctemp;
	int kiter;
	double sig0, sig02, a2;
	double U2, Th2, dv2, dv, sig2, pl;			//仅仅是计算中的中间参数
	double m2s, m3s, Ths, qs, ga;

	if (m0 > 0)
	{
		U = m1 / m0;
		U2 = U * U;
		Th = m2 / m0 - U2;
		q = (m3 - 3. * m0 * U * Th - m0 * U * U2) / m0;
		yta = (m4 - 4. * m3 * U + 6 * m2 * U2 - 3 * m0 * U2 * U2) / m0;
		Th2 = Th * Th;
		if (Th > 0 && yta > Th2 + q * q / Th)		// 将Th=0等边界情形统一归为不可实现！
		{
			if (q != 0.)
			{
				c1 = (yta / Th2 - 3.) / 6.;
				c2 = q * q / 4. / Th / Th2;

				ctemp = c1 * c1 * c1 + c2 * c2;
				if (ctemp >= 1e-2)
				{
					c3 = pow(sqrt(ctemp) + c2, 1. / 3.);
					sig02 = Th * (-c3 + c1 / c3);
				}
				else    // 采用牛顿法，迭代求解根
				{
					kiter = 0;
					a2 = yta - 3. * Th2;
					sig0 = -Th;
					sig02 = sig0 - (2. * sig0 * sig0 * sig0 + a2 * sig0 + q * q) / (6. * sig0 * sig0 + a2);
					while ((2. * sig02 * sig02 * sig02 + a2 * sig02 + q * q < -1e-6 || abs(sig02 - sig0)>1e-6) && kiter < 5000)
					{
						kiter += 1;
						sig0 = sig02;
						sig02 = sig0 - (2. * sig0 * sig0 * sig0 + a2 * sig0 + q * q) / (6. * sig0 * sig0 + a2);
					}
				}
				sig2 = sig02 + Th;

				// Transition Zone：采用Chalons 2017文章中的方法，避免过大的速度节点
				dv = abs(q / sig02);
				if (dv >= dulim)
				{
					pl = dulim + (dumax - dulim) * tanh((dv - dulim) / (dumax - dulim));
					sig2 = Th - abs(q) / pl;
				}

				m2s = m2 - m0 * sig2;
				m3s = m3 - 3. * m1 * sig2;
				Ths = m2s / m0 - U * U;
				qs = m3s / m0 - 3 * U * Ths - U * U2;
				ga = 0.5 * qs / sqrt(qs * qs + 4. * pow(Ths, 3.));
				mpara[0] = m0 * (0.5 + ga);	//w1
				mpara[2] = m0 * (0.5 - ga);	//w2
				mpara[1] = U - sqrt(mpara[2] / mpara[0] * Ths);
				mpara[3] = U + sqrt(mpara[0] / mpara[2] * Ths);
				mpara[4] = sig2;

				// 再加一层处理
				if (mpara[0] < 1e-20)	mpara[1] = 0.;
				if (mpara[2] < 1e-20)	mpara[3] = 0.;
			}
			else if (q == 0. && yta <= 3 * Th2 + 1e-6)  // 意味着程序判定 q==0
			{
				dv2 = 3. * Th2 - yta > 0 ? sqrt(0.5 * (3. * Th2 - yta)) : 0;
				dv = sqrt(dv2);
				mpara[0] = 0.5 * m0;
				mpara[1] = U - dv;
				mpara[2] = 0.5 * m0;
				mpara[3] = U + dv;
				mpara[4] = Th - dv2;
			}
			else
			{
				cerr << "Error: q=0 but yta > 3 * Th^2, not realizable!" << endl;
				exit(1);
			}
		}
		else
		{
			cerr << "Error: m2 to m4 are not realizable !!" << endl;
			exit(1);
		}
	}
	else if (m0 < 0)
	{
		cerr << "Error: m0<0!" << endl;
		exit(1);
	}
	else		// 意味着程序判定m0 = 0
	{
		if (m1 * m2 * m3 * m4)
		{
			cerr << "Error: m0=0 but m1 - m4 are not!" << endl;
			exit(1);
		}
		else
			for (i = 0;i < 5;++i)	mpara[i] = 0.;
	}

	return;
}