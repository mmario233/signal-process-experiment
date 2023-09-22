#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include "mex.h"
#define eps 1e-18
#define PI 3.14159265359
typedef std::complex<double> complex_t;
using namespace std;

//��ת���ӵļ��� 
complex_t W(int k, int n, int N) {
	return complex_t(cos(2 * PI * k * n / N), -sin(2 * PI * k * n / N));
}

//��ʽ�� �� 
complex_t format(complex_t& c) {
	if (fabs(c.real()) < eps) c.real(0);
	if (fabs(c.imag()) < eps) c.imag(0);
	return c;
}

double format(double& c) {
	if (fabs(c) < eps) c = 0;
	return c;
}

//��֤N��2��n����
int bitlen(int N) {
	int n = 0;
	while ((N & 1) == 0) {
		n++;
		N >>= 1;
	}
	return n;
}


int reverse_bit(int n, int len) {//bit��ת 
	int tmp = 0;
	while (len--) {
		tmp += ((n & 1) << len);
		n >>= 1;
	}
	return tmp;

}

//�������� 
void resort(vector<complex_t>& x_n, int N) {
	vector<complex_t> v(x_n);
	int len = bitlen(N);
	for (int i = 0; i < N; ++i) {
		x_n[i] = v[reverse_bit(i, len)];
	}
}

//IFFT,��FFT����һ�� 
void IFFT(vector<complex_t>& x_n) {
	int N = x_n.size();
	int r = bitlen(N);

	vector<complex_t> W(N);

	//Ԥ�ȼ�����ת���� 
	for (int i = 0; i < N; ++i) {
		double angle = i * 2 * PI / N;//IFFT����ת������FFT����ת���Ӳ�һ������ 
		W[i] = complex_t(cos(angle), sin(angle));
	}


	for (int k = 0; k < r; ++k) {//�������� 
		for (int j = 0; j < (1 << k); ++j) {
			int butterfly = 1 << (r - k);
			int p = j * butterfly;
			int s = p + butterfly / 2;
			for (int i = 0; i < butterfly / 2; ++i) {
				complex_t c = x_n[i + p] + x_n[i + s];
				x_n[i + s] = (x_n[i + p] - x_n[i + s]) * W[i * (1 << k)];
				x_n[i + p] = c;
			}
		}
	}

	//�������� 
	resort(x_n, N);
	int cnt = 0;
	cout << "IFFT:" << endl;
	for (int i = 0; i < N; ++i) {
		x_n[i] /= N;//IFFT��FFT����һ��ϵ�� 
		complex_t T = format(x_n[i]);
		cout << T.real() << "+" << T.imag() << "i,";
	}
	cout << endl;
}
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
//     mexPrintf("�Ƿ�Ϊ��%d\n�Ƿ�Ϊ����%d\n",mxIsEmpty(prhs[0]),mxIsComplex(prhs[0]));
//     mexPrintf("��%d\n��%d\n",mxGetM(prhs[0]),mxGetN(prhs[0]));
    double *inDatar;
    double *inDatai;
    inDatar = mxGetPr(prhs[0]);
    inDatai = mxGetPi(prhs[0]);
    vector<complex_t> y_f(mxGetN(prhs[0]),0);
    for (int i = 0;i < mxGetN(prhs[0]);i ++){
        complex_t t {inDatar[i], inDatai[i]};
        y_f[i] = t;
    }
    IFFT(y_f);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),mxCOMPLEX);
    double *zr = mxGetPr(plhs[0]);
    double *zi = mxGetPi(plhs[0]);
    for (int i = 0;i < mxGetN(prhs[0]);i ++){
//         mexPrintf("\%f+%fi\n",y_f[i].real(),y_f[i].imag());
        zr[i] =  y_f[i].real();
        zi[i] =  y_f[i].imag();
    }
}