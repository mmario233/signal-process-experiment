#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#define eps 1e-6 
#define PI 3.14159265354
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

//��ɢ���е�DFT����,ֻ���ʵ������ ,��ȫ����DFT�Ĺ�ʽ������,O(n^2)�ĸ��Ӷ�
void DFT(vector<complex_t>& x_n, vector<complex_t>& X_k) {
	X_k.clear();
	int N = x_n.size();
	for (int k = 0; k < N; ++k) {
		complex_t t(0, 0);
		for (int n = 0; n < N; ++n) {
			t += x_n[n] * W(k, n, N);
		}
		X_k.push_back(format(t));
	}

	int cnt = 0;
	cout << "DFT:" << endl;
	for (int i = 0; i < N; ++i) {
		complex_t T = format(X_k[i]);
		cout << T.real() << "+" << T.imag() << "i,";
	}
	cout << endl;

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


//��2,FFT�㷨ʵ��,O(nlogn)�ĸ��Ӷ�
void FFT(vector<complex_t>& x_n) {
	int N = x_n.size();
	int r = bitlen(N);

	vector<complex_t> W(N);

	//Ԥ�ȼ�����ת���� 
	for (int i = 0; i < N; ++i) {
		double angle = -i * 2 * PI / N;
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
	cout << "FFT:" << endl;
	for (int i = 0; i < N; ++i) {
		complex_t T = format(x_n[i]);
		cout << T.real() << "+" << T.imag() << "i,";
	}
	cout << endl;
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


int main() {
	int n = 0;
	bool choose = 0;
	cout << "Input n, n >=6" << endl;//n������ѡ��һ��ȡ6
	cin >> n;
	int N = pow(2, n);//ȡ����2^n��
	vector<complex_t> x_n(N, 0);
	vector<double> x_n1(N, 0);
	vector<complex_t> X_k(N, 0);
	//����ʱ���е����ز�����ʾ
	complex_t st{ 1 / sqrt(2), 1 / sqrt(2) };
	cout << "ingoing signals:" << endl;
	for (int i = 0; i < N; i++) {
		x_n[i] = st;
		cout << x_n[i].real() << "+" << x_n[i].imag() << "i,";
		if (i == N - 1)
			cout << endl;
	}
	IFFT(x_n);
	//����1ʹ��FFT,0ʹ��DFT
	cout << "input 1 for FFT, 0 for DFT" << endl;
	cin >> choose;
	if (choose)
		FFT(x_n);
	else
		DFT(x_n, X_k);
	system(" pause ");
	return 0;

}