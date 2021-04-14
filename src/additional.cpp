#include "../includes/laba2.h"

double infnorm(double v1, double v2) {
	double max;
	if (fabs(v1) > fabs(v2)) {
        max = fabs(v1);
    } else {
         max = fabs(v2);
    }
	return max;
}

void SweepMethod(double* a, double* b, double* c, double* d, int k, double* res) {
	double* alpha = new double[k];
	double* beta = new double[k];
	alpha[0] = -c[0] / b[0];                    
	beta[0] = d[0] / b[0];
	for (int i = 1; i < k; ++i)
	{
		alpha[i] = -c[i] / (b[i] + a[i] * alpha[i - 1]);
		beta[i] = (d[i] - a[i] * beta[i - 1]) / (b[i] + a[i] * alpha[i - 1]);
	}

	res[k - 1] = beta[k - 1];                
	for (int i = k - 2; i >= 0; --i)
	{
		res[i] = alpha[i] * res[i + 1] + beta[i];
	}
	//cout << endl;
	//OutputMas(res, n + 1);
	delete[] alpha; delete[] beta;
}

double infnorm(double* vec1, double* vec2, int s)
{
	double max=0;
	for (int i = 0; i < s; ++i)
		if (fabs(vec1[i] - vec2[i]) > max) { max = fabs(vec1[i] - vec2[i]); };
	return max;
}

double Solver::a(double x)
{
	double m = 10;
	double step = h / m;
	double res = 0;
	double x_1 = x - h;

	for (int i = 1; i < m; ++i)
	{
		res += (1.0 / K(x_1 + i*step) + 1.0 / K(x_1 + (i-1) * step)) * step / 2;
	}
	res /= h;
	res = 1 / res;
	return res;
};

double Solver::K(double x)
{
	double res;
	if (x < x1) { res = k1; }
	if (x >= x1 && x <= x2) {
		res = k1 * (x - x2) / (x1 - x2) + k2 * (x - x1) / (x2 - x1);
	}
	if (x > x2) { res = k2; }
	//cout << " inside k " << res;
	return res;
};

double Solver::yx(double x)
{
	return 8 * L * L * sin(pi * x / L) / (pi * pi * pi);
}