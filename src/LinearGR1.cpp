#include "../includes/laba2.h"

using namespace std;

void Solver::LinearGR1() //двухслойный метод для решения уравнения при K=K(x) с граничными условиями 1-го рода 
{                                                              // y0 - гран. усл. 1-го рода на левой стенке, yn - гран. усл. 1-го рода на правой стенке, yx - начальное усл.
	double* y_1 = new double[n + 1];
	double* y = new double[n + 1];
	double* amas = new double[n + 1];
	double* a1 = new double[n + 1]; //Вектора для прогонки
	double* b1 = new double[n + 1];
	double* c1 = new double[n + 1];
	double* d1 = new double[n + 1];
	double* u0 = new double[n + 1];
	double norm = 0;

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	for (int i = 0; i <= n; ++i) {
		y_1[i] = yx(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		u0[i] = u1(coord, time);
		coord += h;
	}
	//cout << "time = " << time << " norm  = " << infnorm(y_1, u0, n + 1) << endl;
	if (infnorm(y_1, u0, n + 1) > norm) norm = infnorm(y_1, u0, n + 1);

	for (int i = 1; i <= n; ++i) {
		amas[i] = a(i * h);
	}
	//OutputMas(amas, n);

	if (sigma == 0)                    // явный двухслойный метод (при sigma = 0)
	{
		cout << "sigma = 0" << endl;
		for (int j = 1; j <= k; ++j) {
			y[0] = y0;
			coord = 0;
			time += tau;
			myfile << time << " " << coord << " " << y[0] << endl;
			y[n] = yn;
			for (int i = 1; i <= n - 1; ++i)
			{
				coord += h;
				y[i] = (amas[i + 1] * (y_1[i + 1] - y_1[i]) - amas[i] * (y_1[i] - y_1[i - 1])) * tau / (h * h * rho * c)  + y_1[i];
				u0[i] = u1(coord, time);
				myfile << time << " " << coord << " " << y[i] << endl;
			}
			coord += h;
			myfile << time << " " << coord << " " << y[n] << endl;
			if (infnorm(y, u0, n) > norm) norm = infnorm(y, u0, n);
			swap(y_1, y);
		}
	}
	if (sigma > 0 && sigma < 1)              // неявный шестишаблонный двухслойный метод (при sigma > 0 && < 1)
	{
		cout << "sigma = " << sigma << endl;
		double rcht = rho * c * h / tau;

		a1[0] = 0; //учет граничных условий 
		b1[0] = 1;
		c1[0] = 0;
		d1[0] = y0;

		a1[n] = 0;
		b1[n] = 1;
		c1[n] = 0;
		d1[n] = yn;

		double sh = sigma / h;
		double sh1 = (1 - sigma) / h;
		for (int i = 1; i <= n - 1; ++i)
		{
			a1[i] = -amas[i] * sh;
			b1[i] = rcht + amas[i + 1] * sh + amas[i] * sh;
			c1[i] = -amas[i + 1] * sh;
		}

		for (int j = 1; j <= k; ++j) {
			coord = 0;
			time += tau;
			for (int i = 1; i <= n - 1; ++i)
			{
				d1[i] = sh1 * amas[i] * y_1[i-1] + sh1 * amas[i + 1] * y_1[i + 1] + (rcht - sh1 * (amas[i + 1] + amas[i])) * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			for (int i = 0; i <= n; ++i)
			{
				myfile << time << " " << coord << " " << y[i] << endl;
				u0[i] = u1(coord, time);
				coord += h;
			}
			if (infnorm(y, u0, n) > norm) norm = infnorm(y, u0, n);
			swap(y_1, y);
		}
	}

	if (sigma == 1)               // неявный двухслойный метод (при sigma = 1)
	{
		cout << "sigma = 1" << endl;
		double rcht = rho * c * h/ tau;

		a1[0] = 0; //учет граничных условий 
		b1[0] = 1;
		c1[0] = 0;
		d1[0] = y0;

		a1[n] = 0;
		b1[n] = 1;
		c1[n] = 0;
		d1[n] = yn;

		for (int i = 1; i <= n - 1; ++i)
		{
			a1[i] = -amas[i] / h;
			b1[i] = rcht + amas[i + 1] / h + amas[i] / h;
			c1[i] = -amas[i + 1] / h;
		}

		for (int j = 1; j <= k; ++j) {
			coord = 0;
			time += tau;
			for (int i = 1; i <= n - 1; ++i)
			{
				d1[i] = rcht * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			for (int i = 0; i <= n; ++i)
			{
				myfile << time << " " << coord << " " << y[i] << endl;
				u0[i] = u1(coord, time);
				coord += h;
			}
			//cout << "time = " << time << " norm  = " << infnorm(y, u0, n + 1) << endl;
			if (infnorm(y, u0, n+1) > norm) norm = infnorm(y, u0, n+1);
			//cout << "time = " << time << " norm  = " << norm << endl;
			//OutputMas(y, n + 1); OutputMas(u0, n + 1);
			swap(y_1, y);
		}
	}
	cout << "Norm error: " << norm << endl;
	delete[] y; delete[] y_1; delete[] amas;
	delete[] a1; delete[] b1; delete[] c1; delete[] d1;
}