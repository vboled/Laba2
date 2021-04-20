#include "../includes/Lab2.h"

double const eps = 1e-3;
double const L = 1; //ширина стержня
double const T = 1; //время
int const n = 50; //количество разбиений по пространству
int const k = 50; //количество разбиений по времени
double const u0 = 0.5;
double const Y0 = u0; //гран условия 1-го рода (задана температура) на левой стенке
double const Yn = u0; // на правой стенке  
double ybegin(double x) { return u0 + x*(L - x); }; //начальное условие (пример 1, 2)
double p0(double t) { return 0; }; //гран условия 2-го рода (задан закон Фурье) на левой стенке
double Q = 10;
double t0 = 0.5;
double p1(double t) { if (t < t0) return Q; else return 0; };
double p4(double t) { 
	double res;
	if (t <= 0.5 * t0) res =  2 * Q * t; 
	if (t > 0.5 * t0 && t < t0) res = 2 * Q * (t0 - t); 
	if (t >= t0) res =  0; 
	//cout << res << " ";
	return res;
};
double const h = L / n;
double const tau = T / k;
double const x1 = 0.5;// h* (round(n / 2) - ceil(n / 4));
double const x2 = 0.6;// h* (round(n / 2) + ceil(n / 4));
double const k1 = 0.5;
double const k2 = 1.5;
double const rho = 0.75;
double const c = 1;
double const pi = 3.141592653589793238;
double const e = 2.71828182845904523;

////////////////////////////////////////////////////////////////////////////
//тест 3
double sigma3 = 2;
double kappa3 = 0.5;
double kappa_1 = 1 / kappa3;
double c3 = 5;
double dop = 1 / sigma3;
double u03 = pow(sigma3 * c3 * c3 / kappa3, dop);
//double const alp = 0;
//double const bet = kappa3;
//double const gam = sigma3;

//////////////////////////////////////////////////////////////////////////
//Вариант 5
double const alp = 2;
double const bet = 1.5;
double const gam = 2;
double y5 = u0;

double const constU0 = 0; //граничное условие (3 пример)
double K_u(double u) { return alp + bet * pow(u, gam); }
double y0x_u(double x) { return u0; } //начальные условия
double y0t_u(double t) { return u03 * pow(t, dop); } //граничные условия
double yLt_u(double t) { if (1 < c3 * t) return pow(sigma3 * c3 * 2 * (c3 * t - 1), dop); else return 0; } //граничные условия
double a_u(double yi, double yi_1) { return 0.5 * (K_u(yi) + K_u(yi_1)); } //входит в разностную схему (для решения уравнения)
double u(double x, double t) { if (x < c3 * t) return pow(sigma3 * c3 * 2 * (c3 * t - x), dop); else return 0; };

double K(double x)
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

double a(double x)
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

double u1(double x, double t) //Точное решение теста 1
{
	double k = 2.5;
	return pow(e, -k * pi * pi * t / L / L) * 8 * L * L * sin(pi * x / L) / (pi * pi * pi);
}; 

double u0_1(double x)
{
	return 8 * L * L * sin(pi * x / L) / (pi * pi * pi);
}

int main(int argc, char** argv)
{
	setlocale(LC_ALL, "Russian");
	std::cout.precision(8);
	ClearFile();

	//LinearGR1(Y0, Yn, u0_1, 0.8);
	//LinearGR2(p1,p1,ybegin,1);
	//KvaziImpl(y0x_u, y0t_u, yLt_u, a_u);
	//KvaziImplGR1GR2(y0x_u, y5, p1, a_u);
	//KvaziImplGR2GR1(y0x_u, p4, y5, a_u);
	//LinearGR1GR2(y5, p1, y0x_u, 1);
	LinearGR2GR1(p4, y5, y0x_u, 1);

	cin.get();
	return 0;
}

void ClearFile()
{
	ofstream rewrite("Temperature.dat");
	rewrite.write("", 0);
	rewrite.close();
}

void OutputMas(vector<double>v, int m)
{
	for (int i = 0; i < m; ++i)
		cout << v[i] << " ";
	cout << endl;
}

void multmatrvec(double** mas, vector<double>vec, int s, vector<double>res) //перемножение матрицы и вектора
{
	for (int k = 0; k < s; ++k)
		res[k] = 0;

	for (int i = 0; i < s; ++i)
		for (int j = 0; j < s; ++j)
		{
			res[i] += mas[i][j] * vec[j];
		}
}

double euclnorm(vector<double>v, int m) // Евклидова норма
{
	double res = 0;
	for (int i = 0; i < m; ++i)
	{
		res += v[i] * v[i];
	}
	return sqrt(res);
}

double infnorm(double v1, double v2)
{
	double max;
	if (fabs(v1) > fabs(v2)) { max = fabs(v1); }
	else { max = fabs(v2); }
	return max;
}

void Difference(vector<double>vec1, vector<double>vec2, int s, vector<double>res) //вычисление разности векторов
{
	for (int i = 0; i < s; ++i)
		res[i] = vec1[i] - vec2[i];
}

void freematrix(double** Matrix, int Rows) {
	for (int i = 0; i < Rows; i++) {

		delete[] Matrix[i];
		Matrix[i] = NULL;
	}
	delete[] Matrix;
	Matrix = NULL;
}

void SweepMethod(vector<double>a, vector<double>b, vector<double>c, vector<double>d, int k, vector<double>res)  //Метод правой прогонки
{
	vector<double>alpha(k);
	vector<double>beta(k);
	alpha[0] = -c[0] / b[0];                     // Прямой ход 
	beta[0] = d[0] / b[0];
	for (int i = 1; i < k; ++i)
	{
		alpha[i] = -c[i] / (b[i] + a[i] * alpha[i - 1]);
		beta[i] = (d[i] - a[i] * beta[i - 1]) / (b[i] + a[i] * alpha[i - 1]);
	}

	res[k - 1] = beta[k - 1];                   // Обратный ход
	for (int i = k - 2; i >= 0; --i)
	{
		res[i] = alpha[i] * res[i + 1] + beta[i];
	}
	//cout << endl;
	//OutputMas(res, n + 1);
}

void LinearGR1(double y0, double yn, function yx, double sigma) //двухслойный метод для решения уравнения при K=K(x) с граничными условиями 1-го рода 
{                                                              // y0 - гран. усл. 1-го рода на левой стенке, yn - гран. усл. 1-го рода на правой стенке, yx - начальное усл.
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>amas(n + 1);
	vector<double>a1(n + 1); //Вектора для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);
	vector<double>u0(n + 1);
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
	cout << "норма ошибки = " << norm << endl;
	cout << "результаты записаны! " << endl;
}

void LinearGR2(function y0, function yn, function yx, double sigma) //двухслойный метод для решения уравнения при K=K(x) с граничными условиями 2-го рода 
{                                                                   // y0 - гран. усл. 2-го рода на левой стенке, yn - гран. усл. 2-го рода на правой стенке, yx - начальное усл.
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>amas(n + 1);
	vector<double>a1(n + 1); //Вектора для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);
	vector<double>integral(k + 1);

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	for (int i = 0; i <= n; ++i) {  // начальное условие
		y_1[i] = yx(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	integral[0] = (y_1[0] + y_1[n]) / 2;
	for (int i = 1; i <= n - 1; ++i)
	{
		integral[0] += y_1[i];
	}
	integral[0] *= h;
	cout << "time = " << time << "  int = " << integral[0] << endl;

	for (int i = 1; i <= n; ++i) {
		amas[i] = a(i * h);
	}
	//OutputMas(amas, n);

	if (sigma == 0)                    // явный двухслойный метод (при sigma = 0)
	{
		cout << "sigma = 0" << endl;
		double crht = c * rho * h / (2 * tau);
		for (int j = 1; j <= k; ++j) {
			coord = 0;
			time += tau;
			y[0] = (y0(time-tau) + (crht - amas[1] / h) * y_1[0] + amas[1] / h * y_1[1]) / crht; // с учетом гр. условий на левой границе
			myfile << time << " " << coord << " " << y[0] << endl;
			for (int i = 1; i <= n - 1; ++i)
			{
				coord += h;
				y[i] = (amas[i + 1] * (y_1[i + 1] - y_1[i]) - amas[i] * (y_1[i] - y_1[i - 1])) / (h * h * rho * c) * tau + y_1[i];
				myfile << time << " " << coord << " " << y[i] << endl;
			}
			y[n] = (yn(time-tau) + (crht - amas[n] / h) * y_1[n] + amas[n] / h * y_1[n - 1]) / crht;
			coord += h;
			myfile << time << " " << coord << " " << y[n] << endl;
			swap(y_1, y);
		}
	}
	if (sigma > 0 && sigma < 1)              // неявный шестишаблонный двухслойный метод (при sigma > 0 && < 1)
	{
		cout << "sigma = " << sigma << endl;
		double rcht = rho * c * h / tau;
		double sh = sigma / h;
		double sh1 = (1 - sigma) / h;

		a1[0] = 0; //учет граничных условий 
		b1[0] = rcht / 2 + amas[1] * sh;
		c1[0] = -amas[1] * sh;

		a1[n] = -amas[n] * sh;
		b1[n] = rcht / 2 + amas[n] * sh;
		c1[n] = 0;

		for (int i = 1; i <= n - 1; ++i)
		{
			a1[i] = -amas[i] * sh;
			b1[i] = rcht + amas[i + 1] * sh + amas[i] * sh;
			c1[i] = -amas[i + 1] * sh;
		}

		for (int j = 1; j <= k; ++j) {
			coord = 0;
			time += tau;
			d1[0] = - sigma * y0(time) - (1 - sigma) * y0(time - tau) + (rcht / 2 - amas[1] * sh1) * y_1[0] + amas[1] * y_1[1] * sh1;
			d1[n] = sigma * yn(time) + (1 - sigma) * yn(time - tau) + (rcht / 2 - amas[n] * sh1) * y_1[n] + amas[n] * y_1[n-1] * sh1;
			for (int i = 1; i <= n - 1; ++i)
			{
				d1[i] = sh1 * amas[i] * y_1[i] + sh1 * amas[i + 1] * y_1[i + 1] + (rcht - sh1 * (amas[i + 1] + amas[i])) * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			for (int i = 0; i <= n; ++i)
			{
				myfile << time << " " << coord << " " << y[i] << endl;
				coord += h;
			}
			swap(y_1, y);
		}
	}

	if (sigma == 1)               // неявный двухслойный метод (при sigma = 1)
	{
		cout << "sigma = 1" << endl;
		double rcht = rho * c * h / tau;

		a1[0] = 0; //учет граничных условий 
		b1[0] = rcht / 2 + amas[1] / h;
		c1[0] = -amas[1] / h;

		a1[n] = -amas[n] / h;
		b1[n] = rcht/2+amas[n]/h;
		c1[n] = 0;

		for (int i = 1; i <= n - 1; ++i)
		{
			a1[i] = -amas[i] / h;
			b1[i] = rcht + amas[i + 1] / h + amas[i] / h;
			c1[i] = -amas[i + 1] / h;
		}

		for (int j = 1; j <= k; ++j) {
			coord = 0;
			time += tau;
			d1[0] = y0(time) + rcht / 2 * y_1[0];
			d1[n] = yn(time) + rcht / 2 * y_1[n];
			for (int i = 1; i <= n - 1; ++i)
			{
				d1[i] = rcht * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			for (int i = 0; i <= n; ++i)
			{
				myfile << time << " " << coord << " " << y[i] << endl;
				coord += h;
			}
			swap(y_1, y);
			integral[j] = (y[0] + y[n]) / 2;
			for (int i = 1; i <= n - 1; ++i)
			{
				integral[j] += y[i];
			}
			integral[j] *= h;
			cout << "time = " << time <<" текущий " << integral[j] << "  int = " << fabs(integral[0] - integral[j]) << endl;
		}
		
		
	}

	cout << "результаты записаны! " << endl;
}

void KvaziImpl(function y0_x, function y0_t, function yL_t, function2 a_u, int M) //неявный двухслойный метод для решения уравнения при K=K(x) с граничными условиями 1-го рода 
{
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>a1(n + 1); //Векторы для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);
	vector<double>y_1dop(n + 1);
	vector<double>u0(n + 1);
	double norm = 0;

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	for (int i = 0; i <= n; ++i) { //задаем н.у.: начальный момент времени
		y_1[i] = y0_x(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	double rcht = rho * c * h / tau;

	a1[0] = 0; //учет граничных условий 
	b1[0] = 1;
	c1[0] = 0;

	a1[n] = 0;
	b1[n] = 1;
	c1[n] = 0;

	for (int j = 1; j <= k; ++j) {
		coord = 0;
		time += tau;
		d1[0] = y0_t(time);
		d1[n] = yL_t(time);
		for (int s = 0; s <= M; ++s)
		{
			for (int i = 1; i <= n - 1; ++i)
			{
				a1[i] = -a_u(y_1[i], y_1[i - 1]) / h;
				b1[i] = rcht + a_u(y_1[i + 1], y_1[i]) / h + a_u(y_1[i], y_1[i - 1]) / h;
				c1[i] = -a_u(y_1[i + 1], y_1[i]) / h;
				d1[i] = rcht * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			swap(y_1, y);
		}
		for (int i = 0; i <= n; ++i)
		{
			myfile << time << " " << coord << " " << y[i] << endl;
			u0[i] = u(coord, time);
			coord += h;
		}
		//OutputMas(y,n+1); OutputMas(u0, n + 1);
		if (infnorm(y, u0, n) > norm) norm = infnorm(y, u0, n);
	}

	cout << "норма ошибки = " << norm << endl;
	cout << "результаты записаны! " << endl;
}

void KvaziImpl(function y0_x, function y0_t, function yL_t, function2 a_u) //неявный двухслойный метод для решения уравнения при K=K(x) с граничными условиями 1-го рода 
{
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>a1(n + 1); //Векторы для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);
	vector<double>y_1dop(n + 1);
	vector<double>u0(n + 1);
	int count = 0;
	double norm;

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	ofstream myfile2;
	myfile2.open("Time.dat");
	myfile2 << 0 << " " << 0 << endl;
	//cout << " time = " << time << "  count = " << count << endl;
	for (int i = 0; i <= n; ++i) { //задаем н.у.: начальный момент времени
		y_1[i] = y0_x(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	double rcht = rho * c * h / tau;

	a1[0] = 0; //учет граничных условий 
	b1[0] = 1;
	c1[0] = 0;

	a1[n] = 0;
	b1[n] = 1;
	c1[n] = 0;

	for (int j = 1; j <= k; ++j) {
		coord = 0;
		time += tau;
		count = 0;
		d1[0] = y0_t(time);
		d1[n] = yL_t(time);
		do
		{
			for (int i = 1; i <= n - 1; ++i)
			{
				a1[i] = -a_u(y_1[i], y_1[i - 1]) / h;
				b1[i] = rcht + a_u(y_1[i + 1], y_1[i]) / h + a_u(y_1[i], y_1[i - 1]) / h;
				c1[i] = -a_u(y_1[i + 1], y_1[i]) / h;
				d1[i] = rcht * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			norm = infnorm(y, y_1, n + 1);
			swap(y_1, y);
			count++;
		} while (norm>eps);
		//cout << " time = " << time << "  count = " << count << endl;
		myfile2 << time << " " << count << endl;
		for (int i = 0; i <= n; ++i)
		{
			myfile << time << " " << coord << " " << y[i] << endl;
			u0[i] = u(coord, time);
			coord += h;
		}
		//OutputMas(y,n+1); OutputMas(u0, n + 1);
		if (infnorm(y, u0, n) > norm) norm = infnorm(y, u0, n);
	}

	cout << "норма ошибки = " << norm << endl;
	cout << "результаты записаны! " << endl;
}

double infnorm(vector<double>vec1, vector<double>vec2, int s)
{
	double max=0;
	for (int i = 0; i < s; ++i)
		if (fabs(vec1[i] - vec2[i]) > max) { max = fabs(vec1[i] - vec2[i]); };
	return max;
}

void KvaziImplGR1GR2(function y0_x, double y0, function yn, function2 a_u, int M) //неявный двухслойный метод для решения уравнения при K=K(x) с граничными условиями 1-го рода 
{
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>a1(n + 1); //Векторы для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);
	vector<double>y_1dop(n + 1);
	vector<double>u0(n + 1);
	double norm = 0;

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	for (int i = 0; i <= n; ++i) { //задаем н.у.: начальный момент времени
		y_1[i] = y0_x(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	double rcht = rho * c * h / tau;

	a1[0] = 0; //учет граничных условий 
	b1[0] = 1;
	c1[0] = 0;

	a1[n] = -a_u(y_1[n], y_1[n - 1]) / h;
	b1[n] = rcht + a_u(y_1[n], y_1[n - 1]) / h;
	c1[n] = 0;

	for (int j = 1; j <= k; ++j) {
		coord = 0;
		time += tau;
		d1[0] = y0;
		d1[n] = yn(time);
		for (int s = 0; s <= M; ++s)
		{
			for (int i = 1; i <= n - 1; ++i)
			{
				a1[i] = -a_u(y_1[i], y_1[i - 1]) / h;
				b1[i] = rcht + a_u(y_1[i + 1], y_1[i]) / h + a_u(y_1[i], y_1[i - 1]) / h;
				c1[i] = -a_u(y_1[i + 1], y_1[i]) / h;
				d1[i] = rcht * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			swap(y_1, y);
		}
		for (int i = 0; i <= n; ++i)
		{
			myfile << time << " " << coord << " " << y[i] << endl;
			u0[i] = u(coord, time);
			coord += h;
		}
	}

	cout << "результаты записаны! " << endl;
}

void KvaziImplGR2GR1(function y0_x, function y0, double yn, function2 a_u, int M) //неявный двухслойный метод для решения уравнения при K=K(x) с граничными условиями 1-го рода 
{
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>a1(n + 1); //Векторы для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);
	vector<double>y_1dop(n + 1);
	vector<double>u0(n + 1);
	double norm = 0;

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	for (int i = 0; i <= n; ++i) { //задаем н.у.: начальный момент времени
		y_1[i] = y0_x(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	double rcht = rho * c * h / tau;

	a1[n] = 0; //учет граничных условий 
	b1[n] = 1;
	c1[n] = 0;

	a1[0] = 0;
	b1[0] = rcht + a_u(y_1[1], y_1[0]) / h;
	c1[0] = -a_u(y_1[0], y_1[1]) / h;

	for (int j = 1; j <= k; ++j) {
		coord = 0;
		time += tau;
		d1[0] = y0(time);
		d1[n] = yn;
		for (int s = 0; s <= M; ++s)
		{
			for (int i = 1; i <= n - 1; ++i)
			{
				a1[i] = -a_u(y_1[i], y_1[i - 1]) / h;
				b1[i] = rcht + a_u(y_1[i + 1], y_1[i]) / h + a_u(y_1[i], y_1[i - 1]) / h;
				c1[i] = -a_u(y_1[i + 1], y_1[i]) / h;
				d1[i] = rcht * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			swap(y_1, y);
		}
		for (int i = 0; i <= n; ++i)
		{
			myfile << time << " " << coord << " " << y[i] << endl;
			u0[i] = u(coord, time);
			coord += h;
		}
	}

	cout << "результаты записаны! " << endl;
}

void KvaziImplGR1GR2(function y0_x, double y0, function yn, function2 a_u) //неявный двухслойный метод для решения уравнения при K=K(x) с граничными условиями 1-го рода 
{
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>a1(n + 1); //Векторы для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);
	vector<double>y_1dop(n + 1);
	vector<double>u0(n + 1);
	double norm = 0;
	int count = 0;

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	ofstream myfile2;
	myfile2.open("Time.dat");
	for (int i = 0; i <= n; ++i) { //задаем н.у.: начальный момент времени
		y_1[i] = y0_x(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	double rcht = rho * c * h / tau;

	a1[0] = 0; //учет граничных условий 
	b1[0] = 1;
	c1[0] = 0;

	a1[n] = -a_u(y_1[n], y_1[n - 1]) / h;
	b1[n] = rcht + a_u(y_1[n], y_1[n - 1]) / h;
	c1[n] = 0;

	for (int j = 1; j <= k; ++j) {
		coord = 0;
		time += tau;
		count = 0;
		d1[0] = y0;
		d1[n] = yn(time);
		do
		{
			for (int i = 1; i <= n - 1; ++i)
			{
				a1[i] = -a_u(y_1[i], y_1[i - 1]) / h;
				b1[i] = rcht + a_u(y_1[i + 1], y_1[i]) / h + a_u(y_1[i], y_1[i - 1]) / h;
				c1[i] = -a_u(y_1[i + 1], y_1[i]) / h;
				d1[i] = rcht * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			norm = infnorm(y, y_1, n + 1);
			swap(y_1, y);
			count++;
		} while (norm > eps);
		myfile2 << time << " " << count << endl;
		for (int i = 0; i <= n; ++i)
		{
			myfile << time << " " << coord << " " << y[i] << endl;
			u0[i] = u(coord, time);
			coord += h;
		}
	}

	cout << "результаты записаны! " << endl;
}

void KvaziImplGR2GR1(function y0_x, function y0, double yn, function2 a_u) //неявный двухслойный метод для решения уравнения при K=K(x) с граничными условиями 1-го рода 
{
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>a1(n + 1); //Векторы для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);
	vector<double>y_1dop(n + 1);
	vector<double>u0(n + 1);
	double norm = 0;
	int count = 0;

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	ofstream myfile2;
	myfile2.open("Time.dat");
	for (int i = 0; i <= n; ++i) { //задаем н.у.: начальный момент времени
		y_1[i] = y0_x(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	double rcht = rho * c * h / tau;

	a1[n] = 0; //учет граничных условий 
	b1[n] = 1;
	c1[n] = 0;

	a1[0] = 0;
	b1[0] = rcht + a_u(y_1[1], y_1[0]) / h;
	c1[0] = -a_u(y_1[0], y_1[1]) / h;

	for (int j = 1; j <= k; ++j) {
		coord = 0;
		time += tau;
		d1[0] = y0(time);
		d1[n] = yn;
		do
		{
			for (int i = 1; i <= n - 1; ++i)
			{
				a1[i] = -a_u(y_1[i], y_1[i - 1]) / h;
				b1[i] = rcht + a_u(y_1[i + 1], y_1[i]) / h + a_u(y_1[i], y_1[i - 1]) / h;
				c1[i] = -a_u(y_1[i + 1], y_1[i]) / h;
				d1[i] = rcht * y_1[i];
			}
			SweepMethod(a1, b1, c1, d1, n + 1, y);
			norm = infnorm(y, y_1, n + 1);
			swap(y_1, y);
			count++;
		} while (norm > eps);
		myfile2 << time << " " << count << endl;
		for (int i = 0; i <= n; ++i)
		{
			myfile << time << " " << coord << " " << y[i] << endl;
			u0[i] = u(coord, time);
			coord += h;
		}
	}

	cout << "результаты записаны! " << endl;
}

void LinearGR1GR2(double y0, function yn, function yx, double sigma) //двухслойный метод для решения уравнения при K=K(x) с граничными условиями 2-го рода 
{                                                                   // y0 - гран. усл. 2-го рода на левой стенке, yn - гран. усл. 2-го рода на правой стенке, yx - начальное усл.
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>amas(n + 1);
	vector<double>a1(n + 1); //Вектора для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	for (int i = 0; i <= n; ++i) {  // начальное условие
		y_1[i] = yx(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	for (int i = 1; i <= n; ++i) {
		amas[i] = a(i * h);
	}
	OutputMas(amas, n + 1);

	cout << "sigma = " << sigma << endl;
	double rcht = rho * c * h / tau;
	double sh = sigma / h;
	double sh1 = (1 - sigma) / h;

	a1[0] = 0; //учет граничных условий 
	b1[0] = 1;
	c1[0] = 0;

	a1[n] = -amas[n] * sh;
	b1[n] = rcht / 2 + amas[n] * sh;
	c1[n] = 0;

	for (int i = 1; i <= n - 1; ++i)
	{
		a1[i] = -amas[i] * sh;
		b1[i] = rcht + amas[i + 1] * sh + amas[i] * sh;
		c1[i] = -amas[i + 1] * sh;
	}

	for (int j = 1; j <= k; ++j) {
		coord = 0;
		time += tau;
		d1[0] = y0;
		d1[n] = sigma * yn(time) + (1 - sigma) * yn(time - tau) + (rcht / 2 - amas[n] * sh1) * y_1[n] + amas[n] * y_1[n - 1] * sh1;
		for (int i = 1; i <= n - 1; ++i)
		{
			d1[i] = sh1 * amas[i] * y_1[i] + sh1 * amas[i + 1] * y_1[i + 1] + (rcht - sh1 * (amas[i + 1] + amas[i])) * y_1[i];
		}
		SweepMethod(a1, b1, c1, d1, n + 1, y);
		for (int i = 0; i <= n; ++i)
		{
			myfile << time << " " << coord << " " << y[i] << endl;
			coord += h;
		}
		swap(y_1, y);
	}

	cout << "результаты записаны! " << endl;
}

void LinearGR2GR1(function y0, double yn, function yx, double sigma) //двухслойный метод для решения уравнения при K=K(x) с граничными условиями 2-го рода 
{                                                                   // y0 - гран. усл. 2-го рода на левой стенке, yn - гран. усл. 2-го рода на правой стенке, yx - начальное усл.
	vector<double>y_1(n + 1);
	vector<double>y(n + 1);
	vector<double>amas(n + 1);
	vector<double>a1(n + 1); //Вектора для прогонки
	vector<double>b1(n + 1);
	vector<double>c1(n + 1);
	vector<double>d1(n + 1);

	double time = 0; double coord = 0;
	ofstream myfile;
	myfile.open("Temperature.dat");
	for (int i = 0; i <= n; ++i) {  // начальное условие
		y_1[i] = yx(i * h);
		myfile << time << " " << coord << " " << y_1[i] << endl;
		coord += h;
	}

	for (int i = 1; i <= n; ++i) {
		amas[i] = a(i * h);
	}
	//OutputMas(amas, n + 1);

	cout << "sigma = " << sigma << endl;
	double rcht = rho * c * h / tau;
	double sh = sigma / h;
	double sh1 = (1 - sigma) / h;

	a1[0] = 0; //учет граничных условий 
	b1[0] = rcht / 2 + amas[1] * sh;
	c1[0] = -amas[1] * sh;

	a1[n] = 0;
	b1[n] = 1;
	c1[n] = 0;

	for (int i = 1; i <= n - 1; ++i)
	{
		a1[i] = -amas[i] * sh;
		b1[i] = rcht + amas[i + 1] * sh + amas[i] * sh;
		c1[i] = -amas[i + 1] * sh;
	}

	for (int j = 1; j <= k; ++j) {
		coord = 0;
		time += tau;
		d1[0] = y0(time) + rcht / 2 * y_1[0];
		//d1[0] = -sigma * y0(time) - (1 - sigma) * y0(time - tau) + (rcht / 2 - amas[1] * sh1) * y_1[0] + amas[1] * y_1[1] * sh1;
		d1[n] = yn;
		for (int i = 1; i <= n - 1; ++i)
		{
			d1[i] = sh1 * amas[i] * y_1[i] + sh1 * amas[i + 1] * y_1[i + 1] + (rcht - sh1 * (amas[i + 1] + amas[i])) * y_1[i];
		}
		SweepMethod(a1, b1, c1, d1, n + 1, y);
		for (int i = 0; i <= n; ++i)
		{
			myfile << time << " " << coord << " " << y[i] << endl;
			coord += h;
		}
		swap(y_1, y);
	}

	cout << "результаты записаны! " << endl;
}