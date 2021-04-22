#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

// t0 = 0.5, Q = 10, L = 1;

// Var 3
// p = 1, c = 1, a = 2, b = 0.5, c = 3, u0 = 0.1, k1 = 1, k2 = 0.1, x1 = 1/3, x2 = 2/3;
// Слева P3, Справа u = u0;

// Var 15
// p = 0.5, c = 2, alpha = 5, beta = 0.1, gamma = 4, u0 = 0.04, k1 = 0.1, k2 = 1.5, x1 = 1/4, x2 = 1/2;
// Слева u = u0, Справа P3; 


using namespace std;

const double EPSILON = 1e-8;   //константа сравнения с 0
const double PI = 3.14159265;

typedef double(*func)(const double x, const double u0, const double L);
typedef double(*Fanc)(const double x,const double t);
typedef double(*Func)(const double x);

int MaxT(const int n, const int k, double** u);

void Copy_to_file(const unsigned int n, const unsigned int k, const double h, double** solve, string file_name); //n - кол-во узлов по пространству, h-шаг по пространству 
void copy(int dim, double* a, double* b);
void vivod(const unsigned int DIM, const double* const b);
double norm(const unsigned int DIM, const double* const b, const char flag); //3 нормы вектора с запросом варианта нормы
double* minus_vect(const unsigned int DIM, const double* const b1, const double* const b2); //b1-b2

double** Integro_interpolation(int n, int k, double h, double tao, Fanc ux0, Func u_0t, double c, double p, Func Kx, Func Pt, double sigma);//слева равномерно, справа тепловой поток
double** Integro_interpolation(int n, int k, double h, double tao, Fanc ux0, Func u_0t, Fanc u_Lt, double c, double p, Func Kx, double sigma);//равномерно прогрето с обоих концов
double** Integro_interpolation(int n, int k, double h, double tao, Fanc ux0, double c, double p, Func Kx, Func Pt1, Func Pt2, double sigma);
double** Kvazilin_explicit(int n, int k, double h, double tao, func ux0, func Ku); //n - кол-во узлов по пространству, k - по времени. 
double** Kvazilin_implicit(const unsigned int n, const unsigned int k, double h, double tao, func ux0, func Ku, const unsigned int M); //n - кол-во узлов по пространству, k - по времени, M-число итерация в неявном цикле
double** Kvazilin_implicit_Eps(const unsigned int n, const unsigned int k, double h, double tao, func ux0, func Ku, const double Eps); //n - кол-во узлов по пространству, k - по времени. 
double** Kvazilin_explicit2(int n, int k, double h, double tao, func ux0, func Ku);//n - кол-во узлов по пространству, k - по времени. 
double* progon3d(int DIM, double* a, double* b, double* c, double* d);

double error(int k, int n, double h, double tao, double** a);


double ux0(const double x, const double L)
{   //тест 2
	double u0 = 0.2;
	// if (L-x<EPSILON) { return 0.; }
	// else if (x < EPSILON) { return 0.; }
	// else {return  u0 + x*(L - x); };

	// return u0;//равномерно прогрет
	return  u0 + x*(L - x);//тест 1

	//return sin(x);//свой тест
}

double ux0_0(const double x, const double u0, const double L)
{
	return 0.;
}

double Kx(const double x)
{
	//методичка 1,2 тест
	double k1 = 2;
	double k2 = 0.5;
	double x1 = 0.5;
	double x2 = 2 / 3;

	if (x  <= x1) { return k1; }
	else if (x < x2) { return k1*(x - x2) / (x1 - x2) + k2*(x - x1) / (x2 - x1); }
	else return k2;

}

double Ku(const double u, const double sigma, const double kappa)
{
	return 0.5 + 2 * pow(u, 2);
}

double Pt(const double t)//(тепловой поток)
{  // 5 вар 
	double t0 = 0.5;
	double Q = 10.;

	if (t < t0 && t > EPSILON) {
        return Q;
    }
	else
        return 0.;

	// return 0.;

}

double u_0t(const double t)
{	
	double u0 = 0.2;
	return u0;
	
	// return 0.;//свой тест
}

double u_Lt(const double L,const double t)
{
	double u0 = 0.2;
	return u0;
	
	// return sin(L)*exp(-t);//свой тест
}

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык
	double L = 1.;//длинна стержня
	double T = 0.5;//Время

	//5 вариант
	double c = 2.;
	double p = 0.25;

	double h = 0.01;
	double tao = 0.01;
	int n = L / h + 1;
	int k = T / tao + 1;


	double** u = new double*[k];
	for (int i = 0; i < k; i++) { u[i] = new double[n]; };


	// u = Integro_interpolation(n, k, h, tao, ux0, u_0t, c, p, Kx, Pt, 0.2); // слева равномерно, справа тепловой поток
    string file_name1 = "SOLVE1.txt";
	// Copy_to_file(n, k, h, u, file_name1);
	// cout << MaxT(n, k, u) << endl;

	// u = Integro_interpolation(n, k, h, tao, ux0, u_0t, u_Lt, c, p, Kx, 0.6); // слева и справа равномерно
	// file_name1 = "SOLVE2.txt";
	// Copy_to_file(n, k, h, u, file_name1);
	// cout << MaxT(n, k, u) << endl;

    // u = Integro_interpolation(n, k, h, tao, ux0, c, p, Kx, Pt, Pt, 0.5); // слева и справа тепловой поток
	// file_name1 = "SOLVE3.txt";
	// Copy_to_file(n, k, h, u, file_name1);
	// cout << MaxT(n, k, u) << endl;

    u = Kvazilin_explicit(n, k, h, tao, ux0_0, Ku);
	file_name1 = "SOLVE4.txt";
	Copy_to_file(n, k, h, u, file_name1);
	cout << MaxT(n, k, u) << endl;
	
    // u = Kvazilin_explicit2(n, k, h, tao, ux0_0, Ku);
	// file_name1 = "SOLVE5.txt";
	// Copy_to_file(n, k, h, u, file_name1);
	// cout << MaxT(n, k, u) << endl;

    // u = Kvazilin_implicit(n, k, h, tao, ux0_0, Ku, 3);
	// file_name1 = "SOLVE6.txt";
	// Copy_to_file(n, k, h, u, file_name1);
	// cout << MaxT(n, k, u) << endl;

    // u = Kvazilin_implicit_Eps(n, k, h, tao, ux0_0, Ku, 1e-4);
	// file_name1 = "SOLVE7.txt";
	// Copy_to_file(n, k, h, u, file_name1);
    // cout << MaxT(n, k, u) << endl;

	//vivod(n,u[1]);
	//for (int i = 0; i < n; i++) { cout<< sin(i*h)*exp(-tao)<<" "; }

	// cout <<"погрешность = " <<error(k, n, h, tao, u);

	// char file_name2[] = "SOLVE2.txt";
	// Copy_to_file(n, h, u[10000], file_name2);

	// char file_name3[] = "SOLVE3.txt";
	// Copy_to_file(n, h, u[15000], file_name3);

	//char file_name4[] = "SOLVE4.txt";
	//Copy_to_file(n, h, u[8], file_name4);

	//char file_name5[] = "SOLVE5.txt";
	//Copy_to_file(n, h, u[10], file_name5);

	return 0;
}

int MaxT(const int n, const int k, double** u)
{
	double max = 0;
	int imax = 0;
	double tmp = 0;
	for (int i = 0; i < k; i++)
	{
		tmp = u[i][n / 2];
		if ((tmp - max) > 0.01)
		{
			max = tmp;
			imax = i;
		}
	}
	cout << "max T = " << max << endl;
	return imax;
}
// u[i][n / 2] --> MAX --> i

double** Kvazilin_explicit(int n, int k, double h, double tao, func ux0, func Ku) //n - кол-во узлов по пространству, k - по времени. 
{
	const int N = n - 1;//колличество разбиений

	double** Y = new double*[k];
	Y[0] = new double[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, 0, n*h); }

	double sigma = 2;
	double kappa = 0.5;
	double c = 1;
	double c1 = 5;
	double u0 = sigma*pow(c1, 2) / kappa;

	double u0t;
	double ult;


	double* a = new double[2];

	double* diag1 = new double[N - 1];
	double* diag2 = new double[N - 1];
	double* diag3 = new double[N - 1];
	double* right = new double[N - 1];
	double* vspom = new double[N - 1];
	double A0;
	double BN;


	for (int j = 1; j < k; j++) {

		u0t = pow(u0*j*tao, 1 / sigma);
		if (n*h < c1*j*tao) { ult = pow(u0 / c1*(c1*j*tao - n*h), 1 / sigma); }
		else { ult = 0; };

		Y[j] = new double[n];
		a[1] = 0.5*(Ku(Y[j - 1][1], sigma, kappa) + Ku(Y[j - 1][0], sigma, kappa));
		for (int i = 0; i < N - 1; i++) {
			a[0] = a[1];
			a[1] = 0.5*(Ku(Y[j - 1][i + 1], sigma, kappa) + Ku(Y[j - 1][i], sigma, kappa));
			diag1[i] = -a[0] / pow(h, 2);
			diag3[i] = -a[1] / pow(h, 2);
			diag2[i] = -(diag1[i] + diag3[i] - c / tao);
			right[i] = c / tao*Y[j - 1][i + 1];
		};
		A0 = diag1[0]; BN = diag3[N - 2];
		diag1[0] = 0.; diag3[N - 2] = 0.;

		right[0] -= A0*u0t; right[N - 2] -= BN*ult;

		copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));

		for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

		Y[j][0] = u0t;
		Y[j][N] = ult;
	}

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}

double** Kvazilin_implicit(const unsigned int n, const unsigned int k, double h, double tao, func ux0, func Ku, const unsigned int M) //n - кол-во узлов по пространству, k - по времени. 
{
	const unsigned int N = n - 1;//колличество разбиений

	double** Y = new double*[k];
	Y[0] = new double[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, 0, n*h); }

	double sigma = 2;
	double kappa = 0.5;
	double c = 1;
	double c1 = 5;
	double u0 = sigma*pow(c1, 2) / kappa;

	double u0t;
	double ult;


	double* a = new double[2];

	double* diag1 = new double[N - 1];
	double* diag2 = new double[N - 1];
	double* diag3 = new double[N - 1];
	double* right = new double[N - 1];
	double* vspom1 = new double[N - 1];
	double** vspom2 = new double*[2];//для иттераций неявного решения
	vspom2[0] = new double[n];
	vspom2[1] = new double[n];
	double A0;
	double BN;


	for (int j = 1; j < k; j++) {

		u0t = pow(u0*j*tao, 1 / sigma);
		if (n*h < c1*j*tao) { ult = pow(u0 / c1*(c1*j*tao - n*h), 1 / sigma); }
		else { ult = 0; };

		Y[j] = new double[n];
		copy(n, vspom2[1], Y[j - 1]);
		vspom2[1][0] = u0t;
		vspom2[1][N] = ult;
		for (int m = 1; m <= M; m++) {
			copy(n, vspom2[0], vspom2[1]);
			a[1] = 0.5*(Ku(vspom2[0][1], sigma, kappa) + Ku(vspom2[0][0], sigma, kappa));
			for (int i = 0; i < N - 1; i++) {
				a[0] = a[1];
				a[1] = 0.5*(Ku(vspom2[0][i + 1], sigma, kappa) + Ku(vspom2[0][i], sigma, kappa));
				diag1[i] = -a[0] / pow(h, 2);
				diag3[i] = -a[1] / pow(h, 2);
				diag2[i] = -(diag1[i] + diag3[i] - c / tao);
				right[i] = c / tao*vspom2[0][i + 1];
			};
			A0 = diag1[0]; BN = diag3[N - 2];
			diag1[0] = 0.; diag3[N - 2] = 0.;

			right[0] -= A0*u0t; right[N - 2] -= BN*ult;

			vspom1 = progon3d(N - 1, diag1, diag2, diag3, right);

			for (int i = 1; i < N; i++) vspom2[1][i] = vspom1[i - 1];

			delete[]vspom1;
		}
		copy(n, Y[j], vspom2[1]);
	}

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom2[0];
	delete[] vspom2[1];
	delete[] vspom2;

	return Y;
}

double** Kvazilin_implicit_Eps(const unsigned int n, const unsigned int k, double h, double tao, func ux0, func Ku, const double Eps) //n - кол-во узлов по пространству, k - по времени. 
{
	const unsigned int N = n - 1;//колличество разбиений

	double** Y = new double*[k];
	Y[0] = new double[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, 0, n*h); }

	double sigma = 2;
	double kappa = 0.5;
	double c = 1;
	double c1 = 5;
	double u0 = sigma*pow(c1, 2) / kappa;

	double u0t;
	double ult;


	double* a = new double[2];

	double* diag1 = new double[N - 1];
	double* diag2 = new double[N - 1];
	double* diag3 = new double[N - 1];
	double* right = new double[N - 1];
	double* vspom1 = new double[N - 1];
	double** vspom2 = new double*[2];//для иттераций неявного решения
	vspom2[0] = new double[n];
	vspom2[1] = new double[n];
	double A0;
	double BN;

	int iter;

	for (int j = 1; j < k; j++) {

		u0t = pow(u0*j*tao, 1 / sigma);
		if (n*h < c1*j*tao) { ult = pow(u0 / c1*(c1*j*tao - n*h), 1 / sigma); }
		else { ult = 0; };

		Y[j] = new double[n];
		copy(n, vspom2[1], Y[j - 1]);
		vspom2[1][0] = u0t;
		vspom2[1][N] = ult;
		iter = 0;
		do {
			iter++;
			copy(n, vspom2[0], vspom2[1]);
			a[1] = 0.5*(Ku(vspom2[0][1], sigma, kappa) + Ku(vspom2[0][0], sigma, kappa));
			for (int i = 0; i < N - 1; i++) {
				a[0] = a[1];
				a[1] = 0.5*(Ku(vspom2[0][i + 1], sigma, kappa) + Ku(vspom2[0][i], sigma, kappa));
				diag1[i] = -a[0] / pow(h, 2);
				diag3[i] = -a[1] / pow(h, 2);
				diag2[i] = -(diag1[i] + diag3[i] - c / tao);
				right[i] = c / tao*vspom2[0][i + 1];
			};
			A0 = diag1[0]; BN = diag3[N - 2];
			diag1[0] = 0.; diag3[N - 2] = 0.;

			right[0] -= A0*u0t; right[N - 2] -= BN*ult;

			vspom1 = progon3d(N - 1, diag1, diag2, diag3, right);

			for (int i = 1; i < N; i++) vspom2[1][i] = vspom1[i - 1];

			delete[]vspom1;
		} while (norm(n, minus_vect(n, vspom2[1], vspom2[0]), 'k')>Eps);
		// cout << "итераций на временном слое" << j << " = " << iter << endl;
		copy(n, Y[j], vspom2[1]);
	}

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom2[0];
	delete[] vspom2[1];
	delete[] vspom2;

	return Y;
}

double** Kvazilin_explicit2(int n, int k, double h, double tao, func ux0, func Ku) //n - кол-во узлов по пространству, k - по времени. 
{
	const int N = n ;//колличество разбиений

	double** Y = new double*[k];
	Y[0] = new double[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, 0, n*h); }

	double sigma = 2;
	double kappa = 0.5;
	double c = 1;
	double c1 = 5;
	double u0 = sigma*pow(c1, 2) / kappa;

	double u0t;
	double ult;
	double mu;


	double* a = new double[2];

	double* diag1 = new double[N-1];
	double* diag2 = new double[N-1];
	double* diag3 = new double[N-1];
	double* right = new double[N - 1];
	double* vspom = new double[N - 1];
	double A0;
	double BN;


	for (int j = 1; j < k; j++) {

		u0t = pow(u0*j*tao, 1 / sigma);
		if (n*h < c1*j*tao) {
			Y[j] = new double[n];
			a[1] = 0.5*(Ku(Y[j - 1][1], sigma, kappa) + Ku(Y[j - 1][0], sigma, kappa));
			for (int i = 0; i < N - 1; i++) {
				a[0] = a[1];
				a[1] = 0.5*(Ku(Y[j - 1][i + 1], sigma, kappa) + Ku(Y[j - 1][i], sigma, kappa));
				diag1[i] = -a[0] / pow(h, 2);
				diag3[i] = -a[1] / pow(h, 2);
				diag2[i] = -(diag1[i] + diag3[i] - c / tao);
				right[i] = c / tao*Y[j - 1][i + 1];
			};
			A0 = diag1[0]; BN = diag3[N - 2];
			diag1[0] = 0.; diag3[N - 2] = 0.;
			diag2[N - 2] += 2*BN;
			diag1[N - 2] -= BN;

			right[0] -= A0*u0t; right[N - 2] -= 0.;

			copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));

			for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

			Y[j][0] = u0t;
			//Y[j][N] = Y[j][N - 1] + Y[j - 1][N] - Y[j - 1][N - 1];
		}
		else {
			ult = 0;
			Y[j] = new double[n];
			a[1] = 0.5*(Ku(Y[j - 1][1], sigma, kappa) + Ku(Y[j - 1][0], sigma, kappa));
			for (int i = 0; i < N - 1; i++) {
				a[0] = a[1];
				a[1] = 0.5*(Ku(Y[j - 1][i + 1], sigma, kappa) + Ku(Y[j - 1][i], sigma, kappa));
				diag1[i] = -a[0] / pow(h, 2);
				diag3[i] = -a[1] / pow(h, 2);
				diag2[i] = -(diag1[i] + diag3[i] - c / tao);
				right[i] = c / tao*Y[j - 1][i + 1];
			};
			A0 = diag1[0]; BN = diag3[N - 2];
			diag1[0] = 0.; diag3[N - 2] = 0.;

			right[0] -= A0*u0t; right[N - 2] -= BN*ult;

			copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));

			for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

			Y[j][0] = u0t;
			Y[j][N] = ult;
		};
	}



	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}

double** Integro_interpolation(int n, int k, double h, double tao, Fanc ux0, Func u_0t, double c, double p, Func Kx, Func Pt, double sigma)
{
	const int N = n - 1;//колличество разбиений

	double** Y = new double*[k];
	Y[0] = new double[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h,(n - 1)*h); }

	// double sigma = 0.7;

	double* a = new double[N];



	double* diag1 = new double[N - 1];
	double* diag2 = new double[N - 1];
	double* diag3 = new double[N - 1];
	double* right = new double[N - 1];
	double* vspom = new double[N - 1];
	double A0, BN; double mu, kappa,vspom1;
	double L = N*h; double psi = 0.0;


	for (int i = 0; i < N; i++) { a[i] = Kx((i + 0.5)*h); };

	kappa = (sigma*a[N - 1] / h) / (c*p*h / (2 * tao) + sigma*a[N - 1] / h);

	for (int i = 0; i < N - 1; i++) {
		diag1[i] = sigma / h*a[i];
		diag3[i] = sigma / h*a[i + 1];
		diag2[i] = -(diag1[i] + diag3[i] + c*p*h / tao);
	}

	A0 = diag1[0]; BN = diag3[N - 2];

	diag1[0] = 0.; diag3[N - 2] = 0.;

	diag2[N - 2] += BN*kappa;

	for (int j = 1; j < k; j++) {

		Y[j] = new double[n];
		Y[j][0] = u_0t(j*tao);

		for (int i = 0; i < N - 1; i++) right[i] = -(c*p*h / tao*Y[j - 1][i + 1] + (1 - sigma)*a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);

		mu = (c*p*Y[j - 1][N] * h / (2 * tao) + sigma*Pt(tao*j) + (1 - sigma)*(Pt(tao*(j - 1)) - (Y[j - 1][N] - Y[j - 1][N - 1]) / h)) / (c*p*h / (2 * tao) + sigma*a[N - 1] / h);

		right[0] -= A0*u_0t(j*tao); right[N - 2] -= BN*mu;

		copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));
		

		for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

		Y[j][N] = kappa*Y[j][N - 1] + mu;
	
	}
	
	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}

double** Integro_interpolation(int n, int k, double h, double tao, Fanc ux0, double c, double p, Func Kx, Func Pt1,Func Pt2, double sigma)
{
	const int N = n - 1;//колличество разбиений

	double** Y = new double*[k];
	Y[0] = new double[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, (n - 1)*h); }

	// double sigma = 1.0;

	double* a = new double[N];



	double* diag1 = new double[N - 1];
	double* diag2 = new double[N - 1];
	double* diag3 = new double[N - 1];
	double* right = new double[N - 1];
	double* vspom = new double[N - 1];
	double A0, BN; double mu1, kappa1,mu2,kappa2, vspom1;
	double L = N*h; double I = 0.0;


	for (int i = 0; i < N; i++) { a[i] = Kx((i + 0.5)*h); };

	kappa1 = (sigma*a[0] / h) / (c*p*h / (2 * tao) + sigma*a[0] / h);
	kappa2 = (sigma*a[N - 1] / h) / (c*p*h / (2 * tao) + sigma*a[N - 1] / h);

	for (int i = 0; i < N - 1; i++) {
		diag1[i] = sigma / h*a[i];
		diag3[i] = sigma / h*a[i + 1];
		diag2[i] = -(diag1[i] + diag3[i] + c*p*h / tao);
	}

	A0 = diag1[0]; BN = diag3[N - 2];

	diag1[0] = 0.; diag3[N - 2] = 0.;

	diag2[0] += A0*kappa1;
	diag2[N - 2] += BN*kappa2;

	for (int j = 1; j < k; j++) {

		Y[j] = new double[n];

		for (int i = 0; i < N - 1; i++) right[i] = -(c*p*h / tao*Y[j - 1][i + 1] + (1 - sigma)*a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);

		mu1 = (c*p*Y[j - 1][0] * h / (2 * tao) + sigma*Pt1(tao*j) + (1 - sigma)*(Pt1(tao*(j - 1)) + (Y[j - 1][1] - Y[j - 1][0]) / h)) / (c*p*h / (2 * tao) + sigma*a[0] / h);
		mu2 = (c*p*Y[j - 1][N] * h / (2 * tao) + sigma*Pt2(tao*j) + (1 - sigma)*(Pt2(tao*(j - 1)) - (Y[j - 1][N] - Y[j - 1][N - 1]) / h)) / (c*p*h / (2 * tao) + sigma*a[N - 1] / h);

		right[0] -= A0*mu1; right[N - 2] -= BN*mu2;

		copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));


		for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

		Y[j][0] = kappa1*Y[j][1] + mu1;
		Y[j][N] = kappa2*Y[j][N - 1] + mu2;

		I = 0.0;
		I += 0.5*Y[j][0] * h;
		for (int i = 1; i < N; i++) {
			I+=Y[j][i]*h;
		}
		I += 0.5*Y[j][N] * h;
		// cout << "I_" << j-1 << " = " << I << "; ";
	}

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}

double** Integro_interpolation(int n, int k, double h, double tao, Fanc ux0, Func u_0t, Fanc u_Lt, double c, double p, Func Kx, double sigma)
{
	const int N = n - 1;//колличество разбиений

	double** Y = new double*[k];
	Y[0] = new double[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, (n - 1)*h); }

	//double sigma = 0.2;
	// double sigma = 1.;

	double* a = new double[N];



	double* diag1 = new double[N - 1];
	double* diag2 = new double[N - 1];
	double* diag3 = new double[N - 1];
	double* right = new double[N - 1];
	double* vspom = new double[N - 1];
	double A0, BN; double mu, kappa, vspom1;
	double L = N*h; double psi = 0.0;

	
	for (int i = 0; i < N; i++) { a[i] = Kx((i + 0.5)*h); };


	for (int i = 0; i < N - 1; i++) {
		diag1[i] = sigma / h*a[i];
		diag3[i] = sigma / h*a[i + 1];
		diag2[i] = -(diag1[i] + diag3[i] + c*p*h / tao);
	}

	A0 = diag1[0]; BN = diag3[N - 2];

	diag1[0] = 0.; diag3[N - 2] = 0.;

	for (int j = 1; j < k; j++) {

		Y[j] = new double[n];
		Y[j][0] = u_0t(j*tao);

		for (int i = 0; i < N - 1; i++) right[i] = -(c*p*h / tao*Y[j - 1][i + 1] + (1 - sigma)*a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);


		right[0] -= A0*u_0t(j*tao);
		right[N - 2] -= BN* u_Lt(L,j*tao);

		copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));

		for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

		Y[j][N] = u_Lt(L, j*tao);


		//for (int i = 1; i < N - 1; i++) {
		//vspom1 = (sin((i + 1)*h) - 2 * sin(i*h) + sin((i - 1)*h))*exp(-j*tao) / (h*h) - sin(i*h)*(exp(-(j + 1)*tao) - exp(-j*tao)) / tao;
		//if (fabs(vspom1) > fabs(psi))  psi = vspom1; 
		//};
	
	}

	//cout << "\nпогрешность аппроксимации = " << psi << endl;

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}


double* progon3d(int DIM, double* a, double* b, double* c, double* d)

{
	double* alfa = new double[DIM];
	double* betta = new double[DIM];

	alfa[0] = -c[0] / b[0]; betta[0] = d[0] / b[0];

	for (int i = 1; i < DIM; i++) {
		alfa[i] = -c[i] / (b[i] + a[i] * alfa[i - 1]);
		betta[i] = (-a[i] * betta[i - 1] + d[i]) / (a[i] * alfa[i - 1] + b[i]);
	}

	for (int i = DIM - 2; i > -1; i--) {
		betta[i] += alfa[i] * betta[i + 1];
	}

	delete[] alfa;

	return  betta;
}

double error(int k, int n,double h, double tao, double** a)
{
	double MAXerror = 0.0;
	//double psi = 0.;
	//double vspom1;
	for (int j = 0; j < k; j++)
		for (int i = 0; i < n; i++)
		if (fabs(a[j][i] - sin(i*h)*exp(-j*tao)) > MAXerror) MAXerror = fabs(a[j][i] - sin(i*h)*exp(-j*tao));

	return MAXerror;
}


void Copy_to_file(const unsigned int n, const unsigned int k, const double h, double** solve, string file_name) //n - кол-во узлов по пространству, h-шаг по пространству 
{
	ofstream fout;    // создали переменную для записи в файл
	fout.open(file_name, ios_base::out | ios_base::trunc);
	for (unsigned int i = 0; i < k; i++) {
        for (unsigned int j = 0; j < n; j++) {
		    fout << solve[i][j] << " ";
        }
        fout << endl;
	};
	fout.close();
	fout.clear();
}


void copy(int dim, double* a, double* b)
{
	for (int i = 0; i < dim; i++) a[i] = b[i];
}



void vivod(const unsigned int DIM, const double* const b) //процедура вывода столбца
{
	cout << '(';
	for (int i = 0; i < DIM; i++)
	{
		cout << b[i];
		if (i < DIM - 1) { cout << ','; } //красивая запись
	}
	cout << ")^т\n";
}

double norm(const unsigned int DIM, const double* const b, const char flag) //3 нормы вектора с запросом варианта нормы
{
	double norm = 0;  //получаемая норма вектора
	if (flag == 'k')  //кубическая норма (max)
	{

		for (int i = 0; i < DIM; i++)
		{
			if (norm < abs(b[i])) { norm = abs(b[i]); }  //находим максимум
		}
		return norm;
	}

	if (flag == '1') //октаэдрическая норма (сумма элементов)
	{
		for (int j = 0; j < DIM; j++)
		{
			norm += abs(b[j]);  //производим сложение элементов
		}
		return norm;
	}

	if (flag == '2')  //шаровая норма  (Евклидова)
	{
		for (int i = 0; i<DIM; i++)
		{
			norm += b[i] * b[i];  //суммируем квадраты элементов
		}

		return sqrt(norm); //возвращаем корень квадратный из суммы квадратов
	}
}

double* minus_vect(const unsigned int DIM, const double* const b1, const double* const b2) //b1-b2
{
	double* x = new double[DIM];
	for (int i = 0; i < DIM; i++) {
		x[i] = b1[i] - b2[i];
	}
	return x;
}