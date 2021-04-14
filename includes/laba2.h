#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

typedef double(*function)(double x);
typedef double(*function2)(double x, double y);

class Solver
{
private:
    double const eps = 1e-3;
    double const L = 1; // ширина стержня
    double const T = 1; // время
    int const n = 50; // количество разбиений по пространству
    int const k = 50; // количество разбиений по времени
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

    // Test1
    double const y0 = 0.5, yn = 0.5, sigma = 0.8;

public:
    // methods
    void LinearGR1();

    // additional methods
    double a(double x);
    double K(double x);
    double yx(double x);

    // exactSolutions
    double u1(double x, double t);
};

// additional functions
double  infnorm(double v1, double v2);
void    SweepMethod(double* a, double* b, double* c, double* d,
        int k, double* res);
double  infnorm(double* vec1, double* vec2, int s);