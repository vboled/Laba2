#include "../includes/laba2.h"

double Solver::u1(double x, double t) //Точное решение теста 1
{
	double k = 2.5;
	return pow(e, -k * pi * pi * t / L / L) * 8 * L * L * sin(pi * x / L) / (pi * pi * pi);
};