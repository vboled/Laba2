#pragma once
#include <iostream>
#include <fstream>
#include<string>
#include<utility>
#include <math.h>
#include <iomanip>
#include<vector>
#include<cmath>

using namespace std;
using std::swap;

typedef double(*function)(double x);
typedef double(*function2)(double x, double y);

void OutputMas(vector<double>v, int m);
void multmatrvec(double** mas, vector<double>vec, int s, vector<double>res);
double euclnorm(vector<double>v, int m);
double infnorm(double v1, double v2);
double infnorm(vector<double>vec1, vector<double>vec2, int s);
void ClearFile(); 
void freematrix(double** Matrix, int Rows);
void SweepMethod(vector<double>a, vector<double>b, vector<double>c, vector<double>d, int k, vector<double>res);

void LinearGR1(double y0, double yn, function yx, double sigma);
void LinearGR2(function y0, function yn, function yx, double sigma);
void LinearGR1GR2(double y0, function yn, function yx, double sigma);
void LinearGR2GR1(function y0, double yn, function yx, double sigma);
void KvaziImpl(function y0_x, function y0_t, function yL_t, function2 a_u, int M);
void KvaziImpl(function y0_x, function y0_t, function yL_t, function2 a_u); 
void KvaziImplGR1GR2(function y0_x, double y0, function yn, function2 a_u, int M);
void KvaziImplGR1GR2(function y0_x, double y0, function yn, function2 a_u);
void KvaziImplGR2GR1(function y0_x, function y0, double yn, function2 a_u, int M);
void KvaziImplGR2GR1(function y0_x, function y0, double yn, function2 a_u);