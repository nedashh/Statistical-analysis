#ifndef STAT_H
#define STAT_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

// Data structures
struct ne_simp {
    int n;
    vector<double> p;
    vector<double> x;
    vector<int> r;
    vector<int> nsample;
};

// Declarations of distribution functions
double norm_cdf(double x);
double norm_ppf(double p);
double norm_pdf(double x);
double t_cdf(double x, double f);
double t_ppf(double p, double f);
double chi_cdf(double x, double f);
double chi_ppf(double p, double f);
double f_cdf(double x, double f1, double f2);
double f_ppf(double p, double f1, double f2);
double nct_ppf(double p, double f, double delta);

// Matrix operations
double** InverseMatrix(double** a, int n);
double** MultiplyMatrix(int rowsa, int colsa, int rowsb, int colsb, double** a, double** b);
double** TransMatrix(int m, int n, double** a);

void MLE_Normal(string ff);
void MLE_Weibull(string ff);
void MLS_Normal(string ff);
void MLS_Weibull(string ff);
void GrubbsTest(string ff);
void StudentTest(string ff);
void BartlettTest(string ff);
void ANOVA(string ff);

void cum(int n, double* x, int* r, int& k, double* fcum, double* ycum);
void standart(int n, double* x, double& cp, double& cko);
double NormalMinFunction(vector<double> xsimpl);
double WeibullMinFunction(vector<double> xsimpl);
void CovMatrixMleN(int n, vector<double> x, vector<int> r, double a, double s, double**& v);
void CovMatrixMleW(int n, vector<double> x, vector<int> r, double c, double b, double**& v);
void ordern(int n, double pr, double ps, double& er, double& vrs);
void orderw(int n, double pr, double ps, double& er, double& vrs);
void lmtaprn(double beta, int n, double v11, double v22, double v12, double delta, double& tlow, double& tup);
void nctWeibull(double beta, int n, double v11, double v22, double v12, double delta, double& tlow, double& tup);
void MLS_Base(string ff, int dist_type);
extern ne_simp nesm;

void MleastSquare(int n, int k, double** x, double** y, double**& db, double**& b, double*& yr);
void MleastSquare_weight(int n, int k, double** x, double** y, double** v, double**& db, double**& b, double*& yr);
int neldermead(vector<double>& x0, double eps, double(*func)(vector<double>));

void GenerateChartData(const string& filename, const vector<double>& x_data, const vector<double>& y_data,
    const string& title, const string& xlabel, const string& ylabel);
void GenerateChartData(const string& filename, const vector<double>& x_data, const vector<double>& y_data,
    const string& title, const string& xlabel, const string& ylabel);
void GenerateDistributionPlot(const string& ff, const vector<double>& Xp,
    const vector<double>& Xplow, const vector<double>& Xpup,
    const vector<double>& P, const string& dist_name);
void GenerateQQPlot(const string& filename, const vector<double>& data, const string& title);

void WilcoxonSignedRankTest(string ff);
void WilcoxonRankSumTest(string ff);
bool wilcoxonDistribution(int m, int n, std::vector<double>& cdf, int& fault);
void ShapiroWilkTest(string ff);
void KruskalWallisTest(string ff);

vector<double> CalculateRanks(const vector<double>& data);
#endif