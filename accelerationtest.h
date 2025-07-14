#ifndef ACCELERATIONTEST_H
#define ACCELERATIONTEST_H

#include <iostream>
#include <fstream>
#include<string.h>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include<cmath>
#include<vector>
#include<time.h>
#include<chrono>
#include<windows.h>
#include<algorithm>

using namespace std;

double calculateValue(int row, int col);
double createfilter(vector<double>& matrix,double *pplppower);
double calculateloss(double *lppower,double *pplppower,ofstream& file);

#endif //ACCELERATIONTEST_H