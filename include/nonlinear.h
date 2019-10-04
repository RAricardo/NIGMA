#include <iostream>
#include <vector>
#include <time.h>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <limits>
#include "fparser.hh"
#include "printer.h"

void setF(string);
void setG(string);
void setDervF(string);
double f(double);
double g(double);
double fderv(double);
string incrementalSearch(double, double, int);
string bisection (double, double, double, double, bool);
double regulaFalsiXmAux(double, double);
string regulaFalsi (double, double, double, double, bool);
string fixedPoint(double, double, double, bool);
string newtonMethod(double, double, double, bool);



