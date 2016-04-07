#ifndef Tri_project_H
#define Tri_project_H

#include <iostream>
#include <stdlib.h> // for C qsort 
#include <cmath>
#include <time.h> // for random
#include <fstream>
#include <string>

const int MaxVertices = 500;
const int MaxTriangles = 1000;
const int n_MaxPoints = 10; // for the test programm
const double EPSILON = 0.000001;

struct Triangle{
  int p1, p2, p3;
  double cRadius, xc, yc;
  
};

struct EDGE{
  int p1, p2;
};

struct node{
  double x, y, z;
};

int nodeCompare(const void *v1, const void *v2);
int Triangulate(int nv, node pnode[], Triangle v[], int &ntri);
int CircumCircle(double, double, double, double, double, double, double, 
double, double&, double&, double&);

#endif


