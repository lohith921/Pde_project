#include<stdio.h>
#include<stdlib.h>

struct node{
	int nnum;
	double x,y;
	};
struct triangle{
	int tnum;
	struct node nodes[3];
	};
