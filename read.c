#include<stdio.h>
#include<stdlib.h>

struct triangle{
int tnum;
int n1,n2,n3;
};

int main(){
 int nnodes,i;
 double x, y, dummy, xmin=0, xmax=0, ymin=0, ymax=0, xavg=0, yavg=0;
 double **nodes;
 char* fname;
 FILE* fp;
 fp = fopen("4.node","r");
 fscanf(fp, "%d",&nnodes);
 nodes = (double *)malloc(sizeof(double)*nnodes);
 for(i = 0; i<nnodes; i++){
  nodes[i] = (double *)malloc(sizeof(double)*2);
  fscanf(fp,"%d %f %f",&dummy, &x, &y);
  nodes[i][0] = x;
  nodes[i][1] = y;
 if(x > xmax) 
  xmax = x;
 else if(x < xmin) 
  xmin = x;

 if(y > ymax) 
  ymax = y;
 else if(y < ymin) 
  ymin = y;
 xavg += x;
 yavg += y;
 }
 xavg /= nnodes;
 yavg /= nnodes;
 printf("left is %f , %f\n",(xmin-xavg),(ymin-yavg));
 printf("right is %f, %f\n",(xmax+xavg),(ymin- yavg));
 printf("top is %f , %f\n",(xmin+xmax)/2,(ymax + yavg));
 return 0;
}
