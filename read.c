#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mmap_lib.h"

struct node{
	int nnum;
	double x,y;
	};
struct tree_node {
 	struct triangle *lchild;
	struct triangle *mchild;
	struct triangle *rchild;
	struct triangle *link_to_triangle;
	};
struct triangle{
	int tnum;
	struct node *n1, *n2, *n3;
	//struct triangle *l;
	//struct triangle *c;
	//struct triangle *r;
	struct triangle *nb1, *nb2, *nb3;
	struct tree_node *link_to_tree;
	//struct triangle *pt;
	};
struct mesh{
	struct triangle *T;
	int num_of_nodes, num_of_tris;
	};

void print_tree(struct tree_node *root){
		int i;
		if(root!=NULL){
			for(i=0;i<3;i++){
				printf("N%d.x = %f, N%d.y = %f \n",i,root->nodes[i].x,i, root->nodes[i].y);
			}
			printf("\n");
			print_tree(root->l);
			print_tree(root->c);
			print_tree(root->r);
		}
		else{
			return;
		}
	}
	
void insert(struct tree_node *root, struct triangle *tri_list, struct node p){
		if(point_inside(root,p) && (root->l == NULL && root->c == NULL && root->r == NULL)){
			struct triangle *t1,*t2,*t3;
			t1 = (struct triangle*)malloc(sizeof(struct triangle));
			t2 = (struct triangle*)malloc(sizeof(struct triangle));
			t3 = (struct triangle*)malloc(sizeof(struct triangle));
			t1->nodes[0] = root->nodes[0];
			t1->nodes[1] = root->nodes[1];
			t1->nodes[2] = p; t1->l = NULL; t1->c = NULL; t1->r = NULL;
			t1->pt = root; t1->n1 = t2; t1->n2 = t3;
			root->l = t1;
			t2->nodes[0] = root->nodes[1];
			t2->nodes[1] = root->nodes[2];
			t2->nodes[2] = p; t2->l = NULL; t2->c = NULL; t2->r = NULL;
			t2->pt = root; t2->n1 = t3; t2->n2 = t1;
			root->c = t2;
			t3->nodes[0] = root->nodes[2];
			t3->nodes[1] = root->nodes[0];
			t3->nodes[2] = p; t3->l = NULL; t3->c = NULL; t3->r = NULL;
			t3->pt = root; t3->n1 = t1; t3->n2 = t2;
			root->r = t3;
			
		}
		else if( point_inside(root,p)){
			insert(root->l, p);
			insert(root->c, p);
			insert(root->r, p);
		}
		else{
			return;
		}
	}
double dotp(struct node n1, struct node n2){
	return (n1.x*n2.x + n1.y*n2.y);
}	
int point_inside(struct triangle *t, struct node n){
	struct node v0, v1, v2;
	v0 = t->nodes[0]; v1 = t->nodes[1]; v2 = t->nodes[2];
	double l1,l2,l3,denom;
	denom = (v1.y-v2.y)*(v0.x-v2.x)+(v2.x-v1.x)*(v0.y-v2.y);
	l1 = ((v1.y-v2.y)*(n.x-v2.x)+(v2.x-v1.x)*(n.y-v2.y))/denom;
    	l2 = ((v2.y-v0.y)*(n.x-v2.x)+(v0.x-v2.x)*(n.y-v2.y))/denom;
   	 l3 = 1-l1-l2; 
	//printf("l1 = %f, l2 = %f, l3 = %f\n", l1,l2,l3);
	if( (0.0 < l1 && l1 < 1.0) && (l2 > 0.0 && l2 < 1.0) && (l3 > 0.0 && l3 < 1.0))
		return 1; // point lies inside.
	else
		return 0; // point lies outside.
}

int main(){
    int nnodes,i,*node_num, *t_num,dummy;
    double xavg=0, yavg=0;
    struct node *input_nodes; // to keep track of input nodes.
    struct triangle *T=NULL;  // The big triangle
    struct node temp; // temporarily store x,y values for a node.
    struct mesh *m;
    float x,y,xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0;
    char* fname;
    FILE* fp;
    fp = fopen("4.node","r");
    fscanf(fp, "%d %d %f %f\n",&nnodes,&dummy,&x,&y);
    input_nodes = (struct node* )malloc(sizeof(struct node)*nnodes);
    map_t* node_map;
    node_map = map_create();
    *node_num = 3; *t_num = 0;
    for(i = 0; i<nnodes; i++){
    	fscanf(fp,"%d %f %f\n",&dummy, &x, &y);
    	//printf("x = %d, y = %d\n",x,y);
    	//input_nodes[i].x = (double)x;
    	//input_nodes[i].y = (double)y;
    	map_set(&node_map,node_num,x,y);
    	node_num++;
    	if(x >= xmax) 
  		xmax = x;
    	else if(x <= xmin) 
   		xmin = x;
	if(y >= ymax) 
  		ymax = y;
    	else if(y <= ymin) 
 		ymin = y;
 	}
 	display_map(node_map);
 	/*printf("xmax is %d, xmin is %d, ymax is %d, ymin is %d\n",xmax, xmin,ymax, ymin);	
    printf("left is %f , %f\n",(double)xmin-abs(xmax),(double)ymin-abs(ymax));
    printf("right is %f, %f\n",(double)xmax+abs(xmin),(double)ymin-abs(ymax));
    printf("top is %f , %f\n",(double)(xmin+xmax)/2,(double)ymax+abs(ymin));*/
    // Initializing the big triangle.
    T = (struct triangle*)malloc(sizeof(struct triangle));
    T->l = NULL; T->c = NULL; T->r = NULL; (*T).tnum = ++t_num;
    T->nodes[0].x = (double)xmin-abs(xmax); T->nodes[0].y = (double)ymin-abs(ymax);
    T->nodes[1].x = (double)xmax+abs(xmin); T->nodes[1].y = (double)ymin-abs(ymax);
    T->nodes[2].x = (double)(xmin+xmax)/2; T->nodes[2].y = (double)ymax+abs(ymin);
    temp.x = 0.0; temp.y = 10.0;
    //insert(T,temp);
    //print_tree(T);
    temp.x = 0.0; temp.y = 7.0;
    //insert(T,temp);
    //print_tree(T);
    //printf("checking for inside %d\n", point_inside(T,temp)); 
    return 0;
}
    // *(*(nodes+i)+0) = (double)x;
    // *(*(nodes+i)+1) = (double)y;
    // xavg /= nnodes;
    // yavg /= nnodes;
//	xavg += x;
     //	yavg += y;
//while(  != EOF)
	//printf("dummy= %d x = %f , y = %f \n",dummy,x,y);
// printf("# of nodes is %d %d %d %d\n",nnodes,dummy,x,y);
    //fscanf(fp, " %d %d %d\n",&dummy,&x,&y);
    //printf(" %d %d %d\n", dummy,x,y);
	/*v0.x = t->nodes[3].x-t->nodes[0].x; v0.y = t->nodes[3].y-t->nodes[0].y;
	v1.x = t->nodes[2].x-t->nodes[0].x; v1.y = t->nodes[2].y-t->nodes[0].y;
	v2.x = n.x-t->nodes[0].x; v2.y = n.y-t->nodes[0].y;
	dot00 = dotp(v0, v0);
	dot01 = dotp(v0, v1);
	dot02 = dotp(v0, v2);
	dot11 = dotp(v1, v1);
	dot12 = dotp(v1, v2);
	printf("Dots are, 00=%f, 01=%f,02=%f,11=%f,12=%f\n",dot00,dot01,dot02,dot11,dot12);
	invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	v = (dot00 * dot12 - dot01 * dot02) * invDenom;
	w = 1-u-v;*/
 //insert(T,temp);
    //printf("After inserting first point\n");
    //print_tree(T); 
    //printf("Printing the original tree\n");
    //print_tree(T);
    //temp = (struct node)malloc(sizeof(struct node)*1);
    //temp.x = -11.0; temp.y = -5.0;
