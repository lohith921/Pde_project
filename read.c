#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mmap_lib.h"

struct node{
	int nnum;
	double x,y;
	};
struct tree_node {
 	struct tree_node *lchild = NULL; //left child
	struct tree_node *mchild = NULL;
	struct tree_node *rchild = NULL;
	struct triangle *link_to_triangle = NULL;
	};
struct triangle{
	int tnum;
	struct node *n1, *n2, *n3;
	struct triangle *nb1 = NULL, *nb2 = NULL, *nb3 = NULL;
	struct tree_node *link_to_tree = NULL;
	struct triangle *link_to_next_triangle = NULL;
	//struct triangle *pt;
	//struct triangle *l;
	//struct triangle *c;
	//struct triangle *r;
	};
struct mesh{
	struct triangle *T = NULL;
	int num_of_nodes, num_of_tris;
	};

void print_tree(struct tree_node *root){
		int i;
		if(root != NULL){
			printf("Node nums are %d , %d, %d\n",root->link_to_triangle->n1.num, root->link_to_triangle->n2.num, 
				root->link_to_triangle->n3.num);
			print_tree(root->lchild);
			print_tree(root->mchild);
			print_tree(root->rchild);
		}
		else{
			return;
		}
	}
	
void insert(struct tree_node *root, struct triangle *tri_list, struct node p){
		if(point_inside(root->link_to_triangle,p) && (root->lchild == NULL && root->mchild == NULL && root->rchild == NULL)){
			struct triangle *t1,*t2,*t3;
			t1 = (struct triangle*)malloc(sizeof(struct triangle));
			t2 = (struct triangle*)malloc(sizeof(struct triangle));
			t3 = (struct triangle*)malloc(sizeof(struct triangle));
			t1->n1 = root->link_to_triangle->n1;
			t1->n2 = root->link_to_triangle->n2;
			t1->n3 = p; 
			t1->nb1 = t2; t1->nb2 = t3; t->nb3 = root->link_to_triangle->nb3;
			
			root->lchild = t1;
			t1->link_to_tree = root->lchild;
			root->link_to_triangle->link_to_next_triangle = t1;
			
			t2->n1 = root->n2;
			t2->n2 = root->n3;
			t2->n3 = p; 	t2->nb1 = t3; t2->nb2 = t1; t2->nb3 = root->link_to_triangle->nb1;
			
			root->mchild = t2;
			t2->link_to_tree = root->mchild;
			root->link_to_triangle->link_to_next_triangle = t2;
			
			t3->n1 = root->n3;
			t3->n2 = root->n1;
			t3->n3 = p; 	t3->nb1 = t1; t3->nb2 = 2; t3->nb3 = root->link_to_triangle->nb2;
			
			root->rchild = t3;
			t3->link_to_tree = root->rchild;
			root->link_to_triangle->link_to_next_triangle = t3;
			
		}
		else if( point_inside(root->link_to_triangle,p)){
			insert(root->lchild, p);
			insert(root->mchild, p);
			insert(root->rchild, p);
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
	v0 = t->n1; v1 = t->n2; v2 = t->n3;
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
    struct triangle *T = NULL;  // The big triangle
    struct node temp; // temporarily store x,y values for a node.
    //struct mesh *m;
    struct tree_node *root;
    float x,y,xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0;
    char* fname;
    FILE* fp;
    fp = fopen("4.node","r");
    fscanf(fp, "%d %d %f %f\n",&nnodes,&dummy,&x,&y);
    input_nodes = (struct node*)malloc(sizeof(struct node)*(nnodes+3));
    map_t* node_map;
    node_map = map_create();
    *node_num = 3; *t_num = 0;
    for(i = 0; i < nnodes; i++){
    	fscanf(fp,"%d %f %f\n",&dummy, &x, &y);
    	//printf("x = %d, y = %d\n",x,y);
    	input_nodes[*node_num].num = *node_num;
    	input_nodes[*node_num].x = (double)x;
    	input_nodes[*node_num].y = (double)y;
    	//map_set(&node_map,node_num,x,y); // setting point cloud into the map
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
 	//display_map(node_map);
 	/*printf("xmax is %d, xmin is %d, ymax is %d, ymin is %d\n",xmax, xmin,ymax, ymin);	
    printf("left is %f , %f\n",(double)xmin-abs(xmax),(double)ymin-abs(ymax));
    printf("right is %f, %f\n",(double)xmax+abs(xmin),(double)ymin-abs(ymax));
    printf("top is %f , %f\n",(double)(xmin+xmax)/2,(double)ymax+abs(ymin));*/
    input_nodes[0].num = 0;
    input_nodes[0].x = (double)xmin-abs(xmax);
    input_nodes[0].y = (double)ymin-abs(ymax);
    input_nodes[1].num = 1;
    input_nodes[1].x = (double)xmax+abs(xmin);
    input_nodes[1].y = (double)ymin-abs(ymax);
    input_nodes[2].num = 2;
    input_nodes[2].x = (double)(xmin+xmax)/2;
    input_nodes[2].y = (double)ymax+abs(ymin);
    
    
    // Initializing the big triangle.
    
    T = (struct triangle*)malloc(sizeof(struct triangle));
    T->n1 = input_nodes[0]; T->n2 = input_nodes[1]; T->n3 = input_nodes[2];
    root->link_to_triangle = T;
    T->link_to_tree = root;
    
    for(i = 3; i < nnodes; i++){
    	insert(root,T,input_nodes[i]); // calling insert with root, large triangle T, node
    	}
    print_tree(root);
    
    
    
   // temp.x = 0.0; temp.y = 10.0;
    //insert(T,temp);
    //print_tree(T);
    //temp.x = 0.0; temp.y = 7.0;
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
