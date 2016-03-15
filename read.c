#include<stdio.h>
#include<stdlib.h>
#include<math.h>

struct node{
	int nnum;
	double x,y;
};
struct triangle{
	int tnum;
	struct node nodes[3];
	struct triangle *l;
	struct triangle *c;
	struct triangle *r;
};

void print_tree(struct triangle *root){
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
void insert(struct triangle *root,struct node p){
		if((root->l == NULL) && (root->c == NULL) && (root->r == NULL)){
			struct triangle *t1,*t2,*t3;
			t1 = (struct triangle*)malloc(sizeof(struct triangle));
			t2 = (struct triangle*)malloc(sizeof(struct triangle));
			t3 = (struct triangle*)malloc(sizeof(struct triangle));
			t1->nodes[0] = root->nodes[0];
			t1->nodes[1] = root->nodes[1];
			t1->nodes[2] = p; t1->l = NULL; t1->c = NULL; t1->r = NULL;
			root->l = t1;
			t2->nodes[0] = root->nodes[1];
			t2->nodes[1] = root->nodes[2];
			t2->nodes[2] = p; t2->l = NULL; t2->c = NULL; t2->r = NULL;
			root->c = t2;
			t3->nodes[0] = root->nodes[2];
			t3->nodes[1] = root->nodes[0];
			t3->nodes[2] = p; t3->l = NULL; t3->c = NULL; t3->r = NULL;
			root->r = t3;
		}
	}
double dotp(struct node n1, struct node n2){
	return (n1.x*n2.x + n1.y*n2.y);
}	
int point_inside(struct triangle *t, struct node n){
	struct node v0, v1, v2;
	double u,v,w,invDenom,dot00,dot01,dot02,dot11,dot12;
	v0.x = t->nodes[3].x-t->nodes[0].x; v0.y = t->nodes[3].y-t->nodes[0].y;
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
	w = 1-u-v;
	printf("u = %f, v = %f, w = %f\n", u,v,w);
	if( (u > 0) && (v > 0) && (w > 0))
		return 1;
	else
		return 0;
}
// The following function is used to locate a tree node in which the input node lies
struct triangle* locate(struct triangle *t, struct node n){
	
	
	
}
int main(){
    int nnodes,i;
    double xavg=0, yavg=0;
    struct node *input_nodes; // to keep track of input nodes.
    struct triangle *T=NULL;  // The big triangle
    struct node temp; // temporarily store x,y values for a node.
    int x,y,dummy,xmin=0, xmax=0, ymin=0, ymax=0;
    char* fname;
    FILE* fp;
    fp = fopen("4.node","r");
    fscanf(fp, "%d %d %d %d\n",&nnodes,&dummy,&x,&y);
    input_nodes = (struct node* )malloc(sizeof(struct node)*nnodes);
    for(i = 0; i<nnodes; i++){
    	fscanf(fp,"%d %d %d\n",&dummy, &x, &y);
    	printf("x = %d, y = %d\n",x,y);
    	input_nodes[i].x = (double)x;
    	input_nodes[i].y = (double)y;
    	if(x >= xmax) 
  			xmax = x;
    	else if(x <= xmin) 
   			xmin = x;
	    if(y >= ymax) 
  			ymax = y;
    	else if(y <= ymin) 
 			ymin = y;
 	}
 	/*printf("xmax is %d, xmin is %d, ymax is %d, ymin is %d\n",xmax, xmin,ymax, ymin);	
    printf("left is %f , %f\n",(double)xmin-abs(xmax),(double)ymin-abs(ymax));
    printf("right is %f, %f\n",(double)xmax+abs(xmin),(double)ymin+abs(ymax));
    printf("top is %f , %f\n",(double)(abs(xmin)+xmax)/2,(double)ymax+abs(ymin));*/
    // Initializing the big triangle.
    T = (struct triangle*)malloc(sizeof(struct triangle));
    T->l = NULL; T->c = NULL; T->r = NULL;
    T->nodes[0].x = (double)xmin-abs(xmax); T->nodes[0].y = (double)ymin-abs(ymax);
    T->nodes[1].x = (double)xmax+abs(xmin); T->nodes[1].y = (double)ymin+abs(ymax);
    T->nodes[2].x = (double)(abs(xmin)+xmax)/2; T->nodes[2].y = (double)ymax+abs(ymin);
    //printf("Printing the original tree\n");
    //print_tree(T);
    //temp = (struct node)malloc(sizeof(struct node)*1);
    //temp.x = -11.0; temp.y = -5.0;
    temp.x = 0.0; temp.y = 0.0;
    //insert(T,temp);
    //printf("After inserting first point\n");
    //print_tree(T);  
    printf("checking for inside %d\n", point_inside(T,temp)); 
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
