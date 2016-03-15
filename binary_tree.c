#include<stdio.h>
#include<stdlib.h>

	struct node{
		double x,y;
	};
	struct triangle{
		int tnum;
		struct node nodes[3];
		struct triangle *l;
		struct triangle *c;
		struct triangle *r;
	};
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
	int main(){
		struct triangle *T=NULL;
		struct node p;
		p.x = 10;
		p.y = 10;
		T = (struct triangle*)malloc(sizeof(struct triangle));
		T->nodes[0].x = -20;
		T->nodes[0].y = -14;
		T->nodes[1].x = 20;
		T->nodes[1].y = -14;
		T->nodes[2].x = 10;
		T->nodes[2].y = 14;
		T->l = NULL; T->c = NULL; T->r = NULL;
		print_tree(T);
		insert(T,p);
		print_tree(T);
		return 0;
	}
		
		
		
		
		
	
