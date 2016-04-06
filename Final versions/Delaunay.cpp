#include "Tri_project.h"

#define sendl std::endl
#define scout std::cout
#define scin std::cin

using namespace std; 

/* The following function checks to see if the given point lies inside the circum circle.
 It returns true if the point lies inside, it also returns the computed circum center and circum radius.
*/

int CircumCircle(double xp, double yp, double x1, double y1, double x2, 
		 double y2, double x3, double y3, double &xc, double &yc, double &r){
	double m1, m2, mx1, mx2, my1, my2;
	double deltax, deltay, rsqr, drsqr;

    	/* Checking for points lying on the same line */
    	if (abs(x1 - x2) < EPSILON && abs(x2 - x3) < EPSILON)
      		return(false);
      	/* The following block of code is when x2 and x1 are same */
    	if(abs(x2-x1) < EPSILON){ 
      		m2 = - (x3 - x2) / (y3 - y2);
      		mx2 = (x2 + x3) / 2.0;
      		my2 = (y2 + y3) / 2.0;
      		xc = (x2 + x1) / 2.0;
      		yc = m2 * (xc - mx2) + my2;
    	}
    	/* The following block of code is when x3 and x2 are same */
    	else if(abs(x3 - x2) < EPSILON){ 
        	m1 = - (x2 - x1) / (y2 - y1);
            	mx1 = (x1 + x2) / 2.0;
            	my1 = (y1 + y2) / 2.0;
            	xc = (x3 + x2) / 2.0;
            	yc = m1 * (xc - mx1) + my1;
    	}
    	/* if no two points are on same line */
    	else{
            	m1 = - (x2 - x1) / (y2 - y1); 
            	m2 = - (x3 - x2) / (y3 - y2); 
            	mx1 = (x1 + x2) / 2.0; 
            	mx2 = (x2 + x3) / 2.0;
            	my1 = (y1 + y2) / 2.0;
            	my2 = (y2 + y3) / 2.0;
            	xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2); 
            	yc = m1 * (xc - mx1) + my1; 
    	}
    	deltax = x2 - xc;
    	deltay = y2 - yc;
    	rsqr = dx * dx + dy * dy;
    	r = sqrt(rsqr); 
    	dx = xp - xc;
    	dy = yp - yc;
    	drsqr = dx * dx + dy * dy;
    	return((drsqr < rsqr) ? true : false);
}
///////////////////////////////////////////////////////////////////////////////
// Triangulate() :
//   Triangulation subroutine
//   Takes as input NV vertices in array pnode
//   Returned is a list of ntri triangular faces in the array v
//   These triangles are arranged in a consistent clockwise order.
//   The triangle array 'v' should be malloced to 3 * nv
//   The vertex array pnode must be big enough to hold 3 more points
//   The vertex array must be sorted in increasing x values say
//
//   qsort(p,nv,sizeof(node),nodeCompare);
///////////////////////////////////////////////////////////////////////////////

int Triangulate(int nv, node pnode[], Triangle v[], int &ntri){
	int *complete = NULL;
    	EDGE *edges = NULL; 
    	EDGE *p_EdgeTemp;
    	int trimax, status, inside, i, j, k, nedge, emax;
    	nedge = 0; emax = 200; status = 0;
    	double xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r;
    	double xmin, xmax, ymin, ymax, xmid, ymid;
    	double dx, dy, dmax; 

  	/* Allocate memory for the completeness list, flag for each triangle */
    	trimax = 4 * nv;
    	complete = new int[trimax];
  	/* Allocate memory for the edge list */
   	edges = new EDGE[emax];
 	/* Find the maximum and minimum vertex bounds.
        This is to allow calculation of the bounding triangle */
    	xmin = pnode[0].x;
    	ymin = pnode[0].y;
    	xmax = xmin;
    	ymax = ymin;
    	for(i = 1; i < nv; i++){
    		if (pnode[i].x < xmin) xmin = pnode[i].x;
      		if (pnode[i].x > xmax) xmax = pnode[i].x;
      		if (pnode[i].y < ymin) ymin = pnode[i].y;
      		if (pnode[i].y > ymax) ymax = pnode[i].y;
    	}
    	dx = xmax - xmin;
    	dy = ymax - ymin;
    	dmax = (dx > dy) ? dx : dy;
    	xmid = (xmax + xmin) / 2.0;
    	ymid = (ymax + ymin) / 2.0;
  	/* Set up the supertriangle The supertriangle is the first triangle in the triangle list. */
    	pnode[nv+0].x = xmid - 20 * dmax;
    	pnode[nv+0].y = ymid - dmax;
    	pnode[nv+1].x = xmid;
    	pnode[nv+1].y = ymid + 20 * dmax;
    	pnode[nv+2].x = xmid + 20 * dmax;
    	pnode[nv+2].y = ymid - dmax;
    	v[0].p1 = nv;
    	v[0].p2 = nv+1;
    	v[0].p3 = nv+2;
    	complete[0] = false;
    	ntri = 1;
  	/* Include each point one at a time into the existing mesh */
    	for(i = 0; i < nv; i++){
      		xp = pnode[i].x;
      		yp = pnode[i].y;
      		nedge = 0;
  		/* Set up the edge buffer. If the point (xp,yp) lies inside the circumcircle then the
       	   	   three edges of that triangle are added to the edge buffer and that triangle is removed.
  		*/
      		for(j = 0; j < ntri; j++){
      			if(complete[j])
        			continue;
      			x1 = pnode[v[j].p1].x; 		y1 = pnode[v[j].p1].y;
      			x2 = pnode[v[j].p2].x;      	y2 = pnode[v[j].p2].y;
      			x3 = pnode[v[j].p3].x;      	y3 = pnode[v[j].p3].y;
      			inside = CircumCircle(xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r);
      			if (yc + r < yp) // Suggested: if (xc + r + EPSILON < xp)
        			complete[j] = true;
      			if (inside){
				/* Check that we haven't exceeded the edge list size */
        			if (nedge + 3 >= emax){
        				emax += 100;
        				p_EdgeTemp = new EDGE[emax];
        				for (int i = 0; i < nedge; i++)  // Fix by John Bowman
        					p_EdgeTemp[i] = edges[i];   
               				delete []edges;
          				edges = p_EdgeTemp;
        			}
        			edges[nedge+0].p1 = v[j].p1;
        			edges[nedge+0].p2 = v[j].p2;
       	 			edges[nedge+1].p1 = v[j].p2;
        			edges[nedge+1].p2 = v[j].p3;
        			edges[nedge+2].p1 = v[j].p3;
        			edges[nedge+2].p2 = v[j].p1;
        			nedge += 3;
         			v[j] = v[ntri-1];
        			complete[j] = complete[ntri-1];
        			ntri--;
        			j--;
      			}
    		}
    		/* Tag multiple edges
  		Note: if all triangles are specified anticlockwise then all
  		interior edges are opposite pointing in direction. */
    		for(j = 0; j < nedge - 1; j++){
      			for(k = j + 1; k < nedge; k++){
        			if((edges[j].p1 == edges[k].p2) && (edges[j].p2 == edges[k].p1)){
        				edges[j].p1 = -1;
        				edges[j].p2 = -1;
        				edges[k].p1 = -1;
        				edges[k].p2 = -1;
        			}
         			/* Shouldn't need the following, see note above */
        			if((edges[j].p1 == edges[k].p1) && (edges[j].p2 == edges[k].p2)){
        				edges[j].p1 = -1;
        				edges[j].p2 = -1;
        				edges[k].p1 = -1;
        				edges[k].p2 = -1;
        			}
      			}
    		}
    		/* Form new triangles for the current point, Skipping over any tagged edges.
     	   	All edges are arranged in clockwise order. */
    		for(j = 0; j < nedge; j++) {
      			if(edges[j].p1 < 0 || edges[j].p2 < 0)
      				continue;
      			v[ntri].p1 = edges[j].p1;
      			v[ntri].p2 = edges[j].p2;
      			v[ntri].p3 = i;
      			complete[ntri] = false;
      			ntri++;
    		}
   	}
   	std::cout << "ntri before removing super triangle" << ntri << std::endl;
   	/* Remove triangles with supertriangle vertices. These are triangles which have a vertex number greater than nv */
   	for(i = 0; i < ntri; i++) {
   		if(v[i].p1 >= nv || v[i].p2 >= nv || v[i].p3 >= nv) {
   			v[i] = v[ntri-1];
      			ntri--;
      			i--;
    		}
   	}
  	//std::cout << "ntri before exiting" << ntri << std::endl;
  	delete[] edges;
  	delete[] complete;
  	return 0;
} 

/*void randomize(){
	srand((time_t) time(NULL));  
}

int random(int n){
	return rand()%n; 
}*/

int nodeCompare(const void *v1, const void *v2){
	node *p1, *p2;
	p1 = (node*)v1;
	p2 = (node*)v2;
	if (p1->y < p2->y)
    		return(-1);
  	else if (p1->y > p2->y)
        	return(1);
       	else
        	return(0);
}

void outputtriangle(int &nv, node p[], Triangle v[], int &ntri){
	int X, Y, i = 0;
	int max = 4;
	double x, y;
	//std::cout << "ntri = " << ntri << std::endl;
	for(int i = 0; i < ntri; i++){// replace cout by a compatible lineto to trace
		std::cout << "(" << p[v[i].p1].x << "," << p[v[i].p1].y << ")" << " " 
      		<< "(" << p[v[i].p2].x << "," << p[v[i].p2].y << ")" << " " << 
      		"(" << p[v[i].p3].x << "," << p[v[i].p3].y << ")" << std::endl;
      	}
}

int main(){
	Triangle *v = NULL;
	ifstream inputFile;
	ofstream outFile,nodeOut;
	int X, Y, i, max, ntri, nv, dummy, tmp2, tmp3, tmp4;
	string fname;
	double x, y;
 
	ntri = 0;
  	std::cout << "Please enter the name of file without .node extension" << std::endl;
  	std::cin >> fname;
	inputFile.open(fname + ".node",ios::in);
	inputFile >> max >> tmp2 >> tmp3 >> tmp4;
	nv = max;
	std::cout << "max read is " << max << std::endl;
	node *p = new node[max+3]; 
	for(i=0;i < max; i++)
		inputFile >> dummy >> p[i].x >> p[i].y;
	inputFile.close();
  	v = new Triangle[3 * nv]; 
  	
  	qsort(p, nv, sizeof(node), nodeCompare);
  	
  	clock_t t1,t2;
  	t1 = clock();
  	
  	Triangulate(nv, p, v, ntri);
  	
  	t2 = clock();
  	std::cout << "ntri " << ntri << std::endl;
  	std::cout << "Time taken is " << ((float)t2-(float)t1)/CLOCKS_PER_SEC << std::endl;
  	
  	nodeOut.open(fname + ".d.node", ios::out);
  	nodeOut << max << " " << tmp2 << " " << tmp3 << " " << tmp4  << sendl;
  	for(i = 0; i < max; i++)
  		nodeOut << i+1 << " " << p[i].x << " " << p[i].y << sendl;
  	nodeOut.close();
  	outFile.open(fname + ".d.ele",ios::out);
  	outFile << ntri << std::endl;
  	for (i=0;i < ntri; i++)
  		outFile << (i+1) << " " << (v[i].p1)+1 << " " << (v[i].p2)+1 << " " << (v[i].p3)+1 << sendl;
  	outFile.close();
  	 
  	//outputtriangle(nv, p, v, ntri); // use this fonction to trace the mesh (via 
  	delete []p;// OpenGL, DirectX, ...)
  	delete []v;
  	//p = NULL;
  	//v = NULL;
  	//system("pause"); // remove under Linux
  	return 0;
}
