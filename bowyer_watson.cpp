#include<iostream>
#include<vector>
#include<map>

#define scout std::cout
#define sendl std::endl
struct Vertex{
       double  x,y;
       Vertex(double X,double Y):
       x(X),y(Y){}
       };

struct Edge{
       Vertex *vmin,*vmax;
       Edge(Vertex *v1,Vertex *v2){
       	vmin = std::min(v1,v2);
       	vmax = std::max(v1,v2);
       }      
       bool operator <(const Edge &other) const{
            if (vmin < other.vmin) return true;
            if (vmin > other.vmin) return false;
            if (vmax < other.vmax) return true;
            return false;
       }
};

struct Face{
       Face *F[3];
       Vertex *V[3];
       bool deleted;
       Face(Vertex *v0,Vertex *v1,Vertex *v2){
                   V[0] = v0; V[1] = v1; V[2] = v2;
                   F[0] = F[1] = F[2] = NULL;
                   deleted = false;
                   }
       Edge getEdge(int k){
            return Edge(V[k],V[(k+1)%3]);
       }
       bool inCircle(Vertex *c1){
       		double a, b, c, d, e, f, g, h, i, det;
	 	a = V[0]->x - c1->x; 
	 	b = V[0]->y - c1->y;
 		c = a*a + b*b;
 	
 		d = V[1]->x - c1->x; 
 		e = V[1]->y - c1->y;
 		f = d*d + b*b;
 		
 		g = V[2]->x - c1->x; 
 		h = V[2]->y - c1->y;
 		h = g*g + h*h;
 		
 		det = a*(e*i - h*f) - b*(d*i-g*f) + c*(d*h - e*g);
 		std::cout << det << " "<< (det > 0) << std::endl;
 		return (det >= 0);       
       }
       Vertex centroid(){
       		double c1 = (V[0]->x + V[1]->x + V[2]->x)/3.0;
       		double c2 = (V[0]->y + V[1]->y + V[2]->y)/3.0;
       		return Vertex(c1,c2);
       }
     };

void computeAdjacencies(std::vector<Face*> &cavity){
     std::map<Edge,std::pair<int,Face*> > edgeToFace;
     for(int iFace = 0;iFace < cavity.size();iFace++){
             for(int iEdge = 0;iEdge < 3;iEdge++){
                     Edge edge = cavity[iFace]->getEdge(iEdge);
                     std::map<Edge,std::pair<int,Face*> >::iterator it = edgeToFace.find(edge);
                     if (it == edgeToFace.end()){
                               //edgehasnotyetbeentouched,socreateanentry
                               edgeToFace.insert(std::make_pair(edge,std::make_pair(iEdge,cavity[iFace])));
                     }
                     else{
                               //Connectthetwoneighboringtriangles
                               cavity[iFace]->F[iEdge] = it->second.second;
                               it->second.second->F[it->second.first] = cavity[iFace];
                               //Erase edge from the map
                               edgeToFace.erase(it);
                     }
             }
     }            
}

void delaunayCavity(Face *f,Vertex *v,std::vector<Face*> &cavity,std::vector<Edge> &bnd,std::vector<Face*>
     &otherSide){
     if (f->deleted) return;
        f->deleted = true;//Mark the triangle
        cavity.push_back(f);
        for(int iNeigh = 0;iNeigh < 3;iNeigh++){
                if (f->F[iNeigh] == NULL){
                   bnd.push_back(f->getEdge(iNeigh));
                }
                else if(!f->F[iNeigh]->inCircle(v)){
                   bnd.push_back(f->getEdge(iNeigh));
                   if (!f->F[iNeigh]->deleted){
                       otherSide.push_back(f->F[iNeigh]);
                       f->F[iNeigh]->deleted = true;
                   }
                }
                else delaunayCavity(f->F[iNeigh],v,cavity,bnd,otherSide);
        }
}

double orientationTest(Vertex *a, Vertex *b, Vertex *c){
	double ax = a->x - c->x;
	double ay = a->y - c->y;
	double bx = b->x - c->x;
	double by = b->y - c->y;
	return ax*by - ay*bx;
}

Face* lineSearch(Face *f,Vertex *v){
	while(1){
		if(f == NULL) return NULL;//we should NEVER return here
		if(f->inCircle(v)) return f;
		Vertex c = f->centroid();
		for(int iNeigh = 0;iNeigh < 3;iNeigh++){
			Edge e = f->getEdge(iNeigh);
			if (orientationTest(&c,v,e.vmin) * orientationTest(&c,v,e.vmax) < 0 &&
			   orientationTest(e.vmin,e.vmax,&c) * orientationTest(e.vmin,e.vmax,v) < 0){
				f = f->F[iNeigh];
				break;
			}
		}
	}
}

void delaunayTrgl(std::vector<Vertex*> &S, std::vector<Face*> &T){
	for(int iP = 0;iP < S.size();iP++){
		Face *f = lineSearch(T[0],S[iP]);
		std::vector<Face*> cavity;
		std::vector<Edge> bnd;
		std::vector<Face*> otherSide;
		delaunayCavity(f,S[iP],cavity,bnd,otherSide);
		scout << "error" <<sendl;
		if (bnd.size() != cavity.size()+2) throw;
		for(int i = 0;i < cavity.size();i++){
			//reuse memory slots of invalid elements
			cavity[i]->deleted = false;
			cavity[i]->F[0] = cavity[i]->F[1] = cavity[i]->F[2] = NULL;
			cavity[i]->V[0] = bnd[i].vmin;//V[0];
			cavity[i]->V[1] = bnd[i].vmax; //V[1];
			cavity[i]->V[2] = S[iP];
		}
		unsigned int cSize = cavity.size();
		for(int i = cSize;i < cSize+2;i++){
			Face *newf = new Face(bnd[i].vmin,bnd[i].vmax,S[iP]);
			T.push_back(newf);
			cavity.push_back(newf);
		}
		for(int i = 0;i < otherSide.size();i++)
			if (otherSide[i]) 
				cavity.push_back(otherSide[i]);
		computeAdjacencies(cavity);
	}
}
int main(){
	/*Vertex *ll = new Vertex(-12,-8);
	Vertex *rl = new Vertex(12,-8);
	Vertex *lt = new Vertex(-12,12);
	Vertex *rt = new Vertex(12,12);*/
	Vertex *v1 = new Vertex(-20,-14);
	Vertex *v2 = new Vertex(20,-14);
	Vertex *v3 = new Vertex(0,14);
	Face *f1 = new Face(v1,v2,v3);
	/*Face *f1 = new Face(ll,rl,rt);
	//Face *f2 = new Face(ll,rt,lt);
	f1->F[2] = f2;
	f2->F[0] = f1;*/
	std::vector<Face*> T;
	T.push_back(f1);
	//T.push_back(f2);
	//vertices.push_back(ll);
	//vertices.push_back(rl);
	//vertices.push_back(lt);
	//vertices.push_back(rt);
	//Vertex *v4 = new Vertex(-10,-4);
	Vertex *v8 = new Vertex(0.0,0.0);
	Vertex *v5 = new Vertex(0.0,7.5);
	Vertex *v6 = new Vertex(0.0,10.0);
	//Vertex *v7 = new Vertex(10,-4);
	
	
	std::vector<Vertex*> vertices;
	//vertices.push_back(v4);
	vertices.push_back(v8);
	vertices.push_back(v5);
	vertices.push_back(v6);
	//vertices.push_back(v7);
	
	delaunayTrgl(vertices,T);
	for (int i=0; i < T.size(); i++){
		for(int j=0;j<3; j++){
			scout << "(" << T[i]->V[j]->x << ", " << T[i]->V[j]->y << ") ";
		}
		scout << "deleted " << T[i]->deleted <<sendl;
	}
	return 0;
}
/*
   	scout << V[0]->x << V[0]->y << sendl;
       	scout << V[1]->x << V[1]->y << sendl;
       	scout << V[2]->x << V[2]->y << sendl;
       	scout << c1->x << c1->y << sendl;
       */		
