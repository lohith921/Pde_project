#include<iostream>
#include<vector>

struct Vertex{
       double  x,y,z;
       Vertex(double X,double Y,double Z):
       x(X),y(Y),z(Z){}
       };

struct Edge{
       Vertex *vmin,*vmax;
       Edge(Vertex *v1,Vertex *v2)
       vmin = std::min(v1,v2);
       vmax = std::max(v1,v2);
       }      
       bool operator <(const Edge &other)const{
            if (vmin < other.vmin) return true;
            if (vmin > other.vmin) return false;
            if (vmax < other.vmax) return true;
            return false;
            }
       };

struct Face{
       Face *F[3];
       Vertex*V[3];
       bool deleted;
       Face(Vertex *v0,Vertex *v1,Vertex *v2){
                   V[0]=v0;V[1]=v1;V[2]=v2;
                   F[0]=F[1]=F[2]=NULL;
                   deleted = false;
                   }
       Edge getEdge(int k){
            return Edge(V[k],V[(k+1)%3]);
       }
       bool inCircle(Vertex *c);
       Vertex centroid();
       };

void computeAdjacencies(std::vector<Face*> &cavity){
     std::map<Edge,std::pair<int,Face*>> edgeToFace;
     for(int iFace = 0;iFace < cavity.size();iFace++){
             for(int iEdge = 0;iEdge < 3;iEdge++){
                     Edge edge = cavity[iFace]->getEdge(iEdge);
                     std::map<Edge,std::pair<int,Face*>>::iteratorit = edgeToFace.find(edge);
                     if (it == edgeToFace.end()){
                               //edgehasnotyetbeentouched,socreateanentry
                               edgeToFace.insert(std::make_pair(edge,std::make_pair(iEdge,cavity[iFace])));
                     }
                     else{
                               //Connectthetwoneighboringtriangles
                               cavity[iFace]->F[iEdge] = it->second.second;
                               it->second.second->F[it->second.first] = cavity[iFace];
                               //Eraseedgefromthemap
                               edgeToFace.erase(it);
                     }
             }
     }            
}

void delaunayCavity(Face *f,Vertex *v,std::vector<Face*> &cavity,std::vector<Edge> &bnd,std::vector<Face*>
     &otherSide){
     if (f->deleted) return;
        f->deleted = true;//Markthetriangle
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

Face* lineSearch(Face *f,Vertex *v){
while(1){
if(f==NULL)returnNULL;//weshouldNEVERreturnhere
if(f°>inCircle(v))returnf;
Vertexc=f°>centroid();
for(intiNeigh=0;iNeigh<3;iNeigh++){
Edgee=f°>getEdge(iNeigh);
if(orientationTest(&c,v,e.vmin)*
orientationTest(&c,v,e.vmax)<0&&
orientationTest(e.vmin,e.vmax,&c)*
orientationTest(e.vmin,e.vmax,v)<0){
f=f°>F[iNeigh];
break;
}
}
}

voiddelaunayTrgl(std::vector<Vertex*>&S,std::vector<Face*>&T){
for(intiP=0;iP<S.size();iP++){
Face*f=lineSearch(T[0],S[iP]);
std::vector<Face*>cavity;
std::vector<Edge>bnd;
std::vector<Face*>otherSide;
delaunayCavity(f,S[iP],cavity,bnd,otherSide);
if(bnd.size()!=cavity.size()+2)throw;
for(inti=0;i<cavity.size();i++){
//reusememoryslotsofinvalidelements
cavity[i]°>deleted=false;
cavity[i]°>F[0]=cavity[i]°>F[1]=cavity[i]°>F[2]=NULL;
cavity[i]°>V[0]=bnd[i].V[0];
cavity[i]°>V[1]=bnd[i].V[1];
cavity[i]°>V[2]=S[iP];
}
unsignedintcSize=cavity.size();
for(inti=cSize;i<cSize+2;i++){
Face*newf=newFace(bnd[i].V[0],bnd[i].V[1],S[iP]);
T.push_back[newf];
cavity.push_back(newf);
}
for(inti=0;i<otherSide.size();i++)
if(otherSide[i])cavity.push_back(otherSide[i]);
computeAdjacencies(cavity);
}
}
