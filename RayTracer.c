/*
 * Lightweight Geometric Optics Simulator v1.0. Created for the subject "Fundamentos de Campos y Ondas"
 * 
 * We have written and commented the program in English because we were
 * told that it would not suposse any problem, and besides English is the
 * natural language for computer programming.
 * 
 * Should the teachers encounter any difficulty in the use of the program,
 * or just for any technical inquiry, do not hesitate to contact us.
 * 
 * The Optics project group.
 * 
 * */


//We use "generador.h" library, from the subject "Física Computacional",
//to aid us in the generation of pseudorandom numbers
#include "generador.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI acos(-1.0)

typedef struct Point2D{
	double x;
	double y;
}Point2D;
static Point2D NAP={NAN,NAN};	//Not A Point, with similar functionality as NaN

typedef struct Intersection2D{
	Point2D P;
	double alpha;	//following DIN standard
}Intersection2D;

typedef struct Point2DArray{
	Point2D* data;	//pointer to the first of the nodes, joining them 
					//with straight segments we get the trayectoy of the ray
	int N;	//max number of nodes
}Point2DArray;

typedef struct Segment{
	Point2D P1;
	Point2D P2;
}Segment;

typedef struct Ray{
	Point2DArray  points;
	int c;	//position of last node
	double lambda;	//wavelength
}Ray;

typedef struct Circle{
	Point2D c;	//center
	double r;	//radius
}Circle;

//Note: we will pass by pointer all our rays, even to the functions that
//don't need to modify them, just for the sake of consistency
typedef void (*OpticalSystem)(Ray* r);

//Our OpticalSystems
void oneLens(Ray* r);
void twoLenses(Ray* r);

//Point2D utilities
double distance(Point2D p1, Point2D p2);
int Point2DEquals(Point2D p1,Point2D p2);
void fprintPoint2D(FILE* OUT,Point2D p);

//Ray2D utilities
Ray newRay(Point2D A, Point2D B, double lambda, int N);
void destroyRay(Ray* r);
Point2D getPoint(Ray* r, int n);
void addPoint(Ray* r, Point2D P);
void setPoint(Ray* r, int n, Point2D P);
void fprintRay(FILE* OUT, Ray* r);

//Collision utilities
Intersection2D getRayCircleIntersection(Ray* r,Circle c,int ext);
Intersection2D getRaySegmentIntersection(Ray* r,Segment s,int ext);

//Optics utilites
void calculateRayRefraction(Ray* r, Intersection2D i, double n1, double n2);

//Output utilites with various constants
double XSTART=-20;
double XEND=-15;
double XAXISSTART=-1000;
double XAXISEND=1000;
double YAXISSTART=-1000;
double YAXISEND=1000;
double LAMBDA=0;
double POINTS=10;
void plotDeterministic1DX(FILE* OUT,OpticalSystem os,double s,double e,double n);
void plotDeterministic1DY(FILE* OUT,OpticalSystem os,double x,double s,double e,double n);
void plotDeterministic2DYZ(FILE* OUT,OpticalSystem os, double x, double r, double ny, double nz);
void plotRandom2DYZ(FILE* OUT,OpticalSystem os,double x,double r,double n);

//Theoretical values for the single lens and Clerk lenses
double FOCAL1=14.90797425;
double FOCAL2=24.73037;

int main(int argc, char **argv){
	
	//Usage example
	/*
	double y=3;
	Point2D A={-20,y};
	Point2D B={-15,y};
	Ray r=newRay(A,B,0,10);
	twoLenses(&r);

	setPoint(&r,r.c,(getRaySegmentIntersection(&r,Y)).P);
	fprintRay(OUT,&r);
	*/
	
	FILE* OUT1=fopen("oneLensX.dat","w");
	FILE* OUT2=fopen("twoLensesX.dat","w");

	double h=6.2;
	plotDeterministic1DX(OUT1,oneLens,0,h,100);
	plotDeterministic1DX(OUT2,twoLenses,0,h,100);
	//FILE* OUT=fopen("twoLensesRandomYZ.dat");
	//plotRandom2D(OUT,twoLenses,FOCAL2,0,3,20000);
	
	return 0;
}

//One thick lens
void oneLens(Ray* r){
	Circle c1={{14.25,0},14.25};
	Circle c2={{-3.07,0},-9.47};
	
	double n0=1;
	double n1=1.61765;
	
	Intersection2D i1=getRayCircleIntersection(r,c1,0);
	calculateRayRefraction(r,i1,n0,n1);
	Intersection2D i2=getRayCircleIntersection(r,c2,0);
	calculateRayRefraction(r,i2,n1,n0);
}


//Clark doublet
void twoLenses(Ray* r){
	Circle c1={{14.25,0},14.25};
	Circle c2={{-3.07,0},-9.47};
	Circle c3={{-43.98,0},-52.08};
	
	double n0=1;
	double n1=1.61765;
	double n2=1.68893;
	
	Intersection2D i1=getRayCircleIntersection(r,c1,0);
	calculateRayRefraction(r,i1,n0,n1);

	Intersection2D i2=getRayCircleIntersection(r,c2,0);
	calculateRayRefraction(r,i2,n1,n2);

	Intersection2D i3=getRayCircleIntersection(r,c3,0);
	calculateRayRefraction(r,i3,n2,n0);

}

//Utility function
double distance(Point2D p1, Point2D p2){
	double dx=p2.x-p1.x;
	double dy=p2.y-p1.y;
	
	return sqrt(dx*dx+dy*dy);
}

//Utility function, although the main use case would be to check if a Point2D is NAP
int Point2DEquals(Point2D p1, Point2D p2){
	int x=0;
	int y=0;
	x=(isnan(p1.x) && isnan(p2.x))||(isinff(p1.x)||isnan(p2.x))+(p1.x==p2.x);
	y=(isnan(p1.y) && isnan(p2.y))||(isinff(p1.y)||isnan(p2.y))+(p1.x==p2.y);
	return x & y;
}

//Utility function
void fprintPoint2D(FILE* OUT,Point2D p){
	fprintf(OUT,"%lf %lf",p.x,p.y);
}

//Utility function, allows for a convenient creation of new rays
Ray newRay(Point2D A, Point2D B, double lambda, int N){ 
	//Need at least 2 nodes, there is not point in allowing the user to
	//use just one, he can even modify them later if he needs to
	if(N<2)	N=2;
	
	//Memory allocation, we don't really deallocate them
	Point2D* data=malloc(N*2*sizeof(double));
	data[0]=A;
	data[1]=B;
	int c=1;
	Point2DArray array={data,N};
	Ray ret={array,c,lambda};
	//fprintPoint2D(stdout,getPoint(&ret,1));
	return ret;
}

//When we allocate memory it is always important to have a destructor
void destroyRay(Ray* r){
	free((*r).points.data);
}

//Utility function
Point2D getPoint(Ray* r, int n){
	return *((*r).points.data+n);	
	}

//Utility function
void addPoint(Ray* r, Point2D P){
	int n=++(*r).c;
	*((*r).points.data+n)=P;	
}

//Utility function
void setPoint(Ray* r, int n, Point2D P){
	*((*r).points.data+n)=P;
}

//Utility function, it can print to a file, but if given as FILE* stdout,
//it will print to console
void fprintRay(FILE* OUT, Ray* r){
	int i;
	for(i=0;i<(*r).points.N; i++){
		fprintPoint2D(OUT,*((*r).points.data+i));
		fprintf(OUT,"\n");
	}
}

//We calculate the intersection of a Ray with a Circle. We use the following algorithm:
//http://mathworld.wolfram.com/Circle-LineIntersection.html
Intersection2D getRayCircleIntersection(Ray* r, Circle c, int ext){
	int current=(*r).c;
	Point2D P1=getPoint(r,current-1);
	Point2D P2=getPoint(r,current);
	
	//All of this is directly copied
	double x1=P1.x-c.c.x;
	double y1=P1.y-c.c.y;
	double x2=P2.x-c.c.x;
	double y2=P2.y-c.c.y;
	double dx=x2-x1;
	double dy=y2-y1;
	double dr=sqrt(dx*dx+dy*dy);
	double D=x1*y2-x2*y1;
	
	double delta=c.r*c.r*dr*dr-D*D;
	double x,y;
	
	Point2D P;
	//If delta<0, there is no intersection
	//If delta=0, there is an intersection but it is a tangent,so it is as the Ray didn't interact
	//If delta>0, there is two intersections
	if(delta>0){
		int sgndy=(dy<0)?-1:1;
		double sdelta=sqrt(delta);
		
		x=(D*dy+sgndy*dx*sdelta)/(dr*dr)+c.c.x;
		y=(-D*dx+fabs(dy)*sdelta)/(dr*dr)+c.c.y;
		Point2D A={x,y};
		
		x=(D*dy-sgndy*dx*sdelta)/(dr*dr)+c.c.x;
		y=(-D*dx-fabs(dy)*sdelta)/(dr*dr)+c.c.y;
		Point2D B={x,y};
		
		//This is a method we developed to decide which point the choose
		double oAB=(B.x-A.x)+(B.y-A.y);
		double oRay=(P2.x-P1.x)+(P2.y-P1.y);
		if(oAB*oRay*c.r>0)
			P=A;
		else
			P=B;
		
		//If we have choosen a point behind the Ray, we return NAP (and extended Ray isn't active)
		double oF=(P.x-P1.x)+(P.y-P1.y);
		if(oRay*oF<0 && !ext)
			P=NAP;
		
	}else if(delta==0){
		double d=distance(P1,P2);
		x=P2.x+dx/d;
		y=P2.y+dy/d;
		Point2D T={x,y};
		P=T;
	}else{
		P=NAP;
	}
	//We now calculate the angle, along with a correction in case the radius is negative
	double corr=c.r/fabs(c.r);
	double alpha=atan2(corr*(P.y-c.c.y),(corr*(P.x-c.c.x)));
	
	Intersection2D ret={P,alpha};
	return ret;
}


//We calculate the intersection of a Ray with a Segment. We use the following algorithm:
//https://en.wikipedia.org/wiki/Line-line_intersection
Intersection2D getRaySegmentIntersection(Ray* r, Segment s, int ext){
	int current=(*r).c;
	Point2D P1=getPoint(r,current-1);
	Point2D P2=getPoint(r,current);
	Point2D A=s.P1;
	Point2D B=s.P2;
	
	//All of this is directly copied
	double x1=P1.x;
	double x2=P2.x;
	double x3=A.x;
	double x4=B.x;
	double y1=P1.y;
	double y2=P2.y;
	double y3=A.y;
	double y4=B.y;
	
	double dAB=distance(A,B);	
	double den=((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
	
	//This is for calculating the cross product
	double ux=x1-x3;
	double uy=y1-y3;
	double vx=x4-x3;
	double vy=y4-y3;
	
	Point2D P;
	double alpha;
	if(den==0){	
		if(ux*vy-uy*vx<0)
			P=P1;	//Contained
		else
			P=NAP;	//Parallel
	}else{
		double numx=(x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4);
		double numy=(x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4);
		Point2D T={numx/den,numy/den};
		P=T;
	}
	
	//We calculate the angle and adjust it if the orientation is reversed
	alpha=atan2(y2-y1,x2-x1)+PI/2;
	if(ux*vy-uy*vx<0) alpha+=PI;
	
	double oRay=(P2.x-P1.x)+(P2.y-P1.y);
	double oF=(P.x-P1.x)+(P.y-P1.y);
	
	//If the point is behind the Ray, we return NAP (and extended Ray isn't active)
	if(oRay*oF<0 && !ext){
		P=NAP;
		alpha=NAN;
	}

	//If the point is outside the Segmente, we return NAP
	double dAP=distance(A,P);
	double dBP=distance(B,P);
	if(dAP>dAB || dBP>dAB){
		P=NAP;
		alpha=NAN;
	}
	
	Intersection2D ret={P,alpha};
	return ret;
}




//This function receives the Ray along with the Intersection2D, and calculates
//if it refracts or reflects, adding a new point which determines the new
//trajectory.
void calculateRayRefraction(Ray* r, Intersection2D i, double n1, double n2){
	int c=(*r).c;
	
	Point2D A=getPoint(r,c);
	Point2D B=i.P;
	setPoint(r,c,B);	//We aren't interested in preserving the point that simply indicates direction
		
	//Point
	double Ax=A.x;
	double Ay=A.y;
	//Intersection
	double Bx=B.x;
	double By=B.y;
	
	double alpha=i.alpha;	//Angle between x axis and normal 
	double gamma=atan2(Ay-By,Ax-Bx);	//Angle between x axis and ray
	double e1=gamma-alpha;	//Angle with which the ray intersects the surface
	double e2;	//Angle with which the refracted (or relfected) ray intersects the surface
	
	
	//We calculate the reflection or refraction
	double delta=(n1/n2)*sin(e1);
	if(fabs(delta)<1)	//Refraction
		e2=alpha+PI+asin(delta);
	else	//Reflection
		e2=alpha-e1;
	
	//The next point of the Ray is left at an unit distance of the original
	double Px=Bx+=cos(e2);
	double Py=By+=sin(e2);
	
	Point2D P={Px,Py};
	
	addPoint(r,P);
}

//We plot the intersection of equispaced rays entering parallel to the optical axis with said axis
void plotDeterministic1DX(FILE* OUT,OpticalSystem os, double s, double e, double n){
	double div=(e-s)/n;
	double h;
	Segment X={{XAXISSTART,0},{XAXISEND,0}};
	for(h=s;h<=e;h+=div){
		
		Point2D A={XSTART,h};
		Point2D B={XEND,h};
		Ray r=newRay(A,B,LAMBDA,POINTS);
		os(&r);
		Intersection2D I=getRaySegmentIntersection(&r,X,1);
		
		fprintf(OUT,"%lf\n",I.P.x);
		
		destroyRay(&r);
	}
}

//We plot the intersection of equispaced rays entering parallel to the
//optical axis with a vertical screen placed at an specified x coordinate
void plotDeterministic1DY(FILE* OUT,OpticalSystem os,double x,double s,double e,double n){
	double div=(e-s)/n;
	double h;
	Segment Y={{x,YAXISSTART},{x,YAXISEND}};
	for(h=s;h<=e;h+=div){
		
		Point2D A={XSTART,h};
		Point2D B={XEND,h};
		Ray r=newRay(A,B,LAMBDA,POINTS);
		os(&r);
		Intersection2D I=getRaySegmentIntersection(&r,Y,0);
		
		fprintf(OUT,"%lf\n",I.P.y);
					
		destroyRay(&r);
	}
}

//We plot the simulated 3D intersection of equispaced rays entering parallel
//to the optical axis with a vertical screen placed at an specified x coordinate
void plotDeterministic2DYZ(FILE* OUT,OpticalSystem os, double x, double r, double ny, double nz){
	double divy=2*r/ny;
	double divz=2*r/nz;
	Segment Y={{x,YAXISSTART},{x,YAXISEND}};
	double y,z;
	for(y=-r;y<=r;y+=divy){
		for(z=-r;z<=r;z+=divz){
			
			double d=sqrt(x*x+y*y);
			if(d>r)
				continue;
			
			Point2D A={XSTART,d};
			Point2D B={XEND,d};
			Ray r=newRay(A,B,LAMBDA,POINTS);
			os(&r);
			Intersection2D I=getRaySegmentIntersection(&r,Y,0);		
			double d2=I.P.y;
			x*=d2/d;
			y*=d2/d;
			Point2D P={x,y};
			
			fprintPoint2D(OUT,P);
			fprintf(OUT,"\n");
			
			destroyRay(&r);
		}
	}
}

//We plot the simulated 3D intersection of randomly placed rays entering
//parallel to the optical axis with a vertical screen placed at an specified
//x coordinate
void plotRandom2DYZ(FILE* OUT,OpticalSystem os,double x, double r,double n){
	srand48(seed()); // Set the seed
	Segment Y={{x,YAXISSTART},{x,YAXISEND}};
	int i;
	for(i=0;i<n;i++){
		double x=scalerandom(drand48(),-r,r);
		double y=scalerandom(drand48(),-r,r);
		
		double d=sqrt(x*x+y*y);
		if(d>r){
			i--;
			continue;
		}
		
		Point2D A={XSTART,d};
		Point2D B={XEND,d};
		Ray r=newRay(A,B,LAMBDA,POINTS);
		os(&r);
		Intersection2D I=getRaySegmentIntersection(&r,Y,0);
		double d2=I.P.y;
		x*=d2/d;
		y*=d2/d;
		Point2D P={x,y};
		
		fprintPoint2D(OUT,P);
		fprintf(OUT,"\n");
		
		destroyRay(&r);
	}
}



