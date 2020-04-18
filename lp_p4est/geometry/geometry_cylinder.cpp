#include "geometry_cylinder.h"
#include <iostream>
#include <cmath>


Cylinder::Cylinder():radius(5),length(7),xCen(0),yCen(0),zCen(0){}


bool Cylinder::operator()(double x, double y, double z) const {	
    
    if(z <= (zCen+length/2.) && z >= (zCen-length/2.) && (y*y+x*x)<radius*radius )

        return true;
        
    else

        return false;

}

void Cylinder::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-radius;
	xmax = xCen+radius;
	ymin = yCen-radius;
	ymax = yCen+radius;
	zmin = zCen-length/2;
	zmax = zCen+length/2;
}


