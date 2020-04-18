#include "geometry_disk.h"


Disk::Disk(){
    xcen = 0;
    ycen = 0;
    radius = 5;

}


bool Disk::operator()(double x, double y, double z) const{
    double r = (x-xcen)*(x-xcen) + (y-ycen)*(y-ycen);
    return r<radius*radius;

}

void Disk::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax){
    xmin = xcen - radius;
    xmax = xcen + radius;
    ymin = ycen - radius;
    ymax = ycen + radius;

}


