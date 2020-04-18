#ifndef __GEOMETRY_DISK_H
#define __GEOMETRY_DISK_H

#include "geometry.h"
class Disk: public Geometry{

    public:
        Disk();
        virtual ~Disk(){}
        virtual bool operator()(double x, double y, double z) const;
  	    virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
    
    private:
        double xcen;
        double ycen;
        double radius;
};











#endif
