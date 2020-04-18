#ifndef __GEOMETRY_CYLINDER_H__
#define __GEOMETRY_CYLINDER_H__

#include "geometry.h"

#include <math.h>


class Cylinder :public Geometry{
    
  public:
    Cylinder();

   virtual ~Cylinder() {}
   
   virtual bool operator()(double x, double y, double z) const;
    
   virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
   
   private:
	double radius;
	double length;
	double xCen;
	double yCen;
	double zCen;
   
   };

#endif
