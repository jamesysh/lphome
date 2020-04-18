#ifndef __GEOMETRY_PELLET_H
#define __GEOMETRY_PELLET_H

#include "geometry.h"
#include <math.h>
#include "pellet_solver.h"
class PelletLayer: public Geometry {
public:
	PelletLayer();
	virtual ~PelletLayer() {}
	virtual bool operator()(double x, double y, double z) const;
  	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
	PelletSolver *pellet_solver;		
};



#endif
