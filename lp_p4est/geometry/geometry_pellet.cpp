#include "geometry_pellet.h"

PelletLayer::PelletLayer(){
    
}


bool PelletLayer::operator()(double x, double y, double z) const{
	
	pellet_info_t* pellet;
	size_t pi,pnum = pellet_solver->pelletlist->elem_count;
	double xcen, ycen, zcen;
	double r,r0,r1;
	for(pi = 0;pi<pnum;pi++){
	    pellet = (pellet_info_t*)sc_array_index(pellet_solver->pelletlist,pi);
	    xcen = pellet->xyz[0];	
	    ycen = pellet->xyz[1];	
	    zcen = pellet->xyz[2];
            r0 = pellet->innerradius;
	    r1 = pellet->outerradius;
	    r = sqrt((x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen));
	    if(r<=r1 && r>=r0)
		return true;
	}
	return false;
}


void PelletLayer::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax){
        xmin = ymin = zmin = xmax = ymax = zmax = -1;
	pellet_info_t* pellet;
	size_t pi,pnum = pellet_solver->pelletlist->elem_count;
	double xcen, ycen, zcen;
	double r;
    	pellet = (pellet_info_t*)sc_array_index(pellet_solver->pelletlist,0);
    	xcen = pellet->xyz[0];	
    	ycen = pellet->xyz[1];	
    	zcen = pellet->xyz[2];
    	r = pellet->outerradius;
	xmin = xcen - r;
	ymin = ycen -r;
	zmin = zmin - r;
	xmax = xcen + r;
	ymax = ymax + r;
	zmax = zmax + r;
	for(pi = 1;pi<pnum;pi++){
	    pellet = (pellet_info_t*)sc_array_index(pellet_solver->pelletlist,pi);
	    xcen = pellet->xyz[0];	
	    ycen = pellet->xyz[1];	
	    zcen = pellet->xyz[2];
	    r = pellet->outerradius;
	    xmin = min(xmin,xcen-r);
	    ymin = min(ymin,ycen-r);
	    zmin = min(zmin,zcen-r);
	    xmax = max(xmax,xcen+r);
	    ymax = max(ymax,ycen+r);
	    zmax = max(zmax,zcen+r);
	}

}

