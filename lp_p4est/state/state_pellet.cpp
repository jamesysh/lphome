#include "state_pellet.h"
#include <iostream>

PelletState::PelletState():
m_fDen(1./100), m_fVel(0), m_fPressure(16), xcen(0), ycen(0), zcen(0){}

double PelletState::pressure(double x, double y, double z) {
	return m_fPressure;
}

double PelletState::density(double x, double y, double z){
	return m_fDen;
}

void PelletState::velocity(double x, double y, double z, double& vX, double& vY,double& vZ){
	double dx=x-xcen;
	double dy=y-ycen;
	double dz=z-zcen;
	double r=sqrt(dx*dx+dy*dy+dz*dz);
	vX=dx/r*m_fVel;
	vY=dy/r*m_fVel;
	vZ=dz/r*m_fVel;
    
}
