#ifndef __STATE_PELLET_H__
#define __STATE_PELLET_H__

#include "state.h"
#include <math.h>

class PelletState: public State{
public:
	PelletState();
	virtual ~PelletState(){};
	virtual double pressure(double x, double y, double z);
	virtual double density(double x, double y, double z);
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
	double m_fDen;
	double m_fVel;
	double m_fPressure;
	double xcen;
	double ycen;
	double zcen;
};
#endif
