#ifndef __LP_SOLVER_H__
#define __LP_SOLVER_H__
#include "initializer.h"
#include "particle_data.h"
#include "lp_timer.h"
#include <vector>
#include "ls_solver.h"
#include "octree_manager.h"
#include "particle_viewer.h"
#include "pellet_solver.h"
class LPSolver {

public:
	
    LPSolver(Initializer* init, Global_Data *g, Octree_Manager *octree, ParticleViewer *viewer);
    
    void solve_2d(); 
    void solve_3d();
    ~LPSolver(){}
    
    void moveParticle();
    void solve_upwind(int phase);
    void solve_laxwendroff(); 
    void setInAndOutPointer(pdata_t *pad, double **inpressure, double **outpressure, double **involume, double **outvolume,
        double** invelocity, double **outvelocity, double **insoundspeed, double **outsoundspeed, int dir, int phase);

    void setNeighbourListPointer(pdata_t *pad, sc_array_t** neilist0, sc_array_t **neilist1, int dir );
    
    void computeSpatialDer(int dir,pdata_t *pad, sc_array_t *neighbourlist, const double* inpressure, const double *invelocity,
        double *vel_d, double *vel_dd, double *p_d, double *p_dd);

    void computeSpatialDer(pdata_t *pad, const double* inPressure,  const double* inVolume,
           const double* inVelocityU, const double* inVelocityV, const double* inVelocityW,
             double* Pd, double *volumed, double* Ud, double* Vd, double* Wd);

    void updateLocalSpacing();
    void timeIntegration(
	double realDt, int dir,
	double gravity, double inVolume, double inVelocity, double inPressure, double inSoundSpeed, 
	double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
	double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
	double* outVolume, double* outVelocity, double* outPressure);

    void timeIntegration( double gravity, double inVolume, double inVelocityU, double inVelocityV, double inVelocityW, double inPressure, double inSoundSpeed,
                                      double* Volumed, double* Ud, double* Vd, double *Wd, double *Pd,
                                                        double* outVolume, double* outVelocityU, double* outVelocityV, double* outVelocityW, double* outPressure);
    void computeA3D(double *A ,pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, double distance);
    void computeA2D(double *A ,pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, double distance);


    void computeB3d(double *B, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, const double* indata, indata_t datatype, int dir);
    void computeB2d(double *B, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, const double* indata, indata_t datatype, int dir);

    void computeCFLCondition();
    
    void computeLocalBoundaryAndFluidNum();
    
    void computeTemperature();

    bool adjustDtByWriteTimeInterval();
    private:
    LPTimer *lp_timer; 
    Octree_Manager *octree;
    Global_Data * gdata; 
    ParticleViewer *viewer;
    PelletSolver *pellet_solver;
    int LPFOrder;
    double tstart;
    double tend;
    double currenttime;
    double cfldt;
    double cflcoefficient;
	std::vector<std::vector<int> > m_vDirSplitTable; 
    int splitorder;
    int totalphase;
    int writestep;
    int timestep;
    double writetimeinterval;
    double nextwritetime;
        
};


#endif
