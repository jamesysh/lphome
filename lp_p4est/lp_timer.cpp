#include "lp_timer.h"
#include <iostream>
LPTimer::LPTimer(){

     upwindtimer = 0;
     lwtimer = 0;
     octreetimer = 0;
     neighboursearchertimer = 0;
     heatingmodeltimer = 0;
     boundaryconditiontimer = 0;
     particleimmigrationtimer = 0;
     ghostparticletimer = 0;
}
void LPTimer::updateTimer(double* timer, double dt){

    *timer += dt;

}


void LPTimer::writeTimeAnalysis( int mpirank, double totalt){
    if(mpirank != 0)
        return;
    static double  totalt_old = 0;
    
    printf("---------------------------------TIME ANALYSIS--------------------------------\n");
    printf("%-45s %15s %.2f second %10s %.2f second.\n","Time For Octree Construction:","Current Timestep:",octreedt,"Total:", octreetimer);
    printf("%-45s %15s %.2f second %10s %.2f second.\n","Time For Boundary Condition:","Current Timestep:",boundaryconditiondt,"Total:", boundaryconditiontimer);
    printf("%-45s %15s %.2f second %10s %.2f second.\n","Time For Neighbour Searcher:","Current Timestep:",neighboursearcherdt,"Total:", neighboursearchertimer);
    printf("%-45s %15s %.2f second %10s %.2f second.\n","Time For Heating Model:","Current Timestep:",heatingmodeldt,"Total:", heatingmodeltimer);
    printf("%-45s %15s %.2f second %10s %.2f second.\n","Time For Ghost Particle Construction:","Current Timestep:",ghostparticledt,"Total:", ghostparticletimer);
    printf("%-45s %15s %.2f second %10s %.2f second.\n","Time For Upwind Solver:","Current Timestep:",upwinddt,"Total:", upwindtimer);
    printf("%-45s %15s %.2f second %10s %.2f second.\n","Time For Lax-Wendroff Solver:","Current Timestep:",lwdt,"Total:", lwtimer);
    printf("Runging Time For Current Step is %.2f second.  The Total Running Time is %.2f s.\n",totalt-totalt_old,totalt);
    printf("------------------------------------------------------------------------------\n");
    totalt_old = totalt;
}
