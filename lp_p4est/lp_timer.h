#ifndef __LP_TIMER__
#define __LP_TIMER__

class LPTimer{

    public:
    double upwindtimer, upwinddt;
    double lwtimer, lwdt;
    double octreetimer, octreedt;
    double neighboursearchertimer, neighboursearcherdt;
    double heatingmodeltimer, heatingmodeldt;
    double boundaryconditiontimer, boundaryconditiondt;
    double particleimmigrationtimer, particleimmigrationdt;
    double ghostparticletimer, ghostparticledt;
    
    LPTimer();
    void updateTimer(double *timer, double dt);
    void writeTimeAnalysis(int mpirank, double totalt);

};







#endif
