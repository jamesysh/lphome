#include "boundary_pellet.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include "sc.h"
#include "particle_data.h"
#include "pellet_solver.h"
using namespace std;

PelletInflowBoundary::PelletInflowBoundary():Pinflow(16),Uinflow(0),Vinflow(100){ massflowrate = 0;}

void PelletInflowBoundary::generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx, double dt){
    static bool FIRST = true;
    PelletSolver *p = g->pellet_solver;
    pdata_t*  pad;
    size_t li, lpnum = g->particle_data->elem_count; 
    double qsum_bc = p->qsum_bc;
    double vol_bc = p->vol_bc;
    double pres_bc = p->pres_bc;
    double u_bc = p->u_bc;
    double ss_bc = p->ss_bc;
    double r_bc = p->r_bc;
    double B,C;

    static double mass_fix = dx*dx*dx/Vinflow/sqrt(2);
    double pr = p->pelletradius;
    double x, y, z, d_x, d_y, d_z, dr;
    
    double Ts = 400;
    double gamma = g->gamma;
    double gamma1 = gamma - 1;
    double mu = p->mu;
    double R = 83.1446/mu;
    computeMassFlowRate(g,dx/2);
    
    if(massflowrate == 0){
	B = 0;
	C = 0;
    }
    else{
    	B = (pres_bc+dt*gamma1*(qsum_bc))*vol_bc/ss_bc - u_bc;
    	C = -(massflowrate/4/M_PI/pr/pr)*R*Ts*vol_bc/ss_bc;	
    }
  // double pv  = Vinflow*massflowrate/4/M_PI/pr/pr;
   /* 
    double pv = sqrt(gamma*R*Ts)/2.;
    if(massflowrate != 0){
        Vinflow = 4*M_PI*pr*pr*pv/massflowrate;
        Pinflow = R*Ts/Vinflow;
    }
    */
    
    double pv = (-B+sqrt(B*B-4*C))/2;  //Riemann wave curve

    if(!FIRST && massflowrate!=0 ){
        Vinflow = 4*M_PI*pr*pr*pv/massflowrate;
        Pinflow = R*Ts/Vinflow;
    }

    if(g->mpirank==0){
    	cout<<"pressure on boundary = "<<Pinflow<<endl;
        cout<<"volume on boundary = "<<Vinflow<<endl;
        cout<<"pellet velocity on boundary = "<<pv<<endl;
    	cout<<"sound speed on boundary = "<<g->eos->getSoundSpeed(Pinflow,1./Vinflow)<<endl;
    }
    pelletvelocity = pv;
    
    double xcen = 0;
    double ycen = 0;
    double zcen = 0;
    

   for(li = 0;li<lpnum;li++ ){
       pad = (pdata_t*) sc_array_index(g->particle_data,li);
       x = pad->xyz[0];
       y = pad->xyz[1];
       z = pad->xyz[2];
       d_x = x-xcen;
       d_y = y-ycen;
       d_z = z-zcen;
       dr = sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
       if(dr<pr){
           pad->ifboundary = true;
           pad->flagdelete = g->flagdelete;
       }
   } 
   
    MPI_Barrier(g->mpicomm); 

    ux = px = 0;
    double actualdx = cbrt(Vinflow*mass_fix)*1.2;
    dx = actualdx;
    if(g->mpirank == 0)
        cout<<dx<<endl;

    
    int n = 4.0*M_PI*pr*pr*pr/mass_fix/Vinflow;
    double newpir =0;// pr*4/5;
    
    
    srand(time(NULL)); 
     
    g->gpnum += n;
    for(int i=0;i<n;i++)
    if(g->mpirank == 0){	
    {

        double tx=1,ty=1,tz=1,tr;
		
        while(tx*tx+ty*ty+tz*tz>1)
		
        {
				tx=2*(double)rand()/(double)RAND_MAX-1;
				ty=2*(double)rand()/(double)RAND_MAX-1;
				tz=2*(double)rand()/(double)RAND_MAX-1;
		
    }
        tr=sqrt(tx*tx+ty*ty+tz*tz);
        tx=tx/tr,ty=ty/tr,tz=tz/tr;
        tr=newpir+(double)rand()/(double)RAND_MAX*(pr-newpir);
        x = tx*tr;
        y = ty*tr;
        z = tz*tr;
        d_x = tx*tr;
        d_y = ty*tr;
        d_z = tz*tr;
        dr = sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
        //assert(dr<=pr);
    
   	    double newpv = pv + (dr-pr)*ux;    
   	       
        pad = (pdata_t*)sc_array_push_count(g->particle_data,1);
        pad->v[0] = newpv*d_x/dr;
        pad->v[1] = newpv*d_y/dr;
        pad->v[2] = newpv*d_z/dr;
        
        pad->xyz[0] = x;
        pad->xyz[1] = y;
        pad->xyz[2] = z;

        pad->volume = Vinflow;
        pad->pressure = Pinflow + (dr-pr)*px;
        pad->localspacing = dx;
        pad->mass = mass_fix;//100000000;
        pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure,1./pad->volume);
        if(dr<pr/2){ 
            pad->ifboundary = true;
            pad->flagdelete = !g->flagdelete;
        }
        else{
            pad->ifboundary = false;
            pad->redocount = 0;
            pad->driftv = 0;
            pad->layer = 0;
        }
        }

    }



    FIRST = false;

}
void PelletInflowBoundary::UpdateInflowBoundary(Global_Data* g, EOS* m_pEOS, double dt, double dx){
    
    double massflowrate = 0.6;

    double mass_fix = dx*dx*dx/Vinflow/sqrt(2);

    PelletSolver *p = g->pellet_solver;
    int numberofNewFluid = massflowrate*dt/mass_fix;///10;
    g->gpnum += numberofNewFluid;
    pdata_t*  pad;

     int rnd = rand();

     double pr = p->pelletradius;
     double pv  = Vinflow*massflowrate/4/M_PI/pr/pr;
     double x,y,z;
     double dis = dx/5;
     double offset = 2./numberofNewFluid;
     double increment = M_PI*(3-sqrt(5));
     for(int i=0;i<numberofNewFluid;i++){
           if(g->mpirank == i%g->mpisize){
           double r_random = ((double)rand()/(double)(RAND_MAX))*dis; 
           double  y_tmp =  ((i*offset-1)+offset/2);
           pad = (pdata_t*)sc_array_push_count(g->particle_data,1);
           y = y_tmp * (pr+r_random) + 0;
           
           pad->xyz[1] = y;
           double r = sqrt(1-y_tmp*y_tmp);
           double phi = ((i+rnd)%numberofNewFluid) * increment;

           x = cos(phi)*r*(pr+r_random) + 0;
           
           pad->xyz[0] = x;
           z = sin(phi)*r*(pr+r_random) + 0;
           
           pad->xyz[2] = z;

           double d_x = x;
           double d_y = y;
           double d_z = z;
           double dr = sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
           pad->v[0] = pv*d_x/dr;
           pad->v[1] = pv*d_y/dr;
           pad->v[2] = pv*d_z/dr;
           pad->volume = Vinflow;
           pad->pressure = Pinflow;
           pad->localspacing = dx;
           pad->mass = mass_fix;
           pad->ifboundary = false; 
           pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure,1./pad->volume);
           }
           
           }
     

}


void PelletInflowBoundary::computeMassFlowRate(Global_Data *g,double dx){

    static bool first=true;
   if(first){
       massflowrate = 0;
       first = false;
       return;
       }
   PelletSolver *p = g->pellet_solver;
   pdata_t *pad;
   size_t li, lpnum = g->particle_data->elem_count;
   double x, y, z, dr;
   double pr = p->pelletradius;
   double xcen=0;
   double ycen=0;
   double zcen=0;
   int counter_nei = 0;
   double qsum = 0;
   int counter_g = 0;
   double qsum_g = 0;
   double sublimationenergy = p->sublimationenergy;
   for(li = 0; li<lpnum; li++){
       pad = (pdata_t *)sc_array_index(g->particle_data,li);
       if(pad->ifboundary)
           continue;
       x = pad->xyz[0];
       y = pad->xyz[1];
       z = pad->xyz[2];
       dr = sqrt((x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen))-pr;
       if(dr<dx && dr>0){
            qsum += pad->qplusminus*fabs(x/sqrt(x*x+y*y+z*z));
            counter_nei ++; 
            } 
       }
       
       MPI_Barrier(g->mpicomm);
       
       MPI_Allreduce (&counter_nei, &counter_g, 1, MPI_INT,
                                MPI_SUM, g->mpicomm);
       if(counter_g == 0){
           assert(false);
           }
       
       MPI_Allreduce(&qsum, &qsum_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       
       massflowrate = qsum_g/counter_g*4*M_PI*pr*pr/sublimationenergy;
//       if(g->pellet_solver->heatingmodel == 1)
//           massflowrate *= 2./M_PI;

       P4EST_GLOBAL_ESSENTIALF("The massflowrate is %.16g. Found %d neighbouring particles.\n", massflowrate, counter_g);
}

void PelletInflowBoundary::computeRadialDerivative(Global_Data *g,double dx){

   pdata_t *pad;
   size_t li, lpnum = g->particle_data->elem_count;
   PelletSolver *p = g->pellet_solver;
   double x, y, z, dr;
   double pr = p->pelletradius;
   double xcen=0;
   double ycen=0;
   double zcen=0;
   int counter_nei = 0;
   double psum = 0;
   double usum = 0;
   int counter_g = 0;
   double psum_g = 0;
   double usum_g = 0;
   double vx, vy, vz;
   for(li = 0; li<lpnum; li++){
       pad = (pdata_t *)sc_array_index(g->particle_data,li);
       if(pad->ifboundary)
           continue;
       x = pad->xyz[0];
       y = pad->xyz[1];
       z = pad->xyz[2];
       dr = sqrt((x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen))-pr;
       if(dr>dx && dr < 2*dx){
            vx = pad->v[0];
            vy = pad->v[1];
            vz = pad->v[2];
            psum += pad->pressure;
            usum += (vx*(x-xcen)+vy*(y-ycen)+vz*(z-zcen))/(dr+pr); 
            counter_nei ++; 
           } 
       }
       
       MPI_Barrier(g->mpicomm);
       
       MPI_Allreduce (&counter_nei, &counter_g, 1, MPI_INT,
                                MPI_SUM, g->mpicomm);
       if(counter_g == 0){
           ux = 0;
           px = 0;
           return;
           }

       avg_dis = 1.5*dx;
       MPI_Allreduce(&psum, &psum_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       
       px = (psum_g/counter_g-Pinflow)/avg_dis;

       MPI_Allreduce(&usum, &usum_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       
       ux = (usum_g/counter_g-pelletvelocity)/avg_dis;
       if(px > 0 )
           px = 0;
       if(ux < 0)
           ux = 0;
}





PelletOutflowBoundary::PelletOutflowBoundary():thickness(2.){}


void PelletOutflowBoundary::generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx, double dt){
   size_t li, lpnum = g->particle_data->elem_count;
   double dl = g->domain_len;
   double bl = dl/2.-thickness;
   double x, y, z;
   pdata_t *pad;
   for(li = 0; li<lpnum; li++){
        pad = (pdata_t *) sc_array_index(g->particle_data, li);    
        x = pad->xyz[0];
        y = pad->xyz[1];
        z = pad->xyz[2];
        if(pad->ifboundary){
           if(fabs(x)>bl || fabs(y)>bl || fabs(z)>bl){
           
               pad->flagdelete = !g->flagdelete;
               pad->xyz[0] += dt*pad->v[0];
               pad->xyz[1] += dt*pad->v[1];
               pad->xyz[2] += dt*pad->v[2];
           
           }
        }
        else{
           if(fabs(x)>bl || fabs(y)>bl || fabs(z)>bl){
               pad->ifboundary = true;
               pad->flagdelete = !g->flagdelete;
           }
        }
   }

}

int PelletInflowBoundary::generateInflowCoordinate(vector<double>&xa,vector<double>&ya,vector<double>&za,double dx,Global_Data* g){
    srand(time(NULL));
    PelletSolver *p = g->pellet_solver;
    double pr = p->pelletradius;
    double x, y, z;
    double x0 = -pr/2 + dx/2;;
    double y0 = -pr/2 +dx/2;
    double z0 = -pr/2 + dx/2;
    int num = 2*pr/2/dx+1;
    int counter = 0;
    for(int i=0;i<num;i++){
        for(int j=0;j<num;j++){
            for(int k=0;k<num;k++){
                x = x0+i*dx + (2*(rand()/(double)RAND_MAX)-1)*dx/2;
                y = y0+j*dx + (2*(rand()/(double)RAND_MAX)-1)*dx/2;
                z = z0+k*dx + (2*(rand()/(double)RAND_MAX)-1)*dx/2;
                if((x*x+y*y+z*z)<pr*pr){
                    xa.push_back(x);
                    ya.push_back(y);
                    za.push_back(z);
                    counter++;
                }

            }
        }
    
    }
    return counter;
}






