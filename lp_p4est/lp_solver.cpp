#include <iostream>
#include "lp_solver.h"

using namespace std;



LPSolver::LPSolver(Initializer *init, Global_Data *g, Octree_Manager *o, ParticleViewer *v){
    gdata = g;
    octree = o;
    viewer = v;
    lp_timer = new LPTimer();
    splitorder = 0;
    cfldt = 0;
    cflcoefficient = init->getCFLCoeff();
    tstart = init->getStartTime();
    tend = init->getEndTime();
    LPFOrder = init->getLPFOrder();
    writestep = init->getWriteStep();
    writetimeinterval = init->getWriteTimeInterval();
    currenttime = tstart;
    timestep = 0;
    nextwritetime = writetimeinterval + tstart ;
    if(gdata->dimension == 3){
        m_vDirSplitTable = vector<vector<int> >
        ({{0,1,2},
          {0,2,1},
          {1,0,2},
          {1,2,0},
          {2,0,1},
          {2,1,0}});
        totalphase = 3;
    }
    else if(gdata->dimension == 2){
        m_vDirSplitTable = vector<vector<int> >
        ({{0,1},
          {1,0}});
        totalphase = 2;
    }
    
    if(init->getPelletDistribution()){
        pellet_solver = new PelletSolver(init,gdata);
    }

}

void LPSolver::moveParticle(){
    pdata_t *pad;
    double dt = cfldt;
    size_t li, lpnum = gdata->particle_data->elem_count;
     
    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
       pad->xyz[0] += 0.5*dt*(pad->oldv[0]+pad->v[0]);
       pad->xyz[1] += 0.5*dt*(pad->oldv[1]+pad->v[1]);
       if(gdata->dimension == 3)
           pad->xyz[2] += 0.5*dt*(pad->oldv[2]+pad->v[2]);
       /* 
       if(pad->xyz[0] == 0){
           printf("warningX %f %f %f %d %d\n.",x,pad->oldv[0],pad->v[0],gdata->mpirank,pad->schemeorder);
       }

       if(pad->xyz[1] == 0){
           printf("warningY %f %f %f %d %d\n.",y, pad->oldv[1],pad->v[1],gdata->mpirank, pad->schemeorder);
       }
       if(pad->xyz[2] == 0){
           printf("warningZ %f %f %f %d %d\n.",z,pad->oldv[2],pad->v[2],gdata->mpirank,pad->schemeorder);
       }*/
    }
}

void LPSolver::updateLocalSpacing(){

    pdata_t *pad;
    size_t li, lpnum = gdata->particle_data->elem_count;
         
    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
    
       if(gdata->dimension == 3) 
           pad->localspacing *= cbrt(pad->volume/pad->volumeT1); 
       else if(gdata->dimension == 2)
           pad->localspacing *= sqrt(pad->volume/pad->volumeT2);
      // pad->localspacing = cbrt(pad->mass/pad->sphdensity);
    }

}

//dir: 0 x, 1 y, 2 z
void LPSolver::solve_upwind(int phase){
    
    const int dir = m_vDirSplitTable[splitorder][phase];
    double realdt = gdata->dimension == 3? cfldt/3 : cfldt/2;
    sc_array_t *neighbourlist0;
    sc_array_t *neighbourlist1;
    double *inpressure, *outpressure;
    double *involume, *outvolume;
    double *invelocity, *outvelocity;
    double *insoundspeed, *outsoundspeed;
    double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output
  


  pdata_t * pad;

  size_t li, lpnum = gdata->particle_data->elem_count;
     
  for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
         bool redo = false;
         
         setInAndOutPointer(pad, &inpressure, &outpressure, &involume, &outvolume,
                  &invelocity, &outvelocity, &insoundspeed, &outsoundspeed, dir, phase);
      
         setNeighbourListPointer(pad, &neighbourlist0, &neighbourlist1,dir);
         
         if(*insoundspeed == 0 || *involume == 0){
             pad->schemeorder = 0;
            *outvolume = *involume;
            *outpressure = *inpressure;
            *outvelocity = *invelocity;
            *outsoundspeed = *insoundspeed;
            printf("Detect a particle which has 0 volume or 0 soundspeed!!!.\n"); 
         }
            
         pad->schemeorder = 1;
        
         computeSpatialDer(dir, pad, neighbourlist0, inpressure, invelocity,
        &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0);
     
        computeSpatialDer(dir, pad, neighbourlist1, inpressure, invelocity,
        &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1);
        
        /* 
*/
    timeIntegration( realdt, dir,
	 0,  *involume, *invelocity, *inpressure, *insoundspeed, 
	vel_d_0, vel_dd_0, p_d_0, p_dd_0,
    vel_d_1, vel_dd_1, p_d_1, p_dd_1,
    outvolume, outvelocity, outpressure);
    if(gdata->pelletnumber){
        *outpressure += realdt*(pad->deltaq)*((*insoundspeed)*(*insoundspeed)/(*involume)/(*inpressure)-1);} 
         if(*outpressure < gdata->invalidpressure || 1./(*outvolume)< gdata->invaliddensity)
         {
             /*
                printf("inpressure %.8g outpressure %.8g involume %.8g outvolume %.8g %d\n",*inpressure, *outpressure, *involume, *outvolume ,dir);
                printf("%f %f %f \n",pad->xyz[0],pad->xyz[1],pad->xyz[2]);
            
                pdata_copy_t *padnei;                
                
                
                for(int i=0;i<8;i++){
                    gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist0,i);
                    
                    printf("%f %f %f %f %f \n",padnei->xyz[0],padnei->xyz[1],padnei->xyz[2], padnei->pressure, padnei->volume);
                }
                for(int i=0;i<8;i++){
                    gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist1,i);
                    
                    printf("%f %f %f %f %f \n",padnei->xyz[0],padnei->xyz[1],padnei->xyz[2], padnei->pressure, padnei->volume);
                }
                */
          //  exit(0);
             redo  = true;
         }
        if(isnan(*outvolume) || isnan(*outvelocity) || isnan(*outpressure)) 
        {
                redo = true;}

        if(isinf(*outvolume) || isinf(*outvelocity) || isinf(*outpressure)) 
        {
                redo = true;}


        if(redo == true){

            pad->redocount++;
            
            if(pad->redocount >4){
                pad->schemeorder = 0;
            }
            else{
                li--;
                cout<<"redo upwind for a particle"<<endl;
            }

            
        }
        else{
        
        
            *outsoundspeed = gdata->eos->getSoundSpeed(*outpressure, 1./(*outvolume));
            pad->redocount = 0;        
        }

        if(pad->schemeorder == 0)
        {
            pad->redocount = 0;
            *outvolume = *involume;
            *outpressure = *inpressure;
            *outvelocity = *invelocity;
            *outsoundspeed = *insoundspeed;
        }
        

      
      
      
      }

} //to do

void LPSolver::setInAndOutPointer(pdata_t *pad, double **inpressure, double **outpressure, double **involume, double **outvolume,
        double** invelocity, double **outvelocity, double **insoundspeed, double **outsoundspeed, int dir, int phase){

     if(phase == 0){
        *inpressure = &pad->pressure;
        *outpressure = &pad->pressureT1;
        *involume = &pad->volume;
        *outvolume = &pad->volumeT1;
        *insoundspeed = &pad->soundspeed;
        *outsoundspeed = &pad->soundspeedT1;
     }

     if(phase == 1){
        *inpressure = &pad->pressureT1;
        *outpressure = &pad->pressureT2;
        *involume = &pad->volumeT1;
        *outvolume = &pad->volumeT2;
        *insoundspeed = &pad->soundspeedT1;
        *outsoundspeed = &pad->soundspeedT2;
     }

     if(phase == 2){
        *inpressure = &pad->pressureT2;
        *outpressure = &pad->pressureT1;
        *involume = &pad->volumeT2;
        *outvolume = &pad->volumeT1;
        *insoundspeed = &pad->soundspeedT2;
        *outsoundspeed = &pad->soundspeedT1;
     }
    
     if(dir == 0){
        *invelocity = &pad->v[0];
        *outvelocity = &pad->oldv[0];
     
     }
     if(dir == 1){
     
        *invelocity = &pad->v[1];
        *outvelocity = &pad->oldv[1];
     }

     if(dir == 2){
     
        *invelocity = &pad->v[2];
        *outvelocity = &pad->oldv[2];
     }

}

void LPSolver::setNeighbourListPointer(pdata_t *pad, sc_array_t** neilist0, sc_array_t **neilist1,int dir){
    if(dir == 0){
        *neilist0 = pad->neighbourfrontparticle;
        *neilist1 = pad->neighbourbackparticle;
    }
    else if(dir == 1){
    
        *neilist0 = pad->neighbourrightparticle;
        *neilist1 = pad->neighbourleftparticle;
    }
    else if(dir == 2){
    
        *neilist0 = pad->neighbourupparticle;
        *neilist1 = pad->neighbourdownparticle;
    }
}



void LPSolver::computeSpatialDer(int dir,pdata_t *pad, sc_array_t *neighbourlist, const double* inpressure, const double *invelocity,
        double *vel_d, double *vel_dd, double *p_d, double *p_dd) {
    size_t numrow = pad->redocount + gdata->numrow1st;
    size_t numcol = gdata->numcol1st;
    double distance;
    int info;
    int offset = gdata->dimension == 3? 3:2; // 2 for 2 dimension
    neighbour_info_t *ninfo;
    while( true ){
        if(numrow > neighbourlist->elem_count || numrow > 16){
            *vel_d = 0;
            *vel_dd = 0;
            *p_d = 0;
            *p_dd = 0;
            break;
        }
        ninfo = (neighbour_info_t *)sc_array_index(neighbourlist,0);
        distance = ninfo->distance;
        double A[numrow*numcol];
        
        if(gdata->dimension == 3)
            computeA3D(A, pad, neighbourlist,numrow,distance);
        else if(gdata->dimension == 2)
            computeA2D(A, pad, neighbourlist,numrow,distance);
        
        double B[numrow];
    
        if(gdata->dimension == 3)
            computeB3d(B, pad, neighbourlist, numrow, inpressure ,PRESSURE, dir);
        else if(gdata->dimension == 2)
            computeB2d(B, pad, neighbourlist, numrow, inpressure ,PRESSURE, dir);
        QRSolver qrSolver(numrow,numcol,A);
		
		double result[numcol];
	    info = qrSolver.solve(result,B);
        
        if(info != 0){
            numrow ++;
            printf("Warning: Rank Deficient For a Particle in UPWIND!\n");
        } 
        else{
            
            *p_d = result[dir]/distance; //dir=0 (x), dir=1(y), dir=2(z)
            *p_dd = pad->schemeorder == 2?  result[dir+offset]/distance/distance:0;
           /* 
                cout<<"newone"<<*p_d*distance<<endl;
            if(*p_d > 0){
                cout<<"newone"<<*p_d<<endl;
                
                for(int i=0;i<numrow;i++){
                cout<<B[i]<<endl;
                }
            
            }
            */
             
        if(gdata->dimension == 3)
            computeB3d(B, pad, neighbourlist, numrow, invelocity ,VELOCITY, dir);
        else if(gdata->dimension == 2)
            computeB2d(B, pad, neighbourlist, numrow, invelocity ,VELOCITY, dir);
        
			 qrSolver.solve(result,B);
        
			*vel_d = result[dir]/distance; //dir=0 (x), dir=1(y), dir=2(z)
        
            *vel_dd = pad->schemeorder == 2? result[dir+offset]/distance/distance:0;
        
			if(isnan(*p_d) || isnan(*p_dd) || isnan(*vel_d) || isnan(*vel_dd) ||
			   isinf(*p_d) || isinf(*p_dd) || isinf(*vel_d) || isinf(*vel_dd)) {
                
                numrow ++;
                
            } 
            
            else
            {
                break;
            }
        }
    }

}// to do


void LPSolver::computeA2D(double *A, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, double distance){
    double x, y,  x0, y0;
    double h,k;
    pdata_copy_t * padnei;
    x = pad->xyz[0];
    y = pad->xyz[1];
   
    if(pad->schemeorder == 1){

        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle2d(pad,&padnei,neighbourlist,i);
            x0 = padnei->xyz[0];
            y0 = padnei->xyz[1];
            h = (x0-x)/distance;
            k = (y0-y)/distance;
//           printf("%f %f %f \n",h,k,l);
            //  cout<<k<<endl;
          //  if(k>0)
          //      cout<<"warning"<<endl;
            A[i]            = h;
            A[i + 1*numrow] = k;
   
        }
    
    }
    if(pad->schemeorder == 2){
        for(size_t i=0; i<numrow; i++){
        
            gdata->fetchNeighbourParticle2d(pad,&padnei,neighbourlist,i);
            x0 = padnei->xyz[0];
            y0 = padnei->xyz[1];
            h = (x0-x)/distance;
            k = (y0-y)/distance;
        
            A[i]            = h;
			A[i + 1*numrow] = k;
			A[i + 2*numrow]	= 0.5*h*h;
			A[i + 3*numrow]	= 0.5*k*k;
			A[i + 4*numrow]	= h*k;
	
        } 
    
    
    }


}


void LPSolver::computeA3D(double *A, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, double distance){
    double x, y, z, x0, y0, z0;
    double h,k,l;
    pdata_copy_t * padnei;
    x = pad->xyz[0];
    y = pad->xyz[1];
    z = pad->xyz[2];
   
    if(pad->schemeorder == 1){
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            x0 = padnei->xyz[0];
            y0 = padnei->xyz[1];
            z0 = padnei->xyz[2];
            h = (x0-x)/distance;
            k = (y0-y)/distance;
            l = (z0-z)/distance;
//           printf("%f %f %f \n",h,k,l);
            //  cout<<k<<endl;
          //  if(k>0)
          //      cout<<"warning"<<endl;
            A[i]            = h;
            A[i + 1*numrow] = k;
            A[i + 2*numrow] = l;
   
        }
    
    }
    if(pad->schemeorder == 2){
        for(size_t i=0; i<numrow; i++){
        
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            x0 = padnei->xyz[0];
            y0 = padnei->xyz[1];
            z0 = padnei->xyz[2];
            h = (x0-x)/distance;
            k = (y0-y)/distance;
            l = (z0-z)/distance;
        
			A[i]            = h;
			A[i + 1*numrow] = k;
			A[i + 2*numrow] = l;
			A[i + 3*numrow] = 0.5*h*h;
			A[i + 4*numrow] = 0.5*k*k;
			A[i + 5*numrow] = 0.5*l*l;
			A[i + 6*numrow] = h*k;
			A[i + 7*numrow] = h*l;
			A[i + 8*numrow] = k*l;
        } 
    
    
    }


}

void LPSolver::computeB3d(double *B, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, const double* indata, indata_t datatype, int dir){
    pdata_copy_t *padnei;
    
    
    if(datatype == PRESSURE){
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            B[i] = padnei->pressure-*indata;  
//            cout<<B[i]<<endl;
        }   
    }

    else if(datatype == VELOCITY){
    
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            B[i] = padnei->v[dir] - *indata;  
        }   
    
    }

    else if(datatype == VOLUME){
    
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            B[i] = padnei->volume - *indata;  
        }   
    
    }


}
void LPSolver::computeB2d(double *B, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, const double* indata, indata_t datatype, int dir){
    pdata_copy_t *padnei;
    
    
    if(datatype == PRESSURE){
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle2d(pad,&padnei,neighbourlist,i);
            B[i] = padnei->pressure-*indata;  
//            cout<<B[i]<<endl;
        }   
    }

    else if(datatype == VELOCITY){
    
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle2d(pad,&padnei,neighbourlist,i);
            B[i] = padnei->v[dir] - *indata;  
        }   
    
    }
}


void LPSolver::timeIntegration(
	double realDt, int dir,
	double gravity, double inVolume, double inVelocity, double inPressure, double inSoundSpeed, 
	double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
	double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
	double* outVolume, double* outVelocity, double* outPressure){

    double multiplier1st = gdata->dimension == 3? 3.:2.; // 2 for 2 dimension;
    double multiplier2nd = gdata->dimension == 3? 3./4*cfldt : 1./2*cfldt;


	double K = inSoundSpeed*inSoundSpeed/inVolume/inVolume; 
    
    

	double Pt1st = -0.5*inVolume*K*(vel_d_0+vel_d_1) + 0.5*inVolume*sqrt(K)*(p_d_0-p_d_1) ;
    double Pt2nd = -inVolume*inVolume*pow(K,1.5)*(vel_dd_0-vel_dd_1) + inVolume*inVolume*K*(p_dd_0+p_dd_1);
	double Pt = multiplier1st*Pt1st + multiplier2nd*Pt2nd;
    // Vt
	
    double Vt = -Pt/K;
	// VELt
	double VELt1st = 0.5*inVolume*sqrt(K)*(vel_d_0-vel_d_1) - 0.5*inVolume*(p_d_0+p_d_1);
	
    double VELt2nd = inVolume*inVolume*K*(vel_dd_0+vel_dd_1) - inVolume*inVolume*sqrt(K)*(p_dd_0-p_dd_1);
	double VELt = multiplier1st*VELt1st + multiplier2nd*VELt2nd;
  
	(*outVolume)   = inVolume   + realDt*Vt;
	(*outPressure) = inPressure + realDt*Pt;
	(*outVelocity) = inVelocity + realDt*(VELt+gravity);	
    if(dir == 1 || dir == 2){
      	//pellet_solver->updateStatesByLorentzForce(3.*realDt, inVelocity, inPressure, inVolume, inSoundSpeed, outVelocity, outPressure);
    }
	if(isnan(*outVolume) || isinf(*outVolume) || 
	   isnan(*outPressure) || isinf(*outPressure) ||
	   isnan(*outVelocity) || isinf(*outVelocity)) {
    	printf(" %f %f\n",Pt,VELt);
	    printf("%.16g,%.16g,%.16g\n",*outVolume,*outPressure,*outVelocity);
	    printf("%.16g,%.16g,%.16g %.16g\n",inVolume,inPressure,inVelocity, inSoundSpeed);
        assert(false);   
	}

}

void LPSolver:: computeCFLCondition(){


    static bool FIRST = true;
    if(FIRST){
        cfldt = 0.;
        FIRST = false;
        return;
    }
    pdata_t *pad;
    size_t li, lpnum = gdata->particle_data->elem_count;
     
    double lmindt = 100, gmindt = 0;
    double dt;
    double speed,sc;
    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
         sc = pad->soundspeed;
         speed = sqrt(pad->v[0]*pad->v[0]+pad->v[1]*pad->v[1]+pad->v[2]*pad->v[2]);
         dt = pad->localspacing/max(sc,speed);
         if(dt < lmindt || lmindt == 100){
            lmindt = dt;
        }
    }
    int mpiret = MPI_Allreduce(&lmindt,&gmindt,1,MPI_DOUBLE,MPI_MIN,gdata->mpicomm);

    SC_CHECK_MPI (mpiret);

    cfldt = gmindt * cflcoefficient;

}

void LPSolver:: solve_2d(){

   // computeLocalBoundaryAndFluidNum();
   // viewer->writeResult(0,currenttime);
    
    while(currenttime<tend)
    {
    
        for(size_t id=0; id<gdata->boundarynumber; id++){ 
           gdata->m_vBoundary[id]->generateBoundaryParticle(gdata,gdata->eos,gdata->initlocalspacing, cfldt);
        }
        
        
        gdata->presearch2d();
        gdata->packParticles();
    
        if(gdata->gpnum == 0)
            
        {
          sc_array_destroy_null (&gdata->recevs);
          sc_hash_destroy_null (&gdata->psend);
          sc_array_destroy_null(&gdata->iremain);

          gdata->psend = NULL;
          sc_mempool_destroy (gdata->psmem);
          gdata->psmem = NULL;
            break;
        }
        gdata->communicateParticles();
        gdata->postsearch2d();

        octree->adapt_octree2d(gdata->p4est);
    
        octree->balance_octree2d(NULL,octree->balance_replace2d);
    
        gdata->regroupParticles2d(); 
    
        gdata->partitionParticles2d();
    
        MPI_Barrier(gdata->mpicomm);
        gdata->createViewForOctant2d();
    
        MPI_Barrier(gdata->mpicomm);
    
        gdata->searchNeighbourOctant2d();

        gdata->searchNeighbourParticle2d();
        
        gdata->searchUpwindNeighbourParticle2d(); 
        
  //      gdata->reorderNeighbourList2d();
        
        gdata->computeSPHDensity();
        MPI_Barrier(gdata->mpicomm);
        if(gdata->iffreeboundary) 
            gdata->generateGhostParticle2d();
   
    //gdata->testquad2d();
        computeCFLCondition();
    
        bool iswritestep = adjustDtByWriteTimeInterval(); 
        
        splitorder = (int)rand()%2;
        MPI_Bcast(&splitorder,1,MPI_INT,0,gdata->mpicomm);
       
        P4EST_GLOBAL_ESSENTIALF ("ENTER UPWIND.\n");
        for(int phase = 0;phase < totalphase;phase++){
        solve_upwind(phase);
        MPI_Barrier(gdata->mpicomm);
        gdata->updateViewForOctant2d(phase);
        MPI_Barrier(gdata->mpicomm); 
    }
    
        P4EST_GLOBAL_ESSENTIALF ("FINISH UPWIND.\n");
    if(LPFOrder == 2){
        P4EST_GLOBAL_ESSENTIALF ("ENTER LAX-WENDROFF.\n");
        solve_laxwendroff();
        P4EST_GLOBAL_ESSENTIALF ("FINISH LAX-WENDROFF.\n");
        }
        gdata->updateParticleStates();
   
   
        updateLocalSpacing();
        moveParticle();
        if(iswritestep)
        {
            computeLocalBoundaryAndFluidNum();
            viewer->writeResult(writestep,currenttime);
            if(gdata->printallparticle)
                viewer->writeGhost(writestep);
        }
   
    
        gdata->cleanForTimeStep2d();
    
        gdata->switchFlagDelete();
        
    }

}


void LPSolver::solve_3d(){


    
  //  computeLocalBoundaryAndFluidNum();
    
  //  viewer->writeResult(0,currenttime);
    
    double timer_start = MPI_Wtime();
    
    double timer_end;
    double t1, t2;
    while(currenttime < tend)
    {
        t1 = MPI_Wtime(); 
        for(size_t id=0; id<gdata->boundarynumber; id++){ 
            gdata->m_vBoundary[id]->generateBoundaryParticle(gdata,gdata->eos,gdata->initlocalspacing,cfldt);
         }
    
    MPI_Barrier(gdata->mpicomm); 
    t2 = MPI_Wtime();
    lp_timer->boundaryconditiondt = t2-t1;
    lp_timer->updateTimer(&lp_timer->boundaryconditiontimer,t2-t1);
    


    t1 = MPI_Wtime(); 
    gdata->presearch();
        
    gdata->packParticles();
    
    if(gdata->gpnum == 0)
        
    {
      sc_array_destroy_null (&gdata->recevs);
      sc_hash_destroy_null (&gdata->psend);
      sc_array_destroy_null(&gdata->iremain);

      gdata->psend = NULL;
      sc_mempool_destroy (gdata->psmem);
      gdata->psmem = NULL;
        break;
    } 
    gdata->communicateParticles();
        gdata->postsearch();

        octree->adapt_octree(gdata->p8est); 
    
        octree->balance_octree(NULL,octree->balance_replace);
    
        gdata->regroupParticles(); 
    
        gdata->partitionParticles();
    
        
    P4EST_GLOBAL_ESSENTIALF ("Created view for octants.\n");
        gdata->createViewForOctant();
    
        MPI_Barrier(gdata->mpicomm); 
        
        t2 = MPI_Wtime(); 
        lp_timer->octreedt = t2-t1;
        lp_timer->updateTimer(&lp_timer->octreetimer,t2-t1);


        t1 = MPI_Wtime(); 
        P4EST_GLOBAL_ESSENTIALF ("Search neighbour for octants.\n");
        gdata->searchNeighbourOctant();



    P4EST_GLOBAL_ESSENTIALF ("Search neighbour for particles.\n");
        gdata->searchNeighbourParticle();
        
    P4EST_GLOBAL_ESSENTIALF ("Search upwind neighbour for particles.\n");
        gdata->searchUpwindNeighbourParticle(); 
       
       // gdata->adjustLocalSpacingForParticle();
        MPI_Barrier(gdata->mpicomm); 
    
        t2 = MPI_Wtime(); 
        lp_timer->neighboursearcherdt = t2-t1;
        lp_timer->updateTimer(&lp_timer->neighboursearchertimer,t2-t1);
        P4EST_GLOBAL_ESSENTIALF ("ComputeSPHdensity for particles.\n");
        gdata->computeSPHDensity();



    
        t1 = MPI_Wtime(); 
        if(gdata->iffreeboundary){ 
            
            P4EST_GLOBAL_ESSENTIALF ("generateGhostParticle for particles.\n");
            gdata->generateGhostParticle();
        }
         
            P4EST_GLOBAL_ESSENTIALF ("reorder neighbour for particles.\n");
        gdata->reorderNeighbourList();
        
        MPI_Barrier(gdata->mpicomm); 
        t2 = MPI_Wtime(); 
        lp_timer->ghostparticledt = t2-t1;
        lp_timer->updateTimer(&lp_timer->ghostparticletimer,t2-t1);


        t1 = MPI_Wtime(); 
        if(gdata->pelletnumber){ 
            
            P4EST_GLOBAL_ESSENTIALF ("heatingModell.\n");
            pellet_solver->heatingModel(cfldt);
        }
        
        MPI_Barrier(gdata->mpicomm); 
        t2 = MPI_Wtime(); 
        lp_timer->heatingmodeldt = t2-t1;
        lp_timer->updateTimer(&lp_timer->heatingmodeltimer,t2-t1); 
        computeCFLCondition();
    

        bool iswritestep = adjustDtByWriteTimeInterval(); 
       


        t1 = MPI_Wtime(); 
        MPI_Barrier(gdata->mpicomm); 
        
        splitorder = (int)rand()%6;
        MPI_Bcast(&splitorder,1,MPI_INT,0,gdata->mpicomm);
        
        P4EST_GLOBAL_ESSENTIALF ("ENTER UPWIND.\n");
        for(int phase = 0;phase< totalphase;phase++){
            solve_upwind(phase);
            MPI_Barrier(gdata->mpicomm);
            P4EST_GLOBAL_ESSENTIALF ("Phase %d FINISHED..\n",phase);
            gdata->updateViewForOctant(phase);
            MPI_Barrier(gdata->mpicomm); 
            P4EST_GLOBAL_ESSENTIALF ("Update View Finished.\n");
        }
   
        P4EST_GLOBAL_ESSENTIALF ("FINISH UPWIND.\n");
        t2 = MPI_Wtime(); 
        lp_timer->upwinddt = t2-t1;
        lp_timer->updateTimer(&lp_timer->upwindtimer,t2-t1);
    if(LPFOrder == 2){

        t1 = MPI_Wtime(); 
        P4EST_GLOBAL_ESSENTIALF ("ENTER LAX-WENDROFF.\n");
        solve_laxwendroff();
        P4EST_GLOBAL_ESSENTIALF ("FINISH LAX-WENDROFF.\n");
        MPI_Barrier(gdata->mpicomm); 
        t2 = MPI_Wtime(); 
        lp_timer->lwdt = t2-t1;
        lp_timer->updateTimer(&lp_timer->lwtimer,t2-t1);
    }
    
    pellet_solver->computeBoundaryCondition(gdata,cfldt,gdata->initlocalspacing/5.);
    gdata->updateParticleStates();
    updateLocalSpacing();
//    pellet_solver->computeDriftAcceleration();
    pellet_solver->gradBDrift(cfldt); 
    pellet_solver->updateStatesByLorentzForce(cfldt);
    pellet_solver->neonRadiationCooling(cfldt);
    moveParticle();
    computeTemperature();
        
        if(iswritestep)
        
        {
            pellet_solver->labelLayer();
            computeLocalBoundaryAndFluidNum();
            viewer->writeResult(writestep,currenttime);
       //     viewer->writeComputationNode(writestep);
            if(gdata->printallparticle)
                viewer->writeGhost(writestep);
        //      viewer->writeTXTFile(writestep);
        }
    
        gdata->cleanForTimeStep();
//        gdata->deleteBadStateParticle();
        gdata->switchFlagDelete();
    timer_end = MPI_Wtime();
    lp_timer->writeTimeAnalysis(gdata->mpirank,timer_end-timer_start);
    }

}
void LPSolver:: computeLocalBoundaryAndFluidNum(){

    pdata_t *pad;
    size_t li, lpnum = gdata->particle_data->elem_count;
    size_t boundarycount = 0, fluidcount = 0; 

    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           boundarycount ++;
       else
           fluidcount ++;
    }
    gdata->lfluidnum = fluidcount;
    gdata->lboundarynum = boundarycount;
}


void LPSolver::solve_laxwendroff(){
    
    const double *invelocityu, *invelocityv, *invelocityw;
    const double *inpressure, *involume, *insoundspeed;
    
    double *outvelocityu, *outvelocityv, *outvelocityw;
    double *outpressure, *outvolume, *outsoundspeed;

    int number_of_derivative = gdata->dimension==3? 9:5;
    double Ud[number_of_derivative],Vd[number_of_derivative],Wd[number_of_derivative],
           Pd[number_of_derivative],Volumed[number_of_derivative];
    pdata_t * pad;

    size_t li, lpnum = gdata->particle_data->elem_count;

    for(li = 0; li<lpnum; li++){
    
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
       if(pad->neighbourparticle->elem_count  < (gdata->numrow2nd+pad->redocount))
       {    
           continue;
       }
       
       bool redo = false;
       
       inpressure = &pad->pressure;
       insoundspeed = &pad->soundspeed;
       involume = &pad->volume;
       invelocityu = &pad->v[0];
       invelocityv = &pad->v[1];
       if(gdata->dimension == 3)
           invelocityw = &pad->v[2];
        
       if(gdata->dimension == 3){
          outpressure = &pad->pressureT1;
          outsoundspeed = &pad->soundspeedT1;
          outvolume = &pad->volumeT1;
       }
       else if(gdata->dimension == 2){
           outpressure = &pad->pressureT2;
           outsoundspeed = &pad->soundspeedT2;
           outvolume = &pad->volumeT2;
       }
        
       outvelocityu = &pad->oldv[0];
       outvelocityv = &pad->oldv[1];
       if(gdata->dimension == 3)
           outvelocityw = &pad->oldv[2];
    

       if(*insoundspeed == 0 || *involume == 0){
            pad->schemeorder = 0;
            *outvolume = *involume;
            *outpressure = *inpressure;
            *outsoundspeed = *insoundspeed;
            *outvelocityu = *invelocityu;
            *outvelocityv = *invelocityv;
            if(gdata->dimension == 3)
                *outvelocityw = *invelocityw;
            printf("Detect a particle which has 0 volume or 0 soundspeed!!!.\n"); 
            continue;
       }
      /*
       if((x*x+y*y+z*z)<0.3*0.3)
            continue;
       */
       pad->schemeorder = 2;
        
       computeSpatialDer(pad, inpressure, involume, invelocityu, invelocityv, invelocityw,
          Pd, Volumed, Ud, Vd, Wd);
       double gravity = 0;
       if(gdata->dimension == 3){   
            timeIntegration( gravity, *involume, *invelocityu, *invelocityv, *invelocityw, *inpressure, *insoundspeed,
                            Volumed, Ud, Vd, Wd, Pd,
                            outvolume, outvelocityu, outvelocityv, outvelocityw, outpressure); // output 
       }

       else if(gdata->dimension == 2){
            double temp1 = 0, temp2; 
       
            timeIntegration( gravity, *involume, *invelocityu, *invelocityv, temp1, *inpressure, *insoundspeed,
                            Volumed, Ud, Vd, Wd, Pd,
                            outvolume, outvelocityu, outvelocityv, &temp2, outpressure); // output 
       }
   
         if(gdata->pelletnumber) {
            *outpressure += cfldt*pad->deltaq*((*insoundspeed)*(*insoundspeed)/(*involume)/(*inpressure)-1);
         } 
         if(*outpressure < gdata->invalidpressure || 1./(*outvolume) < gdata->invaliddensity)
         {
            
       /*         
                printf("inpressure %.8g outpressure %.8g involume %.8g outvolume %.8g %s\n",*inpressure, *outpressure, *involume, *outvolume ,"lw");
                printf("%f %f %f %d\n",pad->xyz[0],pad->xyz[1],pad->xyz[2],li);
       
                pdata_copy_t *padnei;                
                
                
                for(int i=0;i<53;i++){
                    gdata->fetchNeighbourParticle(pad,&padnei,pad->neighbourparticle,i);
                    
                    printf("%f %f %f %f %f \n",padnei->xyz[0],padnei->xyz[1],padnei->xyz[2], padnei->pressure, padnei->volume);
                }
          //  exit(0);
        */
             redo  = true;
         }
        if(isnan(*outvolume) || isnan(*outvelocityu) || isnan(*outvelocityv) || isnan(*outpressure)) 
        
        {
            redo = true;
            if(gdata->dimension == 3){
                if(isnan(*outvelocityw))
                {    redo = true;
            
                }
            }
        }

        if(isinf(*outvolume) || isinf(*outvelocityu) || isinf(*outvelocityv) || isinf(*outpressure)) 
        {
            redo = true;
            if(gdata->dimension == 3){
                if(isinf(*outvelocityw))
                {
                    redo = true;
                }
            }
        }
    
        if(redo == true){

            pad->redocount++;
            
            if(pad->redocount > 4){
                pad->schemeorder = 0;
            }
            else{
                li--;
                
                *outvolume = *involume;
                *outpressure = *inpressure;
                *outsoundspeed = *insoundspeed;
                *outvelocityu = *invelocityu;
                *outvelocityv = *invelocityv;
                cout<<"redo laxwendroff for a particle"<<li<<endl;
            }

            
        }
        else{
            
            *outsoundspeed = gdata->eos->getSoundSpeed(*outpressure, 1./(*outvolume));
             pad->redocount = 0;        
        }

        if(pad->schemeorder == 0)
        {

            pad->redocount = 0;
            *outvolume = *involume;
            *outpressure = *inpressure;
            *outsoundspeed = *insoundspeed;
            *outvelocityu = *invelocityu;
            *outvelocityv = *invelocityv;
            if(gdata->dimension == 3)
                *outvelocityw = *invelocityw;
        }
    
    }

}



// for lw scheme
void LPSolver::computeSpatialDer(pdata_t *pad, const double* inPressure,  const double* inVolume,
        const double* inVelocityU, const double* inVelocityV, const double* inVelocityW,
         double* Pd, double *Volumed, double* Ud, double* Vd, double* Wd){

    
    size_t offset = gdata->dimension == 3? 3:2;
    size_t numrow = pad->redocount + gdata->numrow2nd;
    size_t numcol = gdata->numcol2nd;
    double distance;
    int info;
    neighbour_info_t *ninfo;
    sc_array_t *neighbourlist = pad->neighbourparticle;
    
    while(true){
   
        for(size_t i=0;i<numcol;i++)
        {
                Pd[i]=0.0;
                Ud[i]=0.0;
                Vd[i]=0.0;
                Wd[i]=0.0;
                Volumed[i]=0.0;
        }
        if(neighbourlist->elem_count < numrow || numrow > 100){
            pad->schemeorder = 0;
            break;
        }
    
        ninfo = (neighbour_info_t *)sc_array_index(neighbourlist,0);
        distance = ninfo->distance;
        double A[numrow*numcol];
        if(gdata->dimension == 3)
            computeA3D(A, pad, neighbourlist,numrow,distance);
        else if(gdata->dimension == 2)
            computeA2D(A, pad, neighbourlist,numrow,distance);
    
         double B[numrow];
        if(gdata->dimension == 3)
            computeB3d(B, pad, neighbourlist, numrow, inPressure ,PRESSURE, -1); //-1 is meaning less
        else if(gdata->dimension == 2)
            computeB2d(B, pad, neighbourlist, numrow, inPressure ,PRESSURE, -1);
        QRSolver qrSolver(numrow,numcol,A);
		
		double result[numcol];
	    info = qrSolver.solve(result,B);
        
        if(info != 0){
            numrow ++;
            //printf("Warning: Rank Deficient For a Particle in LW!\n");
        } 
        else{
       
            for(size_t i=0;i<numcol;i++){
                
                if(i<offset)
                    Pd[i] = result[i]/distance;
                else
                    Pd[i] = result[i]/distance/distance;
           
            }
        
             
        if(gdata->dimension == 3)
            {computeB3d(B, pad, neighbourlist, numrow, inVelocityU ,VELOCITY, 0);}
        else if(gdata->dimension == 2)
            {computeB2d(B, pad, neighbourlist, numrow, inVelocityU ,VELOCITY, 0);}

		
	    info = qrSolver.solve(result,B);
        for(size_t i=0;i<numcol;i++){
            
            if(i<offset)
                Ud[i] = result[i]/distance;
            else
                Ud[i] = result[i]/distance/distance;
            
            }
        if(gdata->dimension == 3)
            {computeB3d(B, pad, neighbourlist, numrow, inVelocityV ,VELOCITY, 1);}
        else if(gdata->dimension == 2)
            {computeB2d(B, pad, neighbourlist, numrow, inVelocityV ,VELOCITY, 1);}

		
	    info = qrSolver.solve(result,B);
        for(size_t i=0;i<numcol;i++){
            
            if(i<offset)
                Vd[i] = result[i]/distance;
            else
                Vd[i] = result[i]/distance/distance;
            
            }

        if(gdata->dimension == 3){
            computeB3d(B, pad, neighbourlist, numrow, inVelocityW ,VELOCITY, 2);

            info = qrSolver.solve(result,B);
            for(size_t i=0;i<numcol;i++){
                
                if(i<offset)
                    Wd[i] = result[i]/distance;
                else
                    Wd[i] = result[i]/distance/distance;
                
                }
            }

        if(gdata->dimension == 3)
            {computeB3d(B, pad, neighbourlist, numrow, inVolume ,VOLUME, -1);}
        else if(gdata->dimension == 2)
            {computeB2d(B, pad, neighbourlist, numrow, inVolume ,VOLUME, -1);}

        info = qrSolver.solve(result,B);
        for(size_t i=0;i<numcol;i++){
            
            if(i<offset)
                Volumed[i] = result[i]/distance;
            else
                Volumed[i] = result[i]/distance/distance;
            
            }

        bool success = true;
        for(size_t i=0; i<numcol; i++){
            if(isnan(Pd[i]) || isnan(Ud[i]) || isnan(Vd[i]) || isnan(Wd[i]) || isnan(Volumed[i]) ||
               isinf(Pd[i]) || isinf(Ud[i]) || isinf(Vd[i]) || isinf(Wd[i])|| isnan(Volumed[i])) 
            {
                success = false;
            }
        
        }
    
        if(success)
            break;
        else
            numrow ++;
        }
    }

}



//for laxwendroff
void LPSolver::timeIntegration( double gravity, double inVolume, double inVelocityU, double inVelocityV, double inVelocityW, double inPressure, double inSoundSpeed,
                                      double* Volumed, double* Ud, double* Vd, double *Wd, double *Pd,
                                                        double* outVolume, double* outVelocityU, double* outVelocityV, double* outVelocityW, double* outPressure){

      double gamma=inSoundSpeed*inSoundSpeed/inVolume/inPressure;
      double Pinf = gdata->pinf;
      double Dt = cfldt;
      if(gdata->dimension==3)
        {

                double cond = gdata->eos->getElectricConductivity(inPressure,1./inVolume);
                double clight = 2.99792458e7;
                double B = pellet_solver->magneticfield;
                double mu = cond*B*B/clight/clight;
                double div=Ud[0]+Vd[1]+Wd[2];
                double cross=Ud[0]*Ud[0]+Vd[1]*Vd[1]+Wd[2]*Wd[2]+2*Ud[1]*Vd[0]+2*Ud[2]*Wd[0]+2*Vd[2]*Wd[1];
                double Volumet=inVolume*div;
                double VelocityUt=-inVolume*Pd[0];
                double VelocityVt=-inVolume*Pd[1];//-mu*inVelocityV*inVolume;
                double VelocityWt=-inVolume*Pd[2];//-mu*inVelocityW*inVolume;
                double Pt=-gamma*(inPressure+Pinf)*div;// + (gamma-1)*mu*(inVelocityV*inVelocityV+inVelocityW*inVelocityW);
                double Volumett=inVolume*(div*div-Volumed[0]*Pd[0]-Volumed[1]*Pd[1]-Volumed[2]*Pd[2]-inVolume*(Pd[3]+Pd[4]+Pd[5])-cross);
                double VelocityUtt=inVolume*((gamma-1)*Pd[0]*div+gamma*(inPressure+Pinf)*(Ud[3]+Vd[6]+Wd[7])+Ud[0]*Pd[0]+Vd[0]*Pd[1]+Wd[0]*Pd[2]);
                double VelocityVtt=inVolume*((gamma-1)*Pd[1]*div+gamma*(inPressure+Pinf)*(Ud[6]+Vd[4]+Wd[8])+Ud[1]*Pd[0]+Vd[1]*Pd[1]+Wd[1]*Pd[2]);
                double VelocityWtt=inVolume*((gamma-1)*Pd[2]*div+gamma*(inPressure+Pinf)*(Ud[7]+Vd[8]+Wd[5])+Ud[2]*Pd[0]+Vd[2]*Pd[1]+Wd[2]*Pd[2]);
                double Ptt=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+Volumed[2]*Pd[2]+inVolume*(Pd[3]+Pd[4]+Pd[5])+cross);
               
/*
                double extratermV = -mu*inVolume*(Volumed[1]*inVelocityV+Volumed[2]*inVelocityW) - inVolume*inVolume*mu*(Vd[1]+Wd[2]);
                Volumett += extratermV;
                
                double extratermu = -inVolume*(gamma-1)*mu*(2*inVelocityV*Vd[0]+2*inVelocityW*Wd[0]);
                VelocityUtt += extratermu;

                double extratermv = -inVolume*(gamma-1)*mu*(2*inVelocityV*Vd[1]+2*inVelocityW*Wd[1]);// - inVolume*mu*inVelocityV*(div) + inVolume*inVolume*mu*Pd[1]+mu*mu*inVolume*inVolume*inVelocityV;
                VelocityVtt += extratermv;
                
                double extratermw = -inVolume*(gamma-1)*mu*(2*inVelocityV*Vd[2]+2*inVelocityW*Wd[2]);// - inVolume*mu*inVelocityW*(div) + inVolume*inVolume*mu*Pd[2]+mu*mu*inVolume*inVolume*inVelocityW;
                VelocityWtt += extratermw;
                
                double extratermP = -gamma*(gamma-1)*mu*div*(inVelocityW*inVelocityW+inVelocityV*inVelocityV) - 2*mu*((Pd[1]+mu*inVelocityV)*inVolume*inVelocityV+(Pd[2]+mu*inVelocityW)*inVolume*inVelocityW)*(gamma-1)
                       +gamma*inPressure*mu*((Volumed[1]*inVelocityV+Volumed[2]*inVelocityW) + inVolume*(Vd[1]+Wd[2]));
                    
                Ptt += extratermP;
*/
                (*outVolume)=inVolume+Dt*Volumet+0.5*Dt*Dt*Volumett;
                (*outVelocityU)=inVelocityU+Dt*VelocityUt+0.5*Dt*Dt*VelocityUtt;
                (*outVelocityV)=inVelocityV+Dt*VelocityVt+0.5*Dt*Dt*VelocityVtt+Dt*gravity;//TODO: MODIFIED GRAVITY DIRECTION
                (*outVelocityW)=inVelocityW+Dt*VelocityWt+0.5*Dt*Dt*VelocityWtt;
                (*outPressure)=inPressure+Dt*Pt+0.5*Dt*Dt*Ptt;
                //pellet_solver->updateStatesByLorentzForce(Dt, inVelocityV, inVelocityW, inPressure, inVolume, inSoundSpeed, 
                //     outVelocityV, outVelocityW, outPressure);
                
                if(isnan(*outVolume) || isinf(*outVolume) ||
                   isnan(*outPressure) || isinf(*outPressure) ||
                 isnan(*outVelocityU) || isinf(*outVelocityU) ||isnan(*outVelocityV) || isinf(*outVelocityV) ||isnan(*outVelocityW) || isinf(*outVelocityW)) {
                        assert(false);
                }

        }
      else if(gdata->dimension==2)
        {

                double div=Ud[0]+Vd[1];
                double cross=Ud[0]*Ud[0]+Vd[1]*Vd[1]+2*Ud[1]*Vd[0];
                double Volumet=inVolume*div;
                double VelocityUt=-inVolume*Pd[0];
                double VelocityVt=-inVolume*Pd[1];
                double Pt=-gamma*(inPressure+Pinf)*div;
                double Volumett=inVolume*(div*div-Volumed[0]*Pd[0]-Volumed[1]*Pd[1]-inVolume*(Pd[2]+Pd[3])-cross);
                double VelocityUtt=inVolume*((gamma-1)*Pd[0]*div+gamma*(inPressure+Pinf)*(Ud[2]+Vd[4])+Ud[0]*Pd[0]+Vd[0]*Pd[1]);
                double VelocityVtt=inVolume*((gamma-1)*Pd[1]*div+gamma*(inPressure+Pinf)*(Ud[4]+Vd[3])+Ud[1]*Pd[0]+Vd[1]*Pd[1]);
                double Ptt=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+inVolume*(Pd[2]+Pd[3])+cross);
                (*outVolume)=inVolume+Dt*Volumet+0.5*Dt*Dt*Volumett;
                (*outVelocityU)=inVelocityU+Dt*VelocityUt+0.5*Dt*Dt*VelocityUtt;
                (*outVelocityV)=inVelocityV+Dt*VelocityVt+0.5*Dt*Dt*VelocityVtt+Dt*gravity;

                (*outPressure)=inPressure+Dt*Pt+0.5*Dt*Dt*Ptt;
                if(isnan(*outVolume) || isinf(*outVolume) ||
                   isnan(*outPressure) || isinf(*outPressure) ||
                 isnan(*outVelocityU) || isinf(*outVelocityU) ||isnan(*outVelocityV) || isinf(*outVelocityV)) {
                        assert(false);
                }
                Pd[3]=gamma*(inPressure+Pinf)*(inVolume*(Pd[2]+Pd[3]));
                Pd[4]=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+cross);
                Pd[2]=-gamma*(inPressure+Pinf)*div;

        }
}



bool LPSolver::adjustDtByWriteTimeInterval() {
	timestep += 1;
    bool ifadjust;
    if(currenttime+cfldt >= nextwritetime) {
		cfldt = nextwritetime - currenttime;
	    	
        nextwritetime += writetimeinterval;
		writestep ++;
        if(nextwritetime > tend) nextwritetime = tend;
		assert(cfldt >= 0);
		//cout<<"-------TimeController::adjustDtByWriteTimeInterval()-------"<<endl;
		//cout<<"m_fDt="<<m_fDt<<endl;
		//cout<<"-----------------------------------------------------------"<<endl;
		ifadjust = true;
	}
    else if(cfldt == 0.){
        ifadjust = true;
    }
    else{	
        ifadjust = false;
    }
    currenttime += cfldt;
    P4EST_GLOBAL_ESSENTIALF ("MINCFL timestep is %.16g. \n", cfldt);
    P4EST_GLOBAL_ESSENTIALF ("Number of Timesteps: %d     Current Time: %f .\n", timestep, currenttime);
    return ifadjust; // m_fDt did not get adjusted
}


void LPSolver::computeTemperature(){
    pdata_t *pad;
    size_t li, lpnum = gdata->particle_data->elem_count;

    for(li=0;li<lpnum;li++){
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary)
            continue;
        pad->temperature = gdata->eos->getTemperature(pad->pressure,1./pad->volume);
    }
}







