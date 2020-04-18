#include <iostream>
#include "material_lib.h"
#include <map>
#include <math.h>
using namespace std;


Neon::Neon(){
    mu = 20.1797; 
    mass = 3.351e-23;
    Z = 10;
    I = 135.5;
    sublimationEnergy = 1363;
}

double Neon::getOne_Plus_Zstar(double teinf){

    static bool FIRST = true;
    static std::map<double, double> map_one_plus_Zstar;                                                     
     if (FIRST)                                                                                           
	  {                                                                                                    
	    map_one_plus_Zstar[1000] = 4.93908;                                                               
	    map_one_plus_Zstar[2000] = 4.92710;                                                               
	    map_one_plus_Zstar[3000] = 4.93727;                                                               
	    map_one_plus_Zstar[4000] = 4.94743;                                                               
	    map_one_plus_Zstar[6000] = 4.96337;                                                               
	    map_one_plus_Zstar[8000] = 4.97493;                                                               
	    map_one_plus_Zstar[10000] = 4.98376;                                                              
                                                                                                             
	    one_plus_Zstar = map_one_plus_Zstar[teinf];                                          
        cout<<"one plus zstar = "<<one_plus_Zstar<<endl;      
        if (floor(one_plus_Zstar) != 4)                                                                   
	      printf("in %s, line %d, function %s: error in map initialization of 1+Z*\n", __FILE__,__LINE__,\
		     __PRETTY_FUNCTION__);                                                                                        
                                                                                                             
	    FIRST = false;        
      }

     return one_plus_Zstar;

}


Deuterium2::Deuterium2(){
    mu = 4.028;//2.014; 
    mass = 3.34e-24;
    Z = 1;
    I = 19.2;
    sublimationEnergy = 3113.95;
}

double Deuterium2::getOne_Plus_Zstar(double teinf){

    static bool FIRST = true;
    static std::map<double, double> map_one_plus_Zstar;                                                     
     if (FIRST)                                                                                           
	  {                                                                                                    
        map_one_plus_Zstar[500] = 0.951000;
        map_one_plus_Zstar[1000] = 0.952556;                                                               
	    map_one_plus_Zstar[2000] = 1.000534;//0.957404;                                                               
	    map_one_plus_Zstar[3000] = 0.959938;                                                               
	    map_one_plus_Zstar[4000] = 0.961598;                                                               
	    map_one_plus_Zstar[6000] = 0.963751;                                                               
	    map_one_plus_Zstar[8000] = 0.965155;                                                               
	    map_one_plus_Zstar[10000] = 0.966178;                                                              
                                                                                 
	    one_plus_Zstar = map_one_plus_Zstar[teinf];                                          
        cout<<"one plus zstar = "<<one_plus_Zstar<<endl;      
        if (floor(one_plus_Zstar) != 4)                                                                   
	      printf("in %s, line %d, function %s: error in map initialization of 1+Z*\n", __FILE__,__LINE__,\
		     __PRETTY_FUNCTION__);                                                                                        
                                                                                                             
	    FIRST = false;        
      }

     return one_plus_Zstar;

}

Deuterium::Deuterium(){
    mu = 2.014; 
    mass = 3.34e-24;
    Z = 1;
    I = 14.9913;
    sublimationEnergy = 3113.95;
}

double Deuterium::getOne_Plus_Zstar(double teinf){

    static bool FIRST = true;
    static std::map<double, double> map_one_plus_Zstar;                                                     
     if (FIRST)                                                                                           
	  {                                                                                                    
        map_one_plus_Zstar[500] = 0.951000;
        map_one_plus_Zstar[1000] = 0.952556;                                                               
	    map_one_plus_Zstar[2000] = 0.957404;                                                               
	    map_one_plus_Zstar[3000] = 0.959938;                                                               
	    map_one_plus_Zstar[4000] = 0.961598;                                                               
	    map_one_plus_Zstar[6000] = 0.963751;                                                               
	    map_one_plus_Zstar[8000] = 0.965155;                                                               
	    map_one_plus_Zstar[10000] = 0.966178;                                                              
                                                                                 
	    one_plus_Zstar = map_one_plus_Zstar[teinf];                                          
        cout<<"one plus zstar = "<<one_plus_Zstar<<endl;      
        if (floor(one_plus_Zstar) != 4)                                                                   
	      printf("in %s, line %d, function %s: error in map initialization of 1+Z*\n", __FILE__,__LINE__,\
		     __PRETTY_FUNCTION__);                                                                                        
                                                                                                             
	    FIRST = false;        
      }

     return one_plus_Zstar;

}
