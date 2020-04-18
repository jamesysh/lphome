#include <iostream>
#include "eos.h"
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
using namespace std;

#define INTERP_TYPE gsl_interp2d_bicubic // bicubic or bilinear

static const gsl_interp2d_type *T_sound_speed = INTERP_TYPE;
static gsl_spline2d *spline_sound_speed;
static const gsl_interp2d_type *T_temperature = INTERP_TYPE;
static gsl_spline2d *spline_temperature;
static const gsl_interp2d_type *T_conductivity = INTERP_TYPE;
static gsl_spline2d *spline_conductivity;
static const gsl_interp2d_type *T_elec_frac_D = INTERP_TYPE;
static gsl_spline2d *spline_elec_frac_D;
static const gsl_interp2d_type *T_elec_frac_D2 = INTERP_TYPE;
static gsl_spline2d *spline_elec_frac_D2;
static const gsl_interp2d_type *T_diss_frac_D2 = INTERP_TYPE;
static gsl_spline2d *spline_diss_frac_D2;
static const gsl_interp2d_type *T_sound_speed_D2 = INTERP_TYPE;
static gsl_spline2d *spline_sound_speed_D2;


static void adjustStates(double* pres, double* dens,double p0, double p1, double rho0, double rho1);
////////////////////////////////////////////////////////////////////////////////
// Start of PolytropicGasEOS
////////////////////////////////////////////////////////////////////////////////




double PolytropicGasEOS::getEnergy(double pressure, double density) {
	if(((m_fGamma - 1.) * density) != 0) 
		return ( pressure / ((m_fGamma - 1.) * density) );
	else { // divide by zero
		std::cout<<"Error (Divide by zero)! Computing energy by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", density = "<<density<<std::endl;
		assert(false);
	}
}




double PolytropicGasEOS::getSoundSpeed(double pressure, double density) {
	double cs;
	
    if(density != 0)
    cs = m_fGamma * pressure / density;
     
    else {
		std::cout<<"Error (Divide by zero density)! Computing sound speed by EOS: "<<std::endl;
		//std::cout<<"density = "<<density<<std::endl;
		assert(false);
	}
	if(cs > 0) 
		return sqrt(cs);
	else if(cs==0)
		return cs;
	else { // taking square root of a negative number
		std::cout<<"Error (Taking suqre root of a negative number)! Computing sound speed by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", pressure = "<<pressure<<", density = "<<density<<std::endl;
		assert(false);
	}
}

double PolytropicGasEOS::getElectricConductivity(double pressure, double density) {

    double cond = 0.;
    return cond;

}

double PolytropicGasEOS::getTemperature(double pressure, double density) {
  double R,mu;
  
    mu = 20.18;
    R = 83.14;
  return mu*pressure/(R*density)/11604.525;
}


////////////////////////////////////////////////////////////////////////////////
// End of PolytropicGasEOS
////////////////////////////////////////////////////////////////////////////////
void  SahaNeonEOS::getParameters(std::vector<double>& params)
{
  FILE *fp_rho, *fp_pres, *fp_sound_speed, *fp_temperature, *fp_conductivity;
  int dim = 1000;
  double* rho = new double[dim];
  double* pres = new double[dim];
  double* sound_speed = new double[dim*dim];
  double* temperature = new double[dim*dim];
  double* conductivity = new double[dim*dim];
  double rho_entry, pres_entry, sc_entry, t_entry, cond_entry;
  int i,j;

  printf("ENTERED  SahaNeon EOS::getParameters \n");
  char* buffer;
  buffer = getcwd(NULL,0);
  string sub_dir = "/tables_1000/Neon/";
  string main_dir(buffer);
  string dir = main_dir+sub_dir;
  
  if (!(fp_rho = fopen((dir+"rho.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS rho.dat file \n");

  if (!(fp_pres = fopen((dir+"pres.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS pres.dat file \n");

  if (!(fp_sound_speed = fopen((dir+"sound_speed.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS sound_speed.dat file \n");

  if (!(fp_temperature = fopen((dir+"temp.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS temp.dat file \n");

  if (!(fp_conductivity = fopen((dir+"conductivity.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS conductivity.dat file \n");

  printf("Opend files in  SahaEOS::getParameters \n");
 
  i = 0;
  while (fscanf(fp_rho,"%lf",&rho_entry)!=EOF)
    rho[i++] = rho_entry;
  i = 0;
  while (fscanf(fp_pres,"%lf",&pres_entry)!=EOF)
    pres[i++] = pres_entry;

  spline_sound_speed = gsl_spline2d_alloc(T_sound_speed, dim, dim);
  spline_temperature = gsl_spline2d_alloc(T_temperature, dim, dim);
  spline_conductivity = gsl_spline2d_alloc(T_conductivity, dim, dim);

  printf("BUILT 2D splies in  SahaNeon EOS::getParameters \n");

  for (i=0; i<dim; ++i)
    for (j=0; j<dim; ++j)
    {
      fscanf (fp_sound_speed,"%lf",&sc_entry);
      gsl_spline2d_set(spline_sound_speed, sound_speed, i,j, sc_entry);

      fscanf (fp_temperature,"%lf",&t_entry);
      gsl_spline2d_set(spline_temperature, temperature, i,j, t_entry);
      
      fscanf (fp_conductivity,"%lf",&cond_entry);
      gsl_spline2d_set(spline_conductivity, conductivity, i,j, cond_entry);           
    }
  
  fclose(fp_rho);
  fclose(fp_pres);
  fclose(fp_sound_speed);
  fclose(fp_temperature);
  fclose(fp_conductivity);
  
  gsl_spline2d_init(spline_sound_speed, rho, pres, sound_speed, dim, dim);
  gsl_spline2d_init(spline_temperature, rho, pres, temperature, dim, dim);
  gsl_spline2d_init(spline_conductivity, rho, pres, conductivity, dim, dim);
  delete[] rho;
  delete[] pres;
  delete[] sound_speed;
  delete[] temperature;
  delete[] conductivity;
}



double SahaNeonEOS::getEnergy(double pressure, double density) {
  double energy = 0.0;

  return energy;
}


double SahaNeonEOS::getSoundSpeed(double pressure, double density) {
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double sc;

  if (pressure < 1.e-10 || pressure > 200.0 || density < 1.e-12 || density > 2)
    sc = sqrt(1.67*pressure/density);
  else
    sc =  gsl_spline2d_eval(spline_sound_speed, density, pressure, densacc, presacc);

  //  printf("rho = %e  P = %e  sc = %e\n",density,pressure,sc);

  return sc;
}

double SahaNeonEOS::getElectricConductivity(double pressure, double density) {
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double cond;

   if (pressure < 1.e-10 || pressure > 200.0 || density < 1.e-12 || density > 2)
      cond = 0.0;
    else
    cond =  gsl_spline2d_eval(spline_conductivity, density, pressure, densacc, presacc);
 
  return cond;
}

double SahaNeonEOS::getTemperature(double pressure, double density) {
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double T;
  double mu = 20.18;
  double R = 83.14;
  
   if (pressure < 1.e-10 || pressure > 200.0 || density < 1.e-12 || density > 2)
      T = mu*pressure/(R*density)/11604.525;
    else
    T =  gsl_spline2d_eval(spline_temperature, density, pressure, densacc, presacc);
 
  return T;
}

static void make_ionization_table_Deuterium();
static void make_ionization_dissociation_table_Deuterium2();

static double getIonizationFraction_Deuterium(double pressure, double density) 
{
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double fraction = 0;

   adjustStates(&pressure,&density,1.1e-8,199.9,1.1e-10,0.49);

   fraction =  gsl_spline2d_eval(spline_elec_frac_D, density, pressure, densacc, presacc);
  return fraction;

}

static double getIonizationFraction_Deuterium2(double pressure, double density) 
{
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double fraction = 0;

   adjustStates(&pressure,&density,1.1e-8,199.9,1.1e-10,0.49);

   fraction =  gsl_spline2d_eval(spline_elec_frac_D2, density, pressure, densacc, presacc);
  return fraction;

}

 double getDissociationFraction_Deuterium2(double pressure, double density) 
{
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double fraction;

   adjustStates(&pressure,&density,1.1e-8,199.9,1.1e-10,0.49);
   fraction =  gsl_spline2d_eval(spline_diss_frac_D2, density, pressure, densacc, presacc);
 
  return fraction;

}

void SahaDeuterium2EOS::getParameters(std::vector<double>& params){

       params.push_back(m_fGamma);
       printf("ENTERED  SahaDeuterium2 EOS::getParameters \n");
       double pres0 = 1.0e-8;
       double presend = 2.0e2;
       double rho0 =  1.0e-10;
       double rhoend = 0.5;
       int Nrho = 1000, Npres = 1000;
       double *dens_array = new double[Nrho];
       double *pres_array = new double[Npres];
  
       double *sound_speed_D2 = new double[Nrho*Npres];
       double dlogpres = (log10(presend) - log10(pres0))/Npres;
       double dlogrho = (log10(rhoend) - log10(rho0))/Nrho;
       printf("ENTERED SahaDeuterium2  EOS:: Make SoundSpeed Table \n");
       char* buffer;
       
       FILE* fp_sound_speed_D2;
       buffer = getcwd(NULL,0);
  string sub_dir = "/tables_1000/Deuterium2/";
  string main_dir(buffer);
  string dir = main_dir+sub_dir;
      
  if (!(fp_sound_speed_D2 = fopen((dir+"sound_speed_for_D2.dat").c_str(),"r")))
    printf("CANNOT OPEN SahaDeuterium EOS sound_speed_for_D2.dat file \n");
 
  spline_sound_speed_D2 = gsl_spline2d_alloc(T_sound_speed_D2, Nrho, Npres); 
 
    double rho_current, pres_current, ss_current;
     
  for (int i=0; i<Nrho; ++i)
       {
	   rho_current = pow(10,log10(rho0) + i*dlogrho);
	   dens_array[i] = rho_current;
       }
  for (int j=0; j<Npres; ++j)
  
  {		  
      pres_current = pow(10,log10(pres0) + j*dlogpres);
		  
      pres_array[j] = pres_current;
	       
  }
  for(int i=0;i<Nrho; ++i){

      for(int j=0; j<Npres; ++j){
	       
           fscanf(fp_sound_speed_D2, "%lf", &ss_current);
           gsl_spline2d_set(spline_sound_speed_D2, sound_speed_D2, i, j, ss_current);
	       
      }
      
     }
   
  gsl_spline2d_init(spline_sound_speed_D2, dens_array, pres_array, sound_speed_D2, Nrho, Npres);
 
   
  delete[] dens_array;
  delete[] pres_array;
  delete[] sound_speed_D2;
//       make_ionization_table_Deuterium();
       make_ionization_dissociation_table_Deuterium2();
}

double SahaDeuterium2EOS::getElectricConductivity(double pressure, double density) {

    double cond;
    double T = getTemperature(pressure,density);
    double f_d = getDissociationFraction_Deuterium2(pressure,density);
    double m = getIonizationFraction_Deuterium2(pressure,density);
    double nt,fi,lnL,f_one;   
    nt =  density/3.34e-24;
    fi = m*f_d;
    if(fi < 1.e-20)
        cond = 0.;
    else{
        lnL = log(3.6e9*pow(T,1.5)/sqrt(fi*nt));
        cond = 1e11/1.15/(lnL/pow(T,1.5) + 0.054*pow(T,-0.059)*(1.0/fi-1));

        
        }
      
    return cond;
}

double SahaDeuterium2EOS::getEnergy(double pressure, double density) {
  double energy = 0.0;

  return energy;
}


double SahaDeuterium2EOS::getSoundSpeed(double pressure, double density) {
/*  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double sc;


  adjustStates(&pressure,&density,1.e-8,199.9,1.e-10,0.49);

    sc =  gsl_spline2d_eval(spline_sound_speed_D2, density, pressure, densacc, presacc);
*/
  //  printf("rho = %e  P = %e  sc = %e\n",density,pressure,sc);
    
  double sc;

 double gamma = 5.0/3.;

 double alpha_i = 3./2; 
 double alpha_d = 0.327;

 double e_d = 4.48; //eV
 double e_i = 13.6; //eV

 double gamma_eff;    
   
  adjustStates(&pressure,&density,1.e-8,199.9,1.e-10,0.49);
  double f_d = getDissociationFraction_Deuterium2(pressure,density);
  double f_i = getIonizationFraction_Deuterium2(pressure,density);
  f_i = f_i*f_d; 
 double T = getTemperature(pressure,density); 
  double beta_d = alpha_d + e_d/T;
  double beta_i = alpha_i + e_i/T;
   double phi_d = 0.5*f_d*(1-f_d)/(2-f_d);
             
   double phi_i = f_i*(1-f_i)/(2-f_i);
   
   double m = 0.5*(1-f_i);
                                               
   double a = f_d + f_i;

    gamma_eff = (m+a)*pow((1+(phi_d*beta_d + phi_i*beta_i)/(m+a)),2)/((1/(gamma-1))*m+3*a/2+phi_d*pow(beta_d,2)+phi_i*pow(beta_i,2)) - (phi_d+phi_i)/(m+a) + 1;	    
    sc = sqrt(gamma_eff*pressure/density);
   return sc;

}

double SahaDeuterium2EOS::getTemperature(double pressure, double density) {
        
    double mu = 2.014; 
    double R = 83.14;
    double fd = getDissociationFraction_Deuterium2(pressure,density);
    double fi = getIonizationFraction_Deuterium2(pressure,density);
    fi *= fd;
    double T = mu*pressure/((0.5*(1+fd)+fi)*density*R)/11604.525;
    return T;

}
static void adjustStates(double* pres, double* dens,double p0, double p1, double rho0, double rho1){

       if(*pres > p1)
           *pres = p1;
       else if(*pres < p0)
           *pres = p0;

       if(*dens > rho1)
           *dens = rho1;
       else if (*dens < rho0)
           *dens = rho0;


}

static void make_ionization_table_Deuterium(){
       double pres0 = 1.0e-8;
       double presend = 2.0e2;
       double rho0 = 1.0e-10;
       double rhoend = 0.5;
       int Nrho = 1000, Npres = 1000;
       double *dens_array = new double[Nrho];
       double *pres_array = new double[Npres];
       double *elec_frac_array_D = new double[Nrho*Npres];
       double dlogpres = (log10(presend) - log10(pres0))/(Npres-1);
       double dlogrho = (log10(rhoend) - log10(rho0))/(Nrho-1);
       printf("ENTERED SahaDeuterium  EOS:: Make Ionization Table \n");
       char* buffer;
       FILE* fp_elec_frac_D; 
       buffer = getcwd(NULL,0);
  string sub_dir = "/tables_1000/Deuterium/";
  string main_dir(buffer);
  string dir = main_dir+sub_dir;
      
  if (!(fp_elec_frac_D = fopen((dir+"elec_frac_TGAS.dat").c_str(),"r")))
    printf("CANNOT OPEN SahaDeuterium EOS elec_frac_TGAS_D2.dat file \n");
  
  spline_elec_frac_D = gsl_spline2d_alloc(T_elec_frac_D, Nrho, Npres); 
   
     
  double rho_current, pres_current, frac_current_D, frac_current_D2;
     
  for (int i=0; i<Nrho; ++i)
       {
	   rho_current = pow(10,log10(rho0) + i*dlogrho);
	   dens_array[i] = rho_current;
       }
  for (int j=0; j<Npres; ++j)
  
  {		  
      pres_current = pow(10,log10(pres0) + j*dlogpres);
		  
      pres_array[j] = pres_current;
	       
  }
  for(int i=0;i<Nrho; ++i){

      for(int j=0; j<Npres; ++j){
	       
           fscanf(fp_elec_frac_D, "%lf", &frac_current_D);
          gsl_spline2d_set(spline_elec_frac_D, elec_frac_array_D, i, j, frac_current_D);
	       
      }
      
     }
   
 gsl_spline2d_init(spline_elec_frac_D, dens_array, pres_array, elec_frac_array_D, Nrho, Npres);
   
  
  
  fcloseall();
  delete[] dens_array;
  delete[] pres_array;
  delete[] elec_frac_array_D;

}

static void make_ionization_dissociation_table_Deuterium2(){
       double pres0 = 1.0e-8;
       double presend = 2.0e2;
       double rho0 =  1.0e-10;
       double rhoend = 0.5;
       int Nrho = 1000, Npres = 1000;
       double *dens_array = new double[Nrho];
       double *pres_array = new double[Npres];
       double *diss_frac_array_D2 = new double[Nrho*Npres];
  
       double *elec_frac_array_D2 = new double[Nrho*Npres];
       double dlogpres = (log10(presend) - log10(pres0))/(Npres-1);
       double dlogrho = (log10(rhoend) - log10(rho0))/(Nrho-1);
       printf("ENTERED SahaDeuterium2  EOS:: Make Dissociation and Ionization table \n");
       char* buffer;
       FILE* fp_diss_frac_D2;
       
       FILE* fp_elec_frac_D2;
       buffer = getcwd(NULL,0);
  string sub_dir = "/tables_1000/Deuterium2/";
  string main_dir(buffer);
  string dir = main_dir+sub_dir;
      
  if (!(fp_diss_frac_D2 = fopen((dir+"diss_frac_TGAS_D2.dat").c_str(),"r")))
    printf("CANNOT OPEN SahaDeuterium2 EOS diss_frac_TGAS_for_D2.dat file \n");
 
  
  if (!(fp_elec_frac_D2 = fopen((dir+"elec_frac_TGAS_D2.dat").c_str(),"r")))
    printf("CANNOT OPEN SahaDeuterium2 Gas EOS elec_frac_TGAS_for_D2.dat file \n");
 

   
  spline_diss_frac_D2 = gsl_spline2d_alloc(T_diss_frac_D2, Nrho, Npres); 
     
  spline_elec_frac_D2 = gsl_spline2d_alloc(T_elec_frac_D2, Nrho, Npres); 
  
  double rho_current, pres_current, frac_current_D2,ion_frac_current_D2;
     
  for (int i=0; i<Nrho; ++i)
       {
	   rho_current = pow(10,log10(rho0) + i*dlogrho);
	   dens_array[i] = rho_current;
       }
  for (int j=0; j<Npres; ++j)
  
  {		  
      pres_current = pow(10,log10(pres0) + j*dlogpres);
		  
      pres_array[j] = pres_current;
	       
  }
  for(int i=0;i<Nrho; ++i){

      for(int j=0; j<Npres; ++j){
	       
           fscanf(fp_diss_frac_D2, "%lf", &frac_current_D2);
           gsl_spline2d_set(spline_diss_frac_D2, diss_frac_array_D2, i, j, frac_current_D2);
	       
           fscanf(fp_elec_frac_D2, "%lf", &frac_current_D2);
           gsl_spline2d_set(spline_elec_frac_D2, elec_frac_array_D2, i, j, frac_current_D2);
          
      }
      
     }
   
  gsl_spline2d_init(spline_diss_frac_D2, dens_array, pres_array, diss_frac_array_D2, Nrho, Npres);
   
  
  gsl_spline2d_init(spline_elec_frac_D2, dens_array, pres_array, elec_frac_array_D2, Nrho, Npres);
  
  fcloseall();
  delete[] dens_array;
  delete[] pres_array;
  delete[] diss_frac_array_D2;

  delete[] elec_frac_array_D2;
}

/*

static double getTemperature_test(double pressure, double density) {
        
    double mu = 2.014; 
    double R = 83.14;
    double fi = getIonizationFraction_Deuterium2(pressure,density);
    double fd = getDissociationFraction_Deuterium2(pressure,density);
    double T = mu*pressure/((0.5*(1+fd)+fi)*density*R)/11604.525;
    return T;

}
static void makesoundspeedtable() {
       double pres0 = 1.0e-8;
       double presend = 2.0e2;
       double rho0 =  1.0e-10;
       double rhoend = 0.5;
       int Nrho = 1000, Npres = 1000;
       double rho_current,pres_current; 
        double dlogpres = (log10(presend) - log10(pres0))/Npres;
       double dlogrho = (log10(rhoend) - log10(rho0))/Nrho;
 
       double gamma = 5.0/3.;

 double alpha_i = 3./2; 
 double alpha_d = 0.327;

 double e_d = 4.48; //eV
 double e_i = 13.6; //eV

 double gamma_eff;    
 string PT = "sound_speed_for_D2.dat";
 FILE * file1 = fopen(PT.c_str(),"w");

 for (int i=0; i<Nrho; ++i)
       {
	   rho_current = pow(10,log10(rho0) + i*dlogrho);
      
 
       for (int j=0; j<Npres; ++j)
  
  {		  
      pres_current = pow(10,log10(pres0) + j*dlogpres);

      double f_i = getIonizationFraction_Deuterium2(pres_current,rho_current);
      double f_d = getDissociationFraction_Deuterium2(pres_current,rho_current);
      
      double T = getTemperature_test(pres_current,rho_current); 
      double beta_d = alpha_d + e_d/T;
      double beta_i = alpha_i + e_i/T;
   double phi_d = 0.5*f_d*(1-f_d)/(2-f_d);
             
   double phi_i = f_i*(1-f_i)/(2-f_i);
   
   double m = 0.5*(1-f_i);
                                               
   double a = f_d + f_i;

    gamma_eff = (m+a)*pow((1+(phi_d*beta_d + phi_i*beta_i)/(m+a)),2)/((1/(gamma-1))*m+3*a/2+phi_d*pow(beta_d,2)+phi_i*pow(beta_i,2)) - (phi_d+phi_i)/(m+a) + 1;	    
  
    double ss = sqrt(gamma_eff*pres_current/rho_current);
  //  cout<<ss<<endl; 
    fprintf(file1,"%.16g     ",ss);
    }
    

       }
    fclose(file1);          
}
*/
void SahaDeuteriumEOS::getParameters(std::vector<double>& params){

       params.push_back(m_fGamma);
       printf("ENTERED  SahaDeuterium EOS::getParameters \n");
       make_ionization_table_Deuterium();
}

double SahaDeuteriumEOS::getElectricConductivity(double pressure, double density) {

    double cond;
    double T = getTemperature(pressure,density);
    double f_d = 1;
    double m = getIonizationFraction_Deuterium(pressure,density);
    double nt,fi,lnL,f_one;   
    nt =  density/3.34e-24;
    fi = m*f_d;
    if(fi < 1.e-20)
        cond = 0.;
    else{
        lnL = log(3.6e9*pow(T,1.5)/sqrt(fi*nt));
        cond = 1e11/1.15/(lnL/pow(T,1.5) + 0.054*pow(T,-0.059)*(1.0/fi-1));
        
        }
      
    return cond;
}
double SahaDeuteriumEOS::getEnergy(double pressure, double density) {
  double energy = 0.0;

  return energy;
}
double SahaDeuteriumEOS::getTemperature(double pressure, double density) {
        
    double mu = 2.014; 
    double R = 83.14;
    double fd = 1;
    double fi = getIonizationFraction_Deuterium(pressure,density);
    fi *= fd;
    double T = mu*pressure/((0.5*(1+fd)+fi)*density*R)/11604.525;
    return T;

}

double SahaDeuteriumEOS::getSoundSpeed(double pressure, double density) {
  double sc;

 double gamma = 5.0/3.;

 double alpha_i = 3./2; 
 double alpha_d = 0.327;

 double e_d = 4.48; //eV
 double e_i = 13.6; //eV

 double gamma_eff;    
  
  double f_d = 1;
  double f_i = getIonizationFraction_Deuterium(pressure,density);
  f_i = f_i*f_d; 
 double T = getTemperature(pressure,density); 
  double beta_d = alpha_d + e_d/T;
  double beta_i = alpha_i + e_i/T;
   double phi_d = 0.5*f_d*(1-f_d)/(2-f_d);
             
   double phi_i = f_i*(1-f_i)/(2-f_i);
   
   double m = 0.5*(1-f_i);
                                               
   double a = f_d + f_i;

    gamma_eff = (m+a)*pow((1+(phi_d*beta_d + phi_i*beta_i)/(m+a)),2)/((1/(gamma-1))*m+3*a/2+phi_d*pow(beta_d,2)+phi_i*pow(beta_i,2)) - (phi_d+phi_i)/(m+a) + 1;	    
    sc = sqrt(gamma_eff*pressure/density);
   return sc;

}

