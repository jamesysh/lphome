#ifndef __INITIALIZER_H__
#define __INITIALIZER_H__

#include "eos.h"
#include <p4est.h>
#include <p8est.h>

#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p8est_bits.h>
#include <p8est_build.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>

#include "geometry.h"
#include "state.h"
#include <sc.h>
#include <cassert>
#include <string>
using namespace std;
class Global_Data;

class Initializer {
    friend class Global_Data;

    public:
    

	Initializer(const std::string& inputfileName, bool ifDebug=false, const std::string& debugfileName="debug");	
	
	/**
	 * \brief Constructor for restart 
	 * 
	 * 1. Reads the parameter file and data file\n 
	 * 2. Initializes the geometry of fluid objects\n
     * 3. Initializes the state of fluid objects 
	 */	
	Initializer(const std::string& param_filename, const std::string& data_filename, 
	bool ifDebug=false, const std::string& debugfileName="debug");	
	
    void adjustParamsFromDataFile(const string& datafileName);
    
        ~Initializer(){};


	bool getIfDebug() const {return m_bIfDebug;}

	bool getIfRestart() const {return m_bIfRestart;}

    string getDebugfileName() const {return m_sDebugfileName;}
	string getPelletFile() const{return m_sPelletfile};
	double getStartTime() const {return m_fStartTime;}
	
	double getEndTime() const {return m_fEndTime;}

	double getWriteTimeInterval() const {return m_fWriteTimeInterval;}
   
   
	int getWriteStep() const {return m_iWriteStep;}
   
	double getCFLCoeff() const {return m_fCFLCoeff;}
   
	int getDimension() const {return m_iDimension;}	
  
    double getDomainlength() const{return domain_len;}
	bool getRandomDirSplitOrder() const {return m_bRandomDirSplitOrder;}
   
         
	int getLPFOrder() const {return m_iLPFOrder;}	
   
	int getEOSChoice() const {return m_iEOSChoice;} //TODO
   
	double getInitParticleSpacing() const {return m_fInitParticleSpacing;}
   
	double getGravity() const {return m_fGravity;}
   
	bool getUseLimiter() const {return m_iUseLimiter;}
   
    int getInitialPerturbation() const {return m_fInitialPerturbation;}
   
	size_t getNumRow2ndOrder() const {return m_iNumRow2ndOrder;}
   
	size_t getNumRow1stOrder() const {return m_iNumRow1stOrder;}

	size_t getNumCol2ndOrder() const {return m_iNumCol2ndOrder;}

	size_t getNumCol1stOrder() const {return m_iNumCol1stOrder;}

	double getInvalidPressure() const {return m_fInvalidPressure;}
   
	double getInvalidDensity() const {return m_fInvalidDensity;}
   
	double getTimeSearchRadius() const {return m_fTimeSearchRadius;}

    int getInitLevel() const{return initlevel;}
    int getMinLevel() const{return minlevel;}
    int getMaxLevel() const{return maxlevel;}
    int getElemParticle() const{return elem_particle;} 
  
    size_t getBoundaryNumber() const{return m_iBoundaryObjNum;}
    vector<string>& getBoundaryNames() {return m_vBoundaryObjNames;}
    
    string getStateName(){return m_sState;}

    string getGeometryName(){return m_sGeometry;}
    
    double getGamma() const{return m_fGamma;}
    double getPinf() const{return m_fPinf;}
    double getEinf() const{return m_fEinf;}
    bool getIfFreeBoundary() const{return m_bIfFreeBoundary;}


    int getQuadtreeResolution() const {return m_iQuadtreeResolution;}
    int getBinarytreeResolution() const {return m_iBinarytreeResolution;}
    int getHeatingModel() const{return m_iHeatingModel;}
    double getMagneticField() const {return m_iMagneticField;}

	int getPelletDistribution() const {return m_iPelletDistribution;}
    int getPelletMaterial() const {return m_iPelletMaterial;}

    double getPelletRadius() const{return m_fPelletRadius;}
    double getTeinf() const{return m_fTeinf;}
    double getNeinf() const{return m_fNeinf;}

    int getPrintVelocity(){return m_iPrintVelocity;} 

   int getPrintVelocityU(){return m_iPrintVelocityU;}

   int getPrintVelocityV(){return m_iPrintVelocityV;}

   int getPrintVelocityW(){return m_iPrintVelocityW;}

   int getPrintDensity(){return m_iPrintDensity;}

   int getPrintMass(){return m_iPrintMass;}

   int getPrintPressure(){return m_iPrintPressure;}

   int getPrintSoundSpeed(){return m_iPrintSoundSpeed;}

   int getPrintLocalSpacing(){return m_iPrintLocalSpacing;}

   int getPrintTemperature(){return m_iPrintTemperature;}

   int getPrintAllParticle(){return m_iPrintAllParticle;}
    private:


	void readInputfile(const std::string& inputfileName);
    void setParams();

    int m_iPrintAllParticle = 0;
    int m_iPrintVelocity = 0;
    int m_iPrintVelocityU = 0;
    int m_iPrintVelocityV = 0;
    int m_iPrintVelocityW = 0;
    int m_iPrintDensity = 0;
    int m_iPrintMass = 0;
    int m_iPrintPressure = 0;
    int m_iPrintSoundSpeed = 0;
    int m_iPrintLocalSpacing = 0;
    int m_iPrintTemperature = 0;

    double domain_len; 
    int initlevel ; //init level of octree
    int maxlevel ;
    int minlevel;
    int elem_particle; //max number of particles per octant

    bool m_bIfDebug;
    bool m_bIfRestart;
    bool m_bIfFreeBoundary;
	double m_fStartTime;///< simulation start time
	double m_fEndTime; ///< simulation end time	
	double m_fWriteTimeInterval;///< write time interval
	int m_iWriteStep; ///< write step
	double m_fCFLCoeff;///< CFL coeff
	int m_iDimension;///< dimension
	string m_sFilenameSaveInit;
    string m_sDebugfileName;
    string m_sGeometry;
   
	string m_sState;
	string m_sPelletfile;
    size_t m_iBoundaryObjNum;///< number of boundary objects
    
    vector<std::string> m_vBoundaryObjNames; ///< Vector of boundary object names

    bool m_bRandomDirSplitOrder;///< if true then the order of directional splitting is randomly set 1:yes 0:no	

	std::size_t m_iNumRow2ndOrder;///< the smallest number of rows of A to solve 2nd order LPF
	std::size_t m_iNumRow1stOrder;///< the smallest number of rows of A to solve 1st order LPF
	std::size_t m_iNumCol2ndOrder;//TODO///< the number of columns of A when solving 2nd order LPF	
	std::size_t  m_iNumCol1stOrder;//TODO///< the number of columns of A when solving 1st order LPF

    double m_fTimeSearchRadius;///< the radius for neighbour search
	double m_fInvalidPressure;///< if p < invalid pressure => invalid state
	double m_fInvalidDensity;///< volume cannot be negative: if volume < invalid volume => invalid state	

    int m_iLPFOrder;///< the order of Local Polynomial Fitting (LPF)  	
	int m_iEOSChoice;///< choice of eos
	double m_fGamma;///< eos parameter gamma
	double m_fPinf;///< eos parameter pinf (stiffened poly gas) 
	double m_fEinf;///<eos parameter einf (stiffened poly gas) 	
	double m_fInitParticleSpacing;///< the initial particle spacing	
	double m_fGravity;///< gravity
	bool m_iUseLimiter;///< if use limiter or not 1:yes 0:no
	double m_fInitialPerturbation;//<amount of maximal initial perturbation in dx
    int m_iPelletDistribution; 
    int m_iPelletMaterial;
    int m_iQuadtreeResolution;
                        
    int m_iBinarytreeResolution;
    int m_iHeatingModel;
    double m_iMagneticField;
    double m_fPelletRadius;
    double m_fPelletOuterRadius;
    double m_fTeinf;
    double m_fNeinf;


};
#endif // __INITIALIZER_H__
