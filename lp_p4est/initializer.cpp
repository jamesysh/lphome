#include <iostream>
#include <vector>
#include "initializer.h"
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;



Initializer::Initializer(const string& inputfileName, bool ifDebug, const string& debugfileName)
: m_bIfDebug(ifDebug), m_sDebugfileName(debugfileName) {
	
    readInputfile(inputfileName); // read inputfile
	
    setParams();
	
	m_bIfRestart = false;
}

Initializer::Initializer(const string& inputfilename, const string& data_filename,
        bool ifDebug, const string& debugfileName){
        
	    m_bIfRestart = true;
        readInputfile(inputfilename); // read inputfile
        setParams();
        adjustParamsFromDataFile(data_filename);
        
}

void Initializer::readInputfile(const string& inputfileName) {
	
	ifstream ifs(inputfileName); // read the input file
    m_sFilenameSaveInit = "save_init_param"; // filename for saving info for restart

	vector<string> lines;
	string s;
	ofstream save(m_sFilenameSaveInit); // save info for restart
	
	while(getline(ifs, s)) {
		//ofs<<s<<endl; 
		lines.push_back(s);
	}

	istringstream iss;
	size_t i = 0; // input file line number
	
	iss.str(lines[i++]);
	iss>>m_fStartTime;
	save<<"m_fStartTime "<<m_fStartTime<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fEndTime;
	save<<"m_fEndTime "<<m_fEndTime<<endl;

	iss.str(lines[i++]);
	iss>>m_fWriteTimeInterval;
	save<<"m_fWriteTimeInterval "<<m_fWriteTimeInterval<<endl;

	iss.str(lines[i++]);
	iss>>domain_len;
	save<<"domain_len "<<domain_len<<endl;
	
    iss.str(lines[i++]);
	iss>>initlevel;
	save<<"initlevel "<<initlevel<<endl;
	
	iss.str(lines[i++]);
	iss>>minlevel;
	save<<"minlevel "<<minlevel<<endl;
    
	iss.str(lines[i++]);
	iss>>maxlevel;
	save<<"maxlevel "<<maxlevel<<endl;
  
    
	iss.str(lines[i++]);
	iss>>elem_particle;
	save<<"elem_particle "<<elem_particle<<endl;

    iss.str(lines[i++]);
	iss>>m_fCFLCoeff;
	save<<"m_fCFLCoeff "<<m_fCFLCoeff<<endl;

	iss.str(lines[i++]);
	iss>>m_iDimension;
	save<<"m_iDimension "<<m_iDimension<<endl;

	iss.str(lines[i++]);
	iss>>m_sGeometry;
	save<<"m_sGeometry "<<m_sGeometry<<endl;
 
	iss.str(lines[i++]);
	iss>>m_sState;
	save<<"m_sState "<<m_sState<<endl;
	
	
	iss.str(lines[i++]);
	iss>>m_iBoundaryObjNum;
	save<<"m_iBoundaryObjNum "<<m_iBoundaryObjNum<<endl;
		
        for(size_t j=0; j<m_iBoundaryObjNum; j++) {
			string tmpS;
		
			iss.str(lines[i++]);
			iss>>tmpS;
			m_vBoundaryObjNames.push_back(tmpS);
			save<<"m_vBoundaryObjectNames"<<" "<<m_vBoundaryObjNames[j]<<endl;
		}

	iss.str(lines[i++]);
	iss>>m_bIfFreeBoundary;
	save<<"m_bIfFreeBoundary "<<m_bIfFreeBoundary<<endl;
	
    
    iss.str(lines[i++]);
	iss>>m_bRandomDirSplitOrder;
	save<<"m_bRandomDirSplitOrder "<<m_bRandomDirSplitOrder<<endl;

	iss.str(lines[i++]);
	iss>>m_iLPFOrder;
	save<<"m_iLPFOrder "<<m_iLPFOrder<<endl;

	iss.str(lines[i++]);
	iss>>m_iEOSChoice;
	save<<"m_iEOSChoice "<<m_iEOSChoice<<endl;

	iss.str(lines[i++]);
	iss>>m_fGamma;
	save<<"m_fGamma "<<m_fGamma<<endl;

	iss.str(lines[i++]);
	iss>>m_fPinf;
	save<<"m_fPinf "<<m_fPinf<<endl;

	iss.str(lines[i++]);
	iss>>m_fEinf;
	save<<"m_fEinf "<<m_fEinf<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fInitParticleSpacing;
	save<<"m_fInitParticleSpacing "<<m_fInitParticleSpacing<<endl;	

	iss.str(lines[i++]);
	iss>>m_fGravity;
	save<<"m_fGravity "<<m_fGravity<<endl;
    
	iss.str(lines[i++]);
	iss>>m_iUseLimiter;
	save<<"m_iUseLimiter "<<m_iUseLimiter<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fInitialPerturbation;
	save<<"m_fInitialPerturbation "<<m_fInitialPerturbation<<endl;
	
    iss.str(lines[i++]);
	iss>>m_iPelletDistribution;
	save<<"m_iPelletDistribution "<<m_iPelletDistribution<<endl;
      

    if(m_iPelletDistribution){

        iss.str(lines[i++]);
	    iss>>m_iPelletMaterial;
	    save<<"m_iPelletMaterial"<<m_iPelletMaterial<<endl;
       /* 
        iss.str(lines[i++]);
	    iss>>m_fPelletRadius;
	    save<<"m_fPelletRadius"<<m_fPelletRadius<<endl;
        iss.str(lines[i++]);
	    iss>>m_fPelletOuterRadius;
	    save<<"m_fPelletOuterRadius "<<m_fPelletOuterRadius<<endl;
        */

        iss.str(lines[i++]);
	    iss>>m_sPelletfile;
	    save<<"pelletfile "<<m_sPelletfile<<endl;

        iss.str(lines[i++]);
	    iss>>m_fTeinf;
        
        iss.str(lines[i++]);
	    iss>>m_fNeinf;
        
        iss.str(lines[i++]);
	    iss>>m_iHeatingModel;
	    save<<"m_iHeatingModel "<<m_iHeatingModel<<endl;
        
        if(m_iHeatingModel == 1){
        iss.str(lines[i++]);
	    iss>>m_iQuadtreeResolution;
	    save<<"m_iQuadtreeResolution "<<m_iQuadtreeResolution<<endl;
        
        iss.str(lines[i++]);
        iss>>m_iBinarytreeResolution;
	    save<<"m_iBinarytreeResolution "<<m_iBinarytreeResolution<<endl;
        
        iss.str(lines[i++]);
        iss>>m_iMagneticField;
        save<<"m_iMagneticField "<<m_iMagneticField<<endl;
        
        P4EST_GLOBAL_ESSENTIALF ("Magnetic Field is %f Tesla.\n",  m_iMagneticField/10.);
        }


    }
//-------------------SETUP OUTPUT OPTION------------------------
  
    vector<string> output;
    
    output.assign(lines.begin()+i,lines.end());
    vector<string>::iterator result;
    
    result = find(output.begin(),output.end(),"velocity");
    if(result != output.end())
        m_iPrintVelocity = 1;

    result = find(output.begin(),output.end(),"allparticle");
    if(result != output.end())
        m_iPrintAllParticle = 1;

    
    
    result = find(output.begin(),output.end(),"velocityu");
    if(result != output.end())
        m_iPrintVelocityU = 1;
    result = find(output.begin(),output.end(),"velocityv");

    if(result != output.end())
        m_iPrintVelocityV = 1;
   
    result = find(output.begin(),output.end(),"velocityw");

    if(result != output.end())
        m_iPrintVelocityW = 1;
    

    result = find(output.begin(),output.end(),"density");

    if(result != output.end())
        m_iPrintDensity = 1;
 

    result = find(output.begin(),output.end(),"mass");

    if(result != output.end())
        m_iPrintMass = 1;
   
    result = find(output.begin(),output.end(),"pressure");

    if(result != output.end())
        m_iPrintPressure = 1;
    
    result = find(output.begin(),output.end(),"soundspeed");

    if(result != output.end())
        m_iPrintSoundSpeed = 1;
  //pellet case
     
    result = find(output.begin(),output.end(),"localspacing");

    if(result != output.end())
        m_iPrintLocalSpacing = 1;
     
    result = find(output.begin(),output.end(),"temperature");

    if(result != output.end())
        m_iPrintTemperature = 1;
    
      
}

void Initializer::adjustParamsFromDataFile(const string& datafileName) {

	ofstream save(m_sFilenameSaveInit, ofstream::app); // save init info
    
	ifstream ifs(datafileName);
    ifs.precision(16);	
    
    double temp;
    ifs>>temp;
	m_fStartTime = temp;
    cout<<"(RESTART CHANGE): m_fStartTime has been changed to "<<m_fStartTime<<endl;
    m_iWriteStep = (int)(m_fStartTime/m_fWriteTimeInterval);

    }



void Initializer::setParams() {
			
	ofstream save(m_sFilenameSaveInit, ofstream::app); 

	m_iWriteStep = 0; // The default value
	if(m_iDimension==2)
		m_fTimeSearchRadius = 4.;//4
    if(m_iDimension==3)
            m_fTimeSearchRadius = 3.5;//4

	if(m_iDimension==3) {
		m_iNumRow2ndOrder = 53; 
		m_iNumRow1stOrder = 8; 
		m_iNumCol2ndOrder = 9; 
		m_iNumCol1stOrder = 3;

	}
	else if(m_iDimension==2) {
				
		m_iNumRow2ndOrder = 36;//36 
		m_iNumRow1stOrder = 3; 
		m_iNumCol2ndOrder = 5; 
		m_iNumCol1stOrder = 2;	
	}
	// (eos chocie=1:poly 2:spoly)	
	if(m_iEOSChoice==1) 
		m_fInvalidPressure = 0;
	else if(m_iEOSChoice==2) 	
		m_fInvalidPressure = -m_fPinf;
	m_fInvalidDensity = 0;
	
	// save info for restart
	save<<"m_iWriteStep "<<m_iWriteStep<<endl;
	save<<"m_fTimesNeiSearchRadius "<<m_fTimeSearchRadius<<endl;
	save<<"m_fInvalidPressure "<<m_fInvalidPressure<<endl;
	save<<"m_fInvalidDensity "<<m_fInvalidDensity<<endl;
	
}




