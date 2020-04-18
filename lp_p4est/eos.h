
#ifndef __EOS_H__
#define __EOS_H__

#include <cassert>
#include <vector>
class EOS {
protected:
    
    double m_fGamma;	
    int m_iEOSChoice; ///< The eos choice: 1=Polytropic gas; 2=Stiffened Polytropic gas; 3=Saha Eos
public:
    
    /// Destructor
	virtual ~EOS() {};

	/// Getter function of the protected data member m_iEOSChoice
	int getEOSChoice() {return m_iEOSChoice;}
    double getGamma() {return m_fGamma;}	
	virtual void getParameters(std::vector<double>& params) = 0;
	virtual double getEnergy(double pressure, double density) = 0;
 	virtual double getTemperature(double pressure, double density) = 0;
      	virtual double getSoundSpeed(double pressure, double density) = 0;
	virtual double getElectricConductivity(double pressure, double density) = 0;
};





class PolytropicGasEOS : public EOS {
protected:
    
public:
	PolytropicGasEOS(double gamma) {m_iEOSChoice=1; m_fGamma = gamma;}
	
	// Destructor
	virtual ~PolytropicGasEOS() {}
	virtual void getParameters(std::vector<double>& params){params.push_back(m_fGamma);};
	virtual double getEnergy(double pressure, double density);
	virtual double getTemperature(double pressure, double density);
	virtual double getSoundSpeed(double pressure, double density);
	virtual double getElectricConductivity(double pressure, double density);

};

class SahaNeonEOS : public EOS {
protected:

public:
	/// Constructor
	SahaNeonEOS(double gamma)  
    {
        m_iEOSChoice=2;
        m_fGamma = gamma;
    }
	
	// Destructor
	virtual ~SahaNeonEOS() {}	
	virtual void getParameters(std::vector<double>& params);

	virtual double getEnergy(double pressure, double density);
	virtual double getTemperature(double pressure, double density);
	virtual double getSoundSpeed(double pressure, double density);
	virtual double getElectricConductivity(double pressure, double density);

};


class SahaDeuterium2EOS : public EOS {
protected:

public:
	/// Constructor
	SahaDeuterium2EOS(double gamma) {m_iEOSChoice = 3; m_fGamma = gamma;}
	
	// Destructor
	virtual ~SahaDeuterium2EOS() {}	
	virtual void getParameters(std::vector<double>& params);

	virtual double getEnergy(double pressure, double density);
	virtual double getTemperature(double pressure, double density);
	virtual double getSoundSpeed(double pressure, double density);
	virtual double getElectricConductivity(double pressure, double density);

};


class SahaDeuteriumEOS : public EOS {
protected:

public:
	/// Constructor
	SahaDeuteriumEOS(double gamma) {m_iEOSChoice = 4; m_fGamma = gamma;}
	
	// Destructor
	virtual ~SahaDeuteriumEOS() {}	
	virtual void getParameters(std::vector<double>& params);

	virtual double getEnergy(double pressure, double density);
	virtual double getTemperature(double pressure, double density);
	virtual double getSoundSpeed(double pressure, double density);
	virtual double getElectricConductivity(double pressure, double density);

};



#endif
