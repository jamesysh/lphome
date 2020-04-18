#ifndef __PARTICLE_VIEWER_H__
#define __PARTICLE_VIEWER_H__

#include "particle_data.h"

class ParticleViewer {

public:
	/**
	 * \brief                      Constructor 
	 *
	 * \param [in] data            A pointer to the object which holds particle information and data
	 * \param [in] outputfileName  The name of the output file
	 * \param [in] numDigits       The number of digits for the indexing of output file name
	 *                
	 */
	ParticleViewer(Global_Data* data,  const std::string& outputfileName="", int numDigits=7);
	
	/**
	 * \brief virtual destructor
	 *
	 *
*/
     ~ParticleViewer() {}
	
	/**
	 * \brief                Write simulation results to the output file in the .vtk format
	 * \param [in] time      The physical output time
	 * \param [in] writeStep The number of times of output
	 * \return               0 if output success; 1 otherwise
	 */
     void writeResult(int step,double time);
     void writeGhost(int step);
     
     void writeTXTFile(int step);
     
     void writeComputationNode(int timestep);
     Global_Data *gdata;

	std::string  outputfilename;///< The name of the output file 
	int numdigit;///< The number of digits for the indexing of output file name 
    int writestep;
	std::string rightFlush( std::size_t numDigits);

};






#endif
