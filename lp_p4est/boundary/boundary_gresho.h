#ifndef __BOUNDARY_GRESHO_H__
#define __BOUNDARY_GRESHO_H__

#include "boundary.h"
#include "eos.h"
class Gresho2DSolidBoundary: public Boundary {
public:
	/// constructor
	Gresho2DSolidBoundary();

	/// destructor
	virtual ~Gresho2DSolidBoundary() {}
	
	/**
	 * \brief Get a boundary particle based on a fluid particle      
	 * \param [in] x  The x-coordinate of fluid particle
	 * \param [in] y  The y-coordinate of fluid particle
	 * \param [in] z  The z-coordinate of fluid particle
	 * \param [out] xb  The x-coordinate of boundary particle
	 * \param [out] yb  The y-coordinate of boundary particle
	 * \param [out] zb  The z-coordinate of boundary particle		
	 */
	virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz, 
	double& xb, double& yb, double& zb, 
	double& pressureb, double& vxb, double& vyb, double& vzb);
		
    virtual void generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx,double dt);
		
private:
	double radius;	
	double thickness;
	double bo; 
	
};

class Yee3DSolidBoundary: public Boundary {
public:
        /// constructor
        Yee3DSolidBoundary();

        /// destructor
        virtual ~Yee3DSolidBoundary() {}

        /**
         * \brief Get a boundary particle based on a fluid particle      
         * \param [in] x  The x-coordinate of fluid particle
         * \param [in] y  The y-coordinate of fluid particle
         * \param [in] z  The z-coordinate of fluid particle
         * \param [out] xb  The x-coordinate of boundary particle
         * \param [out] yb  The y-coordinate of boundary particle
         * \param [out] zb  The z-coordinate of boundary particle               
         */
        int operator()(double x, double y, double z, double pressure, double volume, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& volumeb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);

        virtual void generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx,double dt);

private:
        double radius;
        double thickness;
        double bo;
	double zlength;
	double zmino,zmaxo,zmin,zmax;
        double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

};

class Yee2DSolidBoundary: public Boundary {
public:
        /// constructor
        Yee2DSolidBoundary();

        /// destructor
        virtual ~Yee2DSolidBoundary() {}

        /**
         * \brief Get a boundary particle based on a fluid particle      
         * \param [in] x  The x-coordinate of fluid particle
         * \param [in] y  The y-coordinate of fluid particle
         * \param [in] z  The z-coordinate of fluid particle
         * \param [out] xb  The x-coordinate of boundary particle
         * \param [out] yb  The y-coordinate of boundary particle
         * \param [out] zb  The z-coordinate of boundary particle               
         */
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        double& xb, double& yb, double& zb,
        double& pressureb, double& vxb, double& vyb, double& vzb);

        virtual void generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx,double dt);

private:
        double radius;
        double thickness;
        double bo;

};

#endif
