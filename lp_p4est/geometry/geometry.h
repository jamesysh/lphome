#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <unordered_map>
#include <string>
class Geometry {
public:
    /// Destructor
	virtual ~Geometry() {}
	virtual bool operator()(double x, double y, double z) const=0; 
  	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax)=0;
};	


class GeometryFactory {
public:
	/**
     * \brief Defines a function pointer pointing to a function which creates objects in the Geometry family 
	 */
	typedef Geometry* (*GeoCreateFunc)();
	
	/**
     * \brief   Returns reference to a Singleton object of class GeometryFactory
	 * 
	 * \param   None 
     *
	 * \return  Reference to a static object of class GeometryFactory 
	 * 
	 *
	 * Example usage: 
	 * \code
	 *          GeometryFactory& factory = GeometryFactory::instance();
	 * \endcode
	 *
	 * \note    This function is implemented based on the lazy Singleton design pattern; 
	 *          therefore, only one GeometryFactory instance is allowed in each program
	 */
	static GeometryFactory& instance(); 
	
	/**
     * \brief      Registers (links) the geometry name \e name 
	 *		       with the function \e func for creating objects in the Geometry family
	 *			   
	 *             After registration, \e name can be used as an argument in the createGeometry member function
	 *             for creating objects of the linked type in the Geometry family
	 *	           
	 *  
	 * \param [in] name the geometry name 
	 * \param [in] func the function pointer pointing to the function that creates objects of a specific type
	 *             in the Geometry family
	 * 
	 * \return     None  
	 *
	 * \note       Instead of using this function directly, consider using the GeometryRegistrar class for
	 *		       the purpose of linking a geometry name and a specific class in the Geometry family. 
	 *             The function is kept public in case one wants to use it directly
	 *		  
	 *
	 */
	void registerGeometry(std::string name, GeoCreateFunc func);
	
	/**
     * \brief      This function creates an object of the class linked to the \name  
	 *
	 * \param [in] name the name linked to a specific class in the Geometry family via the 
	 *				    registrerGeometry member function
	 *  
	 * \return     A Geometry * pointer pointing to an object of a specific class in the Geometry family 
	 * 
	 * Example usage: 
	 * \code
	 *            GeometryFactory& factory = GeometryFactory::instance();
	 *            Geometry* newGeometry = factory.createGeometry(name);
	 * \endcode
	 *
	 */
	Geometry* createGeometry(std::string name);
private:
	std::unordered_map<std::string,GeoCreateFunc> geoTable; ///< hash table for the (name,creatFunction) pair
	GeometryFactory() {} ///< for singleton design pattern
	GeometryFactory(const GeometryFactory& other); ///< for singleton design pattern (Don't implement)
	GeometryFactory& operator=(const GeometryFactory& other); ///< for singleton design pattern (Don't implement)
};


#endif
