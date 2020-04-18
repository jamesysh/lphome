
#ifndef __STATE_H__
#define __STATE_H__

#include <unordered_map>
#include <string>
class State {
public:
	/**
	 * \brief      Calculates pressure based on the Cartesian coordinate (x,y,z) of a particle 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return     The calculated pressure value
	 */
	virtual double pressure(double x, double y, double z)=0;
    /**
	 * \brief      Calculates density based on the Cartesian coordinate (x,y,z) of a particle 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return     The calculated density value
	 */
	virtual double density(double x, double y, double z)=0;

    /**
	 * \brief       Calculates velocity based on the Cartesian coordinate (x,y,z) of a particle 
	 * \param [in]  x  The x-coordinate
	 * \param [in]  y  The y-coordinate
	 * \param [in]  z  The z-coordinate
	 * \param [out] vX  The calculated velocity value in x-coordinate
	 * \param [out] vY  The calculated velocity value in y-coordinate
	 * \param [out] vZ  The calculated velocity value in z-coordinate
	 * \return None
	 */
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ)=0;
	/**
	 * \brief Destructor 
	 */
	virtual ~State() {};
};



class StateFactory {
public:
	/**
     * \brief Defines a function pointer pointing to a function which creates objects in the State family 
	 */
	typedef State* (*StateCreateFunc)();
	
	/**
     * \brief   Returns reference to a Singleton object of class StateFactory
	 * 
	 * \param   None 
     *
	 * \return  Reference to a static object of class StateFactory 
	 * 
	 *
	 * Example usage: 
	 * \code
	 *          StateFactory& factory = StateFactory::instance();
	 * \endcode
	 *
	 * \note    This function is implemented based on the lazy Singleton design pattern; 
	 *          therefore, only one StateFactory instance is allowed in each program
	 */
	static StateFactory& instance(); 
	
	/**
     * \brief      Registers (links) the state name \e name 
	 *		       with the function \e func for creating objects in the State family
	 *			   
	 *             After registration, \e name can be used as an argument in the createState member function
	 *             for creating objects of the linked type in the State family
	 *	           
	 *  
	 * \param [in] name the state name 
	 * \param [in] func the function pointer pointing to the function that creates objects of a specific type
	 *             in the State family
	 * 
	 * \return     None  
	 *
	 * \note       Instead of using this function directly, consider using the StateRegistrar class for
	 *		       the purpose of linking a state name and a specific class in the State family. 
	 *             The function is kept public in case one wants to use it directly
	 *		  
	 *
	 */
	void registerState(std::string name, StateCreateFunc func);

	/**
     * \brief      This function creates an object of the class linked to the \name  
	 *
	 * \param [in] name the name linked to a specific class in the State family via the 
	 *				    registrerState member function
	 *  
	 * \return     A State * pointer pointing to an object of a specific class in the State family 
	 * 
	 * Example usage: 
	 * \code
	 *            StateFactory& factory = StateFactory::instance();
	 *            State* newState = factory.createState(name);
	 * \endcode
	 *
	 */
	State* createState(std::string name);	

private:	
	std::unordered_map<std::string, StateCreateFunc> stateTable; ///< hash table for the (name,creatFunction) pair	
	StateFactory() {}; ///< for singleton design pattern
	StateFactory(const StateFactory& other); ///< for singleton design pattern (Don't implement)
	StateFactory& operator=(const StateFactory& other); ///< for singleton design pattern (Don't implement)

};

#endif // __STATE_H__



