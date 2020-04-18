#include "state.h"
#include <iostream>
#include <cmath>
////////////////////////////////////////////////////////////////////////////////////////
// Start of StateFactory
////////////////////////////////////////////////////////////////////////////////////////

StateFactory& StateFactory::instance() { // singleton
	static StateFactory stateFactory;
	return stateFactory;
}

State* StateFactory::createState(std::string name) {
	const auto result = stateTable.find(name);
	if(result==stateTable.end()) {
		std::cout<<"This state class name is not registered!!!"<<std::endl;
		std::cout<<name<<std::endl;
		return nullptr;
	}
	return (result->second)();
}

void StateFactory::registerState(std::string name, StateCreateFunc func) {
	stateTable.insert({name,func});
}

////////////////////////////////////////////////////////////////////////////////////////
// End of StateFactory
////////////////////////////////////////////////////////////////////////////////////////

