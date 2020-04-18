#include "geometry.h"
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////
// Start of GeometryFactory
////////////////////////////////////////////////////////////////////////////////////////

GeometryFactory& GeometryFactory::instance() { 
	static GeometryFactory geoFactory;
	return geoFactory;
}

void GeometryFactory::registerGeometry(std::string name, GeoCreateFunc func) {
	
	geoTable.insert({name,func});

}

Geometry* GeometryFactory::createGeometry(std::string name) {
	const auto result = geoTable.find(name);
	if(result==geoTable.end()) { // the goemetry class name is not registered
		std::cout<<"This geometry class name is not registered!!!"<<std::endl;
		std::cout<<name<<std::endl;
		return nullptr;
	}
	return (result->second)();
}

////////////////////////////////////////////////////////////////////////////////////////
// End of GeometryFactory
////////////////////////////////////////////////////////////////////////////////////////

