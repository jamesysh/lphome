#ifndef __MATERIAL_LIB__
#define __MATERIAL_LIB__


class Material{

    public:
        double  getMu() const {return mu;} 
        double getMass() const {return mass;}
        double getZ()const {return Z;}
        double getI() const{return I;}
        double getSublimationEnergy() const {return sublimationEnergy;}
        virtual double  getOne_Plus_Zstar(double teinf ) = 0;
    protected: 
        double mu;
        double mass;
        double Z;
        double I;
        double sublimationEnergy;
        double one_plus_Zstar;
};


class Neon: public Material{
    
    public:  
         virtual double  getOne_Plus_Zstar(double teinf ); 
        Neon();
    
        ~Neon();


};

class Deuterium: public Material{

    public:
        virtual double  getOne_Plus_Zstar(double teinf );
        Deuterium();
        ~Deuterium();

};

class Deuterium2: public Material{

    public:
        virtual double  getOne_Plus_Zstar(double teinf );
        Deuterium2();
        ~Deuterium2();

};

#endif
