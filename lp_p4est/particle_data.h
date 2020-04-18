#ifndef __PARTICLE_DATA_H__
#define __PARTICLE_DATA_H__

#include "initializer.h"
#include "boundary.h"



typedef struct pdata{

    double xyz[3]; //coordinates
    double v[3]; //velocity
    double oldv[3];
    double pressure;
    double soundspeed;
    double temperature;
    double volume;
    double sphdensity;
    double pressureT1;
    double volumeT1;
    double soundspeedT1;
    
    double pressureT2;
    double volumeT2;
    double soundspeedT2;
    double mass;
    double localspacing;
    bool ifhasghostneighbour;
    bool ifboundary;   //if a boundary particle
    bool flagdelete;   // delete boundary particle at current timestep: delete if its not the same as flag in gdata
    int8_t schemeorder;
    int8_t redocount;
    int fluidneinum;
    sc_array_t * ghostneighbour;   //pdata_copy_t
    sc_array_t * neighbourparticle;
      
    sc_array_t * neighbourupparticle;
    sc_array_t * neighbourdownparticle;
    sc_array_t * neighbourrightparticle;
    sc_array_t * neighbourleftparticle;
    sc_array_t * neighbourfrontparticle;
    sc_array_t * neighbourbackparticle;
    p4est_gloidx_t      id;
    p4est_gloidx_t      mpirank;
    double leftintegral;
    double rightintegral;
    double deltaq;
    double qplusminus;
    double drifta;
    double driftv;
    int8_t layer;
    double radcool;
} pdata_t;

typedef struct pdata_copy{

    double xyz[3]; //coordinates
    double v[3]; //velocity
    double pressure;
    double soundspeed;
    double temperature;
    double volume;
    double mass;
    double localspacing;
    bool flagboundary;
    bool ifboundary; 

} pdata_copy_t;

/** Data type for payload data inside each quadrant */
typedef struct octant_data
{   
    bool flagboundary;  //if true, octant is at cloud boundary;
    sc_array_t *localneighbourid;
    sc_array_t *ghostneighbourid;
   // sc_array_t* particle_data_view;
    pdata_copy_t localparticle[1000];
    int octantid;
   // int mpirank;
    
    p4est_locidx_t    fluidnum; //number of fluid particle
    p4est_locidx_t    poctant;
  /** Offset into local array of all particles after this quadrant */
    p4est_locidx_t      lpend;
//    p4est_locidx_t ghostneighbourid[300];
  /** counts of local particles remaining on this quadrant and recieved ones */
  p4est_locidx_t      premain, preceive;
}
octant_data_t;

typedef enum pa_mode
{
  PA_MODE_REMAIN,
  PA_MODE_RECEIVE,
  PA_MODE_LOCATE
}
pa_mode_t;

typedef enum indata{
    PRESSURE,
    VELOCITY,
    VOLUME
}
indata_t;
typedef struct neighbour_info{
    
    float distance;
    float distance_weighted;
    size_t quadid;
    size_t parid;
    bool ifremote;        // if the particle is in remote processer
    bool ifghost;        // if the particle is a ghost particle       
    float phi;         //right and left
    float theta;          //up and down
    float sigma; //front and back
} neighbour_info_t;       

typedef struct comm_psend
{
  int                 rank;
  sc_array_t          message;     /** Message data to send */
}
comm_psend_t;

/** Array entry for a process that we send messages to */
typedef struct comm_prank
{
  int                 rank;
  comm_psend_t       *psend;        /**< Points to hash table entry */
}
comm_prank_t;

typedef enum comm_tag
{
  COMM_TAG_PART = P4EST_COMM_TAG_LAST,
  COMM_TAG_FIXED,
  COMM_TAG_CUSTOM,
  COMM_TAG_LAST
}
comm_tag_t;

class PelletSolver;

class Global_Data{

    public:
        friend PelletSolver; 
        Global_Data(Initializer* init);
        ~Global_Data(); 

        void * sc_array_index_begin (sc_array_t * arr);
        
        void sc_array_paste (sc_array_t * dest, sc_array_t * src);
        void adjustCoordByDomain( double xyz[3]);
        void initFluidParticles_distributed();
        void initFluidParticles_hexagonal();
        void initFluidParticles_datafile(string& filename);
        void prerun();
        void presearch();
        void presearch2d();
        void packParticles();
        void postsearch();
        void postsearch2d();
        void resetOctantData();
        void resetOctantData2d();
        void communicateParticles();
        void regroupParticles();
        void regroupParticles2d();
        void partitionParticles();

        void partitionParticles2d();
        void searchNeighbourOctant(); 
        void searchNeighbourOctant2d(); 
        void searchNeighbourParticle();
        void searchNeighbourParticle2d();
        void searchUpwindNeighbourParticle();
        void searchUpwindNeighbourParticle2d();
        void reorderNeighbourList();
        void reorderNeighbourList2d();
        void initParticleNeighbour();
        void setUpwindNeighbourList2d(sc_array_t* nei0, sc_array_t *nei1, sc_array_t *neidest);
        void setUpwindNeighbourList(sc_array_t *list0, sc_array_t *list1, sc_array_t *list2, sc_array_t *list3, sc_array_t *neidest);
        void copyParticle(pdata_copy_t* d, pdata_t *s);
        void generateGhostParticle();
        void generateGhostParticle2d();
        void generateGhostParticleForLW(pdata_t *pad, int dir, int count);
        void fillArrayWithGhostParticle(sc_array_t * array, pdata_t * pad, int count,int dir);
        
        void fillArrayWithGhostParticle2d(sc_array_t * array, pdata_t * pad, int count,int dir);
        
        void computeSPHDensity();
        void updateParticleStates();
        void updateViewForOctant(int phase);
        void updateViewForOctant2d(int phase);
        void addGhostParticle(pdata_copy_t * ghostnei, pdata_t *pad, double dx, double dy, double dz);
        void fetchNeighbourParticle(pdata_t* pad, pdata_copy_t **padnei ,sc_array_t *neighbourlist, size_t index);
        void fetchNeighbourParticle2d(pdata_t* pad, pdata_copy_t **padnei ,sc_array_t *neighbourlist, size_t index);
        void fetchParticle2d(pdata_t* pad, pdata_copy_t **padnei, neighbour_info_t *nei);
        void fetchParticle(pdata_t* pad, pdata_copy_t **padnei, neighbour_info_t *nei);
        void createViewForOctant();
        void createViewForOctant2d();
        void cleanForTimeStep();
        void cleanForTimeStep2d();
        
        void setParticleIDAndRank();
        void computeBoundingBoxLength();        
        void setFlagBoundaryForParticle();
        void testquad();
        void testquad2d();
        void loopquad (p4est_topidx_t tt, p8est_quadrant_t * quad,double lxyz[3], double hxyz[3], double dxyz[3]);
        
        void loopquad2d (p4est_topidx_t tt, p4est_quadrant_t * quad,double lxyz[3], double hxyz[3], double dxyz[3]);

        void split_by_coord ( sc_array_t * in,
        sc_array_t * out[2], pa_mode_t mode, int component,
        const double lxyz[3], const double dxyz[3]);
        void cleanUpArrays();
        void writeVTKFiles();
        void setEOS();
        void switchFlagDelete();
        void adjustLocalSpacingForParticle();
        void deleteBadStateParticle();
        

        Initializer *initializer;
        sc_MPI_Comm mpicomm;
        int dimension;
        int mpisize,mpirank;
        int initlevel;
        double timesearchingradius;
        int maxlevel;
        int minlevel; 
        int elem_particle; //max number of particles per octant
        int eoschoice;
        bool iffreeboundary;
        
        double gamma;
        double pinf;
        double einf;
        double invalidpressure;
        double invaliddensity;

        int LPFOrder;
        double initlocalspacing;
        double initperturbation;
        double gravity;
        bool ifuselimiter;
        size_t numrow1st;
        size_t numrow2nd;
        size_t numcol1st;
        size_t numcol2nd;
        
        double domain_len = 32; 
        double boundingbox_len; //max(y,z)
        double bb[6]; // bounding box of initial fluid particles  

        double lxyz[3],hxyz[3],dxyz[3]; //boundingbox of octant
        Geometry* geometry;
        State* state;        
        EOS* eos;        
	    
        std::vector<Boundary*> m_vBoundary;///< Vector of boundary objects
        size_t boundarynumber;
        bool flagdelete; 
        p4est_locidx_t lpnum; //number of particles on local processor
        p4est_gloidx_t gpnum, gplost; //number of particles on all processor, number of particles on all processers which left domain
        p4est_locidx_t lghostnum;
        p4est_gloidx_t gghostnum;
        p4est_locidx_t lfluidnum;
        p4est_gloidx_t gfluidnum;
        p4est_gloidx_t lboundarynum;
        p4est_locidx_t qremain, qreceive;
        int flagrefine, gflagrefine, flagstartrefine;
        octant_data_t *ghost_data; 
        p8est_ghost_t *ghost;
        p4est_ghost_t *ghost2d;
        sc_array_t *particle_data; //local particle data on process
        
       
       
        int pelletnumber;
        PelletSolver *pellet_solver; 
       
 //-------------------print option----------------------------
        int printallparticle;
        int printvelocity;
        int printdensity;
        int printmass;
        int printpressure;
        int printlocalspacing;
        int printsoundspeed;
        

        
        sc_array_t *pfound; //target process of particle
        
        sc_array_t *iremain; /**< locidx_t Index into padata of stay-local particles */

        sc_array_t *ireceive;/**< Index into particle receive buffer */

        sc_array_t *recevs;   /**< comm_prank_t with one entry per receiver, sorted */

        sc_array_t  *recv_req; /**< sc_MPI_Request for receiving */
         
        sc_array_t *send_req; /**< sc_MPI_Request for sending */

        sc_array_t *prebuf;  /**< pdata_t All received particles */

        sc_array_t *cfound;   /**< char Flag for received particles */

        sc_hash_t  *psend;    /**< comm_psend_t with one entry per receiver */
      
        sc_mempool_t *psmem;    /**< comm_psend_t to use as hash table entries */

      p4est_locidx_t      prevlp;   /**< Count local particles in partition callback */
      p4est_locidx_t      qcount;   /**< Count local quadrants in partition callback */
      sc_array_t         *src_fixed;        /**< int Particle counts per quadrant */
      sc_array_t         *dest_fixed;       /**< int Particle counts per quadrant */

        int octantid;  
        p4est_locidx_t      ireindex, ire2, ireindex2;   /**< Running index into iremain */

        p4est_locidx_t      irvindex, irv2, irvindex2;   /**< Running index into ireceive */
        sc_array_t *irecumu;   //cumulative remain particles count
        sc_array_t *irvcumu;
        sc_array_t *ilh[2],*jlh[2],*klh[2];
      //mesh data
        p8est_connectivity_t *conn;
      
        p8est_t            *p8est;
        
        p4est_connectivity_t *conn2d;
        p4est_t            *p4est;
      


};






#endif
