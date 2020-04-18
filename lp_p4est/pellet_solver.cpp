#include "pellet_solver.h"
#include <algorithm>
#include "sc_notify.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "geometry_pellet.h"
using namespace std;

#define c_light 2.99792458e7
float     Bessel_Kn( int n,float  x);
float     Bessel_I0(float  x);
float     Bessel_K1(float  x);
float     Bessel_I1(float  x);


PelletSolver::PelletSolver(Initializer *init,Global_Data*g){
    gdata = g;
    pelletradius = init->getPelletRadius();
    pelletfile = init->getPelletfile();
    readPelletFile();
    teinf = init->getTeinf();
    neinf = init->getNeinf();
    setPelletMaterial(init->getPelletMaterial());
    heatingmodel = init->getHeatingModel();
    if(heatingmodel){
        elem_particle_box = init->getQuadtreeResolution();
        elem_particle_cell = init->getBinarytreeResolution();
        magneticfield = init->getMagneticField();
    }
    gdata->pellet_solver = this;
    ((PelletLayer*)(gdata->geometry))->pellet_solver = this;    
    
}
void PelletSolver:: readPelletFile(){
	ifstream ifs(pelletfile); // read the pellet file
	if(!ifs){
		cout<<"Pellet File Does not Exist!"<<endl;
		exit(0);
	}

	vector<string> lines;
	string s;
	while(getline(ifs, s)) {
		//ofs<<s<<endl; 
		lines.push_back(s);
	}

	istringstream iss;
	size_t i = 0; // input file line number
	size_t pelletnumber;
	iss.str(lines[i++]);
	iss>>pelletnumber;
        pelletlist = sc_array_new_count(sizeof(pellet_info_t),pelletnumber);
	pellet_info_t *pellet;
	double temp;
	for(size_t j=0;j<pelletnumber;j++){
		iss.str(lines[i++]);
		pellet = (pellet_info_t*) sc_array_index(pelletlist,j);
		iss>>temp;
		pellet->xyz[0] = temp;
		iss>>temp;
		pellet->xyz[1] = temp;
		iss>>temp;
		pellet->xyz[2] = temp;
		iss>>temp;
		pellet->innerradius = temp;
		iss>>temp;
		pellet->outerradius = temp;
	}	
}

static int
comm_node_compare (const void *v1, const void *v2)
{
  computation_node_t *node1 = (computation_node_t *)v1;      

  computation_node_t *node2 = (computation_node_t *)v2;
  int t = 0;
  if(node1->id > node2->id)
      t = 1;
  else if(node1->id < node2->id)
      t = 0;
  else if(node1->xyz[0] > node2->xyz[0])
      t = 1;
  else 
      t = 0;
  return t;
}

static int
comm_prank_compare (const void *v1, const void *v2)
{
  return sc_int_compare (&((const comm_prank_t *) v1)->rank,
                         &((const comm_prank_t *) v2)->rank);
}
static unsigned
psend_hash (const void *v, const void *u)
{
  const comm_psend_t *ps = (const comm_psend_t *) v;


  return ps->rank;
}

static int
psend_equal (const void *v1, const void *v2, const void *u)
{
  const comm_psend_t *ps1 = (const comm_psend_t *) v1;
  const comm_psend_t *ps2 = (const comm_psend_t *) v2;


  return ps1->rank == ps2->rank;
}
static void adjustCoordByDomain( double xyz[3], double dl){
    for(int i=0;i<3;i++){
        xyz[i] *= dl;
        xyz[i] -= dl/2;
    }
    return;
}

static void loopquad2d (PelletSolver* p, p4est_topidx_t tt, p4est_quadrant_t * quad,double lxyz[3], double hxyz[3], double dxyz[3]){

    
  int                 i;
  p4est_qcoord_t      qh;
  qh = P4EST_QUADRANT_LEN (quad->level);
  p4est_qcoord_to_vertex (p->conn, tt, quad->x, quad->y,    lxyz);
  p4est_qcoord_to_vertex (p->conn, tt, quad->x + qh, quad->y + qh,  hxyz);
  adjustCoordByDomain(lxyz,p->gdata->boundingbox_len);
  adjustCoordByDomain(hxyz,p->gdata->boundingbox_len);
  for (i = 0; i < 2; ++i) {

    dxyz[i] = hxyz[i] - lxyz[i];
   
  }

}

static int
psearch_quad2d (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant, int pfirst, int plast,
              p4est_locidx_t local_num, void *point)
{
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data *g = p->gdata;
  /* compute coordinate range of this quadrant */
  loopquad2d (p, which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}

static int
psearch_point2d (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrant, int pfirst, int plast,
               p4est_locidx_t local_num, void *point)
{
  int                 i;
  int                *pfn;
  size_t              zp;
  double             *x;
  PelletSolver *p = (PelletSolver *)p4est->user_pointer;
  Global_Data      *g = p->gdata;
  quadrant_data_t          *qud;
  pdata_t          *pad = (pdata_t *) point;
  /* access location of particle to be searched */
  x = pad->xyz;

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P4EST_DIM; ++i) {
    if (!(g->lxyz[i] <= x[i] && x[i] <= g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  /* convention for entries of pfound:
     -1              particle has not yet been found
     [0 .. mpisize)  particle found on that rank, me or other
   */

  /* find process/quadrant for this particle */
  if (local_num >= 0) {
    /* quadrant is a local leaf */
    zp = sc_array_position (p->particle_data_copy, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* first local match counts (due to roundoff there may be multiple) */
    if (*pfn != g->mpirank) {
      /* particle was either yet unfound, or found on another process */
      /* bump counter of particles in this local quadrant */

      *pfn = g->mpirank;
      *(p4est_locidx_t *) sc_array_push (g->iremain) = (p4est_locidx_t) zp;
      qud = (quadrant_data_t *) quadrant->p.user_data;
      ++qud->premain;
    }
    /* return value will have no effect, but we must return */
    return 0;
  }
  if (pfirst == plast) {
    if (pfirst == g->mpirank) {
      /* continue recursion for local branch quadrant */
      return 1;
    }
    
    /* found particle on a remote process */
    zp = sc_array_position (p->particle_data_copy, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* only count match if it has not been found locally or on lower rank */
    if (*pfn < 0 || (*pfn != g->mpirank && pfirst < *pfn)) {
        *pfn = pfirst;
    }

    /* return value will have no effect, but we must return */
    return 0;
  }

  /* the process for this particle has not yet been found */
  return 1;
}

static int adapt_refine2d (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant)
{
  
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data *g = p->gdata;

  
  
  quadrant_data_t          *oud = (quadrant_data_t *) quadrant->p.user_data;
  
  

  
  
  /* we have set this to -1 in adapt_coarsen */

  if ((double) (oud->premain + oud->preceive ) > p->elem_particle_box){ //|| quadrant->level < 3) {
    /* we are trying to refine, we will possibly go into the replace function */
    g->ire2 = g->ireindex;
    g->ireindex += oud->premain;
    g->irv2 = g->irvindex;
    g->irvindex += oud->preceive;
    
    return 1;
  }
  else {

    /* maintain cumulative particle count for next quadrant */
    g->ireindex += oud->premain;
    g->irvindex += oud->preceive;
    
    p4est_locidx_t *irecumu = (p4est_locidx_t *)sc_array_push(g->irecumu);
    *irecumu = g->ireindex;

    p4est_locidx_t *irvcumu = (p4est_locidx_t *)sc_array_push(g->irvcumu);
    *irvcumu = g->irvindex;
    oud->quadrantid = g->octantid;
    g->octantid ++;
    return 0;
  }

}

void PelletSolver::split_by_coord ( sc_array_t * in,
                sc_array_t * out[2], pa_mode_t mode, int component,
                const double lxyz[3], const double dxyz[3])
{
  p4est_locidx_t      ppos;
  const double       *x;
  size_t              zz, znum;
  pdata_t          *pad;

  sc_array_truncate (out[0]);
  sc_array_truncate (out[1]);

  znum = in->elem_count;
  for (zz = 0; zz < znum; ++zz) {
    ppos = *(p4est_locidx_t *) sc_array_index (in, zz);
    if (mode == PA_MODE_REMAIN) {
      pad = (pdata_t *) sc_array_index (particle_data_copy, ppos);
      x = pad->xyz;
    }
    else if (mode == PA_MODE_RECEIVE) {
      pad = (pdata_t *) sc_array_index (prebuf, ppos);
      x = pad->xyz;
    }
    else {
      P4EST_ASSERT (mode == PA_MODE_LOCATE);
      pad = (pdata_t *) sc_array_index (particle_data_copy, ppos);
      x = pad->xyz;
    }
    if (x[component] <= lxyz[component] + .5 * dxyz[component]) {
      *(p4est_locidx_t *) sc_array_push (out[0]) = ppos;
    }
    else {
      *(p4est_locidx_t *) sc_array_push (out[1]) = ppos;
    }
  }
}
static void  adapt_replace2d (p4est_t * p4est, p4est_topidx_t which_tree,
               int num_outgoing, p4est_quadrant_t * outgoing[],
               int num_incoming, p4est_quadrant_t * incoming[]){

  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg;
  p4est_quadrant_t  **pchild;
  quadrant_data_t          *oud;
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data *g = p->gdata;
  if (num_outgoing == P4EST_CHILDREN) {
    // we are coarsening 
      oud = (quadrant_data_t *) incoming[0]->p.user_data;
    g->ireindex2 += (oud->premain = g->qremain);
    g->irvindex2 += (oud->preceive = g->qreceive);
  }

  else {
    // we are refining 
    // access parent quadrant 

    loopquad2d (p,which_tree, outgoing[0],lxyz,hxyz,dxyz);
    // recover window onto remaining particles for the new family 
    ibeg = g->ire2;
    irem = g->ireindex - ibeg;
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    p->klh[0] = &iview;
    wz = 0;

    p->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      p->split_by_coord (p->jlh[wy], p->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = p->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
      
      }
    }

    // recover window onto received particles for the new family 
    ibeg = g->irv2;
    irem = g->irvindex - ibeg;
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);

    // sort received particles into the children 
    pchild = incoming;
    wz = 0;
    p->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      p->split_by_coord ( p->jlh[wy], p->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = p->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->preceive = (p4est_locidx_t) arr->elem_count);
      }
    }
    p->klh[0] = NULL;
    assert (ibeg == g->irvindex);

  }

}
/*
static void adapt_replace2d (p4est_t * p4est, p4est_topidx_t which_tree,
               int num_outgoing, p4est_quadrant_t * outgoing[],
               int num_incoming, p4est_quadrant_t * incoming[]){
    
  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg;
  p4est_quadrant_t  **pchild;
  quadrant_data_t          *oud;
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data *g = p->gdata;
  
  if (num_outgoing == P4EST_CHILDREN) {
  
  // we are coarsening 
      oud = (quadrant_data_t *) incoming[0]->p.user_data;
    g->ireindex2 += (oud->premain = g->qremain);
    g->irvindex2 += (oud->preceive = g->qreceive);
  }

  else {
    // we are refining 
    // access parent quadrant 
    
    loopquad2d (p,which_tree, outgoing[0],lxyz,hxyz,dxyz);

    // recover window onto remaining particles for the new family 
    ibeg = g->ire2;
    irem = g->ireindex - ibeg;
    
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    p->klh[0] = &iview;
    wz = 0;

    g->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (p->jlh[wy], p->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
      
      }
    }

    p->klh[0] = NULL;

  }

}
*/
void PelletSolver::build_quadtree(){
    pdata_t *pad_s, *pad_d;
    size_t li, lpnum = gdata->particle_data->elem_count;
    particle_data_copy = sc_array_new(sizeof(pdata_t));
    for(li=0;li<lpnum;li++){
        pad_s = (pdata_t*) sc_array_index(gdata->particle_data,li);
        if(pad_s->ifboundary)
            continue;
        pad_d = (pdata_t*)sc_array_push(particle_data_copy);
        memcpy(pad_d,pad_s,sizeof(pdata_t));
    }
//    sc_array_copy(particle_data_copy,gdata->particle_data);
    swapXYZCoordinate();
    conn = p4est_connectivity_new_unitsquare();
    p4est_heating = p4est_new_ext(gdata->mpicomm, conn,1,0,1,sizeof(quadrant_data_t), NULL, this);  
    resetQuadrantData();
    
    }


void PelletSolver::prerun(){
    
   for (int i = 0; i < 2; ++i) {
    ilh[i] = sc_array_new (sizeof (p4est_locidx_t));
    jlh[i] = sc_array_new (sizeof (p4est_locidx_t));
    klh[i] = NULL;
   
   }
    
    }


void PelletSolver::resetQuadrantData(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  quadrant_data_t          *qud;
  for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est_heating->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (quadrant_data_t *) quad->p.user_data;
      qud->premain = 0;
      qud->preceive = 0;
    }
  }

}

void PelletSolver::presearch2d(){
        
  gdata->pfound = sc_array_new_count (sizeof (int), particle_data_copy->elem_count);
  
  gdata->iremain = sc_array_new (sizeof (p4est_locidx_t));
  sc_array_memset (gdata->pfound, -1);

  p4est_search_all (p4est_heating, 0, psearch_quad2d, psearch_point2d, particle_data_copy);



}




void PelletSolver::packParticles(){
  int                 retval;
  int                *pfn;
  size_t              zz, numz;
  void              **hfound;
  p4est_locidx_t      lremain, lsend, llost;
  comm_psend_t       *cps, *there;
  comm_prank_t       *trank;
  pdata_t          *pad;
  
  
  
  gdata->recevs = sc_array_new (sizeof (comm_prank_t));
  
  
  gdata->psmem = sc_mempool_new (sizeof (comm_psend_t));
  numz = gdata->pfound->elem_count;

  gdata->psend = sc_hash_new (psend_hash, psend_equal, NULL, NULL);
  lremain = lsend = llost = 0;
  cps = (comm_psend_t *) sc_mempool_alloc (gdata->psmem);
  cps->rank = -1;

  for (zz = 0; zz < numz; ++zz) {
    pfn = (int *) sc_array_index (gdata->pfound, zz);

    if (*pfn < 0) {
      assert(*pfn == -1);
      ++llost;
      continue;
    
    }
    if (*pfn == gdata->mpirank) {
      ++lremain;
      continue;
    }
    cps->rank = *pfn;
    retval = sc_hash_insert_unique (gdata->psend, cps, &hfound);
  
    there = *((comm_psend_t **) hfound);
  
    if (!retval) {
      assert (there->message.elem_size == sizeof(pdata_t));
      assert (there->message.elem_count > 0);
    }
  
    else {
      assert (there == cps);
      trank = (comm_prank_t *) sc_array_push (gdata->recevs);

      trank->rank = there->rank;
      trank->psend = there;
      sc_array_init (&there->message, sizeof(pdata_t));
      cps = (comm_psend_t *) sc_mempool_alloc (gdata->psmem);
      cps->rank = -1;
    }
  
    pad = (pdata_t *) sc_array_push (&there->message);
    memcpy (pad, sc_array_index (particle_data_copy, zz), sizeof (pdata_t));
  
    ++lsend;
  
  } 
  
  sc_mempool_free (gdata->psmem, cps);
  sc_array_sort (gdata->recevs, comm_prank_compare);

  sc_array_destroy_null (&gdata->pfound);

}

void PelletSolver::communicateParticles(){

  int                 mpiret;
  int                 i;
  int                 num_receivers;
  int                 num_senders;
  int                 count, cucount;
  int                 msglen;
  sc_MPI_Request     *reqs;
  sc_array_t         *notif, *payl;
  sc_array_t         *arr;
  comm_psend_t       *cps;
  comm_prank_t       *trank;
   
  num_receivers = gdata->recevs->elem_count;
 
  
  notif = sc_array_new_count (sizeof (int), num_receivers);
  payl = sc_array_new_count (sizeof (int), num_receivers);   //payload


  for (i = 0; i < num_receivers; ++i) {

    trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);

    *(int *) sc_array_index_int (notif, i) = trank->rank;

    cps = trank->psend;
  
    assert(cps->rank == trank->rank);

    arr = &cps->message;
 
    *(int *) sc_array_index_int (payl, i) = (int) arr->elem_count;
 
  
  }


  sc_notify_ext (notif, NULL, payl, NULL, gdata->mpicomm);

  assert (payl->elem_count == notif->elem_count);

  num_senders = (int) notif->elem_count;

  cucount = 0;
  for (i = 0; i < num_senders; ++i) {
    cucount += *(int *) sc_array_index_int (payl, i);
  }
  prebuf = sc_array_new_count (sizeof(pdata_t), cucount);

  /* post non-blocking receive */
  gdata->recv_req = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);

  cucount = 0;

  for (i = 0; i < num_senders; ++i) {
    count = *(int *) sc_array_index_int (payl, i);
    msglen = count * (int) sizeof(pdata_t);
    mpiret = sc_MPI_Irecv
      (sc_array_index (prebuf, cucount), msglen, sc_MPI_BYTE,
       *(int *) sc_array_index_int (notif, i), 999, gdata->mpicomm,
       (sc_MPI_Request *) sc_array_index_int (gdata->recv_req, i));
    SC_CHECK_MPI (mpiret);
    cucount += count;
  }
  
  assert(cucount == (int) prebuf->elem_count);

  sc_array_destroy_null (&notif);
  sc_array_destroy_null (&payl);

    gdata->send_req = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
    for (i = 0; i < num_receivers; ++i) {
      trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);
      cps = trank->psend;
      arr = &cps->message;
      msglen = (int) (arr->elem_count * arr->elem_size);
      mpiret = sc_MPI_Isend
        (arr->array, msglen, sc_MPI_BYTE, cps->rank, 999,
         gdata->mpicomm, (sc_MPI_Request *) sc_array_index_int (gdata->send_req, i));
      SC_CHECK_MPI (mpiret);
    }

  reqs = (sc_MPI_Request *) sc_array_index (gdata->recv_req,0);
  
  mpiret = sc_MPI_Waitall (num_senders, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&gdata->recv_req);

  num_receivers = (int) gdata->recevs->elem_count;
  reqs = (sc_MPI_Request *) sc_array_index (gdata->send_req,0),
  mpiret = sc_MPI_Waitall (num_receivers, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&gdata->send_req);

  for (i = 0; i < num_receivers; ++i) {
    trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);
    cps = trank->psend;
    sc_array_reset (&cps->message);
  }
  
  sc_array_destroy_null (&gdata->recevs);
  sc_hash_destroy (gdata->psend);

  gdata->psend = NULL;
  sc_mempool_destroy (gdata->psmem);
  gdata->psmem = NULL;
  
}


static int
slocal_quad2d (p4est_t * p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
             void *point)
{
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data      *g = p->gdata;

  /* compute coordinate range of this quadrant */
  loopquad2d (p, which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}


static int
slocal_point2d (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
              void *point)
{
  int                 i;
  char               *cf;
  size_t              zp;
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data      *g = p->gdata;
  quadrant_data_t          *qud;
  double             *x;
  pdata_t          *pad = (pdata_t *) point;

  /* access location of particle to be searched */
  x = pad->xyz;

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P4EST_DIM; ++i) {
    if (!(g->lxyz[i] <= x[i] && x[i] <= g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  if (local_num >= 0) {
    /* quadrant is a local leaf */
    /* first local match counts (due to roundoff there may be multiple) */
    zp = sc_array_position (p->prebuf, point);
    cf = (char *) sc_array_index (g->cfound, zp);
    if (!*cf) {
      /* make sure this particle is not found twice */
      *cf = 1;

      /* count this particle in its target quadrant */
      *(p4est_locidx_t *) sc_array_push (g->ireceive) = (p4est_locidx_t) zp;
      qud = (quadrant_data_t *) quadrant->p.user_data;
      ++qud->preceive;
    }

    /* return value will have no effect */
    return 0;
  }

  /* the leaf for this particle has not yet been found */
  return 1;
}

void PelletSolver::postsearch2d(){
  gdata->ireceive = sc_array_new (sizeof (p4est_locidx_t));
  gdata->cfound = sc_array_new_count (sizeof (char), prebuf->elem_count);
    
  sc_array_memset (gdata->cfound, 0);

  p4est_search_local (p4est_heating, 0, slocal_quad2d, slocal_point2d, prebuf);
  
  sc_array_destroy_null (&gdata->cfound);

}

void PelletSolver::regroupParticles2d(){

  sc_array_t         *newpa;
  p4est_topidx_t      tt;
  p4est_locidx_t      newnum;
  p4est_locidx_t      ppos;
  p4est_locidx_t      lq, prev;
  p4est_locidx_t      qboth, li;
  p4est_locidx_t     *premain, *preceive;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  quadrant_data_t          *qud;
  pdata_t          *pad;

  newnum =
    (p4est_locidx_t) (gdata->iremain->elem_count + gdata->ireceive->elem_count);

  premain = (p4est_locidx_t *) sc_array_index (gdata->iremain,0);
  preceive = (p4est_locidx_t *) sc_array_index (gdata->ireceive,0);
  newpa = sc_array_new_count (sizeof (pdata_t), newnum);
  pad = (pdata_t *) sc_array_index (newpa,0);
  prev = 0;
  
  for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
    
      tree = p4est_tree_array_index (p4est_heating->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (quadrant_data_t *) quad->p.user_data;
      qboth = qud->premain + qud->preceive;
      
      if (qboth == 0) {
        qud->lpend = prev;
        qud->premain = qud->preceive = 0;
        continue;
      }
      prev += qboth;
      for (li = 0; li < qud->premain; ++li) {
        ppos = *premain++;
        memcpy (pad, sc_array_index (particle_data_copy, ppos), sizeof (pdata_t));
        pad ++;
      }
      for (li = 0; li < qud->preceive; ++li) {
        ppos = *preceive++;
        memcpy (pad, sc_array_index (prebuf, ppos), sizeof (pdata_t));
        pad ++;
      }
      qud->lpend = prev;
      qud->premain = qud->preceive = 0;
    }
    
  } 
   
  sc_array_destroy_null (&gdata->iremain);

  sc_array_destroy_null (&prebuf);
  sc_array_destroy_null (&gdata->ireceive);
  sc_array_destroy (particle_data_copy);

  particle_data_copy = newpa;
}

static int  adapt_coarsen2d (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrants[])
{
    int i;
    p4est_locidx_t remain = 0, receive = 0;
    quadrant_data_t *oud;
    PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
    Global_Data      *g = p->gdata;
    
    if (quadrants[1] == NULL ||
      quadrants[0]->level == g->minlevel) {
    
        oud = (quadrant_data_t *) quadrants[0]->p.user_data;
        g->ireindex2 += oud->premain;
        g->irvindex2 += oud->preceive;
    
        return 0;
  }
   
    remain = receive = 0;
  
    for (i = 0; i < P4EST_CHILDREN; ++i) {
        oud = (quadrant_data_t *) quadrants[i]->p.user_data;
        remain += oud->premain;
        receive += oud->preceive;
  }
  if ((double) (remain + receive) <  p->elem_particle_box) {
    /* we will coarsen and adjust ireindex, irvindex in adapt_replace */
    g->qremain = remain;
    g->qreceive = receive;
    return 1;
  }
  else {
    /* we will not coarsen and proceed with next quadrant */
    oud = (quadrant_data_t *) quadrants[0]->p.user_data;
    g->ireindex2 += oud->premain;
    g->irvindex2 += oud->preceive;
    return 0;
  }


}

static void balance_replace2d (p4est_t * p4est, p4est_topidx_t which_tree,
               int num_outgoing, p4est_quadrant_t * outgoing[],
               int num_incoming, p4est_quadrant_t * incoming[]){

  int quadrantid;
  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg, offset, *newoffset;
  p4est_quadrant_t  **pchild, *outoctant;
  quadrant_data_t          *oud, *iud;
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data      *g = p->gdata;
  if (num_outgoing == P4EST_CHILDREN) {
    // we are coarsening 
      oud = (quadrant_data_t *) incoming[0]->p.user_data;
      oud->premain = oud->preceive = 0;
      for(int i=0;i<P4EST_CHILDREN;i++){
      iud = (quadrant_data_t *) outgoing[i]->p.user_data;
      oud->premain += iud->premain;
      oud->preceive += iud->preceive;
      }
      
  }

  else {
    outoctant = outgoing[0];
    iud = (quadrant_data_t *)outoctant->p.user_data;
    quadrantid = iud->quadrantid;
    // we are refining 
    // access parent quadrant 
    loopquad2d (p,which_tree, outgoing[0],lxyz,hxyz,dxyz);
    // recover window onto remaining particles for the new family 
    
    irem = iud->premain ;//sc_array_index(g->irecumu,octantid);
    
    ibeg = *(p4est_locidx_t *) sc_array_index_int(g->irecumu,quadrantid) - irem;
    offset = ibeg;
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    
    p->klh[0] = &iview;
    wz = 0;

    
    p->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      p->split_by_coord (p->jlh[wy], p->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = p->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
        oud->quadrantid = g->octantid;
        g->octantid ++;
        newoffset = (p4est_locidx_t *)sc_array_push(g->irecumu);
        *newoffset = offset + oud->premain;
        offset = offset + oud->premain;

      }
    }
    // recover window onto received particles for the new family 
    
    irem = iud->preceive ;//sc_array_index(g->irecumu,octantid);
    ibeg = *(p4est_locidx_t *) sc_array_index_int(g->irvcumu,quadrantid) - irem;
    offset = ibeg;
    P4EST_ASSERT (irem >= 0);
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);
    P4EST_ASSERT (qod->preceive == irem);

    // sort received particles into the children 
    pchild = incoming;
   
    p->klh[0] = &iview;
    wz = 0;
    
    p->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      p->split_by_coord ( p->jlh[wy], p->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = p->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->preceive = (p4est_locidx_t) arr->elem_count);
        newoffset = (p4est_locidx_t *)sc_array_push(g->irvcumu);
        *newoffset = offset + oud->preceive;
        offset += oud->preceive;
      }
    }
    P4EST_ASSERT (pchild == incoming + P4EST_CHILDREN);

  }
}


void PelletSolver::adaptQuadtree(){
    
    int oldquad = (int)p4est_heating->global_num_quadrants;
    
    while(true){
        gdata->ireindex2 = gdata->irvindex2 = 0;
        p4est_coarsen_ext (p4est_heating, 0, 0, adapt_coarsen2d, NULL, adapt_replace2d);
    
        if((int)p4est_heating->global_num_quadrants == oldquad)
            break;
        else
            oldquad = p4est_heating->global_num_quadrants;
    }
    
    
    
    while(true){
    
        gdata->octantid = 0;
        gdata->irecumu = sc_array_new(sizeof(p4est_locidx_t));
        gdata->irvcumu = sc_array_new(sizeof(p4est_locidx_t));
        gdata->ireindex = gdata->ire2 = 0;
        gdata->irvindex = gdata->irv2 = 0;
        p4est_refine_ext (p4est_heating, 0, 50, adapt_refine2d, NULL, adapt_replace2d);
        if((int)p4est_heating->global_num_quadrants == oldquad)
            break;
        else
            oldquad = p4est_heating->global_num_quadrants;
    
    
        sc_array_destroy(gdata->irecumu);
        sc_array_destroy(gdata->irvcumu);
    }
    
    //p4est_balance_ext(p4est_heating,P4EST_CONNECT_FULL,NULL,balance_replace2d);
    sc_array_destroy(gdata->irecumu);
    sc_array_destroy(gdata->irvcumu);
    }

static int
part_weight2d (p4est_t * p4est,
             p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  p4est_locidx_t      ilem_particles;
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data      *g = p->gdata;
  quadrant_data_t          *qud = (quadrant_data_t *) quadrant->p.user_data;


  ilem_particles = qud->lpend - g->prevlp;

  g->prevlp = qud->lpend;
  *(int *) sc_array_index (g->src_fixed, g->qcount++) =
    (int) (ilem_particles * sizeof (pdata_t));
  return 1+ ilem_particles;
}
void PelletSolver:: partitionParticles2d(){


  sc_array_t         *dest_data;
  p4est_topidx_t      tt;
  p4est_locidx_t      ldatasize, lcount;
  p4est_locidx_t      dest_quads, src_quads;
  p4est_locidx_t      dest_parts;
  p4est_locidx_t      lquad, lq;
  p4est_locidx_t      lpnum;
  p4est_gloidx_t      gshipped;
  p4est_gloidx_t     *src_gfq;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  quadrant_data_t          *qud;
   
  if(gdata->mpisize == 1)
  
      return;

  src_gfq = P4EST_ALLOC (p4est_gloidx_t, gdata->mpisize + 1);

  memcpy (src_gfq, p4est_heating->global_first_quadrant,
          (gdata->mpisize + 1) * sizeof (p4est_gloidx_t));


  src_quads = p4est_heating->local_num_quadrants;

  assert(src_quads == src_gfq[gdata->mpirank+1]-src_gfq[gdata->mpirank]);

  gdata->src_fixed = sc_array_new_count (sizeof (int), src_quads);

  gdata->qcount = 0;
  gdata->prevlp = 0;

  gshipped = p4est_partition_ext (p4est_heating, 1, part_weight2d);
  dest_quads = p4est_heating->local_num_quadrants;

  if (gshipped == 0) {
    sc_array_destroy_null (&gdata->src_fixed);
    P4EST_FREE (src_gfq);
    return;
  }

  gdata->dest_fixed = sc_array_new_count (sizeof (int), dest_quads);

  p4est_transfer_fixed (p4est_heating->global_first_quadrant, src_gfq,
                        gdata->mpicomm, COMM_TAG_FIXED,
                        (int *) gdata->dest_fixed->array,
                        (const int *) gdata->src_fixed->array, sizeof (int));

  ldatasize = (p4est_locidx_t) sizeof (pdata_t);

  dest_parts = 0;

  for (lq = 0; lq < dest_quads; ++lq) {
    dest_parts += *(int *) sc_array_index (gdata->dest_fixed, lq);
  }
  assert(dest_parts % ldatasize == 0); 
  dest_parts /= ldatasize;
  dest_data = sc_array_new_count (sizeof (pdata_t), dest_parts);
  p4est_transfer_custom (p4est_heating->global_first_quadrant, src_gfq,
                         gdata->mpicomm, COMM_TAG_CUSTOM,
                         (pdata_t *) dest_data->array,
                         (const int *) gdata->dest_fixed->array,
                         (const pdata_t *) particle_data_copy->array,
                         (const int *) gdata->src_fixed->array);

  sc_array_destroy_null (&gdata->src_fixed);

  P4EST_FREE (src_gfq);
  sc_array_destroy (particle_data_copy);
  particle_data_copy = dest_data;
  lpnum = 0;
  lquad = 0;
  for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est_heating->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      /* access quadrant */
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (quadrant_data_t *) quad->p.user_data;

      /* back out particle count in quadrant from data size */
      lcount = *(int *) sc_array_index (gdata->dest_fixed, lquad);
      assert (lcount % ldatasize == 0);
      lcount /= ldatasize;
      lpnum += lcount;
      qud->lpend = lpnum;
      ++lquad;
    }
  }
  sc_array_destroy_null (&gdata->dest_fixed);


}
void PelletSolver::destoryQuadtree(){
    
        p4est_destroy (p4est_heating);
      
        p4est_heating = NULL;
      
        p4est_connectivity_destroy (conn);
      
        conn = NULL;
        
        sc_array_destroy_null(&prebuf_integral);  
        sc_array_destroy_null(&particle_data_copy); 
        for(int i=0;i<2;i++){ 
            sc_array_destroy_null (&ilh[i]);
            sc_array_destroy_null (&jlh[i]);
        }
    
    }


static int
compareXCoordinate (const void *p1, const void *p2)
{
  int t = 0;
  pdata_t                * i1 = (pdata_t *) p1;
  pdata_t                * i2 = (pdata_t *) p2;
  
  if((i1->xyz[2] - i2->xyz[2])>0)
      t = 1;
  return t;
}

void PelletSolver::computeLongitudinalLength(){
    size_t li, lpnum = gdata->particle_data->elem_count;
    double max_l = -1000;
    double min_l = 1000;
    double max_lg;
    double min_lg;
    double x;
    pdata_t *pad;
    for(li = 0; li<lpnum; li++){
         pad = (pdata_t *)sc_array_index(gdata->particle_data, li);
         x = pad->xyz[0];
         if(x > max_l)
             max_l = x;
         if(x < min_l)
             min_l = x;
    }
    
    MPI_Barrier(gdata->mpicomm);
    MPI_Allreduce(&max_l, &max_lg, 1, MPI_DOUBLE, MPI_MAX, gdata->mpicomm);
    MPI_Allreduce(&min_l, &min_lg, 1, MPI_DOUBLE, MPI_MIN, gdata->mpicomm);
    longitude_min = min_lg;
    longitude_max = max_lg;
}
void PelletSolver::computeDensityIntegral_interpolation(){

    p4est_topidx_t      tt;
    p4est_locidx_t      lq, lpend, offset = 0;
    p4est_tree_t       *tree;
    p4est_quadrant_t   *quad;
    quadrant_data_t          *qud;
    pdata_t          *pad;
    computation_node_t *node;
    size_t li;
    int  *lnend = (int*)sc_array_index(node_offset,0);
    int length_offset = 0;
    double R = 160.;
    for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
        tree = p4est_tree_array_index (p4est_heating->trees, tt);
        for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
            quad = p4est_quadrant_array_index (&tree->quadrants, lq);
            qud = (quadrant_data_t *) quad->p.user_data;
            lpend = qud->lpend;
            //binary tree construction;            
            if(lpend == offset)
                continue;

            node = (computation_node_t*)sc_array_index_int(nodelist,(*lnend)-1);
            for(li = offset; li< (size_t)lpend; li++){
                pad = (pdata_t *) sc_array_index(particle_data_copy,li);
                if(pad->ifboundary)
                    continue;

                pad->drifta = 2*node->pressure/node->density/R; 
                linearInterpolation(length_offset, *lnend-1, pad);
                
            }
            
            offset = lpend;
            length_offset = *lnend;
            lnend ++;
        }
    }
    sc_array_destroy(nodelist);
    sc_array_destroy(node_offset);
//    assert((size_t)lpend == particle_data_copy->elem_count); 
}

void PelletSolver::linearInterpolation(int startindex, int endindex, pdata_t *pad){
    
    double x0 = pad->xyz[2];
    int start = startindex;
    int end = endindex;
    int mid;
    computation_node_t *node1, *node2,  *nodemid;
    double xmid, x1, x2;
    
    while((end-start)!=1){
        mid = (start+end)/2;
        nodemid = (computation_node_t *)sc_array_index_int(nodelist,mid);
        xmid = nodemid->xyz[0];
        if(xmid >= x0)
            end = mid;
        else
            start = mid;
    }
   
    node1 = (computation_node_t *)sc_array_index_int(nodelist,start);
    node2 = (computation_node_t *)sc_array_index_int(nodelist,end);
    x1 = node1->xyz[0];
    x2 = node2->xyz[0];
    if(!(x0>=x1 && x0<=x2))
        printf("%f %f %f\n",x0,x1,x2);
    assert(x0>=x1 && x0<=x2);
    if(node1->sphdensity<1.e5 && node2->sphdensity<1.e5){ 
        pad->rightintegral = ((x2-x0)*node1->rightintegral + (x0-x1)*node2->rightintegral)/(x2-x1);
        pad->leftintegral = ((x2-x0)*node1->leftintegral + (x0-x1)*node2->leftintegral)/(x2-x1);
    }
    else if(node1->sphdensity<1.e5 && node2->sphdensity>1.e5){

        pad->leftintegral = ((x2-x0)*node1->leftintegral + (x0-x1)*node2->leftintegral)/(x2-x1);
        if(start == startindex)
            pad->rightintegral = node1->rightintegral;
        else{
            node1 = (computation_node_t *)sc_array_index_int(nodelist,start-1);
            node2 = (computation_node_t *)sc_array_index_int(nodelist,end-1);
            x1 = node1->xyz[0];
            x2 = node2->xyz[0];
            assert(x0>=x1 && x0>=x2);
            pad->rightintegral = node2->rightintegral+(node2->rightintegral-node1->rightintegral)/(x2-x1)*(x0-x2);
        }
    }
    else if(node2->sphdensity<1.e5 && node1->sphdensity>1.e5){
       
        pad->rightintegral = ((x2-x0)*node1->rightintegral + (x0-x1)*node2->rightintegral)/(x2-x1);
        if(end == endindex)
            pad->leftintegral = node2->leftintegral;
        else{
            node1 = (computation_node_t *)sc_array_index_int(nodelist,start+1);
            node2 = (computation_node_t *)sc_array_index_int(nodelist,end+1);
            
            x1 = node1->xyz[0];
            x2 = node2->xyz[0];
            assert(x0<=x1 && x0<=x2);
            pad->leftintegral = node2->leftintegral+(node2->leftintegral-node1->leftintegral)/(x2-x1)*(x0-x2);
        }
    }
    else{
    
        pad->rightintegral = ((x2-x0)*node1->rightintegral + (x0-x1)*node2->rightintegral)/(x2-x1);
        pad->leftintegral = ((x2-x0)*node1->leftintegral + (x0-x1)*node2->leftintegral)/(x2-x1);
    }

//    cout<<pad->leftintegral<<" "<<node1->leftintegral<<" "<<node2->leftintegral<<endl;    
}
void PelletSolver::computeDensityIntegral(){

    p4est_topidx_t      tt;
    p4est_locidx_t      lq, lpend, offset = 0;
    p4est_tree_t       *tree;
    p4est_quadrant_t   *quad;
    quadrant_data_t          *qud;
    pdata_t          *pad;
    double dy, dz;
    double lxyz[3],hxyz[3],dxyz[3];
    sc_array_t view;  
    double integral; 
    vector<int> cell_ID;
    vector<double> x_left;
    vector<double> x_right;
    vector<double> x_dividing;
    for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
        tree = p4est_tree_array_index (p4est_heating->trees, tt);
        for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
            quad = p4est_quadrant_array_index (&tree->quadrants, lq);
            qud = (quadrant_data_t *) quad->p.user_data;
            lpend = qud->lpend;
            //binary tree construction;            
            if(lpend == offset)
                continue;
            integral = 0;
            sc_array_init_view(&view,particle_data_copy,offset,lpend-offset);
            sc_array_sort(&view,compareXCoordinate); 
            loopquad2d(this,tt,quad,lxyz,hxyz,dxyz);
            dz = dxyz[0];
            dy = dxyz[1];
            pad = (pdata_t*) sc_array_index_int(&view,0);
            pad = (pdata_t*) sc_array_index_int(&view,lpend-offset-1);
            
            for(int id = 0; id<lpend-offset; id++){
                pad = (pdata_t *) sc_array_index_int(&view,id);
                integral += pad->mass;
                pad->leftintegral = integral/dz/dy;
            }

            integral = 0;

            for(int id = lpend-offset-1; id >= 0; id--){
                pad = (pdata_t *) sc_array_index_int(&view,id);
                integral += pad->mass;
                pad->rightintegral = integral/dz/dy;
            }
            /*
            cell_ID.push_back(0);

            x_left.push_back(x_min);
            x_right.push_back(x_max);

            i_node = 0;
            i_nodeend = 1;
            np_node = lpend - offset;
            while(i_node < i_nodeend){
                node_id = cell_ID[i_node];
                x_lower = x_left[node_id];
                x_upper = x_right[node_id];
                np_node = countNumberinRange(&view,lpend-offset,x_lower,x_upper);
                if(np_node <= elem_particle_cell){ 
                    i_node ++;}
                else{
                    x_center = (x_lower+x_upper)/2;
                    cell_ID.push_back(i_nodeend);
                    cell_ID.push_back(i_nodeend+1);
                    x_left.push_back(x_lower);
                    x_left.push_back(x_center);
                    x_right.push_back(x_center);
                    x_right.push_back(x_upper);
                    x_dividing.push_back(x_center);
                    i_node++;
                    i_nodeend += 2;
                    
                    }
                
                
                
                }
            
            x_dividing.push_back(x_max);
            sort(x_dividing.begin(),x_dividing.end());
            division = x_dividing.size();
            vector<double> left_sum(division,0);
            vector<double> right_sum(division,0);
            vector<double> local_sum(division,0) ;
            
            for(int id_particle=0,id_cell=0;id_particle<lpend-offset;)
                {   
                    pad = (pdata_t *) sc_array_index_int(&view,id_particle);
                    if(pad->xyz[2] <= x_dividing[id_cell]){
                            
                        local_sum[id_cell] += pad->mass;
                        if(pad->mass == 0)
                            assert(false);
                        ++id_particle;
                        
                    }
                    else
                    
                        {    
                             ++id_cell;
                          }            
                }
               
             for(int k = 0; k<division; k++){
                for(int p = 0; p<k+1; p++){
                    left_sum[k] += local_sum[p];
                    
                    }
                    left_sum[k] -= local_sum[k]/2;
                    left_sum[k] = left_sum[k]/dy/dz; 
              
                for(int p = k; p<division; p++){
                    right_sum[k] += local_sum[p];
                    }
                    right_sum[k] -= local_sum[k]/2;
                    right_sum[k] = right_sum[k]/dy/dz; 
                
                }
                
                
            for(int id_particle = 0, id_cell = 0;id_particle<lpend - offset; ){
                
                  pad = (pdata_t *) sc_array_index_int(&view,id_particle);
                  if(pad->xyz[2] <= x_dividing[id_cell])   
                        {
                         pad->leftintegral = left_sum[id_cell];
                         pad->rightintegral = right_sum[id_cell];
                         id_particle++;
                         }
                   else
                         {  
                         id_cell++;
                        }
            }
                
            cell_ID.clear();
            x_left.clear();
            x_right.clear();
            x_dividing.clear();
            */
            offset = lpend;
                
             }   
        
        
        }
    
    
    
    }

void PelletSolver:: generateComputationNode(){
    
            computeLongitudinalLength();
            initComputationNode();
            double *nodesls = new double[nodelist->elem_count];//localspacing
            searchComputationNode(nodelist);
            
            getNodeLocalSpacing(nodelist,nodesls);
            
            int counter_l, counter_r; 
            counter_r = marchComputationNode(nodelist, RIGHT);
            nodelist_right = nodelist_buf;
            counter_l = marchComputationNode(nodelist, LEFT);
            nodelist_left = nodelist_buf;
            
           /* 
    for(size_t li = 0;li<nodelist_left->elem_count;li++){
    
            computation_node_t *node1 = (computation_node_t *) sc_array_index(nodelist_left,li);
            if(node1->xyz[0] != 0)
                cout<<node1->xyz[0]<<endl;
    }  */
                 
            while(true){
                if(counter_r > 0){ 
                    searchComputationNode(nodelist_right);
                    getNodeLocalSpacing(nodelist_right,nodesls);
                }
                if(counter_l > 0){ 
                    searchComputationNode(nodelist_left);
                    getNodeLocalSpacing(nodelist_left,nodesls);
                }
                
                size_t old_count = nodelist->elem_count;
                if((counter_r + counter_l)>0){
                    sc_array_resize(nodelist,old_count+(size_t)counter_r+(size_t)counter_l);
                    if(counter_r>0)
                        sc_array_copy_into(nodelist,old_count,nodelist_right);
                    if(counter_l>0)
                        sc_array_copy_into(nodelist,old_count+counter_r,nodelist_left);
                }
                
                if(counter_r > 0){
                    counter_r = marchComputationNode(nodelist_right, RIGHT);
                    sc_array_destroy(nodelist_right);
                    nodelist_right = nodelist_buf;
                }

                if(counter_l > 0){
                    counter_l = marchComputationNode(nodelist_left, LEFT);
                    sc_array_destroy(nodelist_left);
                    nodelist_left = nodelist_buf;
                }
                
                if((counter_r + counter_l) == 0){
                
                    sc_array_destroy(nodelist_right);
                    sc_array_destroy(nodelist_left);
                    break;
                }
            }
       
         MPI_Barrier(gdata->mpicomm);
        size_t lnnum = nodelist->elem_count;
        
         P4EST_GLOBAL_ESSENTIALF ("Generated %d computational nodes.\n",lnnum);
        delete[] nodesls;    
}
void PelletSolver::reorderComputationNode(){
        
     size_t li, lnnum = nodelist->elem_count;
     computation_node_t *node, *node2;
     int counter, counter_g;
     int id;
     int *lnend;
     nodelist_buf = sc_array_new(sizeof(computation_node_t));
     for(li=0; li<lnnum; li++){
        node = (computation_node_t *)sc_array_index(nodelist,li);
        if(node->mpirank == gdata->mpirank){
            node2 = (computation_node_t *)sc_array_push(nodelist_buf);
            memcpy(node2,node,sizeof(computation_node_t));
        }
     }
    sc_array_destroy(nodelist);
    nodelist = nodelist_buf;
    counter = nodelist->elem_count;
    MPI_Allreduce(&counter, &counter_g, 1, MPI_INT,MPI_SUM,gdata->mpicomm);
    assert(counter_g == (int)lnnum);
    sc_array_sort(nodelist,comm_node_compare); 
    node_offset = sc_array_new(sizeof(int));
    if(counter == 0)
        return;
    node = (computation_node_t *)sc_array_index(nodelist,0);
    id = node->id;
    counter = 1;
    
    for(li=1;li<nodelist->elem_count;li++){
        node = (computation_node_t *)sc_array_index(nodelist,li);
        if(id != node->id){
            lnend = (int *)sc_array_push(node_offset);
            *lnend = counter;  
            
        }
        id = node->id;
        counter ++;
    }

    lnend = (int *)sc_array_push(node_offset);
    *lnend = counter;

    }

void PelletSolver::densityIntegralForNode(){
    
   size_t offset = 0,lnend, li;
   size_t id,idn = node_offset->elem_count;
   double integral_r, integral_l;
   double integral_pressure, integral_rho, integral_pressureold, integral_rhoold;
   computation_node_t *node1, *node2;
   for(id = 0; id<idn; id++){
       integral_r = integral_l = 0.;
       integral_pressure = integral_rho = 0;
       lnend = *(int *)sc_array_index(node_offset,id); 
       
       for(li = offset; li<lnend-1; li++){
            node1 = (computation_node_t *)sc_array_index(nodelist,li);
            node2 = (computation_node_t *)sc_array_index(nodelist,li+1);
            node1->rightintegral = integral_r;
            integral_pressureold = integral_pressure;
            integral_rhoold = integral_rho;
            if(node1->rightintegral == 0)
                node1->rightintegral = 1.e-8;
            integral_r += (node2->xyz[0] - node1->xyz[0]) * (node1->sphdensity+node2->sphdensity)/2.;
            integral_pressure += (node2->xyz[0] - node1->xyz[0]) * max(((node1->pressure+node2->pressure)/2.-0.64),0.);
            integral_rho += (node2->xyz[0] - node1->xyz[0]) * ((node1->density+node2->density)/2.);
            node1->pressure = integral_pressureold;
            node1->density = integral_rhoold;
       }
   
            node1 = (computation_node_t *)sc_array_index(nodelist,lnend-1);
            node1->rightintegral = integral_r;
            node1->pressure = max(integral_pressure,0.);
            node1->density = integral_rho;
            if(node1->rightintegral == 0)
                node1->rightintegral = 1.e-8;
       for(li = lnend-1; li > offset; li--){
            node1 = (computation_node_t *)sc_array_index(nodelist,li);
            node2 = (computation_node_t *)sc_array_index(nodelist,li-1);
            node1->leftintegral = integral_l;
            if(node1->leftintegral == 0)
                node1->leftintegral = 1.e-8;
            integral_l += (node1->xyz[0] - node2->xyz[0]) * (node1->sphdensity+node2->sphdensity)/2.;
       }
   
            node1 = (computation_node_t *)sc_array_index(nodelist,offset);
            node1->leftintegral = integral_l;
            if(node1->leftintegral == 0)
                node1->leftintegral = 1.e-8;
   
      offset = lnend;
   }
}

void PelletSolver:: initComputationNode(){
    
    p4est_topidx_t      tt;
    p4est_locidx_t      lq, lpend, offset = 0;
    p4est_tree_t       *tree;
    p4est_quadrant_t   *quad;
    quadrant_data_t          *qud;
    computation_node_t          *node;
    sc_array_t *nodelist_temp;
    pdata_t *pad;
    size_t li;
    double y_avg, z_avg;
    nodelist_temp = sc_array_new(sizeof(computation_node_t)); 
    for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
        tree = p4est_tree_array_index (p4est_heating->trees, tt);
        for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
            quad = p4est_quadrant_array_index (&tree->quadrants, lq);
            qud = (quadrant_data_t *) quad->p.user_data;
            lpend = qud->lpend;
            //binary tree construction;            
            if(lpend == offset)
                continue;
            
            y_avg = 0;
            z_avg = 0;
            for(li = offset; li<(size_t)lpend ;li++){
                pad = (pdata_t *) sc_array_index(particle_data_copy,li);
                //  if(pad->ifboundary)
                   // continue;
                y_avg += pad->xyz[0];
                z_avg += pad->xyz[1];
            }
            y_avg /= (lpend-offset);
            z_avg /= (lpend-offset);
        //    loopquad2d(this,tt,quad,lxyz,hxyz,dxyz);
            node = (computation_node_t *) sc_array_push(nodelist_temp);
            node->xyz[0] = 0.;
            node->xyz[1] = y_avg;//(lxyz[0]+hxyz[0])/2.;
            node->xyz[2] = z_avg;//(lxyz[1]+hxyz[1])/2.;
            node->mpirank = gdata->mpirank;
            node->id = (int)lq;
            node->localspacing = 0 - gdata->initlocalspacing; 
            node->sphdensity = -1;
            offset = lpend;
        }
    }
    
    size_t lnnum = (int)nodelist_temp->elem_count * (int)sizeof(computation_node_t);
    int lnnum_count[gdata->mpisize];
    int lnnum_offset[gdata->mpisize];
   
    MPI_Allgather(&lnnum,1,MPI_INT,lnnum_count,1,MPI_INT,gdata->mpicomm);
    
    offset = 0;
    for(int i=0; i<(int)gdata->mpisize; i++){
        lnnum_offset[i] = offset;
        
        offset += lnnum_count[i];
        
    }
    offset /= (int)sizeof(computation_node_t);
    nodelist = sc_array_new_count(sizeof(computation_node_t), (size_t) offset);
    MPI_Allgatherv(nodelist_temp->array,lnnum_count[gdata->mpirank],MPI_BYTE,nodelist->array, lnnum_count,
            lnnum_offset,MPI_BYTE,gdata->mpicomm); 
    sc_array_destroy_null(&nodelist_temp);
    
}

static int
psearch_quad (p8est_t * p4est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant, 
              p4est_locidx_t local_num, void *point)
{
  Global_Data      *g = (Global_Data *) p4est->user_pointer;

  /* compute coordinate range of this quadrant */
  g->loopquad ( which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}
static int
compute_node_state (p8est_t * p4est, p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrant,
               p4est_locidx_t local_num, void *point)
{
  
  p8est_tree_t       *tree;
  int                 i;
  double             *loc;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  PelletSolver     *p = g->pellet_solver;
  octant_data_t          *qud, *qudnei;
  computation_node_t          *node = (computation_node_t *) point;
  pdata_copy_t *padnei;
  p8est_quadrant_t *quadnei;
  p4est_locidx_t   *localneiid, *ghostneiid;
  double x, y ,z, x0, y0 ,z0, dx, dy, dz;
  double radius, h ,c, u, dis;
  size_t localsize, ghostsize;
  double count_den = 0.; 
  double b, Ktotal, Ki, den_avg, pre_avg;
  tree = p8est_tree_array_index (p4est->trees, which_tree);
  /* access location of particle to be searched */
  loc = node->xyz;
  /* due to roundoff we call this even for a local leaf */
  
  if(node->sphdensity >= 0)
      return 0;
  for (i = 0; i < P8EST_DIM; ++i) {
    if (!(g->lxyz[i] <= loc[i] && loc[i] < g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }
  if (local_num >= 0) {
        
        assert(node->localspacing > 0);  
        x0 = node->xyz[0];
        y0 = node->xyz[1];
        z0 = node->xyz[2];
        radius = node->localspacing * g->timesearchingradius; 
        qud = (octant_data_t *)quadrant->p.user_data;   
        localsize = qud->localneighbourid->elem_count;
        ghostsize = qud->ghostneighbourid->elem_count;
         h = radius/2.;
         c = 1./M_PI/h/h/h;
         Ktotal = den_avg = pre_avg = 0;
         b = radius/3;
        for(size_t ii=0; ii<localsize; ii++){ //iterate through local neighbour octants
            localneiid = (p4est_locidx_t *)sc_array_index(qud->localneighbourid,ii);
            quadnei = p8est_quadrant_array_index(&tree->quadrants,*localneiid);
            qudnei = (octant_data_t *) quadnei->p.user_data;
            size_t nump = qudnei->poctant;
            for(size_t pid = 0;pid<nump; pid++){
                padnei = &qudnei->localparticle[pid];
                
                if(padnei->ifboundary){
                   // continue;
                }
                x = padnei->xyz[0]; 
                y = padnei->xyz[1]; 
                z = padnei->xyz[2]; 
                dx = x-x0;
                dy = y-y0;
                dz = z-z0;
                dis = sqrt(dx*dx + dy*dy + dz*dz);
                u = dis/h;
//SPH density
                if(u > 2)  continue;
                else if(u <= 1){
                    count_den += c*(1-1.5*u*u+0.75*u*u*u)*padnei->mass; 
                }
                else 
                {
                     
                    count_den += 0.25*c*(2-u)*(2-u)*(2-u)*padnei->mass;
                }
//Weighted avg of density
                Ki = exp(-dis*dis/2/b/b);
                Ktotal += Ki;
                den_avg += Ki * 1./padnei->volume;
                pre_avg += Ki * padnei->pressure;
            }  
        }
    
        for(size_t ii=0; ii<ghostsize; ii++){ //iterate through ghost neighbour octants
            ghostneiid = (p4est_locidx_t *)sc_array_index(qud->ghostneighbourid,ii);
          //  quadnei = p8est_quadrant_array_index(&tree->quadrants,*localneiid);
            qudnei = &g->ghost_data[*ghostneiid];
            size_t nump = qudnei->poctant;
        
            for(size_t pid = 0;pid<nump; pid++){
                padnei = &qudnei->localparticle[pid];
                if(padnei->ifboundary){
                   // continue;
                }
                x = padnei->xyz[0]; 
                y = padnei->xyz[1]; 
                z = padnei->xyz[2]; 
                dx = x-x0;
                dy = y-y0;
                dz = z-z0;
                dis = sqrt(dx*dx + dy*dy + dz*dz);
                u = dis/h;
                if(u > 2)  continue;
//SPH density
                else if(u <= 1){
                    count_den += c*(1-1.5*u*u+0.75*u*u*u)*padnei->mass; 
                }
                else 
                {
                     count_den += 0.25*c*(2-u)*(2-u)*(2-u)*padnei->mass;
                }
            
//Weighted avg of density
                Ki = exp(-dis*dis/2/b/b);
                Ktotal += Ki;
                den_avg += Ki * 1./padnei->volume;
                pre_avg += Ki * padnei->pressure; 
            }  
        
        }

     node->sphdensity = den_avg/Ktotal;//count_den;
     node->pressure = pre_avg/Ktotal;
     if(Ktotal == 0){
      
         node->sphdensity = 0.;
         node->pressure = 0.;
     }
     node->density = node->sphdensity;
     if(p->ifInsidePellet(x0,y0,z0)){
        node->sphdensity  = 1.e8;
    }
     return 0;
  }
  /* the process for this particle has not yet been found */
  return 1;
}

static int
psearch_node (p8est_t * p4est, p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrant,
               p4est_locidx_t local_num, void *point)
{
  
  int                 i;
  double             *loc;
  size_t li, lonum;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  octant_data_t          *qud;
  computation_node_t          *node = (computation_node_t *) point;
  pdata_t *pad;
  double x, y ,z, x0, y0 ,z0;
  double dis, dis_min = -1;
  /* access location of particle to be searched */
  loc = node->xyz;

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P8EST_DIM; ++i) {
    if (!(g->lxyz[i] <= loc[i] && loc[i] < g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  if (local_num >= 0) {
     if(node->localspacing < 0){
        x0 = node->xyz[0];
        y0 = node->xyz[1];
        z0 = node->xyz[2];
        qud = (octant_data_t *)quadrant->p.user_data;   
        lonum = qud->poctant; 
        li = qud->lpend - lonum; 
        for(; li<(size_t )qud->lpend; li++ ){
            pad = (pdata_t *)sc_array_index(g->particle_data,li);
            x = pad->xyz[0];
            y = pad->xyz[1];
            z = pad->xyz[2];
            dis = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
            if(dis_min<0 || dis<dis_min){
                dis_min = dis;
                node->localspacing = pad->localspacing;
            }
        }
        if(dis_min<0)
            {
                /*
                if(p->ifCloseToPellet(1,x0,y0,z0)) 
                    node->localspacing = g->initlocalspacing*10;
                else
                    node->localspacing =  g->initlocalspacing*200;*/

                node->localspacing = 0 - node->localspacing;
                    
//               node->localspacing = 10* g->initlocalspacing;
            }
    
        
     return 0;
    }
  }
  /* the process for this particle has not yet been found */
  return 1;
}

void PelletSolver::searchComputationNode(sc_array_t *nodes){
    
    p8est_search_local(gdata->p8est,0,psearch_quad,psearch_node ,nodes);

}

void PelletSolver::regroupComputationNode(sc_array_t *nodes){
    
    size_t li, lpnum = nodes->elem_count;
    int counter = 0;
    int counter_g;
    int offset;
    int lnnum_count[gdata->mpisize];
    int lnnum_offset[gdata->mpisize];
    sc_array_t *nodelist_temp;
    computation_node_t *node, *node2;
    nodelist_temp = sc_array_new(sizeof(computation_node_t));
    for(li=0;li<lpnum;li++){
        node = (computation_node_t *)sc_array_index(nodes,li);
        if(node->localspacing > 0){
            node2 = (computation_node_t *)sc_array_push(nodelist_temp);
            memcpy(node2,node,sizeof(computation_node_t));
            counter ++;
        }
    }
    MPI_Barrier(gdata->mpicomm);
    MPI_Allreduce(&counter,&counter_g,1,MPI_INT,MPI_SUM,gdata->mpicomm);
    assert((size_t)counter_g == nodes->elem_count);
    counter *= (int)sizeof(computation_node_t);  
    MPI_Allgather(&counter,1,MPI_INT,lnnum_count,1,MPI_INT,gdata->mpicomm);

    offset = 0;
    for(int i=0; i<(int)gdata->mpisize; i++){
        lnnum_offset[i] = offset;
        
        offset += lnnum_count[i];
        
    }
    offset /= (int)sizeof(computation_node_t);
    assert(offset == counter_g);
    nodelist_buf = sc_array_new_count(sizeof(computation_node_t), (size_t) offset);
    MPI_Allgatherv(nodelist_temp->array,lnnum_count[gdata->mpirank],MPI_BYTE,nodelist_buf->array, lnnum_count,
            lnnum_offset,MPI_BYTE,gdata->mpicomm); 
    sc_array_destroy_null(&nodelist_temp);

}

int PelletSolver::marchComputationNode(sc_array_t* nodes, march_direction_t dir ){
        size_t li, lnnum = nodes->elem_count;
        computation_node_t *node1, *node2;
        int counter = 0;
        double march_length;
        if(dir == RIGHT){
           nodelist_buf = sc_array_new(sizeof(computation_node_t));
           for(li=0;li<lnnum;li++){
               node1 = (computation_node_t *)sc_array_index(nodes,li);
            //   march_length = node1->localspacing * 2;
               march_length = min(node1->localspacing,0.2);
               if(node1->xyz[0]  < longitude_max){
                    node2 = (computation_node_t *)sc_array_push(nodelist_buf);
                    node2->xyz[0] = min(node1->xyz[0] + march_length, longitude_max);
                    node2->xyz[1] = node1->xyz[1];
                    node2->xyz[2] = node1->xyz[2];
                    node2->mpirank = node1->mpirank;
                    node2->id = node1->id;
                    node2->localspacing = 0 - node1->localspacing;
                    node2->sphdensity = -1;
                    counter ++;
               }
           }
       } 

        if(dir == LEFT){
           nodelist_buf = sc_array_new(sizeof(computation_node_t));
           for(li=0;li<lnnum;li++){
               node1 = (computation_node_t *)sc_array_index(nodes,li);
           //    march_length = node1->localspacing * 2;
               march_length = min(node1->localspacing,0.2);
               if(node1->xyz[0] > longitude_min){
                    node2 = (computation_node_t *)sc_array_push(nodelist_buf);
                    
                    node2->xyz[0] = max(node1->xyz[0] - march_length, longitude_min);
                    node2->xyz[1] = node1->xyz[1];
                    node2->xyz[2] = node1->xyz[2];
                    node2->mpirank = node1->mpirank;
                    node2->id = node1->id;
                    node2->localspacing = 0 - node1->localspacing;
                    node2->sphdensity = -1;
                    counter ++;
               }
           }
       } 
        
        return counter;
}


void PelletSolver::swapXYZCoordinate(){
    pdata_t *pad;    
    size_t li, lpnum = particle_data_copy->elem_count;
    
    for(li=0; li<lpnum; li++){
        pad = (pdata_t *) sc_array_index(particle_data_copy,li);
        swap(pad->xyz[0],pad->xyz[1]);
        swap(pad->xyz[1],pad->xyz[2]);
        }
    
    } 



int PelletSolver::countNumberinRange(sc_array_t *view, int n, double x, double y) 
    {   
        pdata_t* pad;
        int l = 0, h = n - 1; 
        while (l <= h) { 
            int mid = (l + h) / 2;
            pad = (pdata_t* )sc_array_index_int(view,mid);
            if (pad->xyz[2] >= x) 
                h = mid - 1; 
            else
               l = mid + 1; 
                       } 
        int lower = l;
       l = 0; 
       h = n-1;
         while (l <= h)
         { 
            int mid = (l + h) / 2;
            pad = (pdata_t* )sc_array_index_int(view,mid);
            if (pad->xyz[2] <= y) 
                l = mid + 1; 
            else  
                h = mid - 1;       
                }
        int upper = h;
        return upper-lower+1;
    } 

void PelletSolver::packParticles_phase2(){
    
  int                 retval;
  int                pfn;
  size_t              zz, numz;
  void              **hfound;
  p4est_locidx_t      lremain, lsend;
  comm_psend_t       *cps, *there;
  comm_prank_t       *trank;
  pdata_t          *pad;
  integral_t * integral; 
  
  
  gdata->recevs = sc_array_new (sizeof (comm_prank_t));
  
  gdata->psmem = sc_mempool_new (sizeof (comm_psend_t));
  gdata->psend = sc_hash_new (psend_hash, psend_equal, NULL, NULL);
  numz = particle_data_copy->elem_count;  
    
  lremain = lsend = 0;
  cps = (comm_psend_t *) sc_mempool_alloc (gdata->psmem);
  cps->rank = -1;
  for(zz = 0; zz<numz;zz++){
      pad = (pdata_t *) sc_array_index(particle_data_copy, zz);
      pfn = pad->mpirank;
      if(pfn == gdata->mpirank){
          ++lremain;
          continue;
          }
    
       cps->rank = pfn;

       retval = sc_hash_insert_unique (gdata->psend, cps, &hfound);
            
       there = *((comm_psend_t **) hfound);
      
        if (!retval) {
          assert (there->message.elem_size == sizeof(integral_t));
          assert (there->message.elem_count > 0);
            }
         else{
             assert(there == cps);
             trank = (comm_prank_t *) sc_array_push (gdata->recevs);
             trank->rank = there->rank;
             trank->psend = there;
             sc_array_init (&there->message, sizeof(integral_t));
             cps = (comm_psend_t *) sc_mempool_alloc (gdata->psmem);
             cps->rank = -1;

             }

          integral = (integral_t *) sc_array_push(&there->message);
          integral->leftintegral = pad->leftintegral;
          integral->rightintegral = pad->rightintegral;
          integral->id = pad->id;
          integral->drifta = pad->drifta; 
          ++lsend;
      
      } 
      assert(numz == (size_t)(lsend+lremain));
      sc_mempool_free (gdata->psmem, cps);
      sc_array_sort (gdata->recevs, comm_prank_compare);
    }




void PelletSolver::communicateParticles_phase2(){
    
  int                 mpiret;
  int                 i;
  int                 num_receivers;
  int                 num_senders;
  int                 count, cucount;
  int                 msglen;
  sc_MPI_Request     *reqs;
  sc_array_t         *notif, *payl;
  sc_array_t         *arr;
  comm_psend_t       *cps;
  comm_prank_t       *trank;
    
  num_receivers = gdata->recevs->elem_count;
  notif = sc_array_new_count (sizeof (int), num_receivers);
  payl = sc_array_new_count (sizeof (int), num_receivers);   //payload
    
  for (i = 0; i < num_receivers; ++i) {

    trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);

    *(int *) sc_array_index_int (notif, i) = trank->rank;

    cps = trank->psend;
  
    assert(cps->rank == trank->rank);

    arr = &cps->message;
 
    *(int *) sc_array_index_int (payl, i) = (int) arr->elem_count;
 
  
    }
    

  sc_notify_ext (notif, NULL, payl, NULL, gdata->mpicomm);
    
  assert (payl->elem_count == notif->elem_count);
    
  num_senders = (int) notif->elem_count;


  cucount = 0;
  for (i = 0; i < num_senders; ++i) {
    cucount += *(int *) sc_array_index_int (payl, i);
  }
  prebuf_integral = sc_array_new_count (sizeof(integral_t), cucount);
    
  gdata->recv_req = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);

  cucount = 0;

  for (i = 0; i < num_senders; ++i) {
    count = *(int *) sc_array_index_int (payl, i);
    msglen = count * (int) sizeof(integral_t);
    mpiret = sc_MPI_Irecv
      (sc_array_index (prebuf_integral, cucount), msglen, sc_MPI_BYTE,
       *(int *) sc_array_index_int (notif, i), 888, gdata->mpicomm,
       (sc_MPI_Request *) sc_array_index_int (gdata->recv_req, i));
    SC_CHECK_MPI (mpiret);
    cucount += count;
  }
    
  assert(cucount == (int) prebuf_integral->elem_count);

  sc_array_destroy_null (&notif);
  sc_array_destroy_null (&payl);
    
    gdata->send_req = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
    for (i = 0; i < num_receivers; ++i) {
      trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);
      cps = trank->psend;
      arr = &cps->message;
      msglen = (int) (arr->elem_count * arr->elem_size);
      mpiret = sc_MPI_Isend
        (arr->array, msglen, sc_MPI_BYTE, cps->rank, 888,
         gdata->mpicomm, (sc_MPI_Request *) sc_array_index_int (gdata->send_req, i));
      SC_CHECK_MPI (mpiret);
    }

  reqs = (sc_MPI_Request *) sc_array_index (gdata->recv_req,0);
  
  mpiret = sc_MPI_Waitall (num_senders, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&gdata->recv_req);

  num_receivers = (int) gdata->recevs->elem_count;
  reqs = (sc_MPI_Request *) sc_array_index (gdata->send_req,0),
  mpiret = sc_MPI_Waitall (num_receivers, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&gdata->send_req);
    
  for (i = 0; i < num_receivers; ++i) {
    trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);
    cps = trank->psend;
    sc_array_reset (&cps->message);
  }
    
  sc_array_destroy_null (&gdata->recevs);
  sc_hash_destroy (gdata->psend);

  gdata->psend = NULL;
  sc_mempool_destroy (gdata->psmem);
  gdata->psmem = NULL;
    
    }

void PelletSolver::writeIntegralValue(){
    
    size_t li, num;
    pdata_t *pad, *pad_copy;
    int id;
    integral_t *integral;
    num = prebuf_integral->elem_count;
    size_t counter = 0;
    for(li = 0; li<num ; li++){
        integral = (integral_t *) sc_array_index(prebuf_integral, li);
        id = integral->id;
        pad = (pdata_t *) sc_array_index_int(gdata->particle_data,id);
        pad->leftintegral = integral->leftintegral;
        pad->rightintegral = integral->rightintegral;
        pad->drifta = integral->drifta;
        counter ++;
        }
   
    num = particle_data_copy->elem_count;
    for(li = 0; li<num ;li++){
        pad_copy = (pdata_t *) sc_array_index(particle_data_copy,li);
        if(pad_copy->mpirank != gdata->mpirank)
            continue;
        id = pad_copy->id;
        pad = (pdata_t *) sc_array_index_int(gdata->particle_data,id);
        pad->leftintegral = pad_copy->leftintegral;
        pad->rightintegral = pad_copy->rightintegral;
        pad->drifta = pad_copy->drifta;
        counter ++;
        }
        
 //   assert(counter == gdata->lfluidnum);
    
    }

double getDissociationFraction_Deuterium2(double pressure, double density); 
void PelletSolver::computeHeatDeposition( double dt){
    
    pdata_t *pad;
    
    size_t li, lnump = gdata->particle_data->elem_count;
    double e = heatK*(2.99792458e7)/100;
    double lnLambda = log(2*teinf/I*sqrt(exp(1)/2));
    double tauleft; 
    double tauright;
    double tauinf;
    double taueff;
    double uleft; 
    double uright;
    double qinf;
    double guleft;
    double guright; 
    double nt;
    static double m_fCurrentTime = gdata->initializer->getStartTime(); 
    double k_warmup;

    m_fCurrentTime  += dt;

    if ( m_fCurrentTime > 0.04)
      k_warmup = 1.0;
    else
      k_warmup = m_fCurrentTime/0.04;



    for(li = 0; li<lnump ; li++){
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary)
            continue;
        tauleft = pad->leftintegral/mass*Z;
        tauright = pad->rightintegral/mass*Z;
        tauinf = heatK*heatK*teinf*teinf/(8.0*3.1416*e*e*e*e*lnLambda);
		taueff = tauinf/(0.625+0.55*sqrt(one_plus_Zstar)); //tauinf*sqrt(2./(1.+ZNe));//tauinf/(0.625+0.55*sqrt(one_plus_Zstar));
        uleft = tauleft/taueff;
        uright = tauright/taueff;                               
        qinf=sqrt(2.0/3.1416/masse)*neinf*pow(heatK*teinf,1.5);
        guleft = sqrt(uleft)*Bessel_K1(sqrt(uleft))/4;
        guright = sqrt(uright)*Bessel_K1(sqrt(uright))/4;
        //double fd = getDissociationFraction_Deuterium2(pad->pressure,1./pad->volume);
       // nt = 1./pad->volume/(2*mass*(1-fd) + mass*fd);
        nt=1.0/pad->volume/mass;
        if(heatingmodel == 1){
            pad->deltaq = qinf*nt*Z/tauinf*(guleft+guright)*k_warmup;
            pad->qplusminus = qinf*0.5*(uleft*Bessel_Kn(2,sqrt(uleft))+uright*Bessel_Kn(2,sqrt(uright)))*k_warmup;
        } 
      
        else if(heatingmodel == 0){
            pad->deltaq = qinf*nt*Z/tauinf*guleft*k_warmup;
            pad->qplusminus = qinf*0.5*uleft*Bessel_Kn(2,sqrt(uleft))*k_warmup;
            }
        // if(pad->qplusminus == 0)
         //   cout<<pad->leftintegral<<" "<<pad->rightintegral<<endl;
        
        }

}
float           Bessel_I0(
	float           x)
{
        float   p1 = 1.0;
        float   p2 = 3.5156229;
        float   p3 = 3.0899424;
        float   p4 = 1.2067492;
        float   p5 = 0.2659732;
        float   p6 = 0.360768e-1;
        float   p7 = 0.45813e-2;
	
        float   q1 = 0.39894228;
        float   q2 = 0.1328592e-1;
        float   q3 = 0.225319e-2;
        float   q4 = -0.157565e-2;
        float   q5 = 0.916281e-2;
        float   q6 = -0.2057706e-1;
        float   q7 = 0.2635537e-1;
        float   q8 = -0.1647633e-1;
        float   q9 = 0.392377e-2;
	
	float   ax, y, value;
	
	if (fabs(x) < 3.75)
	  {
	    y = (x/3.75)*(x/3.75);//sqr
	    value = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
	  }
	else
	  {
	    ax = fabs(x);
	    y = 3.75/ax;

	    value = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
	  }

	return value;
}


/* Bessel_I1 returns the modifies Bessel function I1(x) of positive real x  */

float           Bessel_I1(
	float           x)
{
        float   p1 = 0.5;
        float   p2 = 0.87890594;
        float   p3 = 0.51498869;
        float   p4 = 0.15084934;
        float   p5 = 0.2658733e-1;
        float   p6 = 0.301532e-2;
        float   p7 = 0.32411e-3;
	
        float   q1 = 0.39894228;
        float   q2 = -0.3988024e-1;
        float   q3 = -0.362018e-2;
        float   q4 = 0.163801e-2;
        float   q5 = -0.1031555e-1;
        float   q6 = 0.2282967e-1;
        float   q7 = -0.2895312e-1;
        float   q8 = 0.1787654e-1;
        float   q9 = -0.420059e-2;
	
	float   ax, y, value;
	
	if (fabs(x) < 3.75)
	  {
	    y = (x/3.75)*(x/3.75);//sqr
	    value = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
	  }
	else
	  {
	    ax = fabs(x);
	    y = 3.75/ax;

	    value = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
	    if (x < 0)
	      value *= -1.0;
	  }
	return value;
}

/* Bessel_K0 returns the modifies Bessel function K0(x) of positive real x  */

float           Bessel_K0(
	float           x)
{
        float   p1 = -0.57721566;
	float   p2 = 0.4227842;
	float   p3 = 0.23069756;
	float   p4 = 0.348859e-1;
	float   p5 = 0.262698e-2;
	float   p6 = 0.1075e-3;
	float   p7 = 0.74e-5;

	float   q1 = 1.25331414;
	float   q2 = -0.7832358e-1;
	float   q3 = 0.2189568e-1;
	float   q4 = -0.1062446e-1;
	float   q5 = 0.587872e-2;
	float   q6 = -0.25154e-2;
	float   q7 = 0.53208e-3;

	float   y, value;

	if (x <= 2.0)
	  {
	    y = x*x/4.0;
	    value = (-log(x/2.0)*Bessel_I0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
	  }
	else
	  {
	    y = 2.0/x;
	    value = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
	  }
	return value;
}

/* Bessel_K1 returns the modifies Bessel function K1(x) of positive real x  */

float           Bessel_K1(
	float           x)
{
        float   p1 = 1.0;
	float   p2 = 0.15443144;
	float   p3 = -0.67278579;
	float   p4 = -0.18156897;
	float   p5 = -0.01919402;
	float   p6 = -0.110404e-2;
	float   p7 = -0.4686e-4;

	float   q1 = 1.25331414;
	float   q2 = 0.23498619;
	float   q3 = -0.3655620e-1;
	float   q4 = 0.1504268e-1;
	float   q5 = -0.780353e-2;
	float   q6 = 0.325614e-2;
	float   q7 = -0.68245e-3;

	float   y, value;

	if (x <= 2.0)
	  {
	    y = x*x/4.0;
	    value = (log(x/2.0)*Bessel_I1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
	  }
	else
	  {
	    y = 2.0/x;
	    value = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
	  }
	return value;
}
				       		       

/* Bessel_Kn returns the modifies Bessel function Kn(x) of positive real x for n >= 2 */

float           Bessel_Kn(
        int             n,
	float           x)
{
        int    j;
        float  bk, bkm, bkp, tox;

	if (n < 2)
	  {
	    printf("Error in Bessel_Kn(), the order n < 2: n = %d\n",n);
	    assert(false);
	    return 0;
	  }

	tox = 2.0/x;
	bkm = Bessel_K0(x);
	bk = Bessel_K1(x);

	for (j = 1; j < n; j++)
	  {
	    bkp = bkm + j*tox*bk;
	    bkm = bk;
	    bk = bkp;
	  }

	return bk;
}
//Pellet heat deposition calculation




void PelletSolver::heatingModel(double dt){
        
           gdata->setParticleIDAndRank();  
           if(heatingmodel == 1){    
            
            gdata->computeBoundingBoxLength();  //compute max(y,z);
            
            prerun();
            build_quadtree();
          
            presearch2d();
            
            packParticles();
            
            communicateParticles();
            
            postsearch2d();
         
            adaptQuadtree();
            regroupParticles2d();

            partitionParticles2d();
            
            MPI_Barrier(gdata->mpicomm);
            
            P4EST_GLOBAL_ESSENTIALF ("generating computation nodes.\n");
            generateComputationNode(); 
            
            computeNodeStates(); 
            regroupComputationNodeWithStates(nodelist);
            P4EST_GLOBAL_ESSENTIALF ("regroup computation nodes.\n");
            reorderComputationNode();
            P4EST_GLOBAL_ESSENTIALF ("density integral nodes.\n");
            densityIntegralForNode(); 
 //           computeDriftAcceleration(); 
            P4EST_GLOBAL_ESSENTIALF ("density integral particles.\n");
            computeDensityIntegral_interpolation();
           // computeDensityIntegral();
             
            P4EST_GLOBAL_ESSENTIALF ("communicate particles.\n");
            packParticles_phase2();
            
            communicateParticles_phase2();
            
            P4EST_GLOBAL_ESSENTIALF ("write values.\n");
            MPI_Barrier(gdata->mpicomm);
            writeIntegralValue();
            
            MPI_Barrier(gdata->mpicomm);
            
            computeHeatDeposition(dt);
       
            MPI_Barrier(gdata->mpicomm);
            
            destoryQuadtree();
            
            }
        else if(heatingmodel == 0){
            
             computeDensityIntegral1D();
             
             packParticles_phase2();
             communicateParticles_phase2();
             
             writeIntegralValue();
            
            computeHeatDeposition(dt);
            
            sc_array_destroy_null(&prebuf_integral);  
            sc_array_destroy_null(&particle_data_copy); 
            }
        else{
            cout<<"Heating Model does not exsit."<<endl;
            assert(false);
            }
    }


static int
compareRadialDis (const void *p1, const void *p2)
{
  int t = 0;
  pdata_t                * i1 = (pdata_t *) p1;
  pdata_t                * i2 = (pdata_t *) p2;
  
  double d1 = i1->xyz[0]*i1->xyz[0]+i1->xyz[1]*i1->xyz[1]+i1->xyz[2]*i1->xyz[2]; 
  double d2 = i2->xyz[0]*i2->xyz[0]+i2->xyz[1]*i2->xyz[1]+i2->xyz[2]*i2->xyz[2]; 
  if((d1-d2)>0)
      t = 1;
  return t;
}
void PelletSolver::computeDensityIntegral1D(){
    size_t gpnum = gdata->gpnum;
    int lpnum = (int)gdata->particle_data->elem_count;
    
    int lpnum_count[gdata->mpisize];
    int msglen;
    pdata_t *pad;
    
    particle_data_copy = sc_array_new(sizeof(pdata_t)); 
    if(gdata->mpirank != 0){
        MPI_Send(&lpnum,1,MPI_INT,0,666,gdata->mpicomm);        
    }
    else{
        
        for(int i=1;i<(int)gdata->mpisize;i++){
            MPI_Recv(&lpnum_count[i], 1, MPI_INT, i, 666, gdata->mpicomm,MPI_STATUS_IGNORE);
        }
        
        lpnum_count[0] = lpnum;
        
    }
    
    if(gdata->mpirank == 0){
        sc_array_init_count(particle_data_copy,sizeof(pdata_t),gpnum);
        int offset = lpnum_count[0];
        sc_array_copy_into(particle_data_copy,0,gdata->particle_data);
         
        for(int i=1;i<gdata->mpisize;i++){
            msglen = lpnum_count[i]*(int)sizeof(pdata_t);
            MPI_Recv(sc_array_index_int(particle_data_copy,offset),msglen,MPI_BYTE,i,COMM_TAG_PART,gdata->mpicomm,MPI_STATUS_IGNORE);
            offset += lpnum_count[i];
        }
    }
    else{
        msglen = lpnum*(int)sizeof(pdata_t);
        MPI_Send(gdata->particle_data->array,msglen,MPI_BYTE,0,COMM_TAG_PART,gdata->mpicomm);
    
    }
    if(gdata->mpirank == 0){
        
        sc_array_sort(particle_data_copy,compareRadialDis);
        double integral = 0;
        for(int i = gpnum-1; i>=0; i--){
           pad = (pdata_t *)sc_array_index_int(particle_data_copy,i);
           double r = pad->xyz[0]*pad->xyz[0]+pad->xyz[1]*pad->xyz[1]+pad->xyz[2]*pad->xyz[2];
           double temp = pad->mass/4./M_PI/r;
           pad->leftintegral = integral + 0.5*temp; 
           integral += temp; 
        }
    }   

}

void PelletSolver::updateStatesByLorentzForce(double dt,const double invelocity, const double inpressure, const double involume, const double insoundspeed, double* outvelocity, double* outpressure){
    
    if(heatingmodel == 0 || fabs(magneticfield - 0.)<1e-8)
        return;
    double LF, d_v;
    double delta_E;
    double cond;
    double density, sc, pressure;
    double gamma;
    double v = invelocity;
    double vold = v;
    int numtimesteps = 1;
    
    density = 1./involume;
    pressure = inpressure;
    sc = insoundspeed;
    gamma = (sc*sc*density/pressure); 
    for(int i=0; i<numtimesteps; i++){
        cond = gdata->eos->getElectricConductivity(pressure,density);   
       /* 
        LF = -cond*magneticfield*magneticfield/c_light/c_light*v;
        d_v =  LF/density*dt/numtimesteps;
        v += d_v;
        
        if(fabs(d_v) > fabs(v)){
            d_v = -v;
            v = 0;
        }

        *outvelocity = (*outvelocity) + d_v;
        if( v==0 )
            break;
         */   
    

    LF = -cond*magneticfield*magneticfield/c_light/c_light/density*dt/numtimesteps;
    v = (*outvelocity);
    (*outvelocity) = (*outvelocity)*exp(LF);
    delta_E = 0.5*density*(v*v-(*outvelocity)*(*outvelocity));
    (*outpressure) = (*outpressure) + (gamma-1)*delta_E;
    
    } 
//   delta_E = 0.5*density*(vold*vold-v*v); 
//   (*outpressure) = (*outpressure) + (gamma-1)*delta_E;
}
void PelletSolver::updateStatesByLorentzForce(double dt){
    if(heatingmodel == 0 || fabs(magneticfield - 0.)<1e-8)
        return;
    size_t li, lpnum = gdata->particle_data->elem_count;
    pdata_t *pad;
    double vy ,vz;
    double density, sc, pressure;
    double cond;
    double LFy, LFz;
    double delta_E;
    double rad_cool, T;
    int numtimestep = 1;
    for(li = 0; li<lpnum; li++){
        
        pad = (pdata_t *) sc_array_index(gdata->particle_data, li);
        
        if(pad->ifboundary) continue;
        for(int i=0;i<numtimestep;i++){
            vy = pad->v[1]-pad->driftv;
            double vyold = pad->v[1];
            vz = pad->v[2];
            density = 1./pad->volume;
            sc = pad->soundspeed;
            pressure = pad->pressure;
            cond = gdata->eos->getElectricConductivity(pressure,density);
            T = gdata->eos->getTemperature(pressure,density);
            rad_cool = neon_radiation_power_density(density, T);
            LFy = -cond*magneticfield*magneticfield/c_light/c_light/density*dt/numtimestep;
            LFz = -cond*magneticfield*magneticfield/c_light/c_light/density*dt/numtimestep;
            pad->v[1] = vy*(exp(LFy))+pad->driftv;
            pad->v[2] = vz*exp(LFz);
            delta_E = 0.5*density*((vyold*vyold+vz*vz) - (pad->v[1]*pad->v[1]+pad->v[2]*pad->v[2]));
            
//            pad->pressure = pad->pressure - (sc*sc*density/pressure -1 )*rad_cool*dt/numtimestep;
            if(delta_E > 0)
                pad->pressure = pad->pressure + (sc*sc*density/pressure-1)*delta_E;
            
            pad->soundspeed = gdata->eos->getSoundSpeed(pad->pressure, density);
        }
    }

}
/*
void PelletSolver::updateStatesByLorentzForce(double dt){
    if(heatingmodel == 0 || magneticfield == 0)
        return;
    size_t li, lpnum = gdata->particle_data->elem_count;
    pdata_t *pad;
    double vy ,vz;
    double density, sc, pressure;
    double d_vy, d_vz;
    double cond;
    double LFy, LFz;
    double delta_E;
    double rad_cool, T;
    int numtimestep = 500;
    for(li = 0; li<lpnum; li++){
        
        pad = (pdata_t *) sc_array_index(gdata->particle_data, li);
        if(pad->ifboundary) continue;
        vy = pad->oldv[1];
        vz = pad->oldv[2];
        density = 1./pad->volumeT1;
        sc = pad->soundspeedT1;
        pressure = pad->pressureT1;
        
        for(int i=0;i<numtimestep;i++){
            cond = gdata->eos->getElectricConductivity(pressure,density);
            T = gdata->eos->getTemperature(pressure,density);
            rad_cool = neon_radiation_power_density(density, T);
            LFy = -cond*magneticfield*magneticfield/c_light/c_light*vy;
            LFz = -cond*magneticfield*magneticfield/c_light/c_light*vz;
            d_vy = LFy/density*dt/numtimestep; 
            d_vz = LFz/density*dt/numtimestep; 
            pad->v[1] = pad->v[1] + d_vy; //pad->oldv[1]*LFy;// pad->v[1] * exp(LFy);
            pad->v[2] = pad->v[2] + d_vz; // exp(LFz);//pad->oldv[2]*LFz;//pad->v[2] * exp(LFz);
            delta_E = cond*magneticfield*magneticfield/c_light/c_light*(vy*vy+vz*vz)*dt/numtimestep;//0.5*density*((vy*vy+vz*vz) - (pad->v[1]*pad->v[1]+pad->v[2]*pad->v[2]));
            vy += d_vy; 
            vz += d_vz; 
            
           if(fabs(d_vy) > fabs(pad->v[1]))
            {
                pad->v[1] = 0;
                vy = 0;
            }
           if(fabs(d_vz) > fabs(pad->v[2]))
           {
                pad->v[2] = 0;
                vz = 0;
           }

            
            pad->pressure = pad->pressure - (sc*sc*density/pressure -1 )*rad_cool*dt/numtimestep;
            if(delta_E > 0){
                pad->pressure = pad->pressure + (sc*sc*density/pressure-1)*delta_E;
            }
            
            if(pad->pressure<0){
                pad->pressure = pad->pressure + (sc*sc*density/pressure -1 )*rad_cool*dt/numtimestep; 
                pad->soundspeed = gdata->eos->getSoundSpeed(pad->pressure, 1./pad->volume);
            }
            pad->soundspeed = gdata->eos->getSoundSpeed(pad->pressure, 1./pad->volume);
        }

    }

}
*/

void PelletSolver::neonRadiationCooling(double dt){
    
    if(gdata->eoschoice != 2 || heatingmodel == 0 || materialnumber != 0)
        return;
    size_t li, lpnum = gdata->particle_data->elem_count;
    double T, rad_cool;
    double sc, pressure, density;
    int timestep = 1;
    pdata_t *pad;
    for(li = 0; li<lpnum; li++){
        pad = (pdata_t*) sc_array_index(gdata->particle_data, li);
        if(pad->ifboundary) continue;
        density = 1./pad->volume;
        for(int i=0;i<timestep;i++){
            sc = pad->soundspeed;
            pressure = pad->pressure;
            T = gdata->eos->getTemperature(pressure,density);
            rad_cool = neon_radiation_power_density(density, T);
            pad->radcool = rad_cool;
            pad->pressure = pressure - (sc*sc*density/pressure -1 )*rad_cool*dt/timestep;
            if(pad->pressure < 0){
                pad->pressure = pressure;
                break;
            }
            pad->soundspeed = gdata->eos->getSoundSpeed(pad->pressure, density);
        }
    }

}








double PelletSolver::neon_radiation_power_density(
	double rho,
	double T)
{
    double temp,Tl,Tr,Z,Zl,Zr,L,Ll,Lr,k,n_rho;
	int j;
	// radiation power density = n^2 * Z * L * 1.e-15 [W/cm^3]

	n_rho = rho/3.350971e-23; // number density of neon nuclei

	//find lower index in the Neon radiation table
	temp = floor(10*T)/10.0;
	if (temp < 1.0)
	  return 0.0;
	else if (temp < 5.0)
	  j = round(10*(temp - 1));
	else if (temp <= 30.0)
	  j = round((temp - 5)/0.5) + 40;
	else
	  j = 90;
	
	Tl = neon_radiation_data(j,0);
	Tr = neon_radiation_data(j+1,0);
	
	Zl = neon_radiation_data(j,1);
	Zr = neon_radiation_data(j+1,1);
	
	Ll = neon_radiation_data(j,2);
	Lr = neon_radiation_data(j+1,2);

	k = (T - Tl)/(Tr - Tl);
	Z = (1-k)*Zl + k*Zr;
	L = (1-k)*Ll + k*Lr;
	 
	return n_rho*n_rho*Z*L*1.e-17; // power density using ms time unit
}

double PelletSolver::neon_radiation_data(
	int i,
	int j)
{
 static  double rad[91][3];
  bool First = true;
 if(First){
     First = false;
rad[0][0] =   1.000000e+00;
rad[1][0] =   1.100000e+00;
rad[2][0] =   1.200000e+00;
rad[3][0] =   1.300000e+00;
rad[4][0] =   1.400000e+00;
rad[5][0] =   1.500000e+00;
rad[6][0] =   1.600000e+00;
rad[7][0] =   1.700000e+00;
rad[8][0] =   1.800000e+00;
rad[9][0] =   1.900000e+00;
rad[10][0] =  2.000000e+00;
rad[11][0] =  2.100000e+00;
rad[12][0] =  2.200000e+00;
rad[13][0] =  2.300000e+00;
rad[14][0] =  2.400000e+00;
rad[15][0] =  2.500000e+00;
rad[16][0] =  2.600000e+00;
rad[17][0] =  2.700000e+00;
rad[18][0] =  2.800000e+00;
rad[19][0] =  2.900000e+00;
rad[20][0] =  3.000000e+00;
rad[21][0] =  3.100000e+00;
rad[22][0] =  3.200000e+00;
rad[23][0] =  3.300000e+00;
rad[24][0] =  3.400000e+00;
rad[25][0] =  3.500000e+00;
rad[26][0] =  3.600000e+00;
rad[27][0] =  3.700000e+00;
rad[28][0] =  3.800000e+00; 
rad[29][0] =  3.900000e+00; 
rad[30][0] =  4.000000e+00; 
rad[31][0] =  4.100000e+00; 
rad[32][0] =  4.200000e+00; 
rad[33][0] =  4.300000e+00; 
rad[34][0] =  4.400000e+00; 
rad[35][0] =  4.500000e+00; 
rad[36][0] =  4.600000e+00; 
rad[37][0] =  4.700000e+00; 
rad[38][0] =  4.800000e+00; 
rad[39][0] =  4.900000e+00; 
rad[40][0] =  5.000000e+00; 
rad[41][0] =  5.500000e+00; 
rad[42][0] =  6.000000e+00; 
rad[43][0] =  6.500000e+00; 
rad[44][0] =  7.000000e+00; 
rad[45][0] =  7.500000e+00; 
rad[46][0] =  8.000000e+00; 
rad[47][0] =  8.500000e+00; 
rad[48][0] =  9.000000e+00; 
rad[49][0] =  9.500000e+00; 
rad[50][0] =  1.000000e+01; 
rad[51][0] =  1.050000e+01; 
rad[52][0] =  1.100000e+01; 
rad[53][0] =  1.150000e+01; 
rad[54][0] =  1.200000e+01; 
rad[55][0] =  1.250000e+01; 
rad[56][0] =  1.300000e+01; 
rad[57][0] =  1.350000e+01; 
rad[58][0] =  1.400000e+01; 
rad[59][0] =  1.450000e+01; 
rad[60][0] =  1.500000e+01; 
rad[61][0] =  1.550000e+01; 
rad[62][0] =  1.600000e+01; 
rad[63][0] =  1.650000e+01; 
rad[64][0] =  1.700000e+01; 
rad[65][0] =  1.750000e+01; 
rad[66][0] =  1.800000e+01; 
rad[67][0] =  1.850000e+01; 
rad[68][0] =  1.900000e+01; 
rad[69][0] =  1.950000e+01; 
rad[70][0] =  2.000000e+01; 
rad[71][0] =  2.050000e+01; 
rad[72][0] =  2.100000e+01; 
rad[73][0] =  2.150000e+01; 
rad[74][0] =  2.200000e+01; 
rad[75][0] =  2.250000e+01; 
rad[76][0] =  2.300000e+01; 
rad[77][0] =  2.350000e+01; 
rad[78][0] =  2.400000e+01; 
rad[79][0] =  2.450000e+01; 
rad[80][0] =  2.500000e+01; 
rad[81][0] =  2.550000e+01; 
rad[82][0] =  2.600000e+01; 
rad[83][0] =  2.650000e+01; 
rad[84][0] =  2.700000e+01; 
rad[85][0] =  2.750000e+01; 
rad[86][0] =  2.800000e+01; 
rad[87][0] =  2.850000e+01; 
rad[88][0] =  2.900000e+01; 
rad[89][0] =  2.950000e+01; 
rad[90][0] =  3.000000e+01; 
rad[0][1] =   2.855628e-02;
rad[1][1] =   1.463662e-01;
rad[2][1] =   4.330838e-01;
rad[3][1] =   7.347827e-01;
rad[4][1] =   8.955859e-01;
rad[5][1] =   9.587141e-01;
rad[6][1] =   9.825625e-01;
rad[7][1] =   9.920914e-01;
rad[8][1] =   9.962807e-01;
rad[9][1] =   9.985018e-01;
rad[10][1] =  1.000286e+00;
rad[11][1] =  1.002707e+00;
rad[12][1] =  1.006990e+00;
rad[13][1] =  1.014886e+00;
rad[14][1] =  1.028927e+00;
rad[15][1] =  1.052489e+00;
rad[16][1] =  1.089469e+00;
rad[17][1] =  1.143327e+00;
rad[18][1] =  1.215524e+00;
rad[19][1] =  1.303975e+00;
rad[20][1] =  1.402701e+00;
rad[21][1] =  1.503367e+00;
rad[22][1] =  1.597987e+00;
rad[23][1] =  1.681105e+00;
rad[24][1] =  1.750426e+00;
rad[25][1] =  1.806171e+00;
rad[26][1] =  1.849998e+00;
rad[27][1] =  1.884091e+00;
rad[28][1] =  1.910620e+00;
rad[29][1] =  1.931495e+00;
rad[30][1] =  1.948310e+00;
rad[31][1] =  1.962368e+00;
rad[32][1] =  1.974739e+00;
rad[33][1] =  1.986323e+00;
rad[34][1] =  1.997900e+00;
rad[35][1] =  2.010166e+00;
rad[36][1] =  2.023761e+00;
rad[37][1] =  2.039274e+00;
rad[38][1] =  2.057237e+00;
rad[39][1] =  2.078115e+00;
rad[40][1] =  2.102280e+00;
rad[41][1] =  2.275239e+00;
rad[42][1] =  2.498072e+00;
rad[43][1] =  2.690690e+00;
rad[44][1] =  2.819845e+00;
rad[45][1] =  2.900079e+00;
rad[46][1] =  2.954544e+00;
rad[47][1] =  3.000644e+00;
rad[48][1] =  3.049845e+00;
rad[49][1] =  3.109150e+00;
rad[50][1] =  3.181565e+00;
rad[51][1] =  3.266204e+00;
rad[52][1] =  3.358960e+00;
rad[53][1] =  3.454218e+00;
rad[54][1] =  3.546805e+00;
rad[55][1] =  3.633317e+00;
rad[56][1] =  3.712484e+00;
rad[57][1] =  3.784763e+00;
rad[58][1] =  3.851660e+00;
rad[59][1] =  3.915086e+00;
rad[60][1] =  3.976912e+00;
rad[61][1] =  4.038694e+00;
rad[62][1] =  4.101567e+00;
rad[63][1] =  4.166216e+00;
rad[64][1] =  4.232939e+00;
rad[65][1] =  4.301736e+00;
rad[66][1] =  4.372415e+00;
rad[67][1] =  4.444695e+00;
rad[68][1] =  4.518282e+00;
rad[69][1] =  4.592917e+00;
rad[70][1] =  4.668399e+00;
rad[71][1] =  4.744582e+00;
rad[72][1] =  4.821367e+00;
rad[73][1] =  4.898686e+00;
rad[74][1] =  4.976481e+00;
rad[75][1] =  5.054694e+00;
rad[76][1] =  5.133261e+00;
rad[77][1] =  5.212099e+00;
rad[78][1] =  5.291106e+00;
rad[79][1] =  5.370152e+00;
rad[80][1] =  5.449081e+00;
rad[81][1] =  5.527706e+00;
rad[82][1] =  5.605804e+00;
rad[83][1] =  5.683127e+00;
rad[84][1] =  5.759401e+00;
rad[85][1] =  5.834343e+00;
rad[86][1] =  5.907667e+00;
rad[87][1] =  5.979103e+00;
rad[88][1] =  6.048405e+00;
rad[89][1] =  6.115366e+00;
rad[90][1] =  6.179824e+00;
rad[0][2] =   1.294421e-17;
rad[1][2] =   5.032943e-17;
rad[2][2] =   1.175009e-16;
rad[3][2] =   1.618292e-16;
rad[4][2] =   1.633444e-16;
rad[5][2] =   1.487624e-16;
rad[6][2] =   1.340048e-16;
rad[7][2] =   1.251299e-16;
rad[8][2] =   1.267496e-16;
rad[9][2] =   1.459262e-16;
rad[10][2] =  1.946532e-16;
rad[11][2] =  2.921758e-16;
rad[12][2] =  4.680815e-16;
rad[13][2] =  7.667701e-16;
rad[14][2] =  1.254301e-15;
rad[15][2] =  2.028582e-15;
rad[16][2] =  3.232404e-15;
rad[17][2] =  5.064973e-15;
rad[18][2] =  7.782106e-15;
rad[19][2] =  1.167420e-14;
rad[20][2] =  1.702276e-14;
rad[21][2] =  2.405655e-14;
rad[22][2] =  3.293476e-14;
rad[23][2] =  4.376405e-14;
rad[24][2] =  5.663246e-14;
rad[25][2] =  7.163800e-14;
rad[26][2] =  8.890362e-14;
rad[27][2] =  1.085818e-13;
rad[28][2] =  1.308546e-13;
rad[29][2] =  1.559355e-13;
rad[30][2] =  1.840726e-13;
rad[31][2] =  2.155569e-13;
rad[32][2] =  2.507330e-13;
rad[33][2] =  2.900121e-13;
rad[34][2] =  3.338860e-13;
rad[35][2] =  3.829394e-13;
rad[36][2] =  4.378594e-13;
rad[37][2] =  4.993941e-13;
rad[38][2] =  5.685169e-13;
rad[39][2] =  6.461584e-13;
rad[40][2] =  7.333081e-13;
rad[41][2] =  1.340709e-12;
rad[42][2] =  2.244150e-12;
rad[43][2] =  3.337638e-12;
rad[44][2] =  4.502229e-12;
rad[45][2] =  5.694331e-12;
rad[46][2] =  6.914699e-12;
rad[47][2] =  8.175398e-12;
rad[48][2] =  9.488725e-12;
rad[49][2] =  1.086383e-11;
rad[50][2] =  1.230390e-11;
rad[51][2] =  1.380364e-11;
rad[52][2] =  1.534886e-11;
rad[53][2] =  1.691982e-11;
rad[54][2] =  1.849633e-11;
rad[55][2] =  2.006210e-11;
rad[56][2] =  2.160670e-11;
rad[57][2] =  2.312510e-11;
rad[58][2] =  2.461607e-11;
rad[59][2] =  2.608032e-11;
rad[60][2] =  2.751923e-11;
rad[61][2] =  2.893380e-11;
rad[62][2] =  3.032439e-11;
rad[63][2] =  3.169057e-11;
rad[64][2] =  3.303126e-11;
rad[65][2] =  3.434488e-11;
rad[66][2] =  3.562955e-11;
rad[67][2] =  3.688312e-11;
rad[68][2] =  3.810312e-11;
rad[69][2] =  3.928663e-11;
rad[70][2] =  4.043005e-11;
rad[71][2] =  4.152896e-11;
rad[72][2] =  4.257794e-11;
rad[73][2] =  4.357054e-11;
rad[74][2] =  4.449942e-11;
rad[75][2] =  4.535653e-11;
rad[76][2] =  4.613346e-11;
rad[77][2] =  4.682193e-11;
rad[78][2] =  4.741422e-11;
rad[79][2] =  4.790376e-11;
rad[80][2] =  4.828558e-11;
rad[81][2] =  4.855674e-11;
rad[82][2] =  4.871657e-11;
rad[83][2] =  4.876686e-11;
rad[84][2] =  4.871173e-11;
rad[85][2] =  4.855750e-11;
rad[86][2] =  4.831223e-11;
rad[87][2] =  4.798535e-11;
rad[88][2] =  4.758706e-11;
rad[89][2] =  4.712785e-11;
rad[90][2] =  4.661804e-11;
 }
 return rad[i][j];
}

bool PelletSolver::ifInsidePellet(double x, double y, double z){
    bool ifin = false;
    double pr = pelletradius;
    if((x*x+y*y+z*z)<pr*pr)
        ifin = true;

    return ifin;

}

void PelletSolver::computeBoundaryCondition(Global_Data *g, double dt, double dx){
/*   
   static bool first=true;
   if(first){
       first = false;
       return;
       }
*/
   pdata_t *pad;
   size_t li, lpnum = g->particle_data->elem_count;
   double x, y, z, dr;
   double pr = pelletradius;
   double xcen=0;
   double ycen=0;
   double zcen=0;
   
   int counter_nei = 0;
   int counter_g = 0;
   double qsum = 0;
   double qsum_g = 0;
   double vol = 0;
   double vol_g = 0;
   double ur = 0;
   double ur_g = 0;   
   double pres = 0;
   double pres_g = 0;
   double ss = 0;   
   double ss_g = 0;
   double rad = 0;
   double rad_g = 0;

   for(li = 0; li<lpnum; li++){
       pad = (pdata_t *)sc_array_index(g->particle_data,li);
       if(pad->ifboundary)
           continue;
       x = pad->xyz[0];
       y = pad->xyz[1];
       z = pad->xyz[2];

       dr = sqrt((x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen))-pr;
       if(dr>1)  continue;
       double sound = (pad->soundspeed+pad->soundspeedT1)/2;  
       double vx = (pad->v[0]+pad->oldv[0])/2;
       double vy = (pad->v[1]+pad->oldv[1])/2;
       double vz = (pad->v[2]+pad->oldv[2])/2;;
       double vel = sqrt(vx*vx+vy*vy+vz*vz);
       double ssx,ssy,ssz;
       if(vel==0)
       	    ssx=ssy=ssz=0;
       else{
	    ssx = sound*(x-xcen)/(dr+pr);
	    ssy = sound*(y-ycen)/(dr+pr);
	    ssz = sound*(z-zcen)/(dr+pr);
       } 
       double d_x2 = (x-xcen) - (ssx)*dt;
       double d_y2 = (y-ycen) - (ssy)*dt;
       double d_z2 = (z-zcen) - (ssz)*dt;
       double r_shift = (d_x2*d_x2+d_y2*d_y2+d_z2*d_z2);
       if(r_shift<(pr+dx)*(pr+dx) && r_shift>(pr-dx)*(pr-dx) && dr>0){
            qsum += pad->qplusminus/pad->localspacing;
            vol += pad->volume; 
	    pres += pad->pressure;
	    ss += pad->soundspeed;
	    ur += (x*vx+y*vy+z*vz)/(dr+pr); 
	    rad += pad->localspacing;
	    counter_nei ++; 
            } 
       }
       
       MPI_Barrier(g->mpicomm);
       
       MPI_Allreduce (&counter_nei, &counter_g, 1, MPI_INT,
                                MPI_SUM, g->mpicomm);
       
	if(counter_g == 0){
           	cout<<"no neighbours found for bc"<<endl;
		assert(false);
           }
       
       MPI_Allreduce(&qsum, &qsum_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       MPI_Allreduce(&vol, &vol_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       MPI_Allreduce(&pres, &pres_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       MPI_Allreduce(&ur, &ur_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       MPI_Allreduce(&ss, &ss_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       MPI_Allreduce(&rad,&rad_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);

       qsum_bc = qsum_g/counter_g;
       vol_bc = vol_g/counter_g;
       pres_bc = pres_g/counter_g;
       u_bc = ur_g/counter_g;
       ss_bc = ss_g/counter_g;
       r_bc = rad_g/counter_g;
}






void PelletSolver::updateStatesByLorentzForce(double dt, const double invelocityv, const double invelocityw, const double inpressure, const double involume, const double insoundspeed, 
        double* outvelocityv, double* outvelocityw, double* outpressure){
        
    if(heatingmodel == 0 || fabs(magneticfield - 0.)<1e-8)
        return;
    double LFy, d_vy, LFz, d_vz;
    double delta_E;
    double cond;
    double density, sc, pressure;
    double gamma;
    double vy = invelocityv;
    double vz = invelocityw;
    int numtimesteps = 1;
    density = 1./involume;
    pressure = inpressure;
    sc = insoundspeed;
    gamma = (sc*sc*density/pressure); 
    double vzold = vz;
    double vyold = vy;
    for(int i=0; i<numtimesteps; i++){
        cond = gdata->eos->getElectricConductivity(pressure,density);   
       /* 
        LFy = -cond*magneticfield*magneticfield/c_light/c_light*vy;
        d_vy =  LFy/density*dt/numtimesteps;
        LFz = -cond*magneticfield*magneticfield/c_light/c_light*vz;
        d_vz =  LFz/density*dt/numtimesteps;
        vy += d_vy; 
        vz += d_vz;
       if(fabs(d_vy) > fabs(vy))
        {
            d_vy = -vy;
            vy = 0;
        }
       if(fabs(d_vz) > fabs(vz))
       {
            d_vz = -vz;
            vz = 0;
       }
        (*outvelocityv) = (*outvelocityv) + d_vy;
        (*outvelocityw) = (*outvelocityw) + d_vz;
        */
       
        vy = (*outvelocityv);
        vz = (*outvelocityw);
        LFy = -cond*magneticfield*magneticfield/c_light/c_light/density*dt/numtimesteps;
        LFz = -cond*magneticfield*magneticfield/c_light/c_light/density*dt/numtimesteps;
        (*outvelocityv) = (*outvelocityv)*exp(LFy);
        (*outvelocityw) = (*outvelocityw)*exp(LFz);
        delta_E = 0.5*density*((vy*vy+vz*vz)-((*outvelocityv)*(*outvelocityv)+(*outvelocityw)*(*outvelocityw)));
        *outpressure += (gamma-1)*delta_E;
    

    }
   // delta_E = 0.5*density*((vyold*vyold+vzold*vzold)-(vy*vy+vz*vz));
   // *outpressure += (gamma-1)*delta_E;
}



bool PelletSolver:: ifCloseToPellet(double dx,double x, double y, double z){

    double xcen = 0;
    double ycen = 0;
    double zcen = 0;
    double pr = pelletradius;
    if(((x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen)) < (pr+dx)*(pr+dx))
        return true;
    else
        return false;
    
    }

void PelletSolver::setPelletMaterial(int i){
        materialnumber = i;
        if(i == 0){
            material = new Neon();
            }
        else if(i == 1){
            material = new Deuterium2();
            }
        else if(i == 2){
            material = new Deuterium();
            }
        mu = material->getMu();
        mass = material->getMass();
        Z = material->getZ();
        I = material->getI();
        sublimationenergy = material->getSublimationEnergy();
        one_plus_Zstar = material->getOne_Plus_Zstar(teinf);

    }
/*
void PelletSolver::updateStatesByLorentzForce(double dt){
    if(heatingmodel == 0 || fabs(magneticfield - 0.)<1e-8)
        return;
    size_t li, lpnum = gdata->particle_data->elem_count;
    pdata_t *pad;
    double vy ,vz, vyold, vzold;
    double density, sc, pressure;
    double cond;
    double LFy, LFz;
    double delta_E;
    double rad_cool, T;
    int numtimestep = 1;
    for(li = 0; li<lpnum; li++){
        
        pad = (pdata_t *) sc_array_index(gdata->particle_data, li);
        
        if(pad->ifboundary) continue;
        for(int i=0;i<numtimestep;i++){
            vy = pad->v[1];
            vz = pad->v[2];
            vyold = pad->oldv[1];
            vzold = pad->oldv[2];
            density = 1./pad->volume;
            sc = pad->soundspeed;
            pressure = pad->pressure;
            cond = gdata->eos->getElectricConductivity(pressure,density);
            T = gdata->eos->getTemperature(pressure,density);
            //rad_cool = neon_radiation_power_density(density, T);
            LFy = -cond*magneticfield*magneticfield/c_light/c_light/density*dt/numtimestep;
            LFz = -cond*magneticfield*magneticfield/c_light/c_light/density*dt/numtimestep;
            pad->v[1] = vy*exp(LFy);
            pad->v[2] = vz*exp(LFz);
            delta_E = 0.5*density*((vyold*vyold+vzold*vzold) - (pad->v[1]*pad->v[1]+pad->v[2]*pad->v[2]));
            
        //    pad->pressure = pad->pressure - (sc*sc*density/pressure -1 )*rad_cool*dt/numtimestep;
            if(delta_E > 0)
                pad->pressure = pad->pressure + (sc*sc*density/pressure-1)*delta_E;
            
            pad->soundspeed = gdata->eos->getSoundSpeed(pad->pressure, density);
        }
    }

}
*/


void PelletSolver::gradBDrift( double dt) {
    
    if(heatingmodel == 0 || fabs(magneticfield - 0.)<1e-8)
        return;
    size_t li, lpnum = gdata->particle_data->elem_count;
    pdata_t *pad;
    for(li=0; li<lpnum; li++){
        pad = (pdata_t*) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
//        pad->drifta = driftacceleration;
        pad->driftv += pad->drifta * dt;
        pad->v[1] += pad->drifta * dt; 
//        pad->driftv += driftacceleration * dt;
//        pad->v[1] += driftacceleration * dt; 
    }

}

void PelletSolver::computeDriftAcceleration(){
   
   size_t lnend;
   int counter_g;
//   size_t id,idn = node_offset->elem_count;
   double drifta = 0;
   double drifta_g;
   double R = 160.;
  /*
   for(id = 0; id<idn; id++){
       lnend = *(int *)sc_array_index(node_offset,id); 
       node = (computation_node_t*)sc_array_index_int(nodelist,(lnend)-1);
       drifta += 2*node->pressure/node->density/R;    
   }
    int counter = idn;
    MPI_Allreduce(&counter, &counter_g, 1, MPI_INT,MPI_SUM,gdata->mpicomm);
    MPI_Allreduce(&drifta,&drifta_g,1,MPI_DOUBLE,MPI_SUM,gdata->mpicomm);

    drifta_g /= counter_g;
//    cout<<drifta_g<<endl;
    driftacceleration = drifta_g;
*/
    size_t li, lpnum = gdata->particle_data->elem_count;
    drifta = 0;
    int counter = 0;
    pdata_t* pad;
    for(li = 0;li<lpnum;li++){
        pad = (pdata_t*)sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        drifta += pad->drifta;
        counter ++;
    }
    MPI_Barrier(gdata->mpicomm);
    MPI_Allreduce(&counter, &counter_g, 1, MPI_INT,MPI_SUM,gdata->mpicomm);
    MPI_Allreduce(&drifta,&drifta_g,1,MPI_DOUBLE,MPI_SUM,gdata->mpicomm);

    drifta_g /= counter_g;
    driftacceleration = drifta_g;
}

void PelletSolver:: computeNodeStates(){

    p8est_search_local(gdata->p8est,0,psearch_quad,compute_node_state ,nodelist);

}

void PelletSolver::regroupComputationNodeWithStates(sc_array_t *nodes){
    
    size_t li, lpnum = nodes->elem_count;
    int counter = 0;
    int counter_g;
    int offset;
    int lnnum_count[gdata->mpisize];
    int lnnum_offset[gdata->mpisize];
    sc_array_t *nodelist_temp;
    computation_node_t *node, *node2;
    nodelist_temp = sc_array_new(sizeof(computation_node_t));
    for(li=0;li<lpnum;li++){
        node = (computation_node_t *)sc_array_index(nodes,li);
        if(node->sphdensity >= 0){
            node2 = (computation_node_t *)sc_array_push(nodelist_temp);
            memcpy(node2,node,sizeof(computation_node_t));
            counter ++;
        }
    }
    MPI_Barrier(gdata->mpicomm);
    MPI_Allreduce(&counter,&counter_g,1,MPI_INT,MPI_SUM,gdata->mpicomm);
    assert((size_t)counter_g == nodes->elem_count);
    counter *= (int)sizeof(computation_node_t);  
    MPI_Allgather(&counter,1,MPI_INT,lnnum_count,1,MPI_INT,gdata->mpicomm);

    offset = 0;
    for(int i=0; i<(int)gdata->mpisize; i++){
        lnnum_offset[i] = offset;
        
        offset += lnnum_count[i];
        
    }
    offset /= (int)sizeof(computation_node_t);
    assert(offset == counter_g);
    nodelist_buf = sc_array_new_count(sizeof(computation_node_t), (size_t) offset);
    MPI_Allgatherv(nodelist_temp->array,lnnum_count[gdata->mpirank],MPI_BYTE,nodelist_buf->array, lnnum_count,
            lnnum_offset,MPI_BYTE,gdata->mpicomm); 
    sc_array_destroy_null(&nodelist_temp);
    sc_array_destroy(nodelist);
    nodelist = nodelist_buf;
}

void PelletSolver::getNodeLocalSpacing(sc_array_t *nodes, double* ls){
     size_t li, lnnum = nodes->elem_count;    
     computation_node_t *node;
     for(li = 0; li<lnnum; li++){
        node = (computation_node_t *)sc_array_index(nodes,li);
        ls[li] = node->localspacing;
     }
    
     MPI_Barrier(gdata->mpicomm);
      
    
     MPI_Allreduce(MPI_IN_PLACE,ls,(int)lnnum,MPI_DOUBLE,MPI_MAX,gdata->mpicomm);
     for(li = 0; li<lnnum; li++){
        node =(computation_node_t *) sc_array_index(nodes,li);
        node->localspacing = ls[li];
     }

}

void PelletSolver::labelLayer(){
    size_t li, lpnum = gdata->particle_data->elem_count;
    pdata_t *pad;
    double surface = 1.5;
    double y;
    for(li = 0;li<lpnum;li++){
        pad = (pdata_t*)sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary) continue;
        if(pad->layer == 1){
            pad->layer = -1;
            continue;
            }
        y = pad->xyz[1];
        if(y>surface && pad->layer == 0){
            pad->layer = 1;
            }
        }
     
    }
