#include <iostream>
#include "octree_manager.h"
using namespace std;
Octree_Manager:: Octree_Manager(Global_Data *g){

    gdata = g;

}


void Octree_Manager:: build_octree(){

    int mpiret;
    gdata->mpicomm = sc_MPI_COMM_WORLD;

    mpiret = sc_MPI_Comm_size (gdata->mpicomm, &gdata->mpisize);

    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Comm_rank (gdata->mpicomm, &gdata->mpirank);

    SC_CHECK_MPI (mpiret);
  
    
    sc_init (gdata->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  
    p4est_init (NULL, SC_LP_DEFAULT);

    if(gdata->dimension == 3){
//build octree connectivity
    gdata->conn = p8est_connectivity_new_unitcube ();
//build octree
    gdata->p8est = p8est_new_ext (gdata->mpicomm, gdata->conn, 0,
                            0 , 1,
                            sizeof (octant_data_t), NULL, gdata);
    }

    else if(gdata->dimension == 2){
    
//build octree connectivity
    gdata->conn2d = p4est_connectivity_new_unitsquare ();
//build octree
    gdata->p4est = p4est_new_ext (gdata->mpicomm, gdata->conn2d, 0,
                            0 , 1,
                            sizeof (octant_data_t), NULL, gdata);
    
    }
}



void Octree_Manager:: destroy_octree(){

    if(gdata->dimension == 3){
        p8est_destroy (gdata->p8est);
      
        gdata->p8est = NULL;
      
        p8est_connectivity_destroy (gdata->conn);
      
        gdata->conn = NULL;
    }
    else {
            
        p4est_destroy (gdata->p4est);
      
        gdata->p4est = NULL;
      
        p4est_connectivity_destroy (gdata->conn2d);
      
        gdata->conn2d = NULL;
            
    }
    sc_finalize ();

}


void Octree_Manager:: partition_octree(int allow_for_coarsening){
    
   if(gdata->dimension == 3) 
       p8est_partition(gdata->p8est,allow_for_coarsening,NULL); 

   else if(gdata->dimension == 2)

       p4est_partition(gdata->p4est,allow_for_coarsening,NULL); 
}

        
void Octree_Manager:: refine_octree(int recursive, p8est_refine_t refine_fn, p8est_init_t init_fn, p8est_replace_t replace_fn){
   
        p8est_refine_ext(gdata->p8est,recursive,-1,refine_fn,init_fn,replace_fn);
    
}

void Octree_Manager:: refine_octree2d(int recursive, p4est_refine_t refine_fn, p4est_init_t init_fn, p4est_replace_t replace_fn){
   
        p4est_refine_ext(gdata->p4est,recursive,-1,refine_fn,init_fn,replace_fn);
    
}

int Octree_Manager:: adapt_coarsen2d (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrants[])
{
    int i;
    p4est_locidx_t remain = 0, receive = 0;
    octant_data_t *oud;
    Global_Data *g = (Global_Data *)p4est->user_pointer;
    
    if (quadrants[1] == NULL ||
      quadrants[0]->level == g->minlevel) {
    
        oud = (octant_data_t *) quadrants[0]->p.user_data;
        g->ireindex2 += oud->premain;
        g->irvindex2 += oud->preceive;
    
        return 0;
  }
   
    remain = receive = 0;
  
    for (i = 0; i < P4EST_CHILDREN; ++i) {
        oud = (octant_data_t *) quadrants[i]->p.user_data;
        remain += oud->premain;
        receive += oud->preceive;
  }
  if ((double) (remain + receive) < .5 * g->elem_particle) {
    /* we will coarsen and adjust ireindex, irvindex in adapt_replace */
    g->qremain = remain;
    g->qreceive = receive;
    return 1;
  }
  else {
    /* we will not coarsen and proceed with next quadrant */
    oud = (octant_data_t *) quadrants[0]->p.user_data;
    g->ireindex2 += oud->premain;
    g->irvindex2 += oud->preceive;
    return 0;
  }


}

int Octree_Manager:: adapt_coarsen (p8est_t * p8est, p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrants[])
{
    int i;
    p4est_locidx_t remain = 0, receive = 0;
    octant_data_t *oud;
    Global_Data *g = (Global_Data *)p8est->user_pointer;
    
    if (quadrants[1] == NULL ||
      quadrants[0]->level == g->minlevel) {
    
        oud = (octant_data_t *) quadrants[0]->p.user_data;
        g->ireindex2 += oud->premain;
        g->irvindex2 += oud->preceive;
    
        return 0;
  }
   
    remain = receive = 0;
  
    for (i = 0; i < P8EST_CHILDREN; ++i) {
        oud = (octant_data_t *) quadrants[i]->p.user_data;
        remain += oud->premain;
        receive += oud->preceive;
  }
  if ((double) (remain + receive) < .5 * g->elem_particle) {
    /* we will coarsen and adjust ireindex, irvindex in adapt_replace */
    g->qremain = remain;
    g->qreceive = receive;
    return 1;
  }
  else {
    /* we will not coarsen and proceed with next quadrant */
    oud = (octant_data_t *) quadrants[0]->p.user_data;
    g->ireindex2 += oud->premain;
    g->irvindex2 += oud->preceive;
    return 0;
  }


}

void Octree_Manager:: balance_replace2d (p4est_t * p4est, p4est_topidx_t which_tree,
               int num_outgoing, p4est_quadrant_t * outgoing[],
               int num_incoming, p4est_quadrant_t * incoming[]){

  int octantid;
  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg, offset, *newoffset;
  p4est_quadrant_t  **pchild, *outoctant;
  octant_data_t          *oud, *iud;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  if (num_outgoing == P4EST_CHILDREN) {
    // we are coarsening 
      oud = (octant_data_t *) incoming[0]->p.user_data;
      oud->premain = oud->preceive = 0;
      for(int i=0;i<P4EST_CHILDREN;i++){
      iud = (octant_data_t *) outgoing[i]->p.user_data;
      oud->premain += iud->premain;
      oud->preceive += iud->preceive;
      }
      
      oud->poctant = oud->preceive + oud-> premain;
  }

  else {
    outoctant = outgoing[0];
    iud = (octant_data_t *)outoctant->p.user_data;
    octantid = iud->octantid;
    // we are refining 
    // access parent quadrant 
    g->loopquad2d (which_tree, outgoing[0],lxyz,hxyz,dxyz);

    // recover window onto remaining particles for the new family 
    
    irem = iud->premain ;//sc_array_index(g->irecumu,octantid);
    
    ibeg = *(p4est_locidx_t *) sc_array_index(g->irecumu,octantid) - irem;
    offset = ibeg;
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    
    g->klh[0] = &iview;
    wz = 0;

    
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (g->jlh[wy], g->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
        oud->poctant = oud->premain;
        oud->octantid = g->octantid;
        g->octantid ++;
        newoffset = (p4est_locidx_t *)sc_array_push(g->irecumu);
        *newoffset = offset + oud->premain;
        offset = offset + oud->premain;

      }
    }
    // recover window onto received particles for the new family 
    
    irem = iud->preceive ;//sc_array_index(g->irecumu,octantid);
    ibeg = *(p4est_locidx_t *) sc_array_index(g->irvcumu,octantid) - irem;
    offset = ibeg;
    P4EST_ASSERT (irem >= 0);
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);
    P4EST_ASSERT (qod->preceive == irem);

    // sort received particles into the children 
    pchild = incoming;
   
    g->klh[0] = &iview;
    wz = 0;
    
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord ( g->jlh[wy], g->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->preceive = (p4est_locidx_t) arr->elem_count);
        oud->poctant += oud->preceive;
        newoffset = (p4est_locidx_t *)sc_array_push(g->irvcumu);
        *newoffset = offset + oud->preceive;
        offset += oud->preceive;
      }
    }
    P4EST_ASSERT (pchild == incoming + P4EST_CHILDREN);

  }
}


void Octree_Manager:: balance_replace (p8est_t * p8est, p4est_topidx_t which_tree,
               int num_outgoing, p8est_quadrant_t * outgoing[],
               int num_incoming, p8est_quadrant_t * incoming[]){

  int octantid;
  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg, offset, *newoffset;
  p8est_quadrant_t  **pchild, *outoctant;
  octant_data_t          *oud, *iud;
  Global_Data      *g = (Global_Data *) p8est->user_pointer;
  if (num_outgoing == P8EST_CHILDREN) {
    // we are coarsening 
      oud = (octant_data_t *) incoming[0]->p.user_data;
      oud->premain = oud->preceive = 0;
      for(int i=0;i<P8EST_CHILDREN;i++){
      iud = (octant_data_t *) outgoing[i]->p.user_data;
      oud->premain += iud->premain;
      oud->preceive += iud->preceive;
      }
      
      oud->poctant = oud->preceive + oud-> premain;
  }

  else {
    outoctant = outgoing[0];
    iud = (octant_data_t *)outoctant->p.user_data;
    octantid = iud->octantid;
    // we are refining 
    // access parent quadrant 
    g->loopquad (which_tree, outgoing[0],lxyz,hxyz,dxyz);

    // recover window onto remaining particles for the new family 
    
    irem = iud->premain ;//sc_array_index(g->irecumu,octantid);
    
    ibeg = *(p4est_locidx_t *) sc_array_index(g->irecumu,octantid) - irem;
    offset = ibeg;
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    g->split_by_coord ( &iview, g->klh, PA_MODE_REMAIN, 2, lxyz, dxyz);
    
    for (wz = 0; wz < 2; ++wz) {
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (g->jlh[wy], g->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
        oud->poctant = oud->premain;
        oud->octantid = g->octantid;
        g->octantid ++;
        newoffset = (p4est_locidx_t *)sc_array_push(g->irecumu);
        *newoffset = offset + oud->premain;
        offset = offset + oud->premain;

      }
    }
    }
    // recover window onto received particles for the new family 
    
    irem = iud->preceive ;//sc_array_index(g->irecumu,octantid);
    ibeg = *(p4est_locidx_t *) sc_array_index(g->irvcumu,octantid) - irem;
    offset = ibeg;
    P4EST_ASSERT (irem >= 0);
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);
    P4EST_ASSERT (qod->preceive == irem);

    // sort received particles into the children 
    pchild = incoming;
    g->split_by_coord ( &iview, g->klh, PA_MODE_RECEIVE, 2, lxyz, dxyz);
    for (wz = 0; wz < 2; ++wz) {
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord ( g->jlh[wy], g->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->preceive = (p4est_locidx_t) arr->elem_count);
        oud->poctant += oud->preceive;
        newoffset = (p4est_locidx_t *)sc_array_push(g->irvcumu);
        *newoffset = offset + oud->preceive;
        offset += oud->preceive;
      }
    }
    }
    P4EST_ASSERT (pchild == incoming + P4EST_CHILDREN);

  }
}

void Octree_Manager:: adapt_replace2d (p4est_t * p4est, p4est_topidx_t which_tree,
               int num_outgoing, p4est_quadrant_t * outgoing[],
               int num_incoming, p4est_quadrant_t * incoming[]){

  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg;
  p4est_quadrant_t  **pchild;
  octant_data_t          *oud;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  if (num_outgoing == P4EST_CHILDREN) {
    // we are coarsening 
      oud = (octant_data_t *) incoming[0]->p.user_data;
    g->ireindex2 += (oud->premain = g->qremain);
    g->irvindex2 += (oud->preceive = g->qreceive);
    oud->poctant = oud->preceive + oud-> premain;
  }

  else {
    // we are refining 
    // access parent quadrant 
    g->loopquad2d (which_tree, outgoing[0],lxyz,hxyz,dxyz);

    // recover window onto remaining particles for the new family 
    ibeg = g->ire2;
    irem = g->ireindex - ibeg;
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    g->klh[0] = &iview;
    wz = 0;

    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (g->jlh[wy], g->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
        oud->poctant = oud->premain;
      
      }
    }

    // recover window onto received particles for the new family 
    ibeg = g->irv2;
    irem = g->irvindex - ibeg;
    P4EST_ASSERT (irem >= 0);
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);
    P4EST_ASSERT (qod->preceive == irem);

    // sort received particles into the children 
    pchild = incoming;
    wz = 0;
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord ( g->jlh[wy], g->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        P4EST_ASSERT (oud->u.lpend == -1);
        ibeg += (oud->preceive = (p4est_locidx_t) arr->elem_count);
        oud->poctant += oud->preceive;
      }
    }
    g->klh[0] = NULL;
    assert (ibeg == g->irvindex);
    P4EST_ASSERT (pchild == incoming + P4EST_CHILDREN);

  }

}

void Octree_Manager:: adapt_replace (p8est_t * p8est, p4est_topidx_t which_tree,
               int num_outgoing, p8est_quadrant_t * outgoing[],
               int num_incoming, p8est_quadrant_t * incoming[]){

  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg;
  p8est_quadrant_t  **pchild;
  octant_data_t          *oud;
  Global_Data      *g = (Global_Data *) p8est->user_pointer;
  if (num_outgoing == P8EST_CHILDREN) {
    // we are coarsening 
      oud = (octant_data_t *) incoming[0]->p.user_data;
    g->ireindex2 += (oud->premain = g->qremain);
    g->irvindex2 += (oud->preceive = g->qreceive);
    oud->poctant = oud->preceive + oud-> premain;
  }

  else {
    // we are refining 
    // access parent quadrant 
    g->loopquad (which_tree, outgoing[0],lxyz,hxyz,dxyz);

    // recover window onto remaining particles for the new family 
    ibeg = g->ire2;
    irem = g->ireindex - ibeg;
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    g->split_by_coord ( &iview, g->klh, PA_MODE_REMAIN, 2, lxyz, dxyz);
    
    for (wz = 0; wz < 2; ++wz) {
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (g->jlh[wy], g->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
        oud->poctant = oud->premain;
      
      }
    }
    }

    // recover window onto received particles for the new family 
    ibeg = g->irv2;
    irem = g->irvindex - ibeg;
    P4EST_ASSERT (irem >= 0);
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);
    P4EST_ASSERT (qod->preceive == irem);

    // sort received particles into the children 
    pchild = incoming;
    g->split_by_coord ( &iview, g->klh, PA_MODE_RECEIVE, 2, lxyz, dxyz);
    for (wz = 0; wz < 2; ++wz) {
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord ( g->jlh[wy], g->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        P4EST_ASSERT (oud->u.lpend == -1);
        ibeg += (oud->preceive = (p4est_locidx_t) arr->elem_count);
        oud->poctant += oud->preceive;
      }
    }
    }
    assert (ibeg == g->irvindex);
    P4EST_ASSERT (pchild == incoming + P4EST_CHILDREN);

  }

}

int Octree_Manager:: adapt_refine2d (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant)
{
  
  //static p4est_locidx_t maxp = -1;  
  Global_Data      *g = (Global_Data *) p4est->user_pointer;

  p4est_topidx_t      tt;

  p4est_quadrant_t   *quad;
  
  p4est_tree_t       *tree;
  tt = p4est->first_local_tree;
  tree = p4est_tree_array_index (p4est->trees, tt);
  
  quad = p4est_quadrant_array_index (&tree->quadrants, 0);
  octant_data_t          *oud = (octant_data_t *) quadrant->p.user_data;
  
/*  

*/
  
  
  /* we have set this to -1 in adapt_coarsen */

  if ((double) (oud->premain + oud->preceive) > g->elem_particle) {
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
    oud->octantid = g->octantid;
    g->octantid ++;
    
    return 0;
  }

}


int Octree_Manager:: adapt_refine (p8est_t * p8est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant)
{
   
 // static p4est_locidx_t maxp = -1;  
  Global_Data      *g = (Global_Data *) p8est->user_pointer;

  octant_data_t          *oud = (octant_data_t *) quadrant->p.user_data;
/*
  p4est_topidx_t      tt;

  p8est_quadrant_t   *quad;
  
  p8est_tree_t       *tree;
  tt = p8est->first_local_tree;
  tree = p8est_tree_array_index (p8est->trees, tt);
  
  quad = p8est_quadrant_array_index (&tree->quadrants, 0);
  

  if(g->flagstartrefine){
    maxp = -1;
    g->flagstartrefine = 0;
  
  }  
  
  if(quad == quadrant)
  { 
      if(maxp <= g->elem_particle && maxp>-1)
      { g->flagrefine = 0;
      }
      else
      { 
          maxp = -1;
      }
       
      int  mpiret = sc_MPI_Allreduce (&g->flagrefine, &g->gflagrefine, 1, sc_MPI_INT,
                               sc_MPI_MAX, g->mpicomm);
           
      SC_CHECK_MPI (mpiret);
  }
  
  if((oud->premain+oud->preceive) > maxp){
      maxp = oud->premain+oud->preceive;
  }


  */
  
  /* we have set this to -1 in adapt_coarsen */

  if ((double) (oud->premain + oud->preceive) > g->elem_particle) {
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
    oud->octantid = g->octantid;
    g->octantid ++;
    return 0;
  }

}

void Octree_Manager:: adapt_octree2d( p4est_t *p4est){

  
    int oldquad = (int)p4est->global_num_quadrants;

    while(true){
        gdata->ireindex2 = gdata->irvindex2 = 0;
        p4est_coarsen_ext (p4est, 0, 0, adapt_coarsen2d, NULL, adapt_replace2d);
    
        if((int)p4est->global_num_quadrants == oldquad)
            break;
        else
            oldquad = p4est->global_num_quadrants;
    }
    
    gdata->flagrefine = 1;
    gdata->gflagrefine = -1;
    
    gdata->flagstartrefine = 1;

    
    while(true ){   
    
    gdata->irecumu = sc_array_new(sizeof(p4est_locidx_t));
    gdata->irvcumu = sc_array_new(sizeof(p4est_locidx_t));
    
    gdata->octantid = 0;
    gdata->ireindex = gdata->ire2 = 0;
    gdata->irvindex = gdata->irv2 = 0;
    p4est_refine_ext (p4est, 0, gdata->maxlevel, adapt_refine2d, NULL, adapt_replace2d);
    
    if((int)p4est->global_num_quadrants == oldquad)
        break;
    else
        oldquad = p4est->global_num_quadrants;
   /* 
    if(gdata->gflagrefine == 0){
        break;
    }
    */
    sc_array_destroy(gdata->irecumu);
    sc_array_destroy(gdata->irvcumu);
    
   }

   
}

void Octree_Manager:: adapt_octree( p8est_t *p8est){

  
    int oldquad = (int)p8est->global_num_quadrants;
    
    while(true){
        gdata->ireindex2 = gdata->irvindex2 = 0;
        p8est_coarsen_ext (p8est, 0, 1, adapt_coarsen, NULL, adapt_replace);
   
        if((int)p8est->global_num_quadrants == oldquad)
            break;
        else
            oldquad = p8est->global_num_quadrants;
   }
    gdata->flagrefine = 1;
    gdata->gflagrefine = 1;
    
    gdata->flagstartrefine = 1;


    while(true ){   
    
    gdata->irecumu = sc_array_new(sizeof(p4est_locidx_t));
    gdata->irvcumu = sc_array_new(sizeof(p4est_locidx_t));
    
    gdata->octantid = 0;
    gdata->ireindex = gdata->ire2 = 0;
    gdata->irvindex = gdata->irv2 = 0;
    p8est_refine_ext (p8est, 0, gdata->maxlevel, adapt_refine, NULL, adapt_replace);
    
    if((int)p8est->global_num_quadrants == oldquad)
        break;
    else
        oldquad = p8est->global_num_quadrants;
    sc_array_destroy(gdata->irecumu);
    sc_array_destroy(gdata->irvcumu);
    
    }

   
}


void Octree_Manager::balance_octree2d(p4est_init_t init_fn, p4est_replace_t replace_fn){

    p4est_balance_ext(gdata->p4est,P4EST_CONNECT_FULL,init_fn,replace_fn);

    sc_array_destroy(gdata->irecumu);
    sc_array_destroy(gdata->irvcumu);
}

void Octree_Manager::balance_octree(p8est_init_t init_fn, p8est_replace_t replace_fn){

    p8est_balance_ext(gdata->p8est,P8EST_CONNECT_FULL,init_fn,replace_fn);

    sc_array_destroy(gdata->irecumu);
    sc_array_destroy(gdata->irvcumu);
}

static void testfaceside( p8est_iter_face_info_t * info, void *user_data){

    octant_data_t *ghost_data = (octant_data_t *)user_data;
    sc_array_t         *sides = &(info->sides);
    size_t sidescount = sides->elem_count;
    p8est_iter_face_side_t *sidedest, *sidesrc; 
    p8est_quadrant_t *qdest, *qsrc;
    octant_data_t *ouddest, *oudsrc;
    p4est_locidx_t quadid;
    p4est_locidx_t *neighbourid;
    //todo   domain boundary set to falgboundary
      
    for(size_t i = 0;i<sidescount;i++){
       sidedest = p8est_iter_fside_array_index_int(sides,i); 
       if(sidedest->is_hanging){
       
           for(int j=0;j<P8EST_HALF;j++){
                if(sidedest->is.hanging.is_ghost[j])          // is ghost
                    continue;
                qdest = sidedest->is.hanging.quad[j];
                ouddest = (octant_data_t *)qdest->p.user_data;
                if(ouddest->poctant == 0)          //no particles
                    continue;
                
                for(size_t k=0;k<sidescount;k++){
                  //  if(k == i)
                  //      continue;
                    sidesrc = p8est_iter_fside_array_index_int(sides,k); 
                    if(sidesrc->is_hanging){
                        
                        for(int l=0;l<P8EST_HALF;l++){
                            
                            
                            if(sidesrc->is.hanging.is_ghost[l]){
                                oudsrc = &ghost_data[sidesrc->is.hanging.quadid[l]];
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                                *neighbourid = sidesrc->is.hanging.quadid[l];
                            
                            }
                            else{

                                qsrc = sidesrc->is.hanging.quad[l];
                                oudsrc = (octant_data_t *)qsrc->p.user_data;
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                                *neighbourid = sidesrc->is.hanging.quadid[l];
                            
                            }

                        }
                    
                    }
                    else{                        //sidesrc is not hanging

                        if(sidesrc->is.full.is_ghost){
                        
                                oudsrc = &ghost_data[sidesrc->is.full.quadid];
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                                *neighbourid = sidesrc->is.full.quadid;
                        
                        }
                        else{
                        
                                qsrc = sidesrc->is.full.quad;
                                oudsrc = (octant_data_t *)qsrc->p.user_data;
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                                *neighbourid = sidesrc->is.full.quadid;
                        
                        }
                    
                    
                    }
                } 
                
           }
       
       }
       else{   //sidedest is full
       
       
                if(sidedest->is.full.is_ghost)          // is ghost
                    continue;
                qdest = sidedest->is.full.quad;
                ouddest = (octant_data_t *)qdest->p.user_data;
                if(ouddest->poctant == 0)          //no particles
                    continue;
                
                for(size_t k=0;k<sidescount;k++){
                   // if(k == i)
                    //    continue;
                    sidesrc = p8est_iter_fside_array_index_int(sides,k); 
                    if(sidesrc->is_hanging){
                        
                        for(int l=0;l<P8EST_HALF;l++){
                            
                            
                            if(sidesrc->is.hanging.is_ghost[l]){
                                oudsrc = &ghost_data[sidesrc->is.hanging.quadid[l]];
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                                *neighbourid = sidesrc->is.hanging.quadid[l];
                            
                            }
                            else{

                                qsrc = sidesrc->is.hanging.quad[l];
                                oudsrc = (octant_data_t *)qsrc->p.user_data;
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                                *neighbourid = sidesrc->is.hanging.quadid[l];
                           
                            }

                        }
                    
                    }
                    else{                        //sidesrc is not hanging

                        if(sidesrc->is.full.is_ghost){
                        
                                oudsrc = &ghost_data[sidesrc->is.full.quadid];
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                                *neighbourid = sidesrc->is.full.quadid;
                        
                        }
                        else{
                        
                                qsrc = sidesrc->is.full.quad;
                                oudsrc = (octant_data_t *)qsrc->p.user_data;
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                                *neighbourid = sidesrc->is.full.quadid;
                        
                        }
                    
                    
                    }
                } 
                
           
       
       
       }
    
    }

}

static void testedgeside( p8est_iter_edge_info_t * info, void *user_data){


    sc_array_t         *sides = &(info->sides);
    octant_data_t *ghost_data = (octant_data_t *)user_data;
    size_t sidescount = sides->elem_count;
    p8est_iter_edge_side_t *sidedest, *sidesrc; 
    p8est_quadrant_t *qdest, *qsrc;
    octant_data_t *ouddest, *oudsrc;
    p4est_locidx_t quadid;
    p4est_locidx_t *neighbourid;
    
    for(size_t i = 0;i<sidescount;i++){
        sidedest = p8est_iter_eside_array_index_int(sides,i); 
        
        if(sidedest->is_hanging){
  

           for(int j=0;j<2;j++){
                if(sidedest->is.hanging.is_ghost[j])          // is ghost
                    continue;
                qdest = sidedest->is.hanging.quad[j];
                ouddest = (octant_data_t *)qdest->p.user_data;
                if(ouddest->poctant == 0)          //no particles
                    continue;
                
                for(size_t k=0;k<sidescount;k++){
                  //  if(k == i)
                  //      continue;
                    sidesrc = p8est_iter_eside_array_index_int(sides,k); 
                    if(sidesrc->is_hanging){
                        
                        for(int l=0;l<2;l++){
                            
                            
                            if(sidesrc->is.hanging.is_ghost[l]){
                                oudsrc = &ghost_data[sidesrc->is.hanging.quadid[l]];
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                                *neighbourid = sidesrc->is.hanging.quadid[l];
                                  
                            }
                            else{

                                qsrc = sidesrc->is.hanging.quad[l];
                                oudsrc = (octant_data_t *)qsrc->p.user_data;
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                                *neighbourid = sidesrc->is.hanging.quadid[l];
                            
                            }

                        }
                    
                    }
                    else{                        //sidesrc is not hanging

                        if(sidesrc->is.full.is_ghost){
                        
                                oudsrc = &ghost_data[sidesrc->is.full.quadid];
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                                *neighbourid = sidesrc->is.full.quadid;
                        
                        }
                        else{
                        
                                qsrc = sidesrc->is.full.quad;
                                oudsrc = (octant_data_t *)qsrc->p.user_data;
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                                *neighbourid = sidesrc->is.full.quadid;
                        }
                    
                    
                    }
                } 
                
           }

        }
        
       else{   //sidedest is full
       
       
                if(sidedest->is.full.is_ghost)          // is ghost
                    continue;
                qdest = sidedest->is.full.quad;
                ouddest = (octant_data_t *)qdest->p.user_data;
                if(ouddest->poctant == 0)          //no particles
                    continue;
                
                for(size_t k=0;k<sidescount;k++){
                   // if(k == i)
                    //    continue;
                    sidesrc = p8est_iter_eside_array_index_int(sides,k); 
                    if(sidesrc->is_hanging){
                        
                        for(int l=0;l<2;l++){
                            
                            
                            if(sidesrc->is.hanging.is_ghost[l]){
                                oudsrc = &ghost_data[sidesrc->is.hanging.quadid[l]];
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                                *neighbourid = sidesrc->is.hanging.quadid[l];
                            
                            }
                            else{

                                qsrc = sidesrc->is.hanging.quad[l];
                                oudsrc = (octant_data_t *)qsrc->p.user_data;
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                                *neighbourid = sidesrc->is.hanging.quadid[l];
                           
                            }

                        }
                    
                    }
                    else{                        //sidesrc is not hanging

                        if(sidesrc->is.full.is_ghost){
                        
                                oudsrc = &ghost_data[sidesrc->is.full.quadid];
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                                *neighbourid = sidesrc->is.full.quadid;
                        
                        }
                        else{
                        
                                qsrc = sidesrc->is.full.quad;
                                oudsrc = (octant_data_t *)qsrc->p.user_data;
                                if(oudsrc->poctant == 0){
                                    ouddest->flagboundary = 1;
                                    continue;
                                }
                                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                                *neighbourid = sidesrc->is.full.quadid;
                        
                        }
                    
                    
                    }
                } 
       
       }
    
    }
    
    
}


