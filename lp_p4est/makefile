CC = mpic++ 
DEBUG = -g
LAPACK_DIR = /gpfs/home/shyyuan/local/lapack/
P4EST_INC = /home/syuan/p4est/include
P4EST_LIB = /home/syuan/p4est/lib
MAIN_DIR:=${CURDIR}
BOUNDARY_DIR=$(MAIN_DIR)/boundary/
STATE_DIR=$(MAIN_DIR)/state/
GEOMETRY_DIR=$(MAIN_DIR)/geometry/


INCS = -I $(P4EST_INC) -I $(BOUNDARY_DIR) -I $(STATE_DIR) -I $(GEOMETRY_DIR) -I $(MAIN_DIR)
LIBS = -L $(P4EST_LIB) -L $(LAPACK_DIR)
CFLAGS = -Wall -std=c++11 -c  $(DEBUG) $(INCS) $(LIBS) 
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS)
vpath %.h $(GEOMETRY_DIR) $(STATE_DIR) $(BOUNDARY_DIR)

MAIN_OBJS = lp_main.o particle_data.o initializer.o octree_manager.o registrar.o lp_solver.o eos.o particle_viewer.o ls_solver.o hexagonal_packing.o pellet_solver.o material_lib.o lp_timer.o
GEOMETRY_OBJS = geometry.o geometry_pellet.o geometry_disk.o geometry_cylinder.o
STATE_OBJS = state.o state_pellet.o state_gresho.o 
BOUNDARY_OBJS = boundary.o boundary_gresho.o boundary_pellet.o  

B_OBJS = $(foreach OBJ, $(BOUNDARY_OBJS),$(addprefix $(BOUNDARY_DIR),$(OBJ)))
S_OBJS = $(foreach OBJ,$(STATE_OBJS),$(addprefix $(STATE_DIR),$(OBJ)))
G_OBJS = $(foreach OBJ,$(GEOMETRY_OBJS),$(addprefix $(GEOMETRY_DIR),$(OBJ)))



all:  OBJS lp



lp: $(B_OBJS) $(MAIN_OBJS) $(S_OBJS) $(G_OBJS)
	$(CC) $(LFLAGS) $(MAIN_OBJS) $(B_OBJS) $(S_OBJS) $(G_OBJS)  -o lp -lsc -lp4est -llapacke -llapack -lgfortran -lrefblas -lgsl -lgslcblas

OBJS:
	cd $(BOUNDARY_DIR)&&make;
	cd $(STATE_DIR)&&make;
	cd $(GEOMETRY_DIR)&&make


lp_main.o: lp_main.cpp initializer.h particle_data.h
	$(CC) $(CFLAGS) lp_main.cpp
particle_data.o: particle_data.cpp particle_data.h geometry.h initializer.h
	$(CC) $(CFLAGS) particle_data.cpp
initializer.o: initializer.cpp initializer.h
	$(CC) $(CFLAGS) initializer.cpp
octree_manager.o: octree_manager.cpp octree_manager.h particle_data.h
	$(CC) $(CFLAGS) octree_manager.cpp
registrar.o: registrar.h registrar.cpp
	$(CC) $(CFLAGS) registrar.cpp
lp_solver.o: lp_solver.h lp_solver.cpp initializer.h particle_data.h
	$(CC) $(CFLAGS) lp_solver.cpp
eos.o: eos.h eos.cpp
	$(CC) $(CFLAGS) eos.cpp

particle_viewer.o: particle_data.h particle_viewer.cpp particle_viewer.h
	$(CC) $(CFLAGS) particle_viewer.cpp
ls_solver.o: ls_solver.h ls_solver.cpp
	$(CC) $(CFLAGS) ls_solver.cpp
hexagonal_packing.o: hexagonal_packing.h hexagonal_packing.cpp
	$(CC) $(CFLAGS) hexagonal_packing.cpp

pellet_solver.o: pellet_solver.h particle_data.h pellet_solver.cpp
	$(CC) $(CFLAGS) pellet_solver.cpp

material_lib.o: material_lib.cpp material_lib.h
	$(CC) $(CFLAGS) material_lib.cpp
lp_timer.o: lp_timer.cpp lp_timer.h
	$(CC) $(CFLAGS) lp_timer.cpp

clean:
	
	rm *.o
	cd $(GEOMETRY_DIR)&&make clean; 
	cd $(STATE_DIR)&&make clean;
	cd $(BOUNDARY_DIR)&&make clean;
	
	rm lp
