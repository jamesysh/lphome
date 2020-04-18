CC     = g++
DEBUG  = -g
#OMP =
LAPACK_DIR=/gpfs/home/shyyuan/local/lapack-3.8.0
MAIN_DIR:=${CURDIR}/../
BOUNDARY_DIR=$(MAIN_DIR)/boundary/
STATE_DIR=$(MAIN_DIR)/state/
GEOMETRY_DIR=$(MAIN_DIR)/geometry/
INCS   = -I $(BOUNDARY_DIR) -I $(STATE_DIR) -I $(GEOMETRY_DIR) -I $(MAIN_DIR)
LIBS   = -L $(LAPACK_DIR) 
CFLAGS = -Wall -c -std=c++11 $(DEBUG) $(OMP) $(INCS)
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS) $(OMP)
SRC_FILES := $(wildcard *.cpp) 
OBJS := $(SRC_FILES:.cpp=.o)
all: $(OBJS)

geometry.o: geometry.h geometry.cpp
	${CC} ${CFLAGS} geometry.cpp
geometry_pellet.o: geometry.h geometry_pellet.h geometry_pellet.cpp
	${CC} ${CFLAGS} geometry_pellet.cpp
geometry_disk.o: geometry.h geometry_disk.h geometry_disk.cpp
	$(CC) ${CFLAGS} geometry_disk.cpp
geometry_cylinder.o: geometry.h geometry_cylinder.h geometry_cylinder.cpp
	$(CC) $(CFLAGS) geometry_cylinder.cpp


clean:
	rm *.o


