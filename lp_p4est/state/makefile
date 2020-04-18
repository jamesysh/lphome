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


state.o: state.h state.cpp
	${CC} ${CFLAGS} state.cpp
state_pellet.o: state.h state_pellet.h state_pellet.cpp
	$(CC) $(CFLAGS) state_pellet.cpp
state_gresho.o: state.h state_gresho.h state_gresho.cpp
	$(CC) $(CFLAGS) state_gresho.cpp


clean:
	rm *.o


