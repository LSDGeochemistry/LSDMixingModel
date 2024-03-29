# make with make -f mixing_column.make

CC=g++
CFLAGS=-c -O3
OFLAGS = -O3
LDFLAGS=l
SOURCES=mixing_column_old.cpp \
	../VolumeParticleInfo.cpp \
	../CRN_parameters.cpp \
	../LSDParticle.cpp \
	../FT_util.cpp \
	../chronos_particle_info.cpp \
	../flowtube.cpp \
	../CRN_LSDParticle_bins.cpp \
    ../LSDStatsTools.cpp \
    ../TNT/tnt.h
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mixing_column_old.out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
