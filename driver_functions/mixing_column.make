# make with make -f mixing_column.make

CC=g++
CFLAGS=-c -O3 -pg
OFLAGS = -O3 -pg
LDFLAGS=l
SOURCES=CRUNCH_mixing_column.cpp \
	../VolumeParticleInfo.cpp \
	../CRN_parameters.cpp \
	../tParticle.cpp \
	../FT_util.cpp \
	../chronos_particle_info.cpp \
	../flowtube.cpp \
	../CRN_tParticle_bins.cpp \
    ../LSDStatsTools.cpp \
    ../TNT/tnt.h
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mixing_column.exe

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
