# make with make -f CRUNCH_mixing_column.make

CC=g++
CFLAGS=-c -O3
OFLAGS = -O3
LDFLAGS=l
SOURCES=CRUNCH_mixing_column.cpp \
	../CRUNCH_engine.cpp \
	../VolumeParticleInfo.cpp \
	../CRN_parameters.cpp \
	../tParticle.cpp \
	../FT_util.cpp \
	../chronos_particle_info.cpp \
	../flowtube.cpp \
	../CRN_tParticle_bins.cpp \
	../CRUNCH_bins.cpp \
    ../LSDStatsTools.cpp \
    ../TNT/tnt.h
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=CRUNCH_mixing_column_new.exe

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
