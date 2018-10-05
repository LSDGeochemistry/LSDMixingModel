# make with make -f make_CRUNCH_test.make

CC=g++
CFLAGS=-c -O3
OFLAGS = -O3
LDFLAGS=l
SOURCES=CRUNCH_volume_particle_tester.cpp \
	../CRUNCH_engine.cpp \
	../VolumeParticleInfo.cpp \
	../CRN_parameters.cpp \
	../mathutil.cpp \
	../tParticle.cpp \
	../FT_util.cpp \
	../chronos_particle_info.cpp \
	../flowtube.cpp \
	../CRN_tParticle_bins.cpp \
	../CRUNCH_bins.hpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=CRUNCH_vp_test.exe

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
