# make with make -f test_CRUNCH_bins.make

CC=g++
CFLAGS=-c -O3
OFLAGS = -O3
LDFLAGS=l
SOURCES=test_CRUNCH_bins.cpp \
	../CRUNCH_engine.cpp \
	../VolumeParticleInfo.cpp \
	../CRN_parameters.cpp \
	../mathutil.cpp \
	../tParticle.cpp \
	../FT_util.cpp \
	../chronos_particle_info.cpp \
	../flowtube.cpp \
	../CRN_tParticle_bins.cpp \
	../CRUNCH_bins.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test_CRUNCH_bins.exe

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
