# make with make -f test_CRUNCH_engine.make

CC=g++
CFLAGS=-c -O3 -pg
OFLAGS = -O3 -pg
LDFLAGS=l
SOURCES=CRUNCH_eng_tester.cpp \
	../CRUNCH_engine.cpp \
	../VolumeParticleInfo.cpp \
	../CRN_parameters.cpp \
	../mathutil.cpp \
	../tParticle.cpp \
	../FT_util.cpp \
	../chronos_particle_info.cpp \
	../flowtube.cpp \
	../CRN_tParticle_bins.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=CRUNCH_eng_tester.exe

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
