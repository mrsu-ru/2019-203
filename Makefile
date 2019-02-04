CC=g++
CFLAGS=-c -Wall -g
LDFLAGS= -g
SOURCES= \
    main.cpp \
    alekseev.cpp \
    borisovaem.cpp\
    viktorovichev.cpp\
    volkovda.cpp\
    dryginaea.cpp\
    epifanovats.cpp\
    eremkinnv.cpp\
    igaykiniv.cpp\
    kemaykinais.cpp\
    kiselevis.cpp\
    kobyzevkv.cpp\
    makarovaaa.cpp\
    malovava.cpp\
    maslovma.cpp\
    medvedevama.cpp \
    murtazinrr.cpp \
    nachinkinaoa.cpp \
    rusaevon.cpp \
    sirkinaa.cpp \
    sovetnikovakv.cpp \
    tenishevir.cpp \
    tereshkinki.cpp \
    tokarevaaa.cpp \
		shmelevaov.cpp \
		zhalninrv.cpp \
	lab.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=vvm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o vvm
