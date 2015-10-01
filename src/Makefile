CC = gcc 
LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include
#COPT = -O3 -ansi -pedantic -Wall -I$(INCDIR) -L$(LIBDIR)
COPT = -g -ansi -pedantic -Wall -I$(INCDIR) -L$(LIBDIR)
LIBS = -lbiop -lgen -lm -lxml2
EXE  = buildloopdb scanloopdb finddist

BOBJS  = buildloopdb.o
SOBJS  = scanloopdb.o
FOBJS  = finddist.o

all : $(EXE)

buildloopdb : $(BOBJS)
	$(CC) $(COPT) -o $@ $(BOBJS) $(LIBS)

buildloopdb.o : buildloopdb.c distances.h
	$(CC) $(COPT) -c -o $@ $<

scanloopdb : $(SOBJS)
	$(CC) $(COPT) -o $@ $(SOBJS) $(LIBS)

finddist : $(FOBJS)
	$(CC) $(COPT) -o $@ $(FOBJS) $(LIBS)

.c.o : 
	$(CC) $(COPT) -c -o $@ $<

distances.h :
	wget http://www.bioinf.org.uk/abs/abdb/Data/NR_CombinedAb_Chothia.tar.bz2
	wget http://www.bioinf.org.uk/abs/abdb/Data/NR_CombinedHv_Chothia.tar.bz2
	tar xvf NR_CombinedAb_Chothia.tar.bz2
	tar xvf NR_CombinedHv_Chothia.tar.bz2
	mkdir -p abdb
	mv NR_CombinedAb_Chothia/* abdb
	mv NR_CombinedHv_Chothia/* abdb
	./makedistances.pl > distances.h

install : $(EXE)
	mkdir -p ../bin
	cp $(EXE) ../bin

clean :
	\rm -f $(BOBJS) $(SOBJS) $(FOBJS)
	\rm -rf NR_Combined??_Chothia*
	\rm -rf abdb

distclean : clean
	\rm -f $(EXE)



