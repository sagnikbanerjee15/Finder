CC=			gcc
CXX=		g++
#CFLAGS=		-g -Wall -pg #-O2 -m64 -pg
CFLAGS=		-g -Wall -O2 -m64 #-pg
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-DOCC_INTERVAL_DEFAULT=0x80 #-D_FILE_OFFSET_BITS=64
#OBJS=		utils.o bwt.o bwtio.o jigsawaln.o bwtgap.o is.o \
#			bntseq.o bwtmisc.o stdaln.o simple_dp.o \
#			bwaseqio.o bwase.o kstring.o cs2nt.o \
#			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o \
#			bwtsw2_chain.o bamlite.o #bwt_lite.o 

OBJS=		utils.o bwt.o bwtio.o jigsawaln.o bwtgap.o is.o \
			bntseq.o bwtmisc.o stdaln.o simple_dp.o \
			bwaseqio.o bwase.o kstring.o cs2nt.o \
			bamlite.o bwt_lite.o splicesitemap.o splicescore.o


PROG=		olego olegoindex #bwtsw #jigsawse
INCLUDES=	-I./cz	

LIBS=		-lm -lz -lpthread -Lbwt_gen -lbwtgen
LIBS2=		./cz/czlib.a 

SUBDIRS=	. bwt_gen cz

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lib-recur all-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" CXX="$(CXX)" DFLAGS="$(DFLAGS)"  \
				INCLUDES="$(INCLUDES)" $$target || exit 1; \
			cd $$wdir; \
		done;

lib:

olegoindex:lib-recur $(OBJS) jigsawindex.o
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) jigsawindex.o -o $@ $(LIBS) $(LIBS2)

olego:lib-recur $(OBJS) jigsaw.o
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) jigsaw.o -o $@ $(LIBS) $(LIBS2)

bwtsw:lib-recur $(OBJS)  bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwtsw2_chain.o
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwtsw2_chain.o -o $@ $(LIBS) $(LIBS2)

jigsawse:lib-recur $(OBJS) jigsawse.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) jigsawse.o -o $@ $(LIBS)


bwt.o:bwt.h
bwtio.o:bwt.h
jigsaw.o:bwt.h jigsawaln.h kseq.h bwase.h splicesitemap.h splicescore.h 
bwt1away.o:bwt.h jigsawaln.h
bwt2fmv.o:bwt.h
bntseq.o:bntseq.h
bwtgap.o:bwtgap.h jigsawaln.h bwt.h

bwtsw2_core.o:bwtsw2.h bwt.h bwt_lite.h stdaln.h
bwtsw2_aux.o:bwtsw2.h bwt.h bwt_lite.h stdaln.h
bwtsw2_main.o:bwtsw2.h
splicescore.o:bntseq.h bwaseqio.h splicescore.h

cleanlocal:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

clean:cleanlocal-recur
