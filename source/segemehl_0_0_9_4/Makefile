  CC=gcc
  LD=${CC} 
  CFLAGS=-Wall -pedantic -std=c99 -g -O3 -DSHOWALIGN -DDIRECTOUT -DDBGLEVEL=0 -DPROGNFO -Isrc -Ilibs -Ilibs/sufarray  -Lsrc 
  LDFLAGS=-lm -lpthread 
  CTAGS=ctags > tags
  LIBS=-lob -lm -lpthread 

  SEGEMEHLOBJ = libs/md5.o\
  				libs/stringutils.o\
				libs/mathematics.o\
                libs/memman.o\
				libs/fileio.o\
				libs/biofiles.o\
				libs/debug.o\
				libs/info.o\
				libs/karlin.o\
				libs/kdseed.o\
				libs/queue.o\
				libs/sufarray/sufarray.o\
				libs/sufarray/charsequence.o\
				libs/sufarray/multicharseq.o\
				libs/stack.o\
				libs/aluruSort.o\
				libs/sort.o\
				libs/vtprogressbar.o\
				libs/sufarray/mmchar.o\
                libs/bitArray.o\
				libs/bitVector.o\
				libs/bitvectoralg.o\
                libs/manout.o\
				libs/kdchain.o\
				libs/manopt.o\
				libs/container.o\
				libs/vstack.o\
				libs/vqueue.o\
				libs/alignment.o\
				src/kdmatch.o\
				src/segemehl.o


all: segemehl.x  

segemehl.x: ${SEGEMEHLOBJ}
					gcc $(CFLAGS) ${SEGEMEHLOBJ} -o $@ $(LDFLAGS)

clean: 
	rm -f ${SEGEMEHLOBJ}

