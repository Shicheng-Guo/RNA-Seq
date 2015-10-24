  CC=gcc
  LD=${CC} 
  CFLAGS= -Wall -pedantic -std=c99 -O3 -g -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DDBGNFO -DSHOWALIGN -DDBGLEVEL=0 -DPROGNFO -Isrc -Ilibs -Ilibs/sufarray  -Lsrc 
  LDFLAGS=-lm -lpthread -lz -lncurses -lform -lmenu
  CTAGS=ctags > tags
  LIBS=-lob -lm -lpthread -pg 
 
  FILEBINOBJ = libs/debug.o\
			   libs/info.o\
  			   libs/physmem.o\
  			   libs/stringutils.o\
			   libs/fileio.o\
			   libs/memory.o\
			   libs/radixsort.o\
			   libs/fileBins.o\
			   src/filebintest.o\

  GENFASTAOBJ = libs/vtprogressbar.o\
                libs/debug.o\
				libs/info.o\
				libs/seqclip.o\
  				libs/stringutils.o\
				libs/mathematics.o\
				libs/fileio.o\
				libs/zran.o\
  				libs/biofiles.o\
				libs/randseqs.o\
  				libs/memory.o\
				libs/alignment.o\
				libs/sw.o\
				libs/sufarray/charsequence.o\
				libs/manopt.o\
				libs/iupac.o\
				libs/bitVector.o\
				src/genfasta.o

  SEGEMEHLOBJ =	libs/debug.o\
				libs/info.o\
				libs/md5.o\
  				libs/stringutils.o\
				libs/mathematics.o\
                libs/memory.o\
 				libs/zran.o\
				libs/fileio.o\
				libs/biofiles.o\
				libs/karlin.o\
				libs/kdseed.o\
				libs/queue.o\
				libs/sufarray/sufarray.o\
				libs/sufarray/charsequence.o\
				libs/sufarray/multicharseq.o\
			    libs/radixsort.o\
			    libs/fileBins.o\
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
 				libs/sw.o\
				libs/seqclip.o\
				libs/iupac.o\
				libs/hash.o\
				libs/matchfilesfields.o\
				libs/merge.o\
  				src/kdmatch.o\
				src/segemehl.o


  TESTMATOBJ =  libs/vtprogressbar.o\
  				libs/debug.o\
				libs/info.o\
				libs/bitVector.o\
				libs/seqclip.o\
  				libs/stringutils.o\
				libs/mathematics.o\
				libs/fileio.o\
				libs/zran.o\
  				libs/biofiles.o\
				libs/randseqs.o\
				libs/sufarray/charsequence.o\
  				libs/memory.o\
				libs/iupac.o\
				libs/alignment.o\
				libs/sw.o\
				libs/nw.o\
				libs/manopt.o\
				libs/container.o\
				libs/vstack.o\
				libs/vqueue.o\
			    libs/sort.o\
				libs/list.o\
				libs/matepairs.o\
				libs/matchfilesfields.o\
  				libs/matchfiles.o\
				libs/evalmatchfiles.o\
				libs/browsematchfiles.o\
				libs/plotmatchfiles.o\
				libs/ncursesext.o\
				libs/splicesites.o\
				libs/startsites.o\
				libs/matfile.o

   LISTTESTOBJ = libs/memory.o\
   				 libs/vstack.o\
				 libs/vqueue.o\
				 libs/sort.o\
				 libs/debug.o\
				 libs/list.o

    SAXTESTOBJ = libs/memory.o\
				 libs/mathematics.o\
  				 libs/stringutils.o\
				 libs/fileio.o\
				 libs/SAX.o

all: segemehl.x

maribelita.x: ${TESTMATOBJ}
					gcc $(CFLAGS) ${TESTMATOBJ} -o $@ $(LDFLAGS)

listtest.x: ${LISTTESTOBJ}
					gcc ${CFLAGS} ${LISTTESTOBJ} -o $@ $(LDFLAGS)

saxtest.x: ${SAXTESTOBJ}
					gcc ${CFLAGS} ${SAXTESTOBJ} -o $@ $(LDFLAGS)

filebintest.x: ${FILEBINOBJ}
					gcc $(CFLAGS) ${FILEBINOBJ} -o $@ $(LDFLAGS)

segemehl.x: ${SEGEMEHLOBJ}
					echo " " >> src/version.h
					gcc $(CFLAGS) ${SEGEMEHLOBJ} -o $@ $(LDFLAGS)

genfasta.x: ${GENFASTAOBJ}
					gcc $(CFLAGS) ${GENFASTAOBJ} -o $@ $(LDFLAGS)

testgff.x: ${TESTGFFOBJ}
					gcc $(CFLAGS) -DGFFTEST ${TESTGFFOBJ} -o $@ $(LDFLAGS)


clean: 
	rm -f ${TESTMATOBJ} ${TESTMATOBJ}
	rm -f ${SEGEMEHLOBJ} ${GENFASTAOBJ}

