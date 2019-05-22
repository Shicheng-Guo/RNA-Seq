all: blockbuster.c expseq.c
	gcc -g -pedantic blockbuster.c -o blockbuster.x -lm -std=c99 -D_GNU_SOURCE
	gcc -g expseq.c libs/biofilesio.c -o expseq.x

blockbuster: blockbuster.c
	gcc -g -pedantic blockbuster.c -o blockbuster.x -lm -std=c99 -D_GNU_SOURCE

expseq: expseq.c
	gcc -g expseq.c libs/biofilesio.c -o expseq.x -lm -std=c99

clean: 
	rm blockbuster.x
	rm expseq.x
