#ifndef MANOUTFORMATS_H
#define MANOUTFORMATS_H

/*
 * outformats.h
 * definition of used symbols and output formats
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @date Thu Oct  2 09:59:57 CEST 2008
 */

/* definition of used symbols */
#define DESC "%1$s"       /* description of query */
#define QRY_LENGTH "%2$d" /* length of query */
#define SCR "%3$d"        /* score of match */
#define EVALUE "%4$s"     /* evalue of match */
#define QRY_S "%5$d"      /* start position on query */
#define QRY_E "%6$d"      /* end position on query */
#define SEQ_S "%7$d"      /* start position on sequence (absolute) */
#define SEQ_E "%8$d"      /* end position on sequence (absolute) */
#define MAT "%9$d"        /* number of matching symbols */
#define MIS "%10$d"       /* number of mismatching symbols */
#define INS "%11$d"       /* number of inserts */
#define DEL "%12$d"       /* number of deletions */
#define STRAND "%13$c"    /* strand of match */
#define EDIST "%14$d"     /* alignment edist */
#define SEQ "%15$s"       /* sequence of match */
#define SEQ_IDX "%16$s"   /* sequence index in multi fasta file */
#define QUERY "%17$s"     /* the query sequence*/
#define MEOP "%18$s"     /* the meop sequence*/


/* definition of header and format strings */
#define HEAD1 "#colheader\n#descr\tscore\tEvalue\tqstart\tqend\tmatches\tmismatches\tinsertions\tdeletions\tstrand\tsstart\tsend\tsequence\n"
#define FORMAT1 DESC"\t"SCR"\t"EVALUE"\t"QRY_S"\t"QRY_E"\t"MAT"\t"MIS"\t"INS"\t"DEL"\t"STRAND"\t"SEQ_S"\t"SEQ_E"\t"SEQ"\t"SEQ_IDX"\n"
#define HEAD2 "#colheader\n#descr\tscore\tqstart\tqend\tmatches\tmismatches\tinsertions\tdeletions\tstrand\tsstart\tsend\tsequence\n"
#define FORMAT2 DESC"\t"SCR"\t"QRY_S"\t"QRY_E"\t"MAT"\t"MIS"\t"INS"\t"DEL"\t"STRAND"\t"SEQ_S"\t"SEQ_E"\t"SEQ"\n"
#define HEAD3 "#colheader\n#gff-format"
#define FORMAT3 ".\tSEGEMEHL\tmatch\t"SEQ_S"\t"SEQ_E"\t"SCR"\t"STRAND"\t.\t qstart "QRY_S" ; qend "QRY_E" ; mat "MAT" ; mis "MIS" ; ins "INS" ; del "DEL" ; id "DESC" ; qlen "QRY_LENGTH" ; sequence "SEQ"\n"
#define HEAD4 "#colheader\n#descr\tscore\tEvalue\tqstart\tqend\tmatches\tmismatches\tinsertions\tdeletions\tstrand\tsstart\tsend\tsequenceidx\n"
#define FORMAT4 DESC"\t"SCR"\t"EVALUE"\t"QRY_S"\t"QRY_E"\t"MAT"\t"MIS"\t"INS"\t"DEL"\t"STRAND"\t"SEQ_S"\t"SEQ_E"\t"SEQ_IDX"\n"
#define HEAD5 "#\"colheader\"\n#\"descr\"\t\"semi global alignment edist\"\t\"seed score\"\t\"seed Evalue\"\t\"seed qstart\"\t\"seed qend\"\t\"semi global alignment matches\"\t\"semi global alignment mismatches\"\t\"semi global alignment insertions\"\t\"semi global alignment deletions\"\t\"strand\"\t\"start of semi global alignment in subject sequence\"\t\"end of semi global alignment in subject sequence\"\t\"sequenceidx\"\n"
#define FORMAT5 DESC"\t"EDIST"\t"SCR"\t"EVALUE"\t"QRY_S"\t"QRY_E"\t"MAT"\t"MIS"\t"INS"\t"DEL"\t"STRAND"\t"SEQ_S"\t"SEQ_E"\t"SEQ_IDX"\t"MEOP"\n"
#define HEAD6 "#colheader\n#descr\tsstart\tsend\tstrand\tedist\tsequenceidx\n"
#define FORMAT6 DESC"\t"SEQ_S"\t"SEQ_E"\t"STRAND"\t"EDIST"\t"SEQ_IDX"\n"
#define HEAD7 "#colheader\n#descr\tfull alignment edist\tseed score\tseed Evalue\tseed qstart\tseed qend\tmatches\tmismatches\tinsertions\tdeletions\tstrand\tsstart\tsend\tsubject\tquery\tsequenceidx\n"
#define FORMAT7 DESC"\t"EDIST"\t"SCR"\t"EVALUE"\t"QRY_S"\t"QRY_E"\t"MAT"\t"MIS"\t"INS"\t"DEL"\t"STRAND"\t"SEQ_S"\t"SEQ_E"\t"SEQ"\t"QUERY"\t"SEQ_IDX"\n"

/* definition of constant arrays */
const char* HEAD[] = {HEAD1, HEAD2, HEAD3, HEAD4, HEAD5, HEAD6, HEAD7};
const char* FORMAT[] = {FORMAT1, FORMAT2, FORMAT3, FORMAT4, FORMAT5, FORMAT6, FORMAT7};


#endif /* MANOUTFORMATS_H */
