Estimate insert size

We will try to calculate the insert sizes of the two different libraries from the bam files, sometimes when you get data you do not know the insert sizes used, then this can be very helpful. Recall that the columns in the bam-file are: Read name, flag, reference it mapped to, position, mapping quality, cigar, mate reference map, mate position, template length, read sequence, read qualities and the different tags, among here the read group (RG) and the edit distances (NM).

Open the bam-file for library A, the insert size is the number just before the read sequences. You see some are negative other positive, this depends on the direction of the pairs and whether they map on the positive or negative strand.

samtools view HG00418_A.bam | less -S
Extract all insert_sizes - we use cut to only take column 9:

samtools view HG00418_A.bam | cut -f9 > HG00418_A.insertsizes.txt
samtools view HG00418_B.bam | cut -f9 > HG00418_B.insertsizes.txt
Lets calculate it in R (remember to quit R, write q() and press n)

R
a = read.table("HG00418_A.insertsizes.txt")
a.v = a[a[,1]>0,1]
# Lets get rid of outliers and use the 5% and 95% intervals to calculate mean and standard deviation:
mn = quantile(a.v, seq(0,1,0.05))[2]
mx = quantile(a.v, seq(0,1,0.05))[20]
mean(a.v[a.v >= mn & a.v <= mx])    # Mean
sd(a.v[a.v >= mn & a.v <= mx])      # SD

b = read.table("HG00418_B.insertsizes.txt")
b.v = b[b[,1]>0,1]
mn = quantile(b.v, seq(0,1,0.05))[2]
mx = quantile(b.v, seq(0,1,0.05))[20]
mean(b.v[b.v >= mn & b.v <= mx])    # Mean
sd(b.v[b.v >= mn & b.v <= mx])      # SD

