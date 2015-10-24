#########################################################
# R Script to plot read profile of a cluster
# Arguments: <profile file> <png file> <range> <queryid>
#########################################################

ARGV <- commandArgs(TRUE)
d <- read.table(ARGV[1], sep="\t")
plot <- paste(ARGV[2], sep="")
pdf(plot, width=3.5, height=4)

#Get the range for the x and y axis
xrange <- range(1:ARGV[3])
yrange <- range(0,1.2)

#Set up the plot
plot(xrange, yrange, type="n", xlab="alignment coordinate", ylab="normalized read count", cex.lab=1, cex.axis=1)

color <- c("darkred", "darkblue")
linetype <- c(1, 4)
pointer <- c(1, 3)

#Add lines
lines(d[,1], d[,2], type="l", lty=linetype[1], col=color[1], pch=pointer[1])
points(d[,1], d[,2], col=color[1], pch=pointer[1], cex=1)
grid(lwd=1)

#Add title
title(main="Normalized block profile", cex.main=1)

#Add a legend
legend("topright", xrange[1], c(ARGV[4]), cex=0.5, col=color[1], pch=pointer[1], lty=linetype[1], bty="n")

garbage <- dev.off()
q()
