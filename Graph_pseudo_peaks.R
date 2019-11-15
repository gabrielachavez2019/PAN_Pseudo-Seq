#!/usr/bin/Rscript

## Project: Pseudo-Seq
## Purpose: Plot the alredy calculated pseudo-peaks
## Created By: Gabriela Toomer
## Created Date: November 11, 2019
#Read the count data

#var_name=[1]

T0 <- read.table('/home/chris/Pseudo_Bio02/T0.txt', header=FALSE)
#T1 <- read.table('/home/chris/Pseudo_Bio02/T1.txt', header=FALSE)
#T2 <- read.table('/home/chris/Pseudo_Bio02/T2.txt', header=FALSE)
#T3 <- read.table('/home/chris/Pseudo_Bio02/T3.txt', header=FALSE)
#T4 <- read.table('/home/chris/Pseudo_Bio02/T4.txt', header=FALSE)

#How to save it in a png
#png(filename="T1_pp.png", width=960, height=480)
#plot(T1$V2, T1$V1, type="h", xlab="PAN sequence", ylab="Pseudo Ratio", main="T1")
#dev.off()

#png(filename="T1_pp.png", width=960, height=480)
plot(T1$V2, T1$V1, type="h", ylim=c(0,5000), xlim=c(1,1075), xlab="PAN sequence", ylab="Pseudo Ratio", main="T1")
#axis(1, at=c(100,300,500,700,900),labels=c(100,300,500,700,900),las=1)
abline(v=255,lty=2,col="red")
#abline(v=428,lty=2,col="red")
#abline(v=431,lty=2,col="red")
#abline(v=536,lty=2,col="red")
abline(v=650,lty=2,col="red")
abline(v=740,lty=2,col="red")
abline(v=745,lty=2,col="red")
axis(1, at=c(248,660,745),labels=c("248","660","743/745"), las=2, cex=0.02)

#dev.off()
#png(filename="T2_pp.png", width=960, height=480)
plot(T2$V2, T2$V1, type="h", ylim=c(0,5000), xlim=c(1,1075), xlab="PAN sequence", ylab="Pseudo Ratio",main="T2")
#axis(1, at=c(100,300,500,700,900),labels=c(100,300,500,700,900),las=1)
#abline(v=255,lty=2,col="red")
#abline(v=428,lty=2,col="red")
#abline(v=431,lty=2,col="red")
#abline(v=536,lty=2,col="red")
abline(v=650,lty=2,col="red")
abline(v=740,lty=2,col="red")
abline(v=745,lty=2,col="red")
axis(1, at=c(660,745),labels=c("660","743/745"), las=2, cex=0.02)
#dev.off()

#png(filename="T3_pp.png", width=960, height=480)
plot(T3$V2, T3$V1, type="h", ylim=c(0,5000), xlim=c(1,1075), xlab="PAN sequence", ylab="Pseudo Ratio", main="T3")
#axis(1, at=c(100,300,500,700,900),labels=c(100,300,500,700,900),las=1)
#abline(v=255,lty=2,col="red")
#abline(v=428,lty=2,col="red")
#abline(v=431,lty=2,col="red")
#abline(v=536,lty=2,col="red")
abline(v=650,lty=2,col="red")
abline(v=740,lty=2,col="red")
abline(v=745,lty=2,col="red")
axis(1, at=c(660,745),labels=c("660","743/745"), las=2, cex=0.02)

#dev.off()
#png(filename="T4_pp.png", width=960, height=480)
plot(T4$V2, T4$V1, type="h", ylim=c(0,5000), xlim=c(1,1075), xlab="PAN sequence", ylab="Pseudo Ratio", main="T4")
#axis(1, at=c(100,300,500,700,900),labels=c(100,300,500,700,900),las=1)
abline(v=226,lty=2,col="red")
abline(v=650,lty=2,col="red")
abline(v=740,lty=2,col="red")
abline(v=745,lty=2,col="red")
abline(v=762,lty=2,col="red")
axis(1, at=c(660,745),labels=c("660","743/745"), las=2, cex=0.02)
#dev.off()
#png(filename="T0_pp.png", width=960, height=480)
plot(T0$V2, T0$V1, type="h", ylim=c(0,5000), xlim=c(1,1075), xlab="PAN sequence", ylab="Pseudo Ratio", main="T0")
axis(1, at=c(100,300,500,700,900),labels=c(100,300,500,700,900),las=1)
#dev.off()



