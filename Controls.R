

#Read files, the files have two columns one with the posution the second wit the read 
#Controls CLAP for pseudo synthetic RNA

U.CTL <- read.table("SM-16.tex", header = FALSE, sep = ' ')
U.CMC <- read.table("SM-15.tex", header = FALSE, sep = ' ')
Y.CTL <- read.table("SM-14.tex", header = FALSE, sep = ' ')
Y.CMC <- read.table("SM-13.tex", header = FALSE, sep = ' ')

#If you plot that makes no sense, becasue the second column can not be plottted
plot(Y.CMC, type= "l", col=2)
points(Y.CTL, col=3, type = "l")
points(U.CMC, col=4, type = "l")
points(U.CTL, col=5, type = "l")

#I extracted all the positions (column 1) and count the Frequency
Ypos <- as.data.frame(table(Y.CMC$V1))
Yneg <- as.data.frame(table(Y.CTL$V1))
Upos <- as.data.frame(table(U.CMC$V1))
Uneg <- as.data.frame(table(U.CTL$V1))

#Now I have a data frame with two columns one the position the second the Frequency:
#> head(Y.CMC)
#V1                                                                                                   V2
#1  1                                                      GGGAGAGCGAGAACACACCACAACGAAAACGAGCAAAACCCGGACGC
#2  1 GGGAGAGCGAGAACACACCACAACGAAAACGAGCAAAACCCGGACGCAACACAAAAGCGAACAACGCGAAAAAGGACACCGAAGCGGAAGCAAAGACAAC
#3  1                                                        GGGAGAGCGAGAACACACCACAACGAAAACGAGCAAAACCCGGAC
#4  1  GGGAGAGCGAGAACACACCACAACGAAAACGAGCAAAACCCGGACGCAACACAAAAGCGAACAACGCGAAAAAGGACACCGAAGCGGAAGCAAAGACAA
#5  1                                                GGGAGAGCGAGAACACACCACAACGAAAACGAGCAAAACCCGGACGCAACACA
#6  1  GGGAGAGCGAGAACACACCACAACGAAAACGAGCAAAACCCGGACGCAACACAAAAGCGAACAACGCGAAAAAGGACACCGAAGCGGAAGCAAAGACAA
#
#> head(Ypos)
#Var1 Freq
#1    1 1611
#2    2  245
#3    3   60
#4    4  115
#5    5  121
#6    6  549

#Then I can plot the Frequecy for all the data frames
plot(Uneg$Freq, type="l", col="purple", lwd=2,
     xlim=c(7,90),
     ylim=c(0,6000))
lines(Yneg$Freq, col=3, lwd=2)
lines(Upos$Freq, col=4, lwd=2)
lines(Ypos2$Freq, col=1, lwd=2)
#Adds the Pseudo position or where supposed to be
abline(v=43, col=2, cex=3, lwd=1, lty=3)
#Adds the legned o the figure
legend(8,6000, legend = c("U -CMC", "Pseudo -CMC", "U +CMC","Pseudo +CMC"), 
       col = c("purple",3,4,1),
       lty=1, cex = 1)


Ypos
Yneg
summary(Ypos)

#### Only ^T reads and positions
only16 <- read.table("16.only.T.tex", header = FALSE, sep = ' ')
only15 <- read.table("15.only.T.tex", header = FALSE, sep = ' ')
only14 <- read.table("14.only.T.tex", header = FALSE, sep = ' ')
only13 <- read.table("13.only.T.tex", header = FALSE, sep = ' ')

#To generate two columns one with the positio and the other with the frequency
only13.a <- as.data.frame(table(only13$V1))
only16.a <- as.data.frame(table(only16$V1))

#Visualize the Frequency by a dotplot
plot(only16.a$Freq)
points(only13.a$Freq, col=2)
abline(v=43, col="purple")

# Eliminate manually peaks in Y.CMC that are not U.s
write.table(Ypos, "Ypos.tex", sep = ' ')
Ypos2 <- read.table("Ypos.tex", header = TRUE, sep = ',')

#Normalizaton
#Getting info of Frequencies
summary(Uneg)
summary(Upos)
summary(Yneg)
summary(Ypos2)

boxplot(Uneg$Freq, Upos$Freq, Yneg$Freq, Ypos2$Freq)
Ypos2[43,]

#Make libraries even, they all should have same number of total reads
#Pseudo-Freq = log2(+CMC stops at postion/total of reads)/(-CMC stops at position/total of reads)
#Comparisons: 
            # 1  Pseudo+CMC vs Pseudo-CMC
            # 2  U +CMC     vs U -CMC
#How many reads "total read I have per library"?
sum(Ypos2$Freq)
sum(Yneg$Freq)
sum(Upos$Freq)
sum(Uneg$Freq)


U.CTL <- read.table("SM-16.tex", header = FALSE, sep = ' ')
U.CMC <- read.table("SM-15.tex", header = FALSE, sep = ' ')
Y.CTL <- read.table("SM-14.tex", header = FALSE, sep = ' ')
Y.CMC <- read.table("SM-13.tex", header = FALSE, sep = ' ')

#Select random sample to SAME TOTAL NUMBER OF READS this case = 13972

U.CTL.SameTotal <- U.CTL[sample(nrow(U.CTL), 13972), ]
U.CMC.SameTotal <- U.CMC[sample(nrow(U.CMC), 13972), ]
Y.CTL.SameTotal <- Y.CTL[sample(nrow(Y.CTL), 13972), ]
Y.CMC.SameTotal <- Y.CMC[sample(nrow(Y.CMC), 13972), ]

#Check if it works
tail (U.CTL.SameTotal)
tail(Y.CTL.SameTotal)

#I extracted all the positions (column 1) and count the Frequency
Ypos <- as.data.frame(table(Y.CMC.SameTotal$V1))
Yneg <- as.data.frame(table(Y.CTL.SameTotal$V1))
Upos <- as.data.frame(table(U.CMC.SameTotal$V1))
Uneg <- as.data.frame(table(U.CTL.SameTotal$V1))


#Then I can plot the Frequecy for all the data frames
plot(Uneg$Freq, type="l", col="purple", lwd=2,
     xlim=c(17,85),
     ylim=c(0,6000))
lines(Yneg$Freq, col=3, lwd=2)
lines(Upos$Freq, col=4, lwd=2)
lines(Ypos2$Freq, col=1, lwd=2)
#Adds the Pseudo position or where supposed to be
abline(v=43, col=2, cex=3, lwd=2, lty=3)
#Adds the legned o the figure
legend("topleft", inset = 0.05, legend = c("U -CMC", "Pseudo -CMC", "U +CMC","Pseudo +CMC"), 
       col = c("purple",3,4,1),
       lty=1, cex = 1)

sum(Ypos2$Freq)
sum(Yneg$Freq)
sum(Upos$Freq)
sum(Uneg$Freq)
