## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  fig.align = "center"
)

## -----------------------------------------------------------------------------
library(PFExpTools)

#Get gene average across time course
NF54Avg <- as.data.frame(apply(NF54, 1, function(x){mean(x,na.rm=T)}))

#Get log2 transformed data
NF54log2 <- LogTransform(newdata = NF54, refavg = NF54Avg)

#We had Time as a row; it was log transformed. Let's fix that.

NF54log2["Time",] <- NF54["Time",]

#Example plot
plot(NF54log2["Time",which(!is.na(NF54log2["PF3D7_1343700",]))], 
     NF54log2["PF3D7_1343700",which(!is.na(NF54log2["PF3D7_1343700",]))],
     xlab = "hpi", ylab = "log2 Expression", pch = 16, col = "red")



## -----------------------------------------------------------------------------

#impute missing time points, using default number of neighbors
NF54noNAlog2 <- tcImpute(NF54log2)

plot(NF54noNAlog2["Time",], NF54noNAlog2["PF3D7_1343700",],
     xlab = "hpi", ylab = "log2 Expression", pch = 16, col = "red")


## -----------------------------------------------------------------------------

#We need a Time matrix for the function; pull from NF54log2
NF54Time <- data.matrix(NF54noNAlog2["Time",rownames.force = T])

#Smooth data, with default 
NF54smooth <- lSmooth(expdata = NF54noNAlog2, time = NF54Time)

plot(NF54smooth["Time",], NF54smooth["PF3D7_1343700",],
     xlab = "hpi", ylab = "log2 Expression", pch = 16, col = "red")


## -----------------------------------------------------------------------------

#impute time course
NF54imputed <- fftimpute(NF54smooth, NF54Time)

#Again, we operated the Time row. Let's replace it.
NF54imputed["Time",] <- 2:57

plot(NF54imputed["Time",], NF54imputed["PF3D7_1343700",],
     xlab = "hpi", ylab = "log2 Expression", pch = 16, col = "red")


## -----------------------------------------------------------------------------
NF54spline <- splineCourse(NF54log2, NF54Time, fulltimes = 1:48)

plot(NF54spline["Time",], NF54spline["PF3D7_1343700",],
     xlab = "hpi", ylab = "log2 Expression", pch = 16, col = "red")



## -----------------------------------------------------------------------------
#First, let us use the methods demonstrated previously to curate the data

#Get log2 transformed data, using the NF54Avg data allowing for comparisons
NHP4026log2 <- LogTransform(newdata = NHP4026, refavg = NF54Avg)

#We had Time as a row; it was log transformed. Let's fix that.
NHP4026log2["Time",] <- NHP4026["Time",]

#kNN imputation of missing values
NHP4026noNAlog2 <- tcImpute(NHP4026log2)

#We need a Time matrix for the function; pull from NF54log2
NHP4026Time <- data.matrix(NHP4026noNAlog2["Time",rownames.force = T])

#Smooth data, with default 
NHP4026smooth <- lSmooth(expdata = NHP4026noNAlog2, time = NHP4026Time)

#impute to hourly time course
NHP4026imputed <- fftimpute(NHP4026smooth, NHP4026Time)

#Again, we operated the Time row. Let's replace it.
NHP4026imputed["Time",] <- 2:57

plot(NHP4026imputed["Time",], NHP4026imputed["PF3D7_1343700",],
     xlab = "hpi", ylab = "log2 Expression", pch = 16, col = "blue")

## ----fig.height=6, fig.width=6------------------------------------------------
#Stage all NHP4026 time points; generates a heat map by default
NHP4026staging <- stagingByTranscription(NHP4026imputed[!rownames(NHP4026imputed) %in% "Time",],
                                         NF54imputed[!rownames(NF54imputed) %in% "Time",])

#plot NHP4026 stages
plot(seq(1,60), seq(1,60), type = "n", main = "NHP4026 Time Course Correlation to NF54",
     xlab = "Sampling Time (hpi)", ylab = "Peak Correlation to NF54 (hpi)")
abline(a = 0, b = 1, col = "red", lwd = 2)
points(2:57, NHP4026staging[,1], pch = 16, cex = 1)
grid()


## -----------------------------------------------------------------------------


NF54alltimes <- data.matrix(NF54imputed["Time",])
NF54peak <- peakPhase(NF54imputed[!rownames(NF54imputed) %in% "Time",],
                         NF54alltimes, method = 'FFT')

NHP4026alltimes <- data.matrix(NHP4026imputed["Time",])
NHP4026peak <- peakPhase(NHP4026imputed[!rownames(NHP4026imputed) %in% "Time",],
                         NHP4026alltimes, method = 'FFT')

plot(NF54peak[,1], NHP4026peak[,1],
     xlab = "NF54 Peak Time", ylab = "NHP4026 Peak Time", pch = 16, cex = 0.5)


## -----------------------------------------------------------------------------

NF54peak <- peakPhase(NF54imputed[!rownames(NF54imputed) %in% "Time",],
                         NF54alltimes, method = 'poly')


NHP4026peak <- peakPhase(NHP4026imputed[!rownames(NHP4026imputed) %in% "Time",],
                         NHP4026alltimes, method = 'poly')

plot(NF54peak[,1], NHP4026peak[,1],
     xlab = "NF54 Peak Time", ylab = "NHP4026 Peak Time", pch = 16, cex = 0.5)


## ---- eval = FALSE------------------------------------------------------------
#  library(PFExpTools)
#  
#  Mokdata <- fullCurate(GSE = "GSE59097")
#  Mokexp <- Mokdata[[1]]
#  Mokmeta <- Mokdata[[2]]

## ---- eval = FALSE------------------------------------------------------------
#  Mokdata <- fullCurate(GSE = "GSE59097",
#                        method = "blast",
#                        transcripts = "current",
#                        platcols = c(1,2),
#                        aliases = NA,
#                        match = 130,
#                        secmatch = 60,
#                        pullmeta = T,
#                        pct = 0.8)

## ---- eval = FALSE------------------------------------------------------------
#  
#  Gonzalesdata <- fullCurate(GSE = "GSE12515",
#                             platcols = c(1,2),
#                             method = "alias")
#  

## ---- eval = FALSE------------------------------------------------------------
#  Mokdata <- fullCurate(GSE = "GSE59097",
#                        method = "swap",
#                        platcols = c(1,6))

