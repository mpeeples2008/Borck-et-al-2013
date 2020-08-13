http://www.archaeologysouthwest.org/pdf/e-i_index_script.txt

# Feel free to use but please cite:
# Borck, Mills, Peeples, and Clark 2015 Are Social Networks Survival Networks? An Example from the Late Prehispanic U.S. Southwest, Journal of Archaeological Method and Theory 22(1):33-57


####################################################################################################
## Overview ########################################################################################
####################################################################################################

# This script calculates E-I (external-internal) indeces for both similarity matrices or binary 
# networks. The output includes standard and normalized E-I scores by site, standard and normalized 
# E-I scores by group (region), and the actual and expected E-I value for the dataset as a whole 
# based on a permutation test. 
# 
# This script expects two inputs. The first is a symmetric similarity matrix or binary network
# designated as 'sim' in R. The second is a vector or matrix of group (region) designations for 
# every row in the 'sim' matrix defined as 'grp'.
#

####################################################################################################
## Setup ###########################################################################################
####################################################################################################

####### initialize necessary libraries
library(sqldf)
library(statnet)


####################################################################################################
## E-I Index #######################################################################################
####################################################################################################

## Function for calculting the E-I index for all sites in the database
e.i.test <- function(x,y) {
  reg <- as.matrix(y)  # reg is the regional or group designation
  rd <- dim(x)[1]
  results <- matrix(0,rd,3)  
  for (s1 in 1:rd) {  
    for (s2 in 1:rd) {
      x1Temp <- reg[s1,1]
      x2Temp <- reg[s2,1]
      if (x1Temp==x2Temp) {results[s1,2] <- results[s1,2]+x[s1,s2]} 
      if (x1Temp!=x2Temp) {results[s1,1] <- results[s1,1]+x[s1,s2]}}}
  for (j in 1:nrow(x)){  
    if (length(which(x<1 & x>0))>0)
    {results[j,2] <- results[j,2]-1}
    results[j,3] <- (results[j,1]-results[j,2]) / (results[j,1]+results[j,2])}
  row.names(results) <- row.names(x)
  results <- round(results,3)
  out <- as.data.frame(cbind(reg,results)) 
  colnames(out) <- c('Region','E','I','E-I Index')
  return(out)}

## Function for calculting normalized version of the E-I index for all sites in the database
e.i.norm <- function(x,y) {
  reg <- as.matrix(y)
  rd <- dim(x)[1]
  results <- matrix(0,rd,3)
  for (s1 in 1:rd) {
    for (s2 in 1:rd) {
      x1Temp <- reg[s1,1]
      x2Temp <- reg[s2,1]
      if (x1Temp==x2Temp) {results[s1,2] <- results[s1,2]+x[s1,s2]} 
      if (x1Temp!=x2Temp) {results[s1,1] <- results[s1,1]+x[s1,s2]}}}
  if (length(which(x<1 & x>0))>0)
  {results[,2] <- results[,2]-1}
  results[,1] <- results[,1]/max(results[,1])
  results[,2] <- results[,2]/max(results[,2])
  results[,3] <- (results[,1]-results[,2]) / (results[,1]+results[,2])
  row.names(results) <- row.names(x)
  results <- round(results,3)
  out <- as.data.frame(cbind(reg,results))
  colnames(out) <- c('Region','E','I','E-I Index')
  return(out)}

## Function for calculating and standardizing group (region) level E-I index
e.i.group <- function(x,y) {
  EIst <- NULL
  EI <- NULL
  Est <- NULL
  Ist <- NULL
  out <- sqldf('SELECT x.Region, Sum(x.E) AS E, Sum(x.I) AS I
FROM x
GROUP BY x.Region')
  grp2 <- as.matrix(table(y))
  grp2 <- as.matrix(grp2[which(grp2[,1]>0),])
  for (i in 1:nrow(out)) {
    Est[i] <- out[i,2]/(nrow(x)-grp2[i,1])
    Ist[i] <- out[i,3]/(grp2[i,1])}
  EI <- (out[,2]-out[,3])/(out[,2]+out[,3])
  EIst <- (Est-Ist)/(Est+Ist)
  EI <- round(EI,3)
  Est <- round(Est,3)
  Ist <- round(Ist,3)
  EIst <- round(EIst,3)
  out <- cbind(grp2,out,EI,Est,Ist,EIst)
  colnames(out) <- c('NumSites','Region','E','I','E-I Index','E stand','I stand','E-I stand')
  return(out)}


####################################################################################################
## E-I Permutation Test ############################################################################
####################################################################################################


## Function for calculating E-I index for n simulations using randomized versions of original input dataset (non-normalized)
## This may take several seconds to minutes depending on the size of the databse and the number of random runs selected
e.i.perm <- function(x,y) {
  nsim <- 500 ## number of simulations
  E <- NULL
  I <- NULL
  out <- NULL
  EI <- NULL
  final <- matrix(0,1,2)
  for (i in 1:nsim) {
    z <- sample(nrow(x),replace=F)
    x2 <- x
    x <- x[z,z]
    reg <- as.matrix(y)
    rd <- dim(x)[1]
    results <- matrix(0,rd,3)
    for (s1 in 1:rd) {
      for (s2 in 1:rd) {
        x1Temp <- reg[s1,1]
        x2Temp <- reg[s2,1]
        if (x1Temp==x2Temp) {results[s1,2] <- results[s1,2]+x[s1,s2]} 
        if (x1Temp!=x2Temp) {results[s1,1] <- results[s1,1]+x[s1,s2]}}}
    if (length(which(x<1 & x>0))>0) {
      for (j in 1:nrow(x)){
        results[j,2] <- results[j,2]-1}}
    E <- sum(results[,1])
    I <- sum(results[,2])
    out[i] <- (E-I)/(E+I)
    out <- round(out,3)}
  
  results2 <- matrix(0,rd,3)
  for (s1 in 1:rd) {
    for (s2 in 1:rd) {
      x1Temp <- reg[s1,1]
      x2Temp <- reg[s2,1]
      if (x1Temp==x2Temp) {results2[s1,2] <- results2[s1,2]+x2[s1,s2]} 
      if (x1Temp!=x2Temp) {results2[s1,1] <- results2[s1,1]+x2[s1,s2]}}}
  if (length(which(x<1 & x>0))>0) 
  {for (j in 1:nrow(x)){
    results2[j,2] <- results2[j,2]-1}}
  E <- sum(results2[,1])
  I <- sum(results2[,2])
  EI <- (E-I)/(E+I)
  final <- cbind(EI,mean(out))
  colnames(final) <- c('Actual EI','Expected EI')
  return(final)}


####################################################################################################
### Function Calls #################################################################################
####################################################################################################

## sim == symmetric similarity matrix or binary network
## grp == matrix or vector of group (region) designations for every row of sim

# Calculate E-I index by site
site_ei <- e.i.test(sim,grp)

# Calculate normalized E-I index by site
norm_ei <- e.i.norm(sim,grp)

# Calculate and standardize group (region) level E-I index
group_ei <- e.i.group(site_ei,grp)
group_ei_norm <- e.i.group(norm_ei,grp) # use this for normalized scores

# Calculate actual and expected E-I values for dataset based on permutation test
ei_perm <- e.i.perm(sim,grp)


## End of Script

