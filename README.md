# Borck-et-al-2015
Code from Borck et al. 2013 article in Journal of Archaeological Method and Theory

Borck, Lewis, Barbara J. Mills, Matthew A. Peeples, and Jeffery J. Clark. 
2015. Are Social Networks Survival Networks? An Example from the Late Prehispanic Southwest. Journal of Archaeological Method and Theory 22(1):33-57.


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
