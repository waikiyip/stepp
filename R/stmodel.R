#################################################################
#
# stmodel.R
#
#######################
# stepp model         #
#######################
setClassUnion("stmodel",
  c("stmodelGLM", "stmodelKM", "stmodelCI")
)  
