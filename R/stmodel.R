#################################################################
#
# stmodel.R
#
#######################
# stepp model         #
#######################
setClassUnion("stmodel",
  c("stmodelKM", "stmodelCI", "stmodelGLM")
)
