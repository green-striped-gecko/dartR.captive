setOldClass(c("ggplot", "gg"))
# Slot descriptions for output
slotDescriptions <- NULL
slotDescriptions[["InputDf"]] <- "Original genlight input"
slotDescriptions[["SimOutput"]] <- "Genlight object of simulation output"
slotDescriptions[["MergedDf"]] <- "Final dataframe containing results of relatedness analysis"
slotDescriptions[["corOutList"]] <- "Results of correlation analysis"
slotDescriptions[["corVals"]] <- "Output of correlation results bewteen methods"
slotDescriptions[["plotList"]] <- "List of plots"

corrSlotDescriptions <- NULL
corrSlotDescriptions[["rmsePlot"]] <- "Plot of RMSE bewteen actual and estimated related values"
corrSlotDescriptions[["varPlot"]] <- "Plot of variance bewteen actual and estimated related values"

# Base output class (stores original genlight input)
setClass("OutputS4",
         slots = c(BaseInput = "genlight"))

# Display method for OutputS4
setMethod("show", "OutputS4", function(object){
  cat("//// New Object ////\n")
  for(Names in slotNames(object)){
    cat("  @", Names, "\n", sep="")
  }
})

# Base simulation class
setClass("DartSim",
         slots = c(
           type = "character",
           input_data = "genlight",
           table_input = "character",
           sim_input = "character", 
           gen_number = "numeric", 
           number_iterations = "numeric"
         ))

# Generic for simulation method
setGeneric("do_sim", function(object) standardGeneric("do_sim"))


# Method for DartSim
setMethod("do_sim", "DartSim", function(object) {
  ref_table <- gl.sim.WF.table(
    file_var = object@table_input,
    x = object@input_data,
    interactive_vars = FALSE
  )
  
  res_sim <- gl.sim.WF.run(
    file_var = object@sim_input,
    ref_table = ref_table,
    x = object@input_data,
    number_iterations = object@number_iterations,
    every_gen = 1,
    interactive_vars = FALSE,
    sample_percent = 100,
    gen_number_phase2 = object@gen_number
  )
  
  return(res_sim)
})


setClass("finalOutputClass", 
         slots = c(
           "InputDf" = "ANY", 
           "SimOutput" = "ANY", 
           "corOutList" = "ANY", 
           "corVals" = "ANY", 
           "plotList" = "ANY"
         ))

setMethod("show", signature =  "finalOutputClass", 
          definition = function(object){
            nameSlots <- slotNames(object)
            cat("********************\n")
            cat("***   OBJECTS    ***\n")
            cat("********************\n")
            for(i in nameSlots){
              currSlot <- slot(object, i)
              if(is.null(currSlot)){
                cat("  @",i,": NULL","\n", sep = "")
              }else{
                cat("  @",i,": ", slotDescriptions[[i]], "\n", sep = "")
              }
            }
          })





setClass("corOutList", 
         slots = c(
           "rmsePlot" = "ANY", 
           "varPlot" = "ANY"
         ))

setMethod("show", signature =  "corOutList", 
          definition = function(object){
            nameSlots <- slotNames(object)
            for(i in nameSlots){
              currSlot <- slot(object, i)
              if(is.null(currSlot)){
                cat("  @",i,": NULL","\n", sep = "")
              }else{
                cat("  @",i,": ", corrSlotDescriptions[[i]], "\n", sep = "")
              }
            }
          })


printCorVals <- function(corValues){
  
  cat("Number iterations:", length(corValues), "\n")
  for(i in 1:length(corValues)){
    cat("Iteration:", i, "\n")
    cat("Correlation between `related` and:","\n")
    colSelect <- colnames(corValues[[i]])[-which(colnames(corValues[[i]])=="rel")]
    for(j in 1:length(colSelect)){
      colSelect <- colnames(corValues[[i]])[-which(colnames(corValues[[i]])=="rel")]
      cat(colSelect[j],": ", corValues[[i]]["rel", colSelect[j]], "\n", sep = "")
    }
  }
}

setClass("corVals", 
         slots = c(
           "corVals" = "ANY"
         ))

setMethod("show", "corVals", function(object){
  printCorVals(object@corVals)
})

















