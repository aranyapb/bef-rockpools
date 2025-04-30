
# generate the relevant DAGs
get_dag_h1a <- function() {
  
  # create the model using dagitty
  dag_x <- dagitty('
      dag {
      bb="0,0,1,1"
      Biomass [outcome,pos="0.600,0.270"]
      Depth [pos="0.462,0.362"]
      Gamma [pos="0.340,0.273"]
      PC1 [pos="0.469,0.136"]
      PC2 [pos="0.463,0.510"]
      Depth -> Biomass
      Depth -> Gamma
      Gamma -> Biomass
      PC1 -> Biomass
      PC1 -> Gamma
      PC2 -> Biomass
      PC2 -> Gamma
    }
    ')
  return(dag_x)
  
}

get_dag_h1b <- function() {
  # create the model using dagitty
  dag_x <- dagitty('
    dag {
    bb="0,0,1,1"
    Biomass [outcome,pos="0.600,0.270"]
    Depth [pos="0.462,0.362"]
    Gamma [pos="0.340,0.273"]
    PC1 [pos="0.469,0.136"]
    PC2 [pos="0.463,0.510"]
    Depth -> Biomass
    Depth -> Gamma
    Gamma -> Biomass
    PC1 -> Biomass
    PC1 -> Gamma
    PC2 -> Biomass
    PC2 -> Depth
    PC2 -> Gamma
    }')
  return(dag_x)
}





