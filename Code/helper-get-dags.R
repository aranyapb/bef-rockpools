
# generate the relevant DAGs

# hypothesis 1

# causal hypothesis a
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

# causal hypothesis b
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

# hypothesis 2

# causal hypothesis a

# create the model using dagitty
get_dag_h2a <- function() {
  # create model using dagitty
  dag_x <- dagitty('
  dag {
  bb="0,0,1,1"
  Alpha [exposure,pos="0.492,0.246"]
  Biomass [outcome,pos="0.615,0.325"]
  Depth [pos="0.382,0.332"]
  Inselberg [pos="0.486,0.093"]
  Alpha -> Biomass
  Depth -> Alpha
  Depth -> Biomass
  Inselberg -> Alpha
  Inselberg -> Biomass
  Inselberg -> Depth
  }')
  return(dag_x)
}




