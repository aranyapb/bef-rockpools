
# Aranya thesis: DAGs

# tutorial: https://rpubs.com/XC_Prime/dagitty

# load the dagitty library
library(dagitty)

# biomass ~ gamma

# set-up the DAG:
DAG1 <- dagitty(x = 'dag {
  
  BIOMASS [pos = "0.3,0.3"]
  CLIMATE [pos = "0.25,0.35"]
  DEPTH [pos = "0.35,0.35"]
  GAMMA [pos = "0.2,0.3"]
  ALPHA [pos = "0.4,0.3"]
  
  CLIMATE -> BIOMASS
  CLIMATE -> GAMMA
  DEPTH -> BIOMASS
  DEPTH -> ALPHA
  GAMMA -> BIOMASS
  ALPHA -> BIOMASS
  GAMMA -> ALPHA
  
}')
plot(DAG1)

# JOREN'S NOODLING:
DAG1 <- dagitty(x = 'dag {
  
  bb="-2.149,-2.511,2.518,3.303"
  ALPHA [pos="1.464,-0.066"]
  BIOMASS [pos="0.284,-0.075"]
  CLIMATE [pos="-0.458,1.316"]
  DEPTH [pos="0.949,1.281"]
  GAMMA [pos="-0.914,-0.075"]

  ALPHA -> BIOMASS
  CLIMATE -> BIOMASS
  CLIMATE -> GAMMA
  DEPTH -> ALPHA
  DEPTH -> BIOMASS
  GAMMA -> ALPHA [pos="0.329,-2.226"]
  GAMMA -> BIOMASS
  
}')
plot(DAG1)

# get the adjustment set for estimating the effect of gamma on biomass
adjustmentSets(x = DAG1, exposure = "GAMMA", outcome = "BIOMASS", effect = "direct")

# get the adjustment set for estimating the effect of climate on biomass
adjustmentSets(x = DAG1, exposure = "CLIMATE", outcome = "BIOMASS", effect = "direct")

# get the adjustment set for estimating the effect of depth on biomass
adjustmentSets(x = DAG1, exposure = "DEPTH", outcome = "BIOMASS", effect = "direct")

# get the adjustment set for estimating the effect of prominence on gamma
adjustmentSets(x = DAG1, exposure = "PROMINENCE", outcome = "GAMMA", effect = "direct")

# get the adjustment set for estimating the effect of prominence on gamma
adjustmentSets(x = DAG1, exposure = "CLIMATE", outcome = "GAMMA", effect = "direct")

# biomass ~ alpha

# set-up the DAG

DAG2 <- dagitty(x = 'dag {

    ALPHA [pos="0.325,0.5"]
    BIOMASS [pos="0.4,0.4"]
    DEPTH [pos="0.325,0.700"]
    INSELBERG [pos="0.25,0.4"]

    ALPHA -> BIOMASS
    DEPTH -> ALPHA
    DEPTH -> BIOMASS
    INSELBERG -> ALPHA
    INSELBERG -> BIOMASS
    INSELBERG -> DEPTH
    
    }')
plot(DAG2)

# get the adjustment set for estimating the effect of gamma on biomass
adjustmentSets(x = DAG2, exposure = "ALPHA", outcome = "BIOMASS", effect = "direct")

# get the adjustment set for estimating the effect of gamma on biomass
adjustmentSets(x = DAG2, exposure = "INSELBERG", outcome = "BIOMASS", effect = "direct")

# get the adjustment set for estimating the effect of climate on biomass
adjustmentSets(x = DAG2, exposure = "DEPTH", outcome = "BIOMASS", effect = "direct")

### END
