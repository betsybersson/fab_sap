library(readr)
library(dplyr)
library(stringr)
##############################################
## Organize data for use in prediction
##############################################

## read data

df = read_csv("./data/srrs2.dat")
cty = read_csv("./data/cty.dat")

mn = which(df$state %in% "MN")
radon = df$activity[mn]
radon = (radon + sqrt(radon^2+1))/2 # recommended by Price 1993
log.radon = log(radon)

n = length(radon)
y = log.radon

# get the county-level predictor
df.fips = paste0(str_pad(df$stfips,2,"left",pad=0),
                  str_pad(df$cntyfips,3,"left",pad=0),sep="")
usa.fips = paste0(str_pad(cty$stfips,2,"left",pad=0),
                   str_pad(cty$ctfips,3,"left",pad=0),sep="")
unq.fips = unique(df.fips[mn])
usa.rows = match(unq.fips, usa.fips)
uranium = cty$Uppm[usa.rows]
u = as.vector(log (uranium))

# get county index variable
fips.name = as.vector(df.fips[mn])
J = length(unq.fips)
group = rep(1:J,table(fips.name))

# long/lat for adj matrix
loc.df = dplyr::select(cty[usa.rows,],lon,lat)

# create row distance adj mat
W = matrix(data=0,nrow=J,ncol=J)
for ( j in 1:J ){
  
  lon.j = loc.df$lon[j]
  lat.j = loc.df$lat[j]
  
  dist.from.j = sapply(1:J,function(k)sqrt((lon.j-loc.df$lon[k])^2+
                                             (lat.j-loc.df$lat[k])^2) )
  stand.dist = exp(-dist.from.j^2)
  
  W[j,] = stand.dist
  
}
diag(W) = 0
colnames(W) = rownames(W) = unq.fips

## final org for prediction

# remove counties with only 1 data point
small_site = which(table(group)<2)
if(length(small_site)>0){
  loc = which(group%in%small_site)
  y = y[-loc]
  group = group[-loc]
  u = u[-small_site]
  
  W = W[-small_site,-small_site]
}

p = length(unique(group))
n = table(group)
# relabel groups 1: max group
group = rep(1:p,n)
# rename n
names(n) = c(1:p)

W.stand =  diag(1/rowSums(W)) %*% W #row standardize

