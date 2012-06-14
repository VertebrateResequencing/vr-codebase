###source("/lustre/scratch102/user/kw10/work/svpipeline/COHORT/klaudia_pipeline/load_nsta_functions.R")

### Functions ###

### Generate a proper deletion density distribution by adding
### endpoints with probablity zero to the left (min-delta.x) and to
### the right (max+delta.x) of an empirical deletion density and
### normalise the distribution
# dd = empirical density of the expected deletion distribution
get.emp.density <- function(dd)  {
  delta.x <- dd$x[2] - dd$x[1]
  dd$x <- c(min(dd$x)-delta.x, dd$x, max(dd$x)+delta.x)
  dd$y <- c(0, dd$y, 0)
  dd$y <- dd$y/(sum(dd$y)*delta.x)
  dd
}

### Fits a density function to insert size distribution and
### finds the RP separation below density y = 0.001/100
### which is the slope of the tangent of the CDF (=quantiles plot)
### which corresponds to a change of 0.001 new calls within next 100
get.density.thresh <- function(f)  {
  dd <- density(f[f>0 & f < 10000], bw=1, n=1024)  # bandwidth=1, n=1024, default kernal = "gaussian"
  dd <- get.emp.density(dd)
  yt <- 0.001/100
  qu <- round(quantile(f[f>0], probs=0.99))
  indx <- which(dd$x > qu)
  indy <- which(dd$y[indx] < yt)
  x <- round(dd$x[indx][indy][1])
  x
}

########
## read in file with insert sizes in one column

args <- commandArgs(TRUE) 
infile = args[1]
f <- read.delim(infile, stringsAsFactors=F, header=F)  

if (length(f)!=0)  {
    if (any(f>0 & f < 10000))  {
      sde <- round(sd(f[f > 0 & f < 1000]))
      emp.thresh <- get.density.thresh(f)
    }  else  {
      sde <- NA
      emp.thresh <- 0
    }
	write(emp.thresh, stdout())
}
