# ===========================================
# -- computer the pmp and marginal for the 
# -- top 5 variables combined models in the 
# -- study of Fernandez, Ley and Steel(2001)
# ------Written by Shutong Ding--------------
# ------2015-06-03 --------------------------
# ===========================================





#================== Local Functions Begin =================


bin2dec = function(x) {
      sum(x * 2^(rev(seq(along=x)) - 1))
}

dec2bin <- function(x,k){
  tmp <- rev(as.numeric(intToBits(x)))
 # id <- seq_len(match(1,tmp,length(tmp))-1)
 #  tmp[-id]
 tmp[(32-k+1):32]
}

dec2bin_vec =  Vectorize(dec2bin,"x")

combn2model = function(kmax,sub_dim){
  ax = combn(kmax,sub_dim)
  num_mod  = ncol(ax)
  models = matrix(0,num_mod,kmax)
  gix = function(axi,kmax){
      gi = numeric(kmax)
      gi[axi] = 1
      return(gi)
  }
  
  models = t(apply(ax,2,gix,kmax))
  return(models)
  
}



# the log marginal likelihood based on Fernandez, Ley and Steel(2001)
# p. 566 (8)
require("MASS")
log_mml_gprior = function( xymat, midx, kmax, n, g = 1/min(length(y), kmax^2)){
  model = dec2bin(midx,kmax)
  mod_ind = c(1:kmax)*model
  mod_ind = mod_ind[mod_ind!=0] + 1
  kj = length(mod_ind)
  y = xymat[,1]
  if (kmax > (ncol(xymat) - 1) ) {
    stop("The maximum allowed X variables exceeds the limit of data!")
  }
  xmat = cbind(rep(1,n),xymat[,mod_ind])    
  # avoid inverse problem
  # xtx_inv = chol2inv(chol(crossprod(xmat)))   
  # use general inverse 
  xtx_inv = ginv(crossprod(xmat))   
  mxj = diag(1,n) - xmat%*%xtx_inv%*%t(xmat)                                
  lml = kj/2*log(g/(g+1))
  y_demean = y - mean(y)
  lml = lml - (n - 1)/2*log(1/(g+1)*crossprod(y,mxj%*%matrix(y,n,1)) + g/(g+1)*crossprod(y_demean))  
  lml
}

# vectorize the model
log_mml_gprior_vec = Vectorize(log_mml_gprior, "midx" )


# use 1:kmax index
#iks = 1
#ike = 0
#for (kj in 1:kmax) {
#  gc = combn(kmax,kj)
#  ngc = ncol(gc)
#  ike = ike + ngc   
#  models[iks:ike,1:kj] = t(gc)
#  iks = iks + ngc
#}

# use binary index
# vec <- c(0, 1)
# lst <- lapply(numeric(kmax), function(x) vec)
# models = as.matrix(expand.grid(lst))[-1,]
# midxs = apply(models, 1, bin2dec)


# compute the posterior model probabilites based on all combination of models
# using kmax number of variables, 2^kmax in total 
pmp_all_models = function(boot_vec, xymat, kmax , n, keep_cons = FALSE){
  
  # make all the permutations for the model, models
  if ( keep_cons ) {
    num_mod = 2^kmax
    mod_indx = 0:(num_mod-1)
  }
  else {
    num_mod = 2^kmax-1 
    mod_indx = 1:num_mod
  }
  
  lml_all = log_mml_gprior_vec(xymat[boot_vec,], mod_indx , kmax, n)
  lml_max = max(lml_all)
  lml_all = exp(lml_all - lml_max)
  c_lml   =  sum(lml_all)   
  return(lml_all/c_lml)
  
}


bootstrap_pmp_all = function(bootrep, xymat, kmax, keep_cons = FALSE) {
  
  if ( keep_cons ) {
    num_mod = 2^kmax
    mod_indx = 0:(num_mod-1)
  }
  else {
    num_mod = 2^kmax-1 
    mod_indx = 1:num_mod
  }
  
  models = matrix(0,num_mod,kmax)
  models = t(dec2bin_vec(mod_indx , kmax))


  n = nrow(xymat)
  boot_seeds = t(matrix(rep(1:n,bootrep), n, bootrep))
  
   gen_x = function(x){
     sample(x,length(x),replace=TRUE)
   }
  boot_seeds =  apply(boot_seeds,1,gen_x)
  pmp_all = t(apply(boot_seeds, 2, pmp_all_models, xymat, kmax, n, keep_cons))
  
  return(list("pmp"=pmp_all, "model"=models))
  
} 

pmp_sub_models = function(boot_vec, models, xymat, kmax , n, keep_cons = FALSE){
  
  # convert models indexes into indx
  if ( keep_cons ) {
    mod_indx = c(0,apply(models,1,bin2dec))
  }
  else {
    mod_indx = apply(models,1,bin2dec)
  }
  
  lml_all = log_mml_gprior_vec(xymat[boot_vec,], mod_indx , kmax, n)
  lml_max = max(lml_all)
  lml_all = exp(lml_all - lml_max)
  c_lml   =  sum(lml_all)   
  return(lml_all/c_lml)
  
}



bootstrap_pmp_sub = function(bootrep, models, xymat, kmax, keep_cons = FALSE) {
  
  n = nrow(xymat)
  boot_seeds = t(matrix(rep(1:n,bootrep), n, bootrep))
  
  gen_x = function(x){
    sample(x,length(x),replace=TRUE)
  }
  boot_seeds =  apply(boot_seeds,1,gen_x)
  pmp_sub = t(apply(boot_seeds, 2, pmp_sub_models, models, xymat, kmax, n, keep_cons))
  
  return(list("pmp"=pmp_sub, "model"=models))
  
} 












#================== Local Functions End ===================


#--------------- Begin the INPUT --------------

# load the data file 
data_dir = "C:/Users/shudi32/Dropbox/OverOptim/SD_part/Applications/Linear regression_Growth model/R top 5 model"
setwd(data_dir)
datmat   = read.csv(file="gm_steel_xymat.csv",header=T,sep=";")
kmax     = 5
xymat    = as.matrix(datmat[,1:(kmax+1)])
# 
n        = nrow(xymat)
# number of bootstrap replicates 
bootrep  = 100
# whether to keep the constant model in comparsion 
keep_cons= FALSE
# check sub-dimension models 
sub_dim  = 3
#--------------- End  the INPUT --------------- 

# -------------------Set up -------------------
# generate the model indexes, binary
models = combn2model(kmax,sub_dim)




# -------------- excution example ---------------
#log_mml_gprior(xymat, 1, kmax, n)
#log_mml_gprior(xymat, 2, kmax, n)
#log_mml_gprior(xymat, 3, kmax, n)
#log_mml_gprior_vec(xymat,1:3,kmax,n)
#pmp_all_models(1:n, xymat, kmax, n, keep_cons = FALSE )


pmp_boot_sub3 = bootstrap_pmp_sub(bootrep, models, xymat, kmax, keep_cons)

pmp_boot = bootstrap_pmp_all(bootrep, xymat, kmax, keep_cons)
write.table(file="pmp_all_boot.txt",pmp_boot$pmp, row.names = FALSE, col.names=FALSE)
write.table(file="pmp_all_models.txt",pmp_boot$model, row.names = FALSE, col.names=FALSE)


# make some plots 
require(coda)
pmp_sub3 = as.mcmc(pmp_boot_sub3$pmp) 
matplot(pmp_sub3,type="l")





