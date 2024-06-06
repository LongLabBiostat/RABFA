# low-dim
L = 3
n = 50
H=3
p1=30
p2=30
p3=30 
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)
pathway_list=list()
pathway_list[[1]] = rep(p1/3,3)
pathway_list[[2]] = rep(p2/3,3)
pathway_list[[3]] = rep(p3/3,3)

save_data_1 = list()
save_data_2 = list()
save_data_3 = list()
save_data_4 = list()


W = matrix(0,p,L)
for (h in 1:H) {
  pathway_tmp = pathway_list[[h]]
  col_list = list()
  for (i in 1:length(pathway_tmp)) {
    vec_tmp = rep(0,data_dim[h])
    if(i==1){
      vec_tmp[1:sum(pathway_tmp[1:i])] = rnorm(pathway_tmp[i],1.5,0.1)
    }else{
      vec_tmp[(sum(pathway_tmp[1:(i-1)])+1):sum(pathway_tmp[1:i])] = rnorm(pathway_tmp[i],1.5,0.1)
    }
    col_list[[i]] = vec_tmp
  }
  
  w_all =cbind(col_list[[1]],col_list[[2]],col_list[[3]] ) 
  W[ind_s[h]:ind_e[h],] = w_all
}
Z = matrix(rnorm(n*L,0,1),nrow=n,ncol=L)

##
for (mm in 1:100) {
  X_gaussian = W%*%t(Z) + rnorm(p*n,0,2)
  trials = matrix(0,nrow=p,ncol=n)
  data_gaussian = list("U" = Z, "V" = W,"X" = X_gaussian,"trials"=trials)
  save_data_1[[mm]] = data_gaussian
}
saveRDS(save_data_1,paste0('save_data_1.rds'))

##
for (mm in 1:100) {
  X_binary = matrix(0,nrow=p,ncol=n)
  trials = matrix(1,nrow=p,ncol=n)
  mu = W %*% t(Z)
  for (j in 1:p) {
    for (i in 1:n) {
      X_binary[j,i] = rbinom(1,1,1/(1+exp(-mu[j,i])))
    }
  }
  data_binary= list("U" = Z, "V" = W,"X" = X_binary,"trials"=trials)
  save_data_2[[mm]] = data_binary
}
saveRDS(save_data_2,paste0('save_data_2.rds'))

trials = matrix(0,nrow=p,ncol=n)
for (j in 1:p) {
  trials[j,] = sample(seq(1,10,1),1)*rep(1,n)
}
for (mm in 1:100) {
  mu = W %*% t(Z)
  X_binomial = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      X_binomial[j,i] = rbinom(1,trials[j,i],1/(1+exp(-mu[j,i])))
    }
  }
  data_binomial= list("U" = Z, "V" = W,"X" = X_binomial,"trials"=trials)
  save_data_3[[mm]] = data_binomial
}
saveRDS(save_data_3,paste0('save_data_3.rds'))


trials_1 = matrix(0,nrow=p,ncol=n)
trials_2 = matrix(1,nrow=p2,ncol=n)
trials_3 = matrix(0,nrow=p3,ncol=n)
for (j in 1:p3) {
  trials_3[j,] = sample(seq(1,5,1),1)*rep(1,n)
}
for (mm in 1:100) {
  mu = W %*% t(Z)
  mu_1 = mu[1:p1,]
  mu_2 = mu[(1+p1):(p1+p2),]
  mu_3 = mu[(p1+p2+1):p,]
  
  X_gaussian = mu_1 + rnorm(p1*n,0,1)
  
  X_binary = matrix(0,nrow=p2,ncol=n)
  for (j in 1:p2) {
    for (i in 1:n) {
      X_binary[j,i] = rbinom(1,trials_2[j,i],1/(trials_2[j,i]+exp(-mu_2[j,i])))
    }
  }
  
  X_binomial = matrix(0,nrow=p3,ncol=n)
  for (j in 1:p3) {
    for (i in 1:n) {
      X_binomial[j,i] = rbinom(1,trials_3[j,i],1/(1+exp(-mu_3[j,i])))
    }
  }
  
  X_mix =rbind(X_gaussian,X_binary,X_binomial)
  trials =rbind(trials_1,trials_2,trials_3)
  data_mix= list("U" = Z, "V" = W,"X" = X_mix,"trials"=trials)
  save_data_4[[mm]] = data_mix
}
saveRDS(save_data_4,paste0('save_data_4.rds'))


