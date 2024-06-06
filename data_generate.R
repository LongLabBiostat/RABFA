## data generation 
library(Matrix)

############# Generate large ortho W ################
suffix_list = list()
suffix_list[[1]] = "all_ll"
suffix_list[[2]] = "all_l"
suffix_list[[3]] = "all_m"
suffix_list[[4]] = "all_s"
suffix_list[[5]] = "all_ss"
suffix_list[[6]] = "ap_ll"
suffix_list[[7]] = "ap_l"
suffix_list[[8]] = "ap_m"
suffix_list[[9]] = "ap_s"
suffix_list[[10]] = "ap_ss"
suffix_list[[11]] = "full_ll"
suffix_list[[12]] = "full_l"
suffix_list[[13]] = "full_m"
suffix_list[[14]] = "full_s"
suffix_list[[15]] = "full_ss"
suffix_list[[16]] = "ind_ll"
suffix_list[[17]] = "ind_l"
suffix_list[[18]] = "ind_m"
suffix_list[[19]] = "ind_s"
suffix_list[[20]] = "ind_ss"
L = 9
n = 100
H=3

##for low-dim ####
# p1=30
# p2=30
# p3=30
###################

##for high-dim ####
p1=90
p2=60
p3=30
###################
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)
pathway_list=list()
# for low-dim setting
# pathway_list[[1]] = c(10,10,10)
# pathway_list[[2]] = c(20,5,5)
# pathway_list[[3]] = c(15,10,5)

# for high-dim setting
pathway_list[[1]] = rep(15,6)
pathway_list[[2]] = rep(10,6)
pathway_list[[3]] = rep(5,6)

save_data_1 = list()
save_data_2 = list()
save_data_3 = list()
save_data_4 = list()


# ## for high-dim
W = matrix(0,p,L)
for (h in 1:H) {
  pathway_tmp = pathway_list[[h]]
  col_list = list()
  for (i in 1:3) {
    vec_tmp = rep(0,data_dim[h])
    if(i==1){
      vec_tmp[1:sum(pathway_tmp[1:i])] = rnorm(pathway_tmp[i],0,1.5)
    }else{
      vec_tmp[(sum(pathway_tmp[1:(i-1)])+1):sum(pathway_tmp[1:i])] = rnorm(pathway_tmp[i],0,1.5)
    }
    col_list[[i]] = vec_tmp
  }

  w_all =cbind(col_list[[1]],col_list[[2]],col_list[[3]] )

  if(h==1){
    w_ap  =cbind(c(rep(0,3*data_dim[h]/6),rep(1,1*data_dim[h]/6),rep(0,2*data_dim[h]/6))*rnorm(data_dim[h],0,1.5),c(rep(0,4*data_dim[h]/6),rep(1,1*data_dim[h]/6),rep(0,1*data_dim[h]/6))*rnorm(data_dim[h],0,1.5), rep(0,data_dim[h]) )
    w_ind =cbind(c(rep(0,5*data_dim[h]/6),rep(1,1*data_dim[h]/6))*rnorm(data_dim[h],0,1.5),rep(0,data_dim[h]),rep(0,data_dim[h]))
  }else if(h==2){
    w_ap  =cbind(c(rep(0,3*data_dim[h]/6),rep(1,1*data_dim[h]/6),rep(0,2*data_dim[h]/6))*rnorm(data_dim[h],0,1.5), rep(0,data_dim[h]),c(rep(0,4*data_dim[h]/6),rep(1,1*data_dim[h]/6),rep(0,1*data_dim[h]/6))*rnorm(data_dim[h],0,1.5) )
    w_ind =cbind(rep(0,data_dim[h]),c(rep(0,5*data_dim[h]/6),rep(1,1*data_dim[h]/6))*rnorm(data_dim[h],0,1.5),rep(0,data_dim[h]) )
  }else if (h==3){
    w_ap  =cbind(rep(0,data_dim[h]),c(rep(0,3*data_dim[h]/6),rep(1,1*data_dim[h]/6),rep(0,2*data_dim[h]/6))*rnorm(data_dim[h],0,1.5),c(rep(0,4*data_dim[h]/6),rep(1,1*data_dim[h]/6),rep(0,1*data_dim[h]/6))*rnorm(data_dim[h],0,1.5) )
    w_ind =cbind(rep(0,data_dim[h]),rep(0,data_dim[h]),c(rep(0,5*data_dim[h]/6),rep(1,1*data_dim[h]/6))*rnorm(data_dim[h],0,1.5) )
  }
  W[ind_s[h]:ind_e[h],] = cbind(w_all,w_ap,w_ind)
}
Z = matrix(rnorm(n*L,0,1.5),nrow=n,ncol=L)


## gaussian 
noise_level = c(1.5,1,0.5,0.25,0.1) 
for (nn in 1:length(noise_level)) {
  noise = noise_level[nn]
  if(nn==5){
    name_ind = c(5,10,15,20)
  }else{
    name_ind = which((seq(1,20,1)%%5==nn)==TRUE)  
  }
  
  full_ind = name_ind[3]
  all_ind  = name_ind[1]
  ap_ind   = name_ind[2]
  ind_ind  = name_ind[4]
  for (mm in 1:100) {
    X_full = W%*%t(Z) + rnorm(p*n,0,noise)
    X_all  = W[,1:3]%*%t(Z[,1:3]) +rnorm(p*n,0,noise) 
    X_ap   = W[,1:6]%*%t(Z[,1:6]) +rnorm(p*n,0,noise)  
    X_ind  = W[,7:9]%*%t(Z[,7:9]) + rnorm(p*n,0,noise)
    data_full = list("U" = Z, "V" = W,"X" = X_full)
    save_data_1[[mm]] = data_full
    data_ap   = list("U" = Z[,1:6],"V"=W[,1:6],"X" = X_ap)
    save_data_2[[mm]] = data_ap
    data_ind   = list("U" = Z[,7:9],"V"=W[,7:9],"X" = X_ind)
    save_data_3[[mm]] = data_ind
    data_all   = list("U" = Z[,1:3],"V"=W[,1:3],"X" = X_all)
    save_data_4[[mm]] = data_all
  }
  
  saveRDS(save_data_1,paste0('save_data_',suffix_list[[full_ind]],'.rds'))
  saveRDS(save_data_2,paste0('save_data_',suffix_list[[ap_ind]],'.rds'))
  saveRDS(save_data_3,paste0('save_data_',suffix_list[[ind_ind]],'.rds'))
  saveRDS(save_data_4,paste0('save_data_',suffix_list[[all_ind]],'.rds'))
  }











