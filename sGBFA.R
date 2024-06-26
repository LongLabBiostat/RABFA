library(MASS)
library(statmod)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(MCMCpack)
library(nloptr)
library(glmnet)
library(BayesLogit)


##################### working graph ###########################
## x : 1 refers  working graph = true graph
##     2 refers  working graph = adding edges between pathways 
##     3 refers  working graph = random graph
##     4 refers  working graph = removing edges from the true graph 
working_graph <- function(x,pathway_list,H,data_dim,ind_s,ind_e) {
  graph = matrix(0,p,p)
  if(x==1){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = 0.5 * (graph + t(graph))
    return(graph)  
  }else if(x==2){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
        graph_tmp[(node[nn]+1):sum(pathway_tmp[1:nn]),(node[nn]+1):sum(pathway_tmp[1:nn])] = rbinom((pathway_tmp[nn]-1)^2,1,0.3)
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = as.matrix(forceSymmetric(graph,uplo="U"))
    diag(graph)=0
    return(graph)
  }else if(x==3){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind] = rbinom(pathway_tmp[nn]-1,1,0.7)
        
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = 0.5 * (graph + t(graph))
    return(graph)  
  }else if(x==4){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
        # within pathway 
        graph_tmp[(node[nn]+1):sum(pathway_tmp[1:nn]),(node[nn]+1):sum(pathway_tmp[1:nn])] = rbinom((pathway_tmp[nn]-1)^2,1,0.3)
        # across pathway 
        if(sum(pathway_tmp[1:nn])+1 < sum(pathway_tmp)){
          graph_tmp[node[nn]:sum(pathway_tmp[1:nn]),(sum(pathway_tmp[1:nn])+1):sum(pathway_tmp)] = rbinom(pathway_tmp[nn]*(sum(pathway_tmp)-sum(pathway_tmp[1:nn])),1,0.2) 
        }
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = as.matrix(forceSymmetric(graph,uplo="U"))
    diag(graph)=0
    return(graph)
  }
}

# working_graph = function(x,p,L,mask_ind){
#   graph = matrix(0,nrow=p,ncol=p)
#   if(x==1){
#     center_node = c()
#     for (l in 1:L) {
#       tmp_ind = mask_ind[[l]]
#       candi_center = sample(tmp_ind,1)
#       while( candi_center %in% center_node){
#         candi_center = sample(tmp_ind,1)
#       }
#       center_node = c(center_node,candi_center)
#       for (i in 1:length(tmp_ind)) {
#         r_ind = candi_center
#         if(tmp_ind[i]!=candi_center){
#           c_ind = tmp_ind[i]
#           graph[r_ind,c_ind] = 1
#         }
#       }
#     }
#   }else if(x==2){
#     center_node = c()
#     for (l in 1:L) {
#       tmp_ind = mask_ind[[l]]
#       candi_center = sample(tmp_ind,1)
#       while( candi_center %in% center_node){
#         candi_center = sample(tmp_ind,1)
#       }
#       center_node = c(center_node,candi_center)
#       tmp_ind_rmc = tmp_ind[!tmp_ind==candi_center]
#       # center_node edge
#       for (i in 1:length(tmp_ind)) {
#         r_ind = candi_center
#         if(tmp_ind[i]!=candi_center){
#           c_ind = tmp_ind[i]
#           graph[r_ind,c_ind] = 1
#         }
#       }
#       
#       # surrounding node
#       for (i in 1:length(tmp_ind_rmc)) {
#         for (j in i:length(tmp_ind_rmc)) {
#           r_ind = tmp_ind_rmc[i] 
#           c_ind = tmp_ind_rmc[j]
#           graph[r_ind,c_ind] = rbinom(1,1,0.2)
#         }
#       }
#     }
#   }else if(x==3){
#     center_node = c()
#     for (l in 1:L) {
#       tmp_ind = mask_ind[[l]]
#       candi_center = sample(tmp_ind,1)
#       while( candi_center %in% center_node){
#         candi_center = sample(tmp_ind,1)
#       }
#       center_node = c(center_node,candi_center)
#       for (i in 1:length(tmp_ind)) {
#         r_ind = candi_center
#         if(tmp_ind[i]!=candi_center){
#           c_ind = tmp_ind[i]
#           graph[r_ind,c_ind] = rbinom(1,1,0.6)
#         }
#       }
#     }
#   }
#   graph = 0.5 * (graph + t(graph))
#   return(graph)
# }
############### true graph ##########################################
## x = 1 : star-like
## x = 2 : subset of the star-like
## x = 3 : between fully conected 
## x = 4 : adding edges beteen pathways

true_G = function(x,W){
  if(x==1){  
    # star-like pathway
    for (i in 1:L) {
      start_ind = i*n_nz-(n_nz-1)
      end_ind   = start_ind + n_nz-1
      true_graph[start_ind, (start_ind+1):end_ind] = rep(1,n_nz-1)
    }
    true_graph = (forceSymmetric(true_graph,uplo="U"))
    W = W
  }else if(x==2){
    # star-like pathway
    for (i in 1:L) {
      start_ind = i*n_nz-(n_nz-1)
      end_ind   = start_ind + n_nz-1
      true_graph[start_ind, (start_ind+1):end_ind] = rep(1,n_nz-1)
    }
    true_graph = (forceSymmetric(true_graph,uplo="U"))
    # remove the edges 
    ii = which(true_graph!=0,arr.ind = TRUE)    # edges 
    ii = ii[which((ii[,1]-ii[,2])<0),]      # edges in the upper matrix
    for (i in 1:dim(ii)[1]) {
      rm_edge = ii[i,]
      true_graph[rm_edge[1],rm_edge[2]] = rbinom(1,1,0.6)
    }
    true_graph = (forceSymmetric(true_graph,uplo="L"))
    W = W 
  }else if(x==3){
    # between star-like and fully connected graph
    for (i in 1:L) {
      start_ind = i*n_nz-(n_nz-1)
      end_ind   = start_ind + n_nz-1
      true_graph[start_ind, (start_ind+1):end_ind] = rep(1,n_nz-1)
      for (j in (start_ind+1):(end_ind-1)) {
        true_graph[j,(j+1):end_ind] = rbinom(end_ind-j,1,0.25)
      }
    }
    true_graph = (forceSymmetric(true_graph,uplo="U"))
    W = W
  }else if(x==4){
    # based on star-like: adding edges between pathways
    for (i in 1:L) {
      start_ind = i*n_nz-(n_nz-1)
      end_ind   = start_ind + n_nz-1
      true_graph[start_ind, (start_ind+1):end_ind] = rep(1,n_nz-1)
    }
    for ( i in 1:(L-1)) {
      row_l = 2+n_nz*(i-1)
      row_u = n_nz*i
      for (j in (i+1):L) {
        col_l = 2+n_nz*(j-1)
        col_u = n_nz*j
        block = matrix(rbinom((n_nz-1)^2,1,0.1),n_nz-1,n_nz-1)
        true_graph[row_l:row_u,col_l:col_u] = block
      }
    }
    true_graph = (forceSymmetric(true_graph,uplo="U"))
    W_c = W
    to_be_update =which(true_graph!=0,arr.ind = T)
    # to_be_update = to_be_update[which((to_be_update[,2] - to_be_update[,1])>0),]
    for (i in 1:dim(to_be_update)[1]) {
      ss = sum(W[to_be_update[i,][1],]*W[to_be_update[i,][2],])
      if (ss==0){
        W[to_be_update[i,][2],] =  W_c[to_be_update[i,][1],]+W_c[to_be_update[i,][2],]
      }
    }
  }
  ans = list(true_graph,W)
  return(ans)
}

true_G_s_2 = function(x,W){
  if(x==1){  
    # simpliest pathway
    for (i in 1:length(data_dim)) {
      start_ind = ind_s[i]
      end_ind   = p 
      true_graph[start_ind, (start_ind+1):end_ind] = rep(1,end_ind-start_ind)
    }
    true_graph = (forceSymmetric(true_graph,uplo="U"))
    W = W
  }
  ans = list(true_graph,W)
  return(ans)
}

###########################  empty chain ##########################
empty_chain <- function(n,p,T,L,H,rho_ini,phi_ini,w_ini,z_ini,alpha_ini){
  chain_m   = matrix(0,ncol = p, nrow = T+2)
  
  chain_rho   = matrix(0,ncol = p*n, nrow = T+2)
  chain_rho[1,]=as.vector(t(rho_ini))
  
  chain_w   = matrix(0,ncol = p*L, nrow = T+2)
  chain_w[1,]=as.vector(t(w_ini))
  
  chain_phi   = matrix(0,ncol = H*L, nrow = T+2)
  chain_phi[1,]=as.vector(t(phi_ini))
  
  chain_tau   = matrix(0,ncol = p*L, nrow = T+2)
  chain_tau[1,]=as.vector(t(tau_ini))
  
  chain_z   = matrix(0,ncol = L*n, nrow = T+2)
  chain_z[1,]=as.vector(t(z_ini))
  
  chain_alpha =matrix(0,ncol = p*L, nrow = T+2)
  chain_alpha[1,]=as.vector(t(alpha_ini))
  
  list(chain_tau=chain_tau,chain_alpha=chain_alpha,chain_phi=chain_phi,chain_rho=chain_rho,chain_m=chain_m,chain_w=chain_w,chain_z=chain_z)
}



## propto density of alpha_l
# alpha_l : candidate value
# alpha_h : current value
# alpha   : p by 1 alpha vector 
# omega   : the precison matrix 
# tau_l   : tau[j,l]^2
# j : the j-th entry of alpha 


p_density_mmh <- function(alpha_l,alpha_h,tau_l,phi_l,w_l,alpha,omega,j) {
  # alpha_l : candiate
  # alpha_h : old 
  lambda_l = exp(alpha_l)
  lambda_h = exp(alpha_h)
  alp_l_aug = replace(alpha,j,alpha_l)
  omega_j = omega[j,]
  ans = (lambda_l/lambda_h)^(3/2)*exp(-(lambda_l-lambda_h)*(tau_l*phi_l+(w_l^2)*phi_l/tau_l)/2)*exp((-1/(2*nu_2))*(alpha_l-alpha_h)*omega_j%*%(alpha+alp_l_aug-2*nu_1*rep(1,p)))
  return(ans)
}


sample_quantile = function(x){
  probs = c(0.025,0.975)
  ans = quantile(x,probs = probs)
  return(ans)
} 


#### MCMC_algorithm #### 

## arguments of the function

# dt : {0,1,2,3} = {gaussian, binary, binomial,mix}
# T : the number of iterations of MCMC
# L : the number of factors  
# p : the number of features
# n : the number of subjects 
# X : data 


# nu_1  :  mean 
# nu_2  :  variance 
# Sigma :  variance for each entry of m
# Q : covariance matrix for the proposal density
# eta : hyper parameter specifying the prior of omega
# eps : hyper parameter specifying the prior of omega 

# rho_temp : initial value for rho 
# tau_temp : initial value for tau 
# omega_temp : initial value for omega
# inv_omega_temp : initial value for inverse omega
# w_temp : initial value for w
# z_temp : initial value for z
# alpha_temp : initial value for 
# m_temp : initial value for

# th  : threshold to determine the degree of freedom (of W)

##########################################

## values of the function
# w_est : estimation of W
# z_est : estimation of Z
# m_est : estimation of m 
# df    : degree of freedom 
# BIC   : BIC criterion 
# chain_alpha : markov chain of alpha 
# chain_m     : markov chain of m 
# chain_w     : markov chain of W
# chain_z     : markov chain of Z 


sGBFA<- function(dt,T,L,p,n,X,trials,nu_1,nu_2,Sigma,Q,eta,eps,rho_temp,tau_temp,omega_temp,inv_omega_temp,w_temp,z_temp,alpha_temp,m_temp,start,end,H,mu){
  
  # preparation based on the data type 
  
  if(dt==0){
    psi = X
    phi =  as.matrix(rep(1,p))
    kappa = matrix(0,nrow = p,ncol = n)
    rowsum_kappa = rep(0,p)
  }else if(dt==1){
    psi = matrix(0,nrow = p,ncol = n)
    # b = matrix(1,ncol = n,nrow = p) # binary data 
    b = trials # binomial  data 
    kappa = X-b/2
    rowsum_kappa= rowSums(kappa)
  }else if(dt==2){
    psi = matrix(0,nrow = p,ncol = n)
    b = X + trials
    kappa = X-b/2
    rowsum_kappa= rowSums(kappa)
  }else if(dt==3){
    psi =  matrix(0,nrow = p,ncol = n)
    psi[1:n_nz,] = X[1:n_nz,]
    phi = as.matrix(rep(1,n_nz))
    b  = matrix(1,ncol = n,nrow = p)
    b[(2*n_nz+1):p,] = X[(2*n_nz+1):p,] + trials[(2*n_nz+1):p,] 
    kappa = X - b/2
    kappa[1:n_nz,] = 0
    rowsum_kappa= rowSums(kappa)
  }
  
  for (t in 1:T) {
    
    if(dt==0){
      rowsum_rho = rowSums(rho_temp)
      row_rho_psi = rowSums(rho_temp*psi)
    }else if(dt==1 |dt==2){
      rowsum_rho = rowSums(rho_temp)
      row_rho_psi = rep(0,p)
    }else if(dt==3){
      rowsum_rho = rowSums(rho_temp)
      row_rho_psi = rowSums(rho_temp*psi)
    }
    
    ## update m 
    for (j in 1:p) {
      covar_m = (rowsum_rho[j]+Sigma[j]^-1)^-1
      mean_m  = covar_m*(rowsum_kappa[j]+row_rho_psi[j]- (w_temp[j,]%*%z_temp)%*%rho_temp[j,])
      m_temp[j] = rnorm(1,mean_m,sqrt(covar_m))
    }
    
    chain_m[t+1,]=as.vector(m_temp)
    
    
    if(dt==0){
      for (j in 1:p) {
        # rho_j
        r = phi[j,1]+sum((X[j,]-(m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])))^2)
        rho_j = rgamma(1,shape=(phi[j,1]+n)/2 ,rate = r/2 )*rep(1,n)
        rho_temp[j,] = rho_j
      }
    }else if(dt==1 | dt==2){
      for (j in 1:p) {
        mu_j = (m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])) 
        rho_temp[j,] = rpg.devroye(n,b[j,],mu_j)
      }
    }else if(dt==3){
      if(j<=n_nz){
        r = phi[j,1]+sum((X[j,]-(m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])))^2)
        rho_j = rgamma(1,shape=(phi[j,1]+n)/2 ,rate = r/2 )*rep(1,n)
        rho_temp[j,] = rho_j
      }else{
        mu_j = (m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])) 
        rho_temp[j,] = rpg(n,b[j,],mu_j)
      }
    }
    
    
    
    chain_rho[t+1,]=as.vector(t(rho_temp))
    
    # generate alpha_l 
    
    for (l in 1:L) {
      for (j in 1:p) {
        for (len in 1:(length(ind_s)-1)) {
          if(j<ind_s[len+1] & j>=ind_s[len] ){
            tmp_h = len
            break
          }
        }
        candi = rnorm(1,alpha_temp[j,l],1)
        comp =  as.numeric(p_density_mmh(candi,alpha_temp[j,l],tau_temp[j,l],phi_temp[tmp_h,l],w_temp[j,l],alpha_temp[,l],omega_temp,j))   
        if(is.na(comp)){
          comp=0
        }
        # print(comp)
        prob = min(comp,1)
        u = runif(1,0,1)
        #  print(u)
        if(u<prob){
          alpha_temp[j,l]=candi
        }
      }
    }
    
    chain_alpha[t+1,]=as.vector(t(alpha_temp))
    
    
    ## calculate A  
    
    A = eta*(diag(eps,p) + 1) + (alpha_temp-nu_1)%*%t(alpha_temp-nu_1)/nu_2
    
    ## Begin to update Omega
    
    for (j in 1:p) {
      ## j-th diagonal element and j-th column excluding the diag 
      A_jj = A[j,j]     
      a_j = A[-j,j]   
      
      # construct the inverse of Omega_11 via inverse of Omega 
      inv_omega_11 = inv_omega_temp[-j,-j]
      inv_omega_12 = inv_omega_temp[-j,j]
      inv_omega_22 = inv_omega_temp[j,j]
      
      O_11 = inv_omega_11 - inv_omega_12%*%t(inv_omega_12)/inv_omega_22
      
      ## find the index of non-zeros of j-th col 
      nz_ind = which(omega_temp[,j]!=0) # include diagonal
      
      ## remove the diag index: index==j
      nz_ind = nz_ind[nz_ind!=j]        # exclude diagonal
      
      
      ## see if there is any non-zero entry except for the diag
      ## if length(nz_ind)=0: only generate the diag
      if(length(nz_ind)!=0){
        
        ## define the block matrices (mean part)
        ## i) select nz_ind rows and remove the j-th col 
        ## ii) remove the (nz_ind,j) rows and the j-th col 
        ## construct the precision matrix for \omega_12^{(1)}
        
        if(length(nz_ind)>1){
          prec_nz = inv_omega_temp[nz_ind, nz_ind]-inv_omega_temp[nz_ind,j]%*%t(inv_omega_temp[nz_ind,j])/inv_omega_temp[j,j]
          prec_nz_sqrt = chol(prec_nz)  # upper tri matrix 
          mid = rnorm(length(nz_ind),0,1/sqrt(A_jj))
          pt_mean_o1 = (-1/A_jj)*A[nz_ind,j]
          pt_mean_o2 = forwardsolve(t(prec_nz_sqrt),pt_mean_o1)
          mean_nz = backsolve(prec_nz_sqrt,pt_mean_o2)
          non_zero_omega = backsolve(prec_nz_sqrt,mid) + mean_nz 
          
          ## generate the diag entry
          
          xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
          partial = prec_nz_sqrt%*%non_zero_omega
          diag_omega = xi + t(partial)%*%partial
          
        }else if(length(nz_ind)==1){
          prec_nz = inv_omega_temp[nz_ind, nz_ind]-inv_omega_temp[nz_ind,j]%*%t(inv_omega_temp[nz_ind,j])/inv_omega_temp[j,j]
          mean_nz = (-1/(A_jj*prec_nz))*A[nz_ind,j]
          non_zero_omega = rnorm(1,mean_nz, 1/sqrt(prec_nz*A_jj))
          
          ## generate the diag entry
          
          xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
          diag_omega = xi + (non_zero_omega^2)*prec_nz
        }
        
        ## fill in the non-zero entries with the new sample
        ## both col and row d
        omega_temp[nz_ind,j] = non_zero_omega   # col 
        omega_temp[j,nz_ind] = non_zero_omega   # row
        omega_temp[j,j]  = diag_omega
      }else {
        
        ## when each entry = 0 except for the diag  
        ## we only generate the xi with Gamma
        xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
        diag_omega = xi
        omega_temp[j,j]  = diag_omega
      }
      
      # update inv_omega_temp 
      inv_omega_temp[j,j] = 1/xi
      pp = O_11%*%omega_temp[-j,j]
      inv_omega_11 = O_11 + inv_omega_temp[j,j]*pp%*%t(pp)
      inv_omega_12 = -inv_omega_11%*%omega_temp[-j,j]/omega_temp[j,j]
      
      inv_omega_temp[-j,-j] = inv_omega_11
      inv_omega_temp[-j,j] = inv_omega_12
      inv_omega_temp[j,-j] = inv_omega_12 
      
    }
    # generate z_i 
    
    for (i in 1:n) {
      pt_cov = t(w_temp)%*%(w_temp*rho_temp[,i])
      diag(pt_cov) = diag(pt_cov)+ rep(1,L)
      chol_pt_cov = chol(pt_cov) # upper matrix 
      pt_mean_z1 = t(w_temp)%*%(rho_temp[,i]*(psi[,i]-m_temp)+kappa[,i])
      pt_mean_z2 = forwardsolve(t(chol_pt_cov),pt_mean_z1)
      mean_z_i  = backsolve(chol_pt_cov,pt_mean_z2)
      z_temp[,i] = backsolve(chol_pt_cov,rnorm(L,0,1)) + mean_z_i
    }
    chain_z[t+1,]=as.vector(t(z_temp))
    ####################################################
    # update phi
    for (h in 1:H) {
      ind_h_s = ind_s[h]
      ind_h_e = ind_e[h]
      p_h = data_dim[h]
      for (l in 1:L) {
        phi_shape = a_phi + 3*p_h/2
        phi_rate  = b_phi + sum((tau_temp[ind_h_s:ind_h_e,l])*exp(alpha_temp[ind_h_s:ind_h_e,l]))/2+
          exp(alpha_temp[ind_h_s:ind_h_e,l])*(w_temp[ind_h_s:ind_h_e,l])^2/(2*(tau_temp[ind_h_s:ind_h_e,l]))
        phi_temp[h,l] = rgamma(1,shape=phi_shape, rate=phi_rate)
      }
    }
    chain_phi[t+1,]=as.vector(t(phi_temp))
    
    
    # update w_j & tau
    #################################################################
    for (j in 1:p) {
      for (len in 1:(length(ind_s)-1)) {
        if(j<ind_s[len+1] & j>=ind_s[len]){
          tmp_h = len
        }
      }
      cov_pt = z_temp%*%(t(z_temp)*rho_temp[j,])
      diag(cov_pt)= diag(cov_pt)+ (1/tau_temp[j,])*(exp(alpha_temp[j,])*phi_temp[tmp_h,])
      chol_cov_pt = chol(cov_pt)
      
      pt_mean_w1 = t(t(z_temp)*rho_temp[j,])%*%(psi[j,]-m_temp[j]*rep(1,n)+(rho_temp[j,]^-1)*kappa[j,])
      pt_mean_w2 = forwardsolve(t(chol_cov_pt),pt_mean_w1)
      mean_w_j = backsolve(chol_cov_pt,pt_mean_w2)
      w_temp[j,] = backsolve(chol_cov_pt,rnorm(L,0,1)) + mean_w_j
      tau_temp[j,]=1/rinvgauss(L,1/abs(w_temp[j,]),exp(alpha_temp[j,])*phi_temp[tmp_h,])
    }
    chain_w[t+1,]=as.vector(t(w_temp))
    chain_tau[t+1,]=as.vector(t(tau_temp))
    #print(paste0("Task Progress: ", t/T ))
  }
  
  mu_est = matrix(0,nrow=p, ncol=n)
  for (i in start:end) {
    w_t  = matrix(as.numeric(chain_w[i,]),nrow=p ,ncol=L,byrow=T)
    z_t  = matrix(as.numeric(chain_z[i,]),nrow=L ,ncol=n,byrow=T)
    m_t  = matrix(as.numeric(chain_m[i,]),nrow=p ,ncol=1,byrow=T)
    mu_est = mu_est + w_t%*%z_t + as.vector(m_t)
  }
  mu_est = mu_est/(end-start+1)
  # degree of freedom
  error = norm(mu_est-mu,type = "F")
  interval = apply(chain_w[start:end,],2,sample_quantile)
  df = sum(interval[1,]*interval[2,]>0)
  
  
  
  list(error=error,interval=interval,df=df,like=likelihood_g,chain_m=chain_m,chain_w=chain_w,chain_z=chain_z,mu_est=mu_est)
  #  list(error=error,interval=interval,df=df,prec_mat=prec_mat,bic_1=bic_1,bic_2=bic_2,chain_phi=chain_phi,chain_alpha=chain_alpha,chain_m=chain_m,chain_w=chain_w,chain_z=chain_z,mu_est=mu_est)
}





