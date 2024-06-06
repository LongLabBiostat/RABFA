# source("/home/qiyiwen/function_package/s-GBFA.R")

d_ind = 1
suffix_ind = 6
suffix_list = list()
suffix_graph = list()
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

suffix_graph[[1]] = "g_1"
suffix_graph[[2]] = "g_2"
suffix_graph[[3]] = "g_3"
suffix_graph[[4]] = "g_4"


for (mm in (2*d_ind-1):(2*d_ind)) {
  for (gg in 1:4) {
    arg = mm
    set.seed(19402+arg)
    n = 100
    H=3
    T = 1000
    # n_nz = p/L  # number of non-zero elements of each col (W)
    p1 = 30
    p2 = 30
    p3 = 30
    p = p1+p2+p3 
    data_dim = c(p1,p2,p3)
    ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
    ind_e = c(p1,p1+p2,p1+p2+p3)
    
    graph_ind = gg
    ## Data Generating or import
    # save_data_mask = readRDS("../save_data.rds")
    data_1 = save_data_mask[[mm]]
    W = data_1$V
    Z = t(data_1$U)
    X = data_1$X
    L = dim(W)[2]
    grid_L  = c(L-1,L,L+1)
    # X: p x n matrix 
    mu = W %*% Z
    trials = matrix(0,nrow = p, ncol = n)
    data_tbs = list('true_W'= W,'true_Z'=Z,'true_X'=X)
    #####################################################
    
    Sigma = rep(0.01,p)  # variance for each entry of m
    
    Q = diag(4,nrow = p,ncol = p) # proposal density
    
    ######## working graph####################################
    pathway_list=list()
    pathway_list[[1]] = c(10,10,10)
    pathway_list[[2]] = c(20,5,5)
    pathway_list[[3]] = c(15,10,5)
    
    graph = working_graph(graph_ind,pathway_list=pathway_list,H=H,data_dim =data_dim,ind_s=ind_s,ind_e=ind_e ) 
    data_tbs[['G']] = graph 
    ########################################################
    
    w_ini_l = list()
    z_ini_l = list()
    L_seed = 1113
    for (l in 1:length(grid_L)) {
      set.seed(l%%3+L_seed)
      LL= grid_L[l]
      w_ini_l[[l]]=matrix(rnorm(LL*p,0,1),nrow = p,ncol = LL)
    }
    
    for (l in 1:length(grid_L)) {
      set.seed(l%%3+L_seed)
      LL= grid_L[l]
      z_ini_l[[l]] = matrix(rnorm(LL*n,0,1),nrow=LL,ncol=n)
    }
    
    ## initialization 
    
    grid_nu_1 = c(-1,-0.5,0,0.5,1)   # mean 
    #  grid_nu_1 = c(1)   # mean 
    
    grid_nu_2 = c(0.5,1) # variance 
    #  grid_nu_2 = c(0.5) # variance 
    
    a_phi = 1
    b_phi = 1
    eta = 10
    eps = 0.2
    tun_arg = 1
    
    
    for (l in 1:length(grid_L)) {
      L = grid_L[l]
      rho_ini = matrix(1,nrow=p,ncol=n )
      
      tau_temp = matrix(1,nrow=p,ncol=L )
      tau_ini  = tau_temp
      
      w_ini = w_ini_l[[l]] 
      z_ini = z_ini_l[[l]] 
      
      omega_ini = diag(1,nrow=p,ncol=p)  # compatible with graph
      
      ind = which(graph!=0,arr.ind = T)  # ind of the nonzero elements
      for (i in 1:dim(ind)[1]) {
        if(ind[i,1]>ind[i,2]){
          r_ind = ind[i,1]
          c_ind = ind[i,2]
          omega_ini[r_ind,c_ind] = 0.05
        }
      }
      
      omega_temp = as.matrix(forceSymmetric(omega_ini,uplo = 'L'))
      
      inv_omega_temp = solve(omega_temp) 
      
      for (j in 1:length(grid_nu_1)) {
        for (i in 1:length(grid_nu_2)) {
          nu_1 = grid_nu_1[j]
          nu_2 = grid_nu_2[i]
          alpha_ini = matrix(nu_1,nrow=p,ncol=L )
          phi_ini = matrix(1,nrow=H,ncol=L)
          mcmc_box = empty_chain(n,p,T,L=L, H=H,rho_ini,phi_ini,w_ini,z_ini,alpha_ini)
          chain_rho = mcmc_box$chain_rho
          chain_w = mcmc_box$chain_w
          chain_z = mcmc_box$chain_z
          chain_alpha = mcmc_box$chain_alpha
          chain_m = mcmc_box$chain_m
          chain_tau = mcmc_box$chain_tau
          chain_phi = mcmc_box$chain_phi
          w_temp   = matrix(as.numeric(chain_w[1,]),nrow=p ,ncol=L,byrow=T)
          z_temp   = matrix(as.numeric(chain_z[1,]),nrow=L ,ncol=n,byrow=T)
          rho_temp   = matrix(as.numeric(chain_rho[1,]),nrow=p ,ncol=n,byrow=T)
          alpha_temp = matrix(as.numeric(chain_alpha[1,]),nrow=p ,ncol=L,byrow=T)
          phi_temp = matrix(as.numeric(chain_phi[1,]),nrow=H ,ncol=L,byrow=T)
          m_temp = rep(0,p)
          trial_result = sGBFA(0,T=T,L=L,p=p,n=n,X,trials=trials,nu_1=nu_1,nu_2=nu_2,Sigma,Q,eta,eps,rho_temp,tau_temp,omega_temp,inv_omega_temp,w_temp,z_temp,alpha_temp,m_temp,start=400,end=T,H=3,mu=mu)       
          trial_name = paste0('mcmc_',tun_arg)
          data_tbs[[trial_name]] = trial_result
          tun_arg = tun_arg+1 
        }
      }
    }
    
  }
}



