get.GMod_Free_OC_2<-function(N,pE,pT,ntrial,phiT,C_T,phi_E,delta){
  library(Iso)
  
  d_opt = rep(0,ntrial)
  ndose<-length(pE)
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = 1 + yytox
    bt = 1 + nn - yytox
    Tox_prob = 1 - pbeta(phiT, at, bt)
    AT_naive = which(Tox_prob < C_T)
    if (length(AT_naive)==0){
      AT=AT_naive
    }
    else{
      full_seq = seq(min(AT_naive),max(AT_naive),1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
    }
    return(AT)
  }
  
  #### mod
  est_MA<-function(n,yeff,num_dose){
    AIC=rep(0,num_dose)
    est = matrix(0,nrow = num_dose,ncol = num_dose)
    for (j_dagger in 1:num_dose) {
      if (j_dagger==1){
        yE1 = sum(yeff)
        n1 = sum(n)
        est1 = yE1/n1
        est[j_dagger,] = est1
        L = prod((est1^yE1)*(1-est1)^(n1-yE1))
        AIC[(j_dagger)] = -2*log(L)+2*length(unique(est1))
      }
      else{
        yE1_b_star =  yeff[1:(j_dagger-1)]
        yE1_j_dagger = sum(yeff[j_dagger:num_dose])
        yE1 = c(yE1_b_star,yE1_j_dagger)
        n1 = c(n[1:(j_dagger-1)],sum(n[j_dagger:num_dose]))
        est1 = pava(yE1/n1,w=n1)
        est[j_dagger,] = c(est1[1:(j_dagger-1)],
                           rep(est1[j_dagger],length(yeff[j_dagger:num_dose])))
        L = prod((est1^yE1)*(1-est1)^(n1-yE1))
        AIC[(j_dagger)] = -2*log(L)+2*length(unique(est1))
      }
    }
    wegt_ind = exp(-0.5*AIC)
    wegt = wegt_ind/sum(wegt_ind)
    est_final = c(wegt %*% est) 
    return(est_final)
  }
  
  
  
  for(trial in 1:ntrial){
    yE_ob = matrix(NA,nrow = ndose,ncol = max(N))
    for (j in 1:ndose) {
      yE_ob[j,1:N[j]]=rbinom(N[j],1,pE[j])
    }
    
    yE = apply(yE_ob,1,sum,na.rm=T)
    yE_all = c(yE_ob)
    yE_all = yE_all[!is.na(yE_all)]
    
    yT<-rbinom(ndose,N,pT)
    n<-N
    #########
    S=adm_tox(n=n,yT)
    final_est = est_MA(n=n,yeff = yE,ndose)[S]
    A1 = which(final_est>=delta*final_est[length(final_est)])
    A2 = which(final_est>phi_E)
    A = intersect(A1,A2)
    if(length(A)==0){
      d_opt[trial]<-0
    } else{
      d_opt[trial]<-min(A)
    }
  }
  
  selpercent = rep(0,ndose+1)
  ###Summarize results
  for (k in 0:ndose) {
    selpercent[(k+1)] = sum(d_opt == k)/ntrial*100
  }
  print("selection probability")
  cat(formatC(selpercent, digits = 1, format = "f"), sep = " ","\n")
}

get.GMod_Free_OC_2(N=c(10,10,10),pE=c(0.5,0.5,0.5),pT=c(0.1,0.15,0.18),
             ntrial = 10000,
             phiT = 0.3,C_T = 0.9,phi_E = 0.3,delta = 0.9)
##########
GMod_Free_decision2<-function(N,yE,yT,phi_E,delta,phiT,C_T){
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = 1 + yytox
    bt = 1 + nn - yytox
    Tox_prob = 1 - pbeta(phiT, at, bt)
    AT_naive = which(Tox_prob < C_T)
    if (length(AT_naive)==0){
      AT=AT_naive
    }
    else{
      full_seq = seq(min(AT_naive),max(AT_naive),1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
    }
    return(AT)
  }
  
  #### mod
  est_MA<-function(n,yeff,num_dose){
    AIC=rep(0,num_dose)
    est = matrix(0,nrow = num_dose,ncol = num_dose)
    for (j_dagger in 1:num_dose) {
      if (j_dagger==1){
        yE1 = sum(yeff)
        n1 = sum(n)
        est1 = yE1/n1
        est[j_dagger,] = est1
        L = prod((est1^yE1)*(1-est1)^(n1-yE1))
        AIC[(j_dagger)] = -2*log(L)+2*length(unique(est1))
      }
      else{
        yE1_b_star =  yeff[1:(j_dagger-1)]
        yE1_j_dagger = sum(yeff[j_dagger:num_dose])
        yE1 = c(yE1_b_star,yE1_j_dagger)
        n1 = c(n[1:(j_dagger-1)],sum(n[j_dagger:num_dose]))
        est1 = pava(yE1/n1,w=n1)
        est[j_dagger,] = c(est1[1:(j_dagger-1)],
                           rep(est1[j_dagger],length(yeff[j_dagger:num_dose])))
        L = prod((est1^yE1)*(1-est1)^(n1-yE1))
        AIC[(j_dagger)] = -2*log(L)+2*length(unique(est1))
      }
    }
    wegt_ind = exp(-0.5*AIC)
    wegt = wegt_ind/sum(wegt_ind)
    est_final = c(wegt %*% est) 
    return(est_final)
  }
  
  
  ndose<-length(yT)
  n<-N
  S=adm_tox(n=n,yT)
  final_est = est_MA(n=n,yeff = yE,ndose)[S]
  A1 = which(final_est>=delta*final_est[length(final_est)])
  A2 = which(final_est>phi_E)
  A = intersect(A1,A2)
  
  if(length(A)==0){
    d_opt<-0
  } else{
    d_opt<-min(A)
  }
  
  print( paste("Dose selection:", d_opt)) 
  print( paste("Toxicity admissible set:",paste(S, collapse = ",")))
  print( paste("Efficacy estimates:", paste(final_est , collapse = ",")))
}

GMod_Free_decision2(c(15,15,15),c(0,6,7),c(0,6,5),0.3,0.9,0.3,0.9)
