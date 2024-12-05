########### Simulation code for GMod-Free ############
get.GMod.Free.OC<-function(N,pE,pT,ntrial,alpha,boot_method,phiT,C_T,Delta,delta){
  
  LRT_Iso<-function(yE,n){
    p_null_orig<-sum(yE)/sum(n)
    L_null_orig<-(p_null_orig^sum(yE))*(1-p_null_orig)^(sum(n)-sum(yE))
    p_alter_orig<-pava(y=yE/n,w=n)
    L_alter_orig<-prod((p_alter_orig^yE)*(1-p_alter_orig)^(n-yE))
    test_stat_orig = -2*log(L_null_orig/L_alter_orig)
    test_stat = rep(0,10000)
    num_dose = length(n)
    for (i in 1:10000) {
      y_star = matrix(NA,nrow = num_dose,ncol = max(n))
      for (k in 1:num_dose) {
        y_star[k,1:n[k]] = rbinom(n[k],1,p_null_orig)
      }
      y_E = apply(y_star,1,sum,na.rm=T)
      p_null<-sum(y_E)/sum(n)
      L_null<-(p_null^sum(y_E))*(1-p_null)^(sum(n)-sum(y_E))
      p_alter<-pava(y=y_E/n,w=n)
      L_alter<-prod((p_alter^y_E)*(1-p_alter)^(n-y_E))
      test_stat[i] = -2*log(L_null/L_alter)
    }
    return(sum(test_stat>test_stat_orig)/10000)
  }
  
  LRT_Iso2<-function(yE,n,y){
    p_null_orig<-sum(yE)/sum(n)
    L_null_orig<-(p_null_orig^sum(yE))*(1-p_null_orig)^(sum(n)-sum(yE))
    p_alter_orig<-pava(y=yE/n,w=n)
    L_alter_orig<-prod((p_alter_orig^yE)*(1-p_alter_orig)^(n-yE))
    test_stat_orig = -2*log(L_null_orig/L_alter_orig)
    test_stat = rep(0,10000)
    num_dose = length(n)
    for (i in 1:10000) {
      nw_sample = sample(y,replace = T)
      y_star = matrix(NA,nrow = num_dose,ncol = max(n))
      y_star[1,1:n[1]] = nw_sample[1:n[1]]
      for (k in 2:num_dose) {
        y_star[k,1:n[k]] = nw_sample[(cumsum(n)[k-1]+1):cumsum(n)[k]]
      }
      y_E = apply(y_star,1,sum,na.rm=T)
      p_null<-sum(y_E)/sum(n)
      L_null<-(p_null^sum(y_E))*(1-p_null)^(sum(n)-sum(y_E))
      p_alter<-pava(y=y_E/n,w=n)
      L_alter<-prod((p_alter^y_E)*(1-p_alter)^(n-y_E))
      test_stat[i] = -2*log(L_null/L_alter)
    }
    return(sum(test_stat>test_stat_orig)/10000)
  }
  
  poc=rep(0,ntrial)
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
    AIC=rep(0,num_dose-1)
    est = matrix(0,nrow = (num_dose-1),ncol = num_dose)
    for (j_dagger in 2:num_dose) {
      yE1_b_star =  yeff[1:(j_dagger-1)]
      yE1_j_dagger = sum(yeff[j_dagger:num_dose])
      yE1 = c(yE1_b_star,yE1_j_dagger)
      n1 = c(n[1:(j_dagger-1)],sum(n[j_dagger:num_dose]))
      est1 = pava(yE1/n1,w=n1)
      est[j_dagger-1,] = c(est1[1:(j_dagger-1)],
                           rep(est1[j_dagger],length(yeff[j_dagger:num_dose])))
      L = prod((est1^yE1)*(1-est1)^(n1-yE1))
      AIC[(j_dagger-1)] = -2*log(L)+2*length(unique(est1))
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
    
    if(boot_method==1){
      pvalue<-LRT_Iso(yE,n)
    }else
    {pvalue<-LRT_Iso2(yE,n,yE_all)
    }
    if(pvalue>alpha){
      poc[trial]=0
      d_opt[trial]<-0
    }else{
      poc[trial]=1
      #########
      S=adm_tox(n=n,yT)
      final_est = est_MA(n=n,yeff = yE,ndose)[S]
    
      A1 = which(final_est>=delta*final_est[length(final_est)])
      A2 = which(final_est-final_est[1]>=Delta)

      A = intersect(A1,A2)
      if(length(A)==0){
        d_opt[trial]<-0
      } else{
        d_opt[trial]<-min(A)
      }
    }
  }
  poc_percentage = sum(poc)/ntrial*100
  print("PoC verified percentage")
  cat(formatC(poc_percentage, digits = 1, format = "f"), sep = " ","\n")
  selpercent = rep(0,ndose+1)
  ###Summarize results
  for (k in 0:ndose) {
    selpercent[(k+1)] = sum(d_opt == k)/ntrial*100
  }
  print("selection probability")
  cat(formatC(selpercent, digits = 1, format = "f"), sep = " ","\n")
}

GMod_Free_OC(N=c(20,20,20,20),pE=c(0.12,0.45,0.45,0.45),
             pT=c(0.05,0.1,0.15,0.2),ntrial=100,alpha=0.1,
             boot_method=1,phiT=0.3,C_T=0.9,Delta=0.2,delta=0.9)


######### Real trial implementation code ########
GMod_Free_decision<-function(N,yE,yE_ind,yT,alpha,boot_method,phiT,C_T,Delta,delta){
  
  LRT_Iso<-function(yE,n){
    p_null_orig<-sum(yE)/sum(n)
    L_null_orig<-(p_null_orig^sum(yE))*(1-p_null_orig)^(sum(n)-sum(yE))
    p_alter_orig<-pava(y=yE/n,w=n)
    L_alter_orig<-prod((p_alter_orig^yE)*(1-p_alter_orig)^(n-yE))
    test_stat_orig = -2*log(L_null_orig/L_alter_orig)
    test_stat = rep(0,10000)
    num_dose = length(n)
    for (i in 1:10000) {
      y_star = matrix(NA,nrow = num_dose,ncol = max(n))
      for (k in 1:num_dose) {
        y_star[k,1:n[k]] = rbinom(n[k],1,p_null_orig)
      }
      y_E = apply(y_star,1,sum,na.rm=T)
      p_null<-sum(y_E)/sum(n)
      L_null<-(p_null^sum(y_E))*(1-p_null)^(sum(n)-sum(y_E))
      p_alter<-pava(y=y_E/n,w=n)
      L_alter<-prod((p_alter^y_E)*(1-p_alter)^(n-y_E))
      test_stat[i] = -2*log(L_null/L_alter)
    }
    return(sum(test_stat>test_stat_orig)/10000)
  }
  
  LRT_Iso2<-function(yE,n,y){
    p_null_orig<-sum(yE)/sum(n)
    L_null_orig<-(p_null_orig^sum(yE))*(1-p_null_orig)^(sum(n)-sum(yE))
    p_alter_orig<-pava(y=yE/n,w=n)
    L_alter_orig<-prod((p_alter_orig^yE)*(1-p_alter_orig)^(n-yE))
    test_stat_orig = -2*log(L_null_orig/L_alter_orig)
    test_stat = rep(0,10000)
    num_dose = length(n)
    for (i in 1:10000) {
      nw_sample = sample(y,replace = T)
      y_star = matrix(NA,nrow = num_dose,ncol = max(n))
      y_star[1,1:n[1]] = nw_sample[1:n[1]]
      for (k in 2:num_dose) {
        y_star[k,1:n[k]] = nw_sample[(cumsum(n)[k-1]+1):cumsum(n)[k]]
      }
      y_E = apply(y_star,1,sum,na.rm=T)
      p_null<-sum(y_E)/sum(n)
      L_null<-(p_null^sum(y_E))*(1-p_null)^(sum(n)-sum(y_E))
      p_alter<-pava(y=y_E/n,w=n)
      L_alter<-prod((p_alter^y_E)*(1-p_alter)^(n-y_E))
      test_stat[i] = -2*log(L_null/L_alter)
    }
    return(sum(test_stat>test_stat_orig)/10000)
  }
  
  
  
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
    AIC=rep(0,num_dose-1)
    est = matrix(0,nrow = (num_dose-1),ncol = num_dose)
    for (j_dagger in 2:num_dose) {
      yE1_b_star =  yeff[1:(j_dagger-1)]
      yE1_j_dagger = sum(yeff[j_dagger:num_dose])
      yE1 = c(yE1_b_star,yE1_j_dagger)
      n1 = c(n[1:(j_dagger-1)],sum(n[j_dagger:num_dose]))
      est1 = pava(yE1/n1,w=n1)
      est[j_dagger-1,] = c(est1[1:(j_dagger-1)],
                           rep(est1[j_dagger],length(yeff[j_dagger:num_dose])))
      L = prod((est1^yE1)*(1-est1)^(n1-yE1))
      AIC[(j_dagger-1)] = -2*log(L)+2*length(unique(est1))
    }
    wegt_ind = exp(-0.5*AIC)
    wegt = wegt_ind/sum(wegt_ind)
    est_final = c(wegt %*% est) 
    return(est_final)
  }
  
  
  ndose<-length(yT)
  n<-N
  yE_all = yE_ind
  
  if(boot_method==1){
    pvalue<-LRT_Iso(yE,n)
  }else
  {pvalue<-LRT_Iso2(yE,n,yE_all)
  }
  
  if(pvalue>alpha){
    poc=0
    d_opt<-0
  }else{
    poc=1
    #########
    S=adm_tox(n=n,yT)
    efficacy_est = est_MA(n=n,yeff = yE,ndose)
    final_est = est_MA(n=n,yeff = yE,ndose)[S]
    #A = which(final_est>=0.8*final_est[3] & final_est>=phiE)
    A1 = which(final_est>=delta*final_est[length(final_est)])
    A2 = which(final_est-final_est[1]>=Delta)
    #A2 = which(final_est>=phiE)
    A = intersect(A1,A2)
    if(length(A)==0){
      d_opt<-0
    } else{
      d_opt = min(A)-1
    }
  }
  print( paste("PoC p-value", pvalue)) 
  print( paste("PoC verified?", ifelse(poc==1,"Yes","No"))) 
  print( paste("Dose selection:", d_opt)) 
  print( paste("Toxicity admissible set:",paste(S-1, collapse = ",")))
  print( paste("Efficacy estimates:", paste(efficacy_est, collapse = ",")))
  
}

GMod_Free_decision(c(15,15,15),c(0,6,7),c(rep(0,15),0,1,0,1,1,0,0,1,0,0,0,1,0,0,1,1,1,0,0,1,0,1,0,1,0,0,0,1,1,0),
                   c(0,6,5),0.1,1,0.3,0.9,0.2,0.9)

GMod_Free_decision(c(15,15,15),c(0,6,7),c(rep(0,15),0,1,0,1,1,0,0,1,0,0,0,1,0,0,1,1,1,0,0,1,0,1,0,1,0,0,0,1,1,0),
                   c(0,6,5),0.1,2,0.3,0.9,0.2,0.9)
