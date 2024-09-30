library(intsurv)
library(survminer)
library(survival)
library(readr)
library(randomForestSRC)
library(ranger)
library(alabama)
library(dplyr)

# CHANGE THE TAU VALUE!!!!
# CHECK LIN PRED!
# CHECK PATH.CC
iter=1
num_iter=1
tau = 0.0005
gamma.weibull = 0.7
# path.cc <-"/Users/fanyiyang/Desktop/JAS\ revision/"
path.cc <- "../"

c_index=function(time,status,risk_score,beta){
  permissible <- 0 # comparable pairs
  concord <- 0 # completely concordance
  n <- length(time)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Exclude incomparable pairs
      if ((time[i] < time[j] & status[i] == 1 ) ) {
        permissible <- permissible + 1
        concord = concord + 1/(1+exp(-beta*(risk_score[i]-risk_score[j])))
      }
      
      if((time[j] < time[i] & status[j] == 1)){
        permissible <- permissible + 1
        concord = concord + 1/(1+exp(-beta*(risk_score[j]-risk_score[i])))
      }
      
    }
  }
  C_index <- concord / permissible
  names(C_index) <- "C index"
  return(round(C_index, 6))
}

beta=20
num_dataset=500/num_iter
initial_loss=rep(0,num_dataset)
super_loss=rep(0,num_dataset)
forest_loss=rep(0,num_dataset)
initial_ibs=rep(0,num_dataset)
super_ibs=rep(0,num_dataset)
forest_ibs=rep(0,num_dataset)

super_surv_mean_one_thrid=rep(0,num_dataset)
super_surv_median_one_thrid=rep(0,num_dataset)
initial_surv_mean_one_thrid=rep(0,num_dataset)
initial_surv_median_one_thrid=rep(0,num_dataset)
forest_surv_mean_one_thrid=rep(0,num_dataset)
forest_surv_median_one_thrid=rep(0,num_dataset)
super_haz_mean_one_thrid=rep(0,num_dataset)
super_haz_median_one_thrid=rep(0,num_dataset)
initial_haz_mean_one_thrid=rep(0,num_dataset)
initial_haz_median_one_thrid=rep(0,num_dataset)
forest_haz_mean_one_thrid=rep(0,num_dataset)
forest_haz_median_one_thrid=rep(0,num_dataset)

super_surv_mean_two_thrid=rep(0,num_dataset)
super_surv_median_two_thrid=rep(0,num_dataset)
initial_surv_mean_two_thrid=rep(0,num_dataset)
initial_surv_median_two_thrid=rep(0,num_dataset)
forest_surv_mean_two_thrid=rep(0,num_dataset)
forest_surv_median_two_thrid=rep(0,num_dataset)
super_haz_mean_two_thrid=rep(0,num_dataset)
super_haz_median_two_thrid=rep(0,num_dataset)
initial_haz_mean_two_thrid=rep(0,num_dataset)
initial_haz_median_two_thrid=rep(0,num_dataset)
forest_haz_mean_two_thrid=rep(0,num_dataset)
forest_haz_median_two_thrid=rep(0,num_dataset)

for (dataset in 1:num_dataset){
  try(
    {
    train = read_csv(paste0(path.cc,"data01/train",(num_dataset*(iter-1)+dataset),".csv"))
    test = read_csv(paste0(path.cc,"data01/test",(num_dataset*(iter-1)+dataset),".csv"))
    structure=TRUE
    model1=ranger(Surv(time, status) ~ .,data=train,num.tree=500,mtry=ncol(train)-2,splitrule = "logrank", min.node.size = 5)
    
    set.seed(300) # seed 109, num_leaf=4, min leaf size=13
    # seed 7305, min leaf size=11
    
    model2=ranger(Surv(time, status) ~ .,data=train,num.tree=1,mtry=ncol(train)-2,splitrule = "logrank", min.node.size = 5) # sample without rep
    # model2=rfsrc(Surv(time, status) ~ ., data=train, mtry=ncol(train)-2,sampsize=nrow(train),nsplit=0,ntree=1) # sample without rep
    
    predict2=predict(model2, data=train)
    chf_time=predict2$unique.death.times
    chf_value=predict2$chf
    
    leaf <- list() #create a list for the first four CHF curves
    leaf[[1]] <- chf_value[1,]
    leaf_allocation <- list()
    leaf_allocation[[1]] <- c(1)
    group=rep(NA,nrow(train)) #create a dictionary to reflect which group each observation is in
    group[1]=1
    
    #categorize the initial chf function for each leaf
    for (i in 2:nrow(train)) {
      found <- 0
      for (j in 1:length(leaf)) {
        if (all(chf_value[i,] == leaf[[j]])) {
          leaf_allocation[[j]] <- c(leaf_allocation[[j]], i)
          found <- 1
          group[i]=j
        }
      }
      if (found == 0) {
        leaf[[length(leaf) + 1]] <- chf_value[i,]
        leaf_allocation[[length(leaf_allocation) + 1]] <- i
        group[i]=length(leaf_allocation)
      }
    }
    
    # find the leaf with the least observation
    # print(min(table(group)))
    
    # treeInfo(model2,1)
    leaf_size=as.vector(table(group))
    single_leaf=as.vector(which(leaf_size<=2))
    normal_leaf=as.vector(which(leaf_size>2))
    
    reassign_dict=as.vector(1:length(leaf))
    for (i in 1:length(single_leaf)){
      new_assignment=sample(normal_leaf,1)
      group[which(group==single_leaf[i])]=new_assignment
      reassign_dict[single_leaf[i]]=new_assignment
    }
    
    # print(min(table(group)))
    # print(table(group))
    num_leaf=length(table(group))
    loose_group=as.numeric(names(table(group)))
    for (i in 1:length(group)){
      group[i]=which(loose_group==group[i])
    }
    for (i in 1:length(leaf)){
      reassign_dict[i]=which(loose_group==reassign_dict[i])
    }
    # print(table(group))
    
    num_leaf=length(table(group))
    # print(leaf_allocation)
    # print(group)
    
    # keep the first group allocation for phase diagram
    first_group=group
    
    #create overall lists of chf function x and y values
    chf_value_dict=list(leaf) # a list of lists of chf values
    time=list()
    for (i in 1:num_leaf){
      time[[i]]=chf_time
    }
    time_dict=list(time) # a list of lists of time points
    
    if (min(table(group))>=2){
      fold=10
      sampleframe = rep(1:fold, ceiling( nrow(train)/fold ) )
      CV_index=sample(sampleframe,nrow(train) ,  replace=FALSE )
      # CV_index=rep(1:10,nrow(train)/10)
      CV_losses=list()
    }else{
      fold=nrow(train)
      CV_index=1:nrow(train)
      CV_losses=list()
    }
    step_count=nrow(train)/fold
    
    row_time=rep(NA,nrow(train))
    row_event=rep(NA,nrow(train))
    risk_score=list()
    for (i in 1:num_leaf){
      risk_score[[i]]=rep(NA,nrow(train))
    }
    
    
    for (step in 1:fold){
      CV_test=train[CV_index==step,]
      CV_train=train[CV_index!=step,]
      group=first_group[CV_index!=step] # if 10-fold, drop 10 elements in group
      if (length(table(group))!=num_leaf){
        structure=FALSE
        break
      }
      
      predict2=predict(model2, data=CV_train)
      # chf_time=predict2$unique.death.times
      # chf_value=predict2$chf
      
      train_set=cbind(CV_train,group)
      # Perform the log-rank test
      
      result=tryCatch(
        expr={
          pairwise_survdiff(Surv(time, status)~group, data=train_set)
        },
        error = function(e) {
          #what should be done in case of exception?
          
          message("There was an error message.") # prints structure of exception
          return(list(NULL))
        }
      )
      
      if (length(result)==1){
        structure=FALSE
        break
      }
      result=pairwise_survdiff(Surv(time, status)~group, data=train_set)
      
      pvalue=result$p.value
      index=which(pvalue == max(pvalue,na.rm=T), arr.ind = TRUE)
      row=as.integer(rownames(pvalue)[index[1,1]])
      col=as.integer(colnames(pvalue)[index[1,2]])
      
      result_combined_2=survfit(Surv(time,status)~1, type = "fleming-harrington",data=subset(train_set,group%in%c(row,col)))
      
      chf_value_combined_2=result_combined_2$cumhaz
      chf_time_combined_2=result_combined_2$time
      
      dict=as.vector(1:num_leaf)
      group_dict=list(dict)
      
      if (num_leaf>1){
        for (i in 2:(num_leaf)){
          #record the current chf function value and time
          chf=chf_value_dict[[i-1]]
          time=time_dict[[i-1]]
          for (j in 1:num_leaf){
            if (group_dict[[i-1]][j]%in%c(row,col)){
              chf[[j]]=result_combined_2$cumhaz
              time[[j]]=result_combined_2$time
            }
          }
          
          chf_value_dict[[i]]=chf
          time_dict[[i]]=time
          
          #update the group assignment and decide the next combination
          
          group=ifelse(group%in%c(row,col), min(row,col), group)
          #print(paste0("Combine group ",row," and ",col," into group ",min(row,col)))
          
          for (ii in 1:num_leaf){
            if (dict[ii]==max(row,col)){
              dict[ii]=min(row,col)
            }
          }
          group_dict[[i]]=dict
          
          if (i != num_leaf){
            train_2=cbind(CV_train,group)
            result=pairwise_survdiff(Surv(time, status)~group, data=train_2)
            pvalue=result$p.value
            index=which(pvalue == max(pvalue,na.rm=T), arr.ind = TRUE)
            row=as.integer(rownames(pvalue)[index[1,1]])
            col=as.integer(colnames(pvalue)[index[1,2]])
            
            result_combined_2=survfit(Surv(time,status)~1, type = "fleming-harrington",data=subset(train_2,group%in%c(row,col)))
          }
        }
      } # end of fitting the trees
      # print(nrow(subset(train_2,group%in%c(row,col))))
      # print(row)
      # print(col)
      # print("Next")
      
      test_predict=predict(model2, data=CV_test)
      chf_predict=test_predict$chf
      test_group=rep(NA,nrow(CV_test))
      
      for (i in 1:nrow(CV_test)) {
        for (j in 1:length(leaf)) {
          if (all(chf_predict[i,] == leaf[[j]])) {
            test_group[i]=reassign_dict[j]
          }
        }
      }
      
      for (j in 1:nrow(CV_test)){
        row_time[(step-1)*step_count+j]=CV_test$time[j]
        row_event[(step-1)*step_count+j]=CV_test$status[j]
        for (i in 1:num_leaf){
          value=chf_value_dict[[i]][[test_group[j]]]
          time=time_dict[[i]][[test_group[j]]]
          
          index=length(time) # find the index of the largest time value smaller than t
          while(tau<time[index] && index>1){
            index=index-1
          }
          risk_score[[i]][(step-1)*step_count+j]=1-exp(-value[index])
          
        }
        
      }
      
    }
    
    # cv_c_index=rep(0,num_leaf)
    # for (i in 1:num_leaf){
    #   cv_c_index[i]=cIndex(row_time,row_event,as.vector(risk_score[[i]]))[1]
    # }
    # # cv_c_index
    
    if (structure==FALSE){
      next
    }
    
    
    fn=function(x){
      result=0
      for (i in 1:num_leaf){
        result=result+x[i]*risk_score[[i]]
      }
      fn=-c_index(row_time,row_event,as.vector(result),beta)[1]
      
      fn
    }
    
    heq=function(x){
      h=rep(0,1)
      for (i in 1:num_leaf){
        h[1]=h[1]+x[i]
      }
      h[1]=h[1]-1
      h
    }
    
    hin=function(x){
      h=rep(NA,1)
      for (i in 1:num_leaf){
        h[i]=x[i]
      }
      h
    }
    
    set.seed(1111)
    p0=runif(num_leaf)
    ans=constrOptim.nl(par=p0, fn=fn, heq=heq, hin=hin)
    ans$par
    
    parameter=ans$par/sum(ans$par)
    
    
    group=first_group # if 10-fold, drop 10 elements in group
    
    predict2=predict(model2, data=train)
    chf_time=predict2$unique.death.times
    chf_value=predict2$chf
    
    chf_value_dict=list(leaf) # a list of lists of chf values
    time=list()
    for (i in 1:num_leaf){
      time[[i]]=chf_time
    }
    time_dict=list(time) # a list of lists of time points
    
    
    test_risk_score=list()
    for (i in 1:num_leaf){
      test_risk_score[[i]]=rep(NA,nrow(test))
    }
    
    train_set=cbind(train,group)
    
    # Perform the log-rank test
    result <- pairwise_survdiff(Surv(time, status)~group, data=train_set)
    
    pvalue=result$p.value
    index=which(pvalue == max(pvalue,na.rm=T), arr.ind = TRUE)
    row=as.integer(rownames(pvalue)[index[1,1]])
    col=as.integer(colnames(pvalue)[index[1,2]])
    
    result_combined_2=survfit(Surv(time,status)~1, type = "fleming-harrington",data=subset(train_set,group%in%c(row,col)))
    
    chf_value_combined_2=result_combined_2$cumhaz
    chf_time_combined_2=result_combined_2$time
    
    dict=as.vector(1:num_leaf)
    group_dict=list(dict)
    
    if (num_leaf>1){
      for (i in 2:(num_leaf)){
        #record the current chf function value and time
        chf=chf_value_dict[[i-1]]
        time=time_dict[[i-1]]
        for (j in 1:num_leaf){
          if (group_dict[[i-1]][j]%in%c(row,col)){
            chf[[j]]=result_combined_2$cumhaz
            time[[j]]=result_combined_2$time
          }
        }
        
        chf_value_dict[[i]]=chf
        time_dict[[i]]=time
        
        #update the group assignment and decide the next combination
        
        group=ifelse(group%in%c(row,col), min(row,col), group)
        #print(paste0("Combine group ",row," and ",col," into group ",min(row,col)))
        
        for (ii in 1:num_leaf){
          if (dict[ii]==max(row,col)){
            dict[ii]=min(row,col)
          }
        }
        group_dict[[i]]=dict
        
        if (i != num_leaf){
          train_2=cbind(train,group)
          result=pairwise_survdiff(Surv(time, status)~group, data=train_2)
          pvalue=result$p.value
          index=which(pvalue == max(pvalue,na.rm=T), arr.ind = TRUE)
          row=as.integer(rownames(pvalue)[index[1,1]])
          col=as.integer(colnames(pvalue)[index[1,2]])
          
          result_combined_2=survfit(Surv(time,status)~1, type = "fleming-harrington",data=subset(train_2,group%in%c(row,col)))
        }
      }
    } # end of fitting the trees
    
    # print(nrow(subset(train_2,group%in%c(row,col))))
    # print(row)
    # print(col)
    # print("Next")
    
    test_predict=predict(model2, data=test)
    chf_predict=test_predict$chf
    test_group=rep(NA,nrow(test))
    
    for (i in 1:nrow(test)) {
      for (j in 1:length(leaf)) {
        if (all(chf_predict[i,] == leaf[[j]])) {
          test_group[i]=reassign_dict[j]
        }
      }
    }
    
    for (j in 1:nrow(test)){
      for (i in 1:num_leaf){
        value=chf_value_dict[[i]][[test_group[j]]]
        time=time_dict[[i]][[test_group[j]]]
        
        index=length(time) # find the index of the largest time value smaller than t
        while(tau<time[index] && index>1){
          index=index-1
        }
        test_risk_score[[i]][j]=1-exp(-value[index])
        
      }
      
    }
    
    test_c_index=rep(0,num_leaf)
    for (i in 1:num_leaf){
      test_c_index[i]=cIndex(test$time,test$status,as.vector(test_risk_score[[i]]))[1]
    }
    
    #test_c_index
    
    # plot(cv_c_index,test_c_index)
    # text(cv_c_index,test_c_index, labels = 1:length(cv_c_index), pos = 4)
    
    
    result=0
    for (i in 1:num_leaf){
      result=result+parameter[i]*test_risk_score[[i]]
    }
    super_loss[dataset]=cIndex(test$time,test$status,as.vector(result))[1]
    initial_loss[dataset]=test_c_index[1]
    
    test_predict=predict(model1, data=test)
    chf_predict=test_predict$chf
    chf_predict_time=test_predict$unique.death.times
    test_risk_score2=rep(0,nrow(test))
    
    for (j in 1:nrow(test)){
      index=length(chf_predict_time) # find the index of the largest time value smaller than t
      while(tau<chf_predict_time[index] && index>1){
        index=index-1
      }
      test_risk_score2[j]=1-exp(-chf_predict[j,index])
    }
    
    forest_loss[dataset]=cIndex(test$time,test$status,as.vector(test_risk_score2))[1]
    
    
    X=test
   
    true.lin.pred <- as.data.frame(10*(0.3*X[,1] -0.4*X[,2]+0.3*X[,4] -0.1*X[,6]-0.4*X[,7]+0.25*X[,3]+0.15*X[,8]+0.3*X[,9]-0.2*X[,10]))
    
    survival.true <- function(t){
      survival =  exp(-exp(true.lin.pred)*t^gamma.weibull)
      return(survival)
    }
    
    hazard.true <- function(t){
      hazard =  exp(true.lin.pred)*t^gamma.weibull
      return(hazard)
    }
    
  
    t=tau/3*1
    true=survival.true(t)
    hazard=hazard.true(t)
    for (j in 1:nrow(test)){
      for (i in 1:num_leaf){
        value=chf_value_dict[[i]][[test_group[j]]]
        time=time_dict[[i]][[test_group[j]]]
  
        index=length(time) # find the index of the largest time value smaller than t
        while(t<time[index] && index>1){
          index=index-1
        }
        test_risk_score[[i]][j]=1-exp(-value[index])
  
      }
  
    }
  
    predicted=1-test_risk_score[[1]]
    initial_surv_mean_one_thrid[dataset]=mean(abs(true-predicted)$X1)
    initial_surv_median_one_thrid[dataset]=median(abs(true-predicted)$X1)
    
    initial_haz_median_one_thrid[dataset]=median(abs(hazard+log(predicted))$X1)
    initial_haz_mean_one_thrid[dataset]=mean(abs(hazard+log(predicted))$X1)
    
  
    for (j in 1:nrow(test)){
      index=length(chf_predict_time) # find the index of the largest time value smaller than t
      while(t<chf_predict_time[index] && index>1){
        index=index-1
      }
      test_risk_score2[j]=1-exp(-chf_predict[j,index])
    }
  
    predicted=1-test_risk_score2
    forest_surv_mean_one_thrid[dataset]=mean(abs(true-predicted)$X1)
    forest_surv_median_one_thrid[dataset]=median(abs(true-predicted)$X1)
    
    forest_haz_median_one_thrid[dataset]=median(abs(hazard+log(predicted))$X1)
    forest_haz_mean_one_thrid[dataset]=mean(abs(hazard+log(predicted))$X1)
  
    result=0
    for (i in 1:num_leaf){
      result=result+parameter[i]*test_risk_score[[i]]
    }
    predicted=1-as.vector(result)
    super_surv_mean_one_thrid[dataset]=mean(abs(true-predicted)$X1)
    super_surv_median_one_thrid[dataset]=median(abs(true-predicted)$X1)
    
    super_haz_median_one_thrid[dataset]=median(abs(hazard+log(predicted))$X1)
    super_haz_mean_one_thrid[dataset]=mean(abs(hazard+log(predicted))$X1)
    
    
    t=tau/3*2
    true=survival.true(t)
    hazard=hazard.true(t)
    for (j in 1:nrow(test)){
      for (i in 1:num_leaf){
        value=chf_value_dict[[i]][[test_group[j]]]
        time=time_dict[[i]][[test_group[j]]]
        
        index=length(time) # find the index of the largest time value smaller than t
        while(t<time[index] && index>1){
          index=index-1
        }
        test_risk_score[[i]][j]=1-exp(-value[index])
        
      }
      
    }
    
    predicted=1-test_risk_score[[1]]
    initial_surv_mean_two_thrid[dataset]=mean(abs(true-predicted)$X1)
    initial_surv_median_two_thrid[dataset]=median(abs(true-predicted)$X1)
    
    initial_haz_median_two_thrid[dataset]=median(abs(hazard+log(predicted))$X1)
    initial_haz_mean_two_thrid[dataset]=mean(abs(hazard+log(predicted))$X1)
    
    
    for (j in 1:nrow(test)){
      index=length(chf_predict_time) # find the index of the largest time value smaller than t
      while(t<chf_predict_time[index] && index>1){
        index=index-1
      }
      test_risk_score2[j]=1-exp(-chf_predict[j,index])
    }
    
    predicted=1-test_risk_score2
    forest_surv_mean_two_thrid[dataset]=mean(abs(true-predicted)$X1)
    forest_surv_median_two_thrid[dataset]=median(abs(true-predicted)$X1)
    
    forest_haz_median_two_thrid[dataset]=median(abs(hazard+log(predicted))$X1)
    forest_haz_mean_two_thrid[dataset]=mean(abs(hazard+log(predicted))$X1)
    
    result=0
    for (i in 1:num_leaf){
      result=result+parameter[i]*test_risk_score[[i]]
    }
    predicted=1-as.vector(result)
    super_surv_mean_two_thrid[dataset]=mean(abs(true-predicted)$X1)
    super_surv_median_two_thrid[dataset]=median(abs(true-predicted)$X1)
    
    super_haz_median_two_thrid[dataset]=median(abs(hazard+log(predicted))$X1)
    super_haz_mean_two_thrid[dataset]=mean(abs(hazard+log(predicted))$X1)
    
    initial=0
    super=0
    forest=0

    for (index in 1:20){
      t=tau/20*index
      true=survival.true(t)
    for (j in 1:nrow(test)){
      for (i in 1:num_leaf){
        value=chf_value_dict[[i]][[test_group[j]]]
        time=time_dict[[i]][[test_group[j]]]

        index=length(time) # find the index of the largest time value smaller than t
        while(t<time[index] && index>1){
          index=index-1
        }
        test_risk_score[[i]][j]=1-exp(-value[index])

      }

    }

    predicted=1-test_risk_score[[1]]
    initial=initial+sum((true-predicted)^2)/nrow(test)/20
    for (j in 1:nrow(test)){
      index=length(chf_predict_time) # find the index of the largest time value smaller than t
      while(t<chf_predict_time[index] && index>1){
        index=index-1
      }
      test_risk_score2[j]=1-exp(-chf_predict[j,index])
    }

    predicted=1-test_risk_score2
    forest=forest+sum((true-predicted)^2)/nrow(test)/20

    result=0
    for (i in 1:num_leaf){
      result=result+parameter[i]*test_risk_score[[i]]
    }
    predicted=1-as.vector(result)
    super=super+sum((true-predicted)^2)/nrow(test)/20
    }

    initial_ibs[dataset]=initial
    super_ibs[dataset]=super
    forest_ibs[dataset]=forest
    print(dataset)
  }
  )
}

summaries=data.frame(forest_loss,initial_loss,super_loss,
                      forest_ibs,initial_ibs,super_ibs,
                   super_surv_mean_one_thrid,initial_surv_mean_one_thrid,forest_surv_mean_one_thrid,
                   super_surv_median_one_thrid,initial_surv_median_one_thrid,forest_surv_median_one_thrid,
                   super_haz_mean_one_thrid,initial_haz_mean_one_thrid,forest_haz_mean_one_thrid,
                   super_haz_median_one_thrid,initial_haz_median_one_thrid,forest_haz_median_one_thrid,
                   super_surv_mean_two_thrid,initial_surv_mean_two_thrid,forest_surv_mean_two_thrid,
                   super_surv_median_two_thrid,initial_surv_median_two_thrid,forest_surv_median_two_thrid,
                   super_haz_mean_two_thrid,initial_haz_mean_two_thrid,forest_haz_mean_two_thrid,
                   super_haz_median_two_thrid,initial_haz_median_two_thrid,forest_haz_median_two_thrid
                      )

summary_no_zero=summaries[apply(summaries!=0, 1, all),]
# Q1=apply(summary_no_zero, 2, quantile, probs = 0.25,type = 1)
# Q2=apply(summary_no_zero, 2, quantile, probs = 0.5,type = 1)
# Q3=apply(summary_no_zero, 2, quantile, probs = 0.75,type = 1)
# summaries=rbind(Q1,Q2,Q3,summary_no_zero)

write.csv(summary_no_zero,file=paste0(path.cc, 'results/summary0',iter,'.csv'), row.names = F)


# write.csv(as.matrix(forest_loss),file=paste0(path.cc, 'results/forest_loss.csv'), row.names = F)
# write.csv(as.matrix(initial_loss),file=paste0(path.cc, 'results/single_loss.csv'), row.names = F)
# write.csv(as.matrix(super_loss),file=paste0(path.cc, 'results/super_loss.csv'), row.names = F)
# 
# a=cbind(super_abs_haz_one_thrid,initial_abs_haz_one_thrid,forest_abs_haz_one_thrid,super_abs_haz_two_thrid,initial_abs_haz_two_thrid,forest_abs_haz_two_thrid)
# boxplot(a)

