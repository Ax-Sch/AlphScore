# Precision Recall F1 section
calculateTrueFalsePosNeg<-function(cutoff, score){
  predicted_val=unlist(testSet[,score])>cutoff
  true_val=testSet$outcome
  correct_prediction=predicted_val==true_val
  TPFP_tibble<-tibble(
  TruePositives=sum(correct_prediction & predicted_val, na.rm=TRUE),
  FalsePositives=sum(!correct_prediction & predicted_val, na.rm=TRUE),
  FalseNegatives=sum(!correct_prediction & !predicted_val, na.rm=TRUE),
  TrueNegatives=sum(correct_prediction & !predicted_val, na.rm=TRUE),
  )
  return(TPFP_tibble)
}

Precision<-function(TPFP_tibble){
  Precision = TPFP_tibble$TruePositives / (TPFP_tibble$TruePositives + TPFP_tibble$FalsePositives)
  return(Precision)
}

Recall<-function(TPFP_tibble){
  Recall = TPFP_tibble$TruePositives / (TPFP_tibble$TruePositives + TPFP_tibble$FalseNegatives)
  return(Recall)
}

Spec<-function(TPFP_tibble){
  Specifity=TPFP_tibble$TrueNegatives/(TPFP_tibble$TrueNegatives + TPFP_tibble$FalsePositives)
  return(Specifity)
}

Fscore<-function(Precision, Recall){
  Fscore = (2*Precision*Recall) / sum(Precision, Recall)
  return(Fscore)
}

calculatePrecision<-function(cutoff, score){
  TPFP_tibble<-calculateTrueFalsePosNeg(cutoff, score)
  return(Precision(TPFP_tibble))
}

calculateRecall<-function(cutoff, score){
  TPFP_tibble<-calculateTrueFalsePosNeg(cutoff, score)
  return(Recall(TPFP_tibble))
}

calculateSpec<-function(cutoff, score){
  TPFP_tibble<-calculateTrueFalsePosNeg(cutoff, score)
  return(Spec(TPFP_tibble))
}

calculateFscore<-function(cutoff, score){
  Precisionval=calculatePrecision(cutoff, score)
  Recallval=calculateRecall(cutoff, score)
  Fscoreval<-Fscore(Precisionval,Recallval)
  return(Fscoreval)
}

generate_table<-function(scorename){
    score_vals<-unlist(testSet[,scorename])
  fscores=tibble(cutoff=seq(from=round(min(score_vals, na.rm=TRUE),2), 
                            to=max(score_vals, na.rm=TRUE), by=0.01))
  fscores$Precision<-sapply(fscores$cutoff,function(x) 
    {calculatePrecision(x, scorename)})
  fscores$Recall<-sapply(fscores$cutoff,function(x) 
  {calculateRecall(x, scorename)})
  fscores$Specifity<-sapply(fscores$cutoff,function(x) 
  {calculateSpec(x, scorename)})
  fscores$Fscore<-sapply(fscores$cutoff, function(x) 
  {calculateFscore(x, scorename)})
return(fscores)
}
