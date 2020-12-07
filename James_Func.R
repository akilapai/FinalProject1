#Part 2 function
James_Func<-function(DATA=NA, MOD=NA, SIZE=1){
  
  betas = NULL
  
  # iterate through unique protein numbers
  for (protein.number in unique(DATA[['ProteinNumber']])){
    
    # create a data subset for each protein
    protein.sub <- subset(DATA, DATA[['ProteinNumber']] == protein.number)
    
    # specifying the number of boostrapped resamples to generate
    for(draw in 1:SIZE) {
      
      # create a resample of the protein data subset, with the same size as the subset
      resample = protein.sub[sample(1:nrow(protein.sub),size=nrow(protein.sub),replace=T),]
      
      # fit a linear regression model of the expression data in the resample, using MOD as predictor(s)
      fit<-lm(paste0('Expression~', MOD), data=resample)
      
      # bind the model coefficients into a dataframe
      betas <- as.data.frame(rbind(betas, coef(fit)))
      
      # add the protien number to the coefficients
      betas$ProteinNumber[nrow(betas)]<-protein.number
    }
  }
  RES<-as.data.frame(betas)
  return(RES)
}

