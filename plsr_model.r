
library(tidyverse)
library(magrittr)
library(ggpubr)

library(pls)
library(caTools)
library(vip)


# This is a function to perform a single PLSR iteration on a dataset.
# The "dataset" variable should be a dataframe or tibble containing the Mordred descriptors
# and the dependent variable, with molecule names as rownames.
# "Validation" sets the validation to be used in the model, which in the case of 
# this paper is "LOOCV".
# "Trainsplit" is the proportion to be put into the training set, i.e. 0.8 or 0.9.
# "DepVarName" is the name of the dependent variable in the "dataset" dataframe.
# "ncomp" is the number of components to be used in the model. In this paper, a PLS
# model was ran beforehand in order to ascertain the optimal number of components.
# It can also be set to "min" to find the number of components with the lowest

PLSRSplit_v2 <- function(dataset, Validation, TrainSplit, DepVarName,
                         ncomp = "min")
{
  # Take the dataset, change the dependent variable name for ease of manipulation
  dataset <- dataset %>%
    rename("DepVarName" = DepVarName) %>%
    rownames_to_column("molecule")%>% # Take molecule names out of rownames and add as column
    rownames_to_column("numbID") # Take number location of each compound. This is used to identify which compounds are in the test and train sets.
  
  TestSet <- dataset %>%
    sample_frac((1-TrainSplit)) # Take samples of the dataset to be used as the test set
  
  # Extract compound list
  TestSplitCmpd <- TestSet$molecule
  
  # Exclude the test set from the trainset
  TrainSet <- dataset %>%
    anti_join(TestSet)
  # Extract compound list
  TrainSplitCmpd <- TrainSet$molecule

  # Make some dataframes containing just the dependent variable and the molecule identifiers
  TestSet_depvar <- TestSet %>% select(DepVarName, molecule, numbID)
  TrainSet_depvar <- TrainSet %>% select(DepVarName, molecule, numbID)
  TrainSet_depvar <- TrainSet %>%
    select(DepVarName)
  TestSet_nodepvar <- TestSet %>%
    column_to_rownames("molecule")%>%
    select(-c(numbID, DepVarName))
  
  # Remove identifying columns for the purpose of actually running the model
  TrainSet1 <- TrainSet %>%
    select(-c(molecule, numbID))
  
  # Run PLSR on the train set
  PLSR_Train <- plsr(DepVarName~., data = TrainSet1, validation = Validation, ncomp = ncomp)
  
  # Extract the predicted values for the train set as a tibble
  PLSR_predicted_train <- as_tibble(PLSR_Train[["fitted.values"]])
  
  # Take the column that uses the number of components specified
  PLSR_predicted_train <- PLSR_predicted_train %>%
    select(ncomp) %>%
    rename(Y_train_pred= 1)
  
  # Extract RMSEP from trainset
  TrainsetRMSE <- pls::RMSEP(PLSR_Train, estimate =  "train")
  TrainsetRMSE <- TrainsetRMSE[["val"]]
  # ncomp+1 because the first column is the intercept
  TrainsetRMSE <- TrainsetRMSE[ncomp+1]
  
  # calculate R2 for the trainset
  r_squared_train <- 1 - sum((TrainSet_depvar$DepVarName - PLSR_predicted_train$Y_train_pred)^2) / sum((TrainSet_depvar$DepVarName - mean(TrainSet_depvar$DepVarName))^2)
  
  # Test the model on the other data and repeat the previous extractions
  PLSR_Test <- predict(PLSR_Train, TestSet_nodepvar, ncomp = ncomp, scale = TRUE)
  PLSR_testset_df <- as.tibble(PLSR_Test)
  TestsetRMSE <- pls::RMSEP(PLSR_Train, newdata=TestSet, estimate =  "test")
  TestsetRMSE <- TestsetRMSE[["val"]]
  
  TestsetRMSE <- TestsetRMSE[ncomp+1]
  
  
  # Calculate Q squared for the test set
  ss_tot <- sum((TestSet_depvar$DepVarName - mean(TestSet_depvar$DepVarName))^2)
  
  ss_res <- sum((TestSet_depvar$DepVarName - PLSR_Test)^2)
  
  q_squared <- 1 - ss_res / ss_tot
  

  # Extract the predicted values and bring back the compound names from earlier (lost in the modelling function)
    PLSR_testset_df <- as.tibble(PLSR_Test)
    rownames(PLSR_testset_df) <- c(TestSplitCmpd)
    rownames(PLSR_predicted_train) <- c(TrainSplitCmpd)
    
    # Now we need dataframes of measured and predicted variables for both test and train set
    TrainSet_depvar <- TrainSet %>%
      select(DepVarName, molecule) %>%
      rename("Measured" = DepVarName)
    
    
    PLSR_predicted_train <- PLSR_predicted_train %>%
      rownames_to_column("molecule") %>%
      full_join(TrainSet_depvar)%>%
      rename("Predicted" = Y_train_pred) %>%
      mutate("TrainOrTest" = "Train")%>%
      mutate("RMSE" = TrainsetRMSE)
    
    # Return a dataframe with Q2, R2, predicted value, actual value, test/train flag
    
    PLSR_testset_df <- PLSR_testset_df %>%
      rownames_to_column("molecule")%>%
      full_join(TestSet_depvar) %>%
      mutate("Q2" = q_squared)%>%
      mutate("R2" = r_squared_train) %>%
      rename("Predicted" = 2) %>%
      rename("Measured" = 3) %>%
      mutate("TrainOrTest" = "Test")%>%
      mutate("RMSE" = TestsetRMSE) %>%
      full_join(PLSR_predicted_train)
    

  return(PLSR_testset_df)

}



# This is a wrapper function to perform many repeats of the previous function
# The input is the same as before, but with the addition of "Repeats", which is
# an integer indicating how any repeats are desired.


PLSRSplit_repeats_v2 <- function(dataset, Validation, TrainSplit, DepVarName,
                                 ncomp = "min", Repeats)
{
  # Take the input number and assign to a new variable, because sometimes
  # the input argument doesn't work as-is.
  repeats <- Repeats
  
  # Perform the first iteration outside of a for loop, to populate a dataframe.
  # This can probably be done quicker by creating an empty dataframe with the appropriate
  # column names and nrows = Repeats.
  PLSRRepeat_1 <- PLSRSplit_v2(dataset = dataset,
                               Validation = Validation,
                               TrainSplit = TrainSplit, 
                               DepVarName =DepVarName,
                               ncomp = ncomp)
# Add a repeat ID column 
    Q2_df <- PLSRRepeat_1 %>%
    mutate("Repeat" = 1) %>%
    na.omit() %>%
    select(-1)%>%
    distinct()
  
    # Perform the other repeats in a for loop
  for(R in 2:repeats)
  {
    print(paste0("Repeat number ", R)) # Keep track of progress via print()
    PLSRRepeat_R <- PLSRSplit_v2(dataset = dataset,
                                 Validation = Validation,
                                 TrainSplit = TrainSplit, 
                                 DepVarName = DepVarName,
                                 ncomp = ncomp)
    # Same as before
    PLSRRepeat_R <- PLSRRepeat_R %>%
      mutate("Repeat" = R) %>%
      #   mutate(as.numeric(Repeat)) %>%
      select(-1) %>%
      na.omit()
    # Join to df outside of loop
    Q2_df <- Q2_df %>%
      full_join(PLSRRepeat_R) %>%
      na.omit() %>%
      distinct()
  }
    # return final dataframe
  return(Q2_df)
}  