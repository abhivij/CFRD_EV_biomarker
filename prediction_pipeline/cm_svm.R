# label.train <- output_labels.train
# label.test <- output_labels.test
# kernel <- "sigmoid"

svm_model <- function(data.train, label.train, data.test, label.test, 
                      classes, kernel = "sigmoid", random_seed = 1000,
                      ...){
  
  kernel_name <- paste(toupper(substring(kernel, 1, 1)), substring(kernel, 2), sep = "")
  model_name <- paste(kernel_name, "Kernel SVM")

  
  try({
    set.seed(random_seed)
    model <- e1071::svm(data.train, factor(label.train$Label, levels = classes), 
                        probability = TRUE, kernel = kernel)
    
    pred.train <- predict(model, data.train, probability = TRUE)
    pred_prob.train <- data.frame(attr(pred.train, 'probabilities'))[classes[2]]
    
    pred <- predict(model, data.test, probability = TRUE)
    pred_prob <- data.frame(attr(pred, 'probabilities'))[classes[2]]

    
    result_df.train <- data.frame("TrueLabel" = label.train$Label,
                                  "Pred_prob" = pred_prob.train[,1],
                                  "PredictedLabel" = pred.train,
                                  "Type" = "train")
    
    result_df.test <- data.frame("TrueLabel" = label.test$Label,
                                 "Pred_prob" = pred_prob[,1],
                                 "PredictedLabel" = pred,
                                 "Type" = "test")
    
    result_df <- rbind(result_df.train %>%
                         rownames_to_column("sample"), 
                       result_df.test %>%
                         rownames_to_column("sample"))
    
    result_df <- result_df %>%
      mutate(Pred_prob = as.double(Pred_prob))
    
    
    return (result_df)
    
  })
  
  
}
