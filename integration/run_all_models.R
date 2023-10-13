run_all_models <- function(x.train, y.train, x.test, y.test, classes){
  
  result_df <- rbind(
    log_reg_model(x.train, y.train, x.test, y.test, classes = classes),
    log_reg_model(x.train, y.train, x.test, y.test, classes = classes, regularize = 'l1'),
    log_reg_model(x.train, y.train, x.test, y.test, classes = classes, regularize = 'l2'),
    log_reg_model(x.train, y.train, x.test, y.test, classes = classes, regularize = 'el'),
    svm_model(x.train, y.train, x.test, y.test, classes = classes),
    svm_model(x.train, y.train, x.test, y.test, classes = classes, kernel = 'radial'),
    rf_model(x.train, y.train, x.test, y.test, classes = classes)
  )
  
  return (result_df)
}
