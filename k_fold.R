partition.by.fold <- function(fold.count, holdout.fold, full.data) {
    # create a vector maping each row to a fold
    row.to.fold <- sample(cut(seq(nrow(full.data)), fold.count, label = F))
    train.data <- full.data[row.to.fold != holdout.fold, , drop = F]
    holdout.data <- full.data[row.to.fold == holdout.fold, , drop = F]
    
    list(train.data = train.data, holdout.data = holdout.data)
}

# full.data is a data.frame with the variable to predict named 'value'.  All
# other columns are assumed to be inputs.  this method returns the R^2 and the
# RMSE of the out of sample values
compute.kfold.performance <- function(seed, fold.count, holdout.fold, full.data, 
    model_func, prediction_func) {
    library(caret)
    
    set.seed(seed)
    d <- partition.by.fold(fold.count, holdout.fold, full.data)
    
    model <- model_func(d$train.data)
    predictions <- prediction_func(model, d$holdout.data)
    
    postResample(predictions, d$holdout.data$value)
}

# full.data is a data.frame with the variable to predict named 'value'.  All
# other columns are assumed to be inputs.  this method returns the R^2 and the
# RMSE of the out of sample values
compute.sampled.performance <- function(seed, holdout.fraction, n_folds, holdout.i, 
    full.data, model_func, prediction_func) {
    library(caret)
    
    set.seed(seed)
    # generate a sequence of possible seeds, and pick the one corresponding to which
    # pass we're executing
    seed_seq <- as.integer(runif(n_folds, 0, 1e+08))
    set.seed(seed_seq[holdout.i])
    
    all.idx <- seq(nrow(full.data))
    train.idx <- sample(all.idx, (1 - holdout.fraction) * nrow(full.data))
    validate.idx <- all.idx[-train.idx]
    
    model <- model_func(full.data[train.idx, ])
    predictions <- prediction_func(model, full.data[validate.idx, ])
    
    postResample(predictions, full.data$value[validate.idx])
}
