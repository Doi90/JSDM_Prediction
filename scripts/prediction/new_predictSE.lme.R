new_predictSE.lme <- function (mod, newdata, se.fit = TRUE, print.matrix = FALSE, 
          level = 0, ...) 
{
  if (!identical(level, 0)) 
    stop("\nThis function does not support computation of predicted values\n", 
         "or standard errors for higher levels of nesting\n")
  fixed <- formula(mod$terms)[-2]    ### Only line that has been changed
  tt <- terms(mod)
  TT <- delete.response(tt)
  newdata <- as.data.frame(newdata)
  mfArgs <- list(formula = fixed, data = newdata)
  dataMix <- do.call("model.frame", mfArgs)
  contr <- mod$contrasts
  for (i in names(dataMix)) {
    if (inherits(dataMix[, i], "factor") && !is.null(contr[[i]])) {
      levs <- levels(dataMix[, i])
      levsC <- dimnames(contr[[i]])[[1]]
      if (any(wch <- is.na(match(levs, levsC)))) {
        stop(paste("Levels", paste(levs[wch], collapse = ","), 
                   "not allowed for", i))
      }
      attr(dataMix[, i], "contrasts") <- contr[[i]][levs, 
                                                    , drop = FALSE]
    }
  }
  m <- model.frame(TT, data = dataMix)
  des.matrix <- model.matrix(TT, m)
  newdata <- des.matrix
  fix.coef <- fixef(mod)
  ncoefs <- length(fix.coef)
  names.coef <- labels(fix.coef)
  nvals <- dim(newdata)[1]
  int.yes <- any(names.coef == "(Intercept)")
  if (!int.yes) 
    stop("\nThis function does not work with models excluding the intercept terms\n")
  formula <- character(length = ncoefs)
  nbetas <- ncoefs - 1
  if (int.yes & nbetas >= 1) {
    formula <- paste("Beta", 1:nbetas, sep = "")
    formula <- c("Beta0", formula)
  }
  else {
    if (int.yes & nbetas == 0) {
      formula <- "Beta0"
    }
  }
  inters <- rep(NA, ncoefs)
  for (m in 1:ncoefs) {
    inters[m] <- attr(regexpr(pattern = ":", text = names.coef[m]), 
                      "match.length")
  }
  names.cov <- paste("cov", 1:ncoefs - 1, sep = "")
  if (!int.yes) {
    names.cov <- paste("cov", 1:ncoefs, sep = "")
  }
  id <- which(inters == 1)
  for (k in 1:length(id)) {
    names.cov[id[k]] <- paste("inter", k, sep = "")
  }
  formula2 <- character(length = ncoefs)
  for (b in 1:ncoefs) {
    formula2[b] <- paste(formula[b], names.cov[b], sep = "*")
  }
  if (int.yes) {
    formula2[1] <- "Beta0"
  }
  eq.space <- parse(text = as.expression(paste(formula2, collapse = "+")), 
                    srcfile = NULL)
  no.space <- gsub("[[:space:]]+", "", as.character(eq.space))
  equation <- parse(text = as.expression(no.space))
  if (identical(se.fit, TRUE)) {
    part.devs <- list()
    for (j in 1:ncoefs) {
      part.devs[[j]] <- D(equation, formula[j])
    }
  }
  ncovs <- ncoefs - length(id)
  cov.values <- list()
  if (int.yes && ncovs == 1) {
    cov.values[[1]] <- 1
  }
  if (int.yes && ncovs > 1) {
    cov.values[[1]] <- rep(1, nvals)
    for (q in 2:ncoefs) {
      cov.values[[q]] <- newdata[, labels(fix.coef)[q]]
    }
  }
  else {
    for (q in 1:ncoefs) {
      cov.values[[q]] <- newdata[, labels(fix.coef)[q]]
    }
  }
  names(cov.values) <- names.cov
  cov.values.mat <- matrix(data = unlist(cov.values), nrow = nvals, 
                           ncol = ncoefs)
  if (identical(se.fit, TRUE)) {
    predicted.SE <- matrix(NA, nrow = nvals, ncol = 2)
    colnames(predicted.SE) <- c("Pred.value", "SE")
    rownames(predicted.SE) <- 1:nvals
    part.devs.eval <- list()
    part.devs.eval[[1]] <- 1
    for (w in 1:nvals) {
      if (int.yes && ncovs > 1) {
        for (p in 2:ncoefs) {
          part.devs.eval[[p]] <- cov.values[names.cov[p]][[1]][w]
        }
      }
      part.devs.solved <- unlist(part.devs.eval)
      vcmat <- vcov(mod)
      mat_partialdevs <- as.matrix(part.devs.solved)
      mat_tpartialdevs <- t(part.devs.solved)
      var_hat <- mat_tpartialdevs %*% vcmat %*% mat_partialdevs
      SE <- sqrt(var_hat)
      predicted.vals <- fix.coef %*% cov.values.mat[w, 
                                                    ]
      predicted.SE[w, 1] <- predicted.vals
      predicted.SE[w, 2] <- SE
    }
    out.fit.SE <- list(fit = predicted.SE[, "Pred.value"], 
                       se.fit = predicted.SE[, "SE"])
  }
  else {
    predicted.SE <- matrix(NA, nrow = nvals, ncol = 1)
    colnames(predicted.SE) <- c("Pred.value")
    rownames(predicted.SE) <- 1:nvals
    for (w in 1:nvals) {
      predicted.vals <- fix.coef %*% cov.values.mat[w, 
                                                    ]
      predicted.SE[w, 1] <- predicted.vals
    }
    out.fit.SE <- predicted.SE
    colnames(out.fit.SE) <- "fit"
  }
  if (identical(print.matrix, TRUE)) {
    out.fit.SE <- predicted.SE
    if (identical(se.fit, TRUE)) {
      colnames(out.fit.SE) <- c("fit", "se.fit")
    }
    else {
      colnames(out.fit.SE) <- c("fit")
    }
  }
  return(out.fit.SE)
}
