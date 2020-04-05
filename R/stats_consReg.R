#' @export
#' @rdname logLik.ConsReg

logLik.ConsReg <- function(object, REML = F,...){

  if(object$family.name == 'gaussian'){
    if (inherits(object, "mlm"))
      stop("'logLik.lm' does not support multiple responses")
    res <- object$residuals
    p <- object$rank
    N <- length(res)
    w <- rep.int(1, N)
    N0 <- N
    if (REML)
      N <- N - p
    val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) +
                                       log(sum(w * res^2))))
    if (REML)
      val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
    attr(val, "nall") <- N0
    attr(val, "nobs") <- N
    attr(val, "df") <- p + 1
    class(val) <- "logLik"
    val
  }else{
    if (!missing(...))
      warning("extra arguments discarded")
    fam <- object$family.name
    p <- object$rank
    if (fam %in% c("gaussian", "Gamma", "inverse.gaussian"))
      p <- p + 1
    val <- object$metrics$LogLik
    attr(val, "nobs") <- sum(!is.na(object$residuals))
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
  }
}


#' @export
#' @rdname AIC.ConsReg

AIC.ConsReg <- function(object, ..., k = 2){

  ll = stats::logLik
  if (!missing(...)) {
    lls <- lapply(list(object, ...), ll)
    vals <- sapply(lls, function(el) {
      no <- attr(el, "nobs")
      c(as.numeric(el), attr(el, "df"), if (is.null(no)) NA_integer_ else no)
    })
    val <- data.frame(df = vals[2L, ], ll = vals[1L, ])
    nos <- stats::na.omit(vals[3L, ])
    if (length(nos) && any(nos != nos[1L]))
      warning("models are not all fitted to the same number of observations")
    val <- data.frame(df = val$df, AIC = -2 * val$ll + k *
                        val$df)
    Call <- match.call()
    Call$k <- NULL
    row.names(val) <- as.character(Call[-1L])
    val
  }
  lls <- ll(object)
  -2 * as.numeric(lls) + k * attr(lls, "df")

}

#' @export
#' @rdname AIC.ConsRegArima
AIC.ConsRegArima <- function(object, ..., k = 2){
  object$aic
}




