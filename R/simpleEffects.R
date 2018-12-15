#' Simple Effects Analysis for Factorial ANOVAs
#' 
#' This function performs simple effects analysis on ANOVA models, allowing
#' easy splitting of models by a given factor. Additionaly, this function
#' supports using the pooled error term and degrees of freedom from the 
#' unsplit model for increased statistical power. 
#' 
#' @usage 
#' simpleEffects(mod, split, pooled = FALSE, returns = 'summary',
#'   correction = 'GG')
#' 
#'
#' @param mod An \code{afex_aov} ANOVA model returned by \code{aov_car}, \code{aov_4},
#'   or \code{aov_ez}.
#' @param split A \code{character} vector (of length 1) indicating the name
#'   of the factor the in the ANOVA model to split the analysis on.
#' @param pooled \code{logical} indicating whether to use the pooled error term
#'   and degrees of freedom from the unsplit model for increased statistical
#'   power when splitting a model. Defaults to \code{FALSE}.
#' @param returns \code{character} vector indicating the format in which to 
#'   return the results. Can be 'summary', 'table', 'nice', or 'afex_aov'.
#'   'table' and 'nice' are concise and easy to read, but omit sphericity
#'   test information for repeated-measures designs. 'afex_aov' returns the
#'   full afex anova object. Defaults to 'summary'.
#' @param correction \code{character} vector specifying which method of
#'   sphericity correction to use for repeated-measures factors, can be 'GG'
#'   (Greenhouse-Geisser), 'HF' (Huynh-Feldt), or 'none'. Defaults to 'GG'.
#'
#'
#' @return A \code{\link{list}} of ANOVA tables. The format of the tables
#'   depends on the value of the \code{returns} parameter.
#'
#' @author Austin Hurst
#'
#'   
#' @name simpleEffects
#' @aliases simpleEffects
#' @export simpleEffects
#' @importFrom afex aov_ez
#' @importFrom car Anova
#' 
#' 
#' @encoding UTF-8
#'

simpleEffects <- function(mod, split, returns='summary', pooled=FALSE, correction='GG') {
    
    dat <- mod$data$long
    bt <- names(attr(mod, 'between'))
    wn <- names(attr(mod, 'within'))
    cov <- NULL
    
    if (!is.null(bt)) {
        for (var in bt) {
            if (is.null(var)) { cov <- c(cov, var) }
        }
        bt <- bt[!bt %in% c(split, cov)]
        if (length(bt) == 0) { bt <- NULL }
    }
    if (!is.null(wn)) {
        if (split %in% wn & pooled) {
            stop(
                paste0("Pooling error for repeated-measures factors is likely to ",
                    "yield strongly biased results and is not supported."),
                call. = FALSE
            )
        }
        wn <- wn[!wn %in% split]
        if (length(wn) == 0) { wn <- NULL }
    }
    
    if (!is.null(wn) & pooled & (correction != 'none' | returns == 'summary')) {
        warning(
            paste0("Using pooled error alters the results of sphericity tests and ",
            "corrections.\nYou may want to try without pooled error for comparison."),
            call. = FALSE
        )
    }
    
    levels <- unique(dat[, split])
    effs <- lapply(levels, function(level) {
        sub.dat <- dat[dat[, split] == level, ]
        sub.mod <- aov_ez(
            id = attr(mod, 'id'),
            dv = attr(mod, 'dv'),
            data = sub.dat,
            between = bt,
            within = wn,
            covariate = cov,
            observed = attr(mod$anova_table, 'observed'),
            type = attr(mod, 'type'),
            anova_table = list(correction = correction)
        )
        if (pooled) {
            sm <- sub.mod$lm # split model lm()
            if (is.null(wn)) {
                sub.mod$Anova <- Anova(sm, error = mod$lm, type = attr(mod, 'type'))
            } else {
                SSPE <- wcrossprod(
                    residuals(mod$lm), 
                    w = rep(1, nrow(model.matrix(sm)))
                )
                sub.mod$Anova <- Anova(
                    sm,
                    idata = sub.mod$data$idata,
                    idesign = as.formula(paste0('~', paste0(wn, collapse='*'))),
                    SSPE = SSPE,
                    error.df = mod$Anova$error.df,
                    type = attr(mod, 'type')
                )
            }
            es <- attr(mod$anova_table, 'es')
            old_table <- sub.mod$anova_table
            sub.mod$anova_table <- do.call('anova', c(
                object = list(sub.mod),
                observed = list(attr(mod$anova_table, 'observed')),
                list(correction = correction)
            ))
            if (es != 'none') {
                sub.mod$anova_table[, es] <- old_table[, es]
            }
        }
        if (returns == 'table') {
            sub.mod$anova_table
        } else if (returns == 'summary') {
            summary(sub.mod, multivariate=FALSE)
        } else if (returns == 'nice'){
            nice(sub.mod, correction=correction)
        } else if (returns == 'afex_aov') {
            sub.mod
        } else {
            stop(paste0('Unrecognized return type. Must be one of ',
                '"summary", "table", "nice", or "afex_aov".'), call. = FALSE)
        }
    })
    names(effs) <- as.character(levels)
    effs
}

