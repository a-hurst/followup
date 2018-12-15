#' Linear Contrasts and Polynomial Trends for ANOVAs
#' 
#' These functions perform linear contrast and polynomial trend analyses on 
#' ANOVAs and ANCOVA models produced by the \pkg{afex} package, providing 
#' mean-squared errors, F-statistics, and (partial/contrast) eta-squared 
#' effect sizes in a manner similar to SPSS. They serve as a wrapper around
#' the \code{emmeans::contrast} method from the \pkg{emmeans} package, which
#' only reports t.scores for contrasts and does not provide effect sizes or MSE.
#' Splitting the contrasts across levels of a factor is supported, as is using
#' the pooled error term from the omnibus model when splitting on a
#' between-subjects factor.
#' 
#' @usage 
#' getContrasts(mod, factor, contr, split = NULL, pooled = FALSE,
#'              adj = 'none')
#' 
#' getTrends(mod, factor, split = NULL, pooled = FALSE, adj = 'none')
#' 
#'
#' @param mod An \code{afex_aov} ANOVA model to use for contrast analysis.
#' @param factor \code{character} vector (of length 1) indicating the name of
#'   the factor to perform contrasts on in the provided ANOVA model.
#' @param contr \code{list} of contrasts to perform on the model.
#' @param split \code{character} vector indicating the name of the factor (if
#'    any) to split the analysis on. Defaults to \code{NULL} (no split).
#' @param pooled \code{logical} indicating whether to use the pooled
#'   error term and degrees of freedom from the unsplit model for increased
#'   statistical power when splitting a contrast. Defaults to \code{FALSE}.
#' @param adj \code{character} vector indicating which method of p-value
#'   adjustment for multiple comparisons should be used when generating the 
#'   output table (e.g. "tukey", "bonferroni"). Defaults to "none".
#'
#'
#' @return A \code{\link[emmeans]{summary_emm}} and \code{\link{data.frame}}
#'   containing the contrast or trend results.
#'
#' @author Austin Hurst
#'
#'   
#' @name getContrasts
#' @aliases getContrasts getTrends
#' @export getContrasts getTrends
#' @importFrom emmeans emmeans contrast
#' 
#' 
#' @encoding UTF-8
#'

getContrasts <- function(mod, factor, contr, split=NULL, pooled=FALSE, adj='none') {

    # Perform contrasts on ANOVA model
    
    contrasts <- contrast(emmeans(mod, factor, by = split), contr, adjust = adj)
    contrasts.table <- summary(contrasts)
    has.within <- length(attr(mod, 'within')) > 0
    factor.is.within <- has.within & factor %in% names(attr(mod, 'within'))

    # Calculate contrast sums of squares
    
    n.per.group <- nrow(mod$data$long)/nrow(unique(mod$data$long[factor]))
    ests <- summary(contrasts)$estimate
    SScontr <- ests^2 / rowSums(attr(contrasts, 'misc')$con.coef^2 / n.per.group)
    
    # Add MSE, F-ratios, and effect sizes for all contrasts to table
    
    pvals <- contrasts.table$p.value # remove p values and later put at end of table
    contrasts.table$p.value <- NULL
    contrasts.table$MSE <- SScontr/(contrasts.table$t.ratio^2)
    contrasts.table$'F' <- contrasts.table$t.ratio^2
    
    if ( (has.within & factor.is.within) | is.ANCOVA(mod) ) {
        SSE <- contrasts.table$MSE * contrasts.table$df
        contrasts.table$pes <- SScontr / (SScontr + SSE)
    } else if ( has.within & !factor.is.within ) {
        ANOVA.SSa <- summary(mod$Anova)$univariate.tests[factor, 'Sum Sq']
        contrasts.table$etaSq <- SScontr / ANOVA.SSa
    } else {
        ANOVA.SSa <- mod$Anova[factor, 'Sum Sq']
        contrasts.table$etaSq <- SScontr / ANOVA.SSa
    }
    
    contrasts.table$SE <- NULL
    contrasts.table$t.ratio <- NULL
    contrasts.table$p.value <- pvals
    
    if (!is.null(split)) {
        simple <- suppressWarnings(simpleEffects(mod, split, returns='afex_aov'))
        unpooled <- do.call('rbind', lapply(simple, function(s) {
            getContrasts(s, factor, contr, NULL, FALSE, adj)
        }))
        if ('pes' %in% names(unpooled)) {
            contrasts.table$pes <- unpooled$pes
        } else {
            contrasts.table$etaSq <- unpooled$etaSq
        }
        if (!pooled) {
            contrasts.table$df <- unpooled$df
            contrasts.table$MSE <- unpooled$MSE
            contrasts.table$'F' <- unpooled$'F'
            contrasts.table$p.value <- unpooled$p.value
        } else {
            unsplit <- getContrasts(mod, factor, contr, NULL, FALSE, adj)
            contrasts.table$MSE <- unsplit$MSE
        }
    }
    
    contrasts.table
}


getTrends <- function(mod, factor, split=NULL, pooled=FALSE, adj="none") {
    getContrasts(mod, factor, 'poly', split, pooled, adj)
}


is.ANCOVA <- function(mod) {
    any(unlist(lapply(attr(mod, 'between'), is.null)))
}
