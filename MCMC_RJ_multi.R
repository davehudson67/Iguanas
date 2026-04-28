## sampler for more than one variable
sampler_RJ_indicator_multi <- nimbleFunction(
    name = 'sampler_RJ_indicator_multi',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## note: target is the indicator variable,
        ## control$targetNodes are the variables conditionally in the model
        ## control list extraction
        coefNodes     <- control$targetNodes
        nCoefs        <- length(coefNodes)
        if(nCoefs <= 1) stop("Need at least two targets for RJ_indicator_multi sampler.")
        proposalScale <- control$scale
        proposalMean  <- control$mean
        ## node list generation
        calcNodes <- model$getDependencies(c(coefNodes, target))
        calcNodesReduced <- model$getDependencies(target)
    },
    run = function() {
        currentIndicator <- model[[target]]
        if(currentIndicator == 0) {   ## propose addition of coefNodes
            currentLogProb <- model$getLogProb(calcNodesReduced)
            proposalCoefs <- rnorm(nCoefs, proposalMean, proposalScale)
            logProbForwardProposal <- sum(dnorm(proposalCoefs, proposalMean, proposalScale, log = TRUE))
            model[[target]] <<- 1
            values(model, coefNodes) <<- proposalCoefs
            proposalLogProb <- model$calculate(calcNodes)
            logAcceptanceProb <- proposalLogProb - currentLogProb - logProbForwardProposal
        } else {                      ## propose removal of coefNodes
            currentLogProb <- model$getLogProb(calcNodes)
            currentCoefs <- values(model, coefNodes)
            logProbReverseProposal <- sum(dnorm(currentCoefs, proposalMean, proposalScale, log = TRUE))
            values(model, coefNodes) <<- rep(0, nCoefs)
            model[[target]] <<- 0
            model$calculate(calcNodes)
            logAcceptanceProb <- model$getLogProb(calcNodesReduced) - currentLogProb + logProbReverseProposal
        }
        accept <- decide(logAcceptanceProb)
        if(accept) { copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else     { copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE) }
    },
    methods = list(
        reset = function() { }
    )
)

## configure for multiple targets
configureRJ_multi <- function(conf, targetNodes, indicatorNodes = NULL, priorProb = NULL, control = list(mean = NULL, scale = NULL, fixedValue = NULL)) {
    model <- conf$model
    nNodes <- length(targetNodes)
    fixedValue <- if(!is.null(control$fixedValue)) control$fixedValue else 0
    mean       <- if(!is.null(control$mean))       control$mean       else 0
    scale      <- if(!is.null(control$scale))      control$scale      else 1
    ## repeat values for multiple nodes if only one value is provided
    if(length(fixedValue) != nNodes)
        if(length(fixedValue) == 1) fixedValue <- rep(fixedValue, nNodes) else stop("configureRJ: inconsistent length of 'fixedValue' argument and specified number of 'targetNodes'.")
    if(length(mean) != nNodes)
        if(length(mean) == 1) mean <- rep(mean, nNodes) else stop("configureRJ: inconsistent length of 'mean' argument and specified number of 'targetNodes'.")
    if(length(scale) != nNodes)
        if(length(scale) == 1) scale <- rep(scale, nNodes) else stop("configureRJ: inconsistent length of 'scale' argument and specified number of 'targetNodes'.")
    ## flag for indicators and prior
    indicatorFlag <- !is.null(indicatorNodes)
    priorFlag     <- !is.null(priorProb)
    ## must provide either indicatorNodes or priorProb
    if(indicatorFlag == priorFlag) stop("configureRJ: Provide 'indicatorNodes' or 'priorProb' vector")
    ## fixedValue can be used only with priorProb
    if(indicatorFlag && any(fixedValue != 0)) warning("configureRJ: 'fixedValue' can be provided only when using 'priorProb'; it will be ignored.")
    ##
    if(priorFlag) {    ## no indicator variables; use RJ_fixed_prior sampler
        ## check that priorProb values are in [0,1]
        if(any(priorProb < 0 | priorProb > 1)) stop("configureRJ: elements in priorProb must be probabilities in [0,1].")
        ## if one value for prior is given, it is used for each variable
        if(length(priorProb) != nNodes)
            if(length(priorProb) == 1) priorProb <- rep(priorProb, nNodes) else stop("configureRJ: Length of 'priorProb' vector must match 'targetNodes' length.")
        for(i in 1:nNodes) {
            nodeAsScalar <- model$expandNodeNames(targetNodes[i], returnScalarComponents = TRUE)
            ## if the node is multivariate throw an error
            if(any(model$isMultivariate(targetNodes[i])))
                stop(paste0("configureRJ: '", targetNodes[i], "' is multivariate; only univariate priors can be used with reversible jump sampling."))
            ## Create RJ control list for the node
            nodeControl <- list(priorProb = priorProb[i], mean = mean[i],
                                scale = scale[i], fixedValue = fixedValue[i])
            for(j in 1:length(nodeAsScalar)) {
                currentConf <- conf$getSamplers(nodeAsScalar[j])
                ## check on node configuration
                if(length(currentConf) == 0) {
                    warning(paste0("configureRJ: There are no samplers for '", nodeAsScalar[j],"'. Skipping it."))
                    next
                }
                else if(any(sapply(currentConf,'[[','name') == 'RJ_fixed_prior' |
                                sapply(currentConf,'[[','name') == 'RJ_indicator' |
                                    sapply(currentConf,'[[','name') == 'RJ_toggled')) {
                    stop(paste0("configureRJ: node '", nodeAsScalar[j],"' is already configured for reversible jump."))
                }
                else if(length(currentConf) > 1)
                    warning(paste0("configureRJ: There is more than one sampler for '", nodeAsScalar[j], "'. Only the first will be used by RJ_toggled sampler, and others will be removed."))
                ## substitute node sampler
                conf$removeSamplers(nodeAsScalar[j])
                conf$addSampler(type = sampler_RJ_fixed_prior,
                                    target = nodeAsScalar[j],
                                    control = nodeControl)
                conf$addSampler(type = sampler_RJ_toggled,
                                    target = nodeAsScalar[j],
                                    control = list(samplerType = currentConf[[1]], fixedValue = nodeControl$fixedValue))
            }
        }
    }
    ##
    if(indicatorFlag) {   ## indicator variables; use RJ_indicator sampler
        if(length(indicatorNodes) != nNodes)
            stop("configureRJ: Length of 'indicatorNodes' vector must match 'targetNodes' length.")
        for(i in 1:nNodes) {
            nodeAsScalar <- model$expandNodeNames(targetNodes[i], returnScalarComponents = TRUE)
            indicatorsAsScalar <- model$expandNodeNames(indicatorNodes[i], returnScalarComponents = TRUE)
            ## if the node is multivariate throw an error
            if(any(model$isMultivariate(targetNodes[i]))) 
                stop(paste0("configureRJ: '", targetNodes[i], "' is multivariate; only univariate nodes can be used with reversible jump sampling."))
            ## check that length of indicatorNodes matches targetNodes
            if(length(nodeAsScalar) != length(indicatorsAsScalar))
                stop(paste0("configureRJ: indicatorNodes node '", indicatorNodes[i] ,"' does not match '", targetNodes[i], "' size."))
            nodeControl  = list(mean = mean[i], scale = scale[i])
            for(j in 1:length(nodeAsScalar)) {
                nodeControl$targetNode <- nodeAsScalar[j]
                currentConf <- conf$getSamplers(nodeAsScalar[j])
                ## check on node configuration
                if(length(currentConf) == 0) {
                    warning(paste0("configureRJ: There are no samplers for '", nodeAsScalar[j],"'. Skipping it."))
                } else if(any(sapply(currentConf,'[[','name') == 'RJ_fixed_prior' |
                                  sapply(currentConf,'[[','name') == 'RJ_indicator' |
                                      sapply(currentConf,'[[','name') == 'RJ_toggled')) {
                    stop(paste0("configureRJ: Node '", nodeAsScalar[j],"' is already configured for reversible jump."))
                } else if(length(currentConf) > 1) {
                    warning(paste0("configureRJ: There is more than one sampler for '", nodeAsScalar[j], "'. Only the first will be used by RJ_toggled sampler, and others will be removed."))
                }
                ## check whether indicator nodes have been used before
                ## and if so then link indicators to multiple nodes
                currentConfInd <- conf$getSamplers(indicatorsAsScalar[j])
                if(length(currentConfInd) > 1) stop("Error with currentConfInd")
                if(any(sapply(currentConfInd, '[[', 'name') == 'RJ_indicator' | 
                        sapply(currentConfInd, '[[', 'name') == 'RJ_indicator_multi')) {
                    currentRJInd <- currentConfInd[sapply(currentConfInd, '[[', 'name') == 'RJ_indicator']
                    currentRJ <- currentConf[sapply(currentConfInd, '[[', 'name') == 'RJ_indicator']
                    if(length(currentRJInd) > 0) {
                        ## extract current coefficients
                        coefs <- sapply(lapply(currentRJInd, '[[', 'control'), '[[', 'targetNode')
                        ## append to new node
                        coefs <- c(coefs, nodeAsScalar[j])
                        nodeControl1 <- nodeControl
                        nodeControl1$targetNode <- NULL
                        nodeControl1$targetNodes <- coefs
                        ## Add multi reversible jump sampler for the indicatorNodes variable
                        conf$removeSamplers(indicatorsAsScalar[j])
                        conf$addSampler(type = sampler_RJ_indicator_multi,
                                            target = indicatorsAsScalar[j],
                                            control = nodeControl1)
                        ## Add sampler for the new coefficient variable (when is in the model)
                        conf$removeSamplers(nodeAsScalar[j])
                        conf$addSampler(type = sampler_RJ_toggled,
                                            target = nodeAsScalar[j],
                                            control = list(samplerType = currentRJ[[1]]))
                    }
                    currentRJInd <- currentConfInd[sapply(currentConfInd, '[[', 'name') == 'RJ_indicator_multi']
                    currentRJ <- currentConf[sapply(currentConfInd, '[[', 'name') == 'RJ_indicator_multi']
                    if(length(currentRJInd) > 0) {
                        ## extract current coefficients
                        coefs <- sapply(lapply(currentRJInd, '[[', 'control'), '[[', 'targetNodes')
                        ## append to new node
                        coefs <- c(coefs, nodeAsScalar[j])
                        nodeControl1 <- nodeControl
                        nodeControl1$targetNode <- NULL
                        nodeControl1$targetNodes <- coefs
                        ## Add multi reversible jump sampler for the indicatorNodes variable
                        conf$removeSamplers(indicatorsAsScalar[j])
                        conf$addSampler(type = sampler_RJ_indicator_multi,
                                            target = indicatorsAsScalar[j],
                                            control = nodeControl1)
                        ## Add sampler for the new coefficient variable (when is in the model)
                        conf$removeSamplers(nodeAsScalar[j])
                        conf$addSampler(type = sampler_RJ_toggled,
                                            target = nodeAsScalar[j],
                                            control = list(samplerType = currentRJ[[1]]))
                    }
                } else {
                    ## Add reversible jump sampler for the indicatorNodes variable
                    conf$removeSamplers(indicatorsAsScalar[j])
                    conf$addSampler(type = sampler_RJ_indicator,
                                        target = indicatorsAsScalar[j],
                                        control = nodeControl)
                    ## Add sampler for the coefficient variable (when is in the model)
                    conf$removeSamplers(nodeAsScalar[j])
                    conf$addSampler(type = sampler_RJ_toggled,
                                        target = nodeAsScalar[j],
                                        control = list(samplerType = currentConf[[1]]))
                }
            }
        }
    }
    return(invisible(NULL))
}

