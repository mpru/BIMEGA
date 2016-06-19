# Fix for NOTEs from R CMD check (no visible binding for global variable)
if(getRversion() >= "2.15.1")  {
    utils::globalVariables(c('SNPprobes', 'ProbeAnnotation', 'BatchData'))
}

#' BIMEGA: Bivariate Gaussian Mixture model for DNA methylation and gene expression data in cancer.
#' 
#' BIMEGA identifies DNA methylation driven genes by jointly modeling DNA methylation and gene expression data in cancer vs. normal and looking 
#' for homogeneous subpopulations. Matched gene expression data (e.g. from microarray technology or RNA sequencing) 
#' is also used to identify functional DNA methylation events by requiring a significant association between methylation 
#' and gene expression of a particular gene.
#' @param METcancer	Matrix with the methylation data of cancer tissue with genes in rows and samples in columns.
#' @param METnormal	Matrix with the normal methylation data of the same genes as in METcancer. Again genes in rows and samples in columns. The samples do not have to match with the cancer data.
#' @param MAcancer Gene expression data for cancer tissue with genes in rows and samples in columns.
#' @param MAnormal Gene expression data for normal tissue with genes in rows and samples in columns (optional data set).
#' @param listOfGenes Vector with genes names to be evaluated, names must coincide with the names of the rows of METcancer.
#' @param filter Logical indicating if the polynomial regression to select genes with significative relation between methylation and gene expression should be performed (default: TRUE).
#' @param NoNormalMode Logical indicating if the methylation states found in the cancer samples should be compared to the normal samples (default: TRUE).
#' @param OutputRoot Path to store the BIMEGA results object.
#' @return A list with the following components:
#' \item{MethylationStates}{Matrix with for all genes the Methylation states 
#'     using DM-value (i.e. Differential methylation values) that are defined as 
#'     the methylation value with respect to the average normal methylation for a gene.}
#' \item{NrComponents}{The number of methylation states for each gene.}
#' \item{Models}{Bivariate Gaussian mixture model parameters for each gene.}
#' \item{MethylationDrivers}{Genes identified as transcriptionally predictive and differentially methylated by BIMEGA.}
#' \item{MixtureStates}{A list with the DM-values for each gene that is functional and differential.}
#' \item{Classifications}{Matrix with integers indicating to which mixture component each sample was assigned to for each gene.}
#' \item{FunctionalGenesResults}{Matrix with information on the polynomial regression fit for each driver gene.}
#' @export
#' @examples 
#' # load the data sets needed for BIMEGA
#' data(METcancer)
#' data(METnormal)
#' data(MAcancer)
#' data(MAnormal)
#' 
#' # run BIMEGA on a small set of example data
#' BIMEGAresults <- BIMEGA(METcancer, METnormal, MAcancer, MAnormal)
#' 
BIMEGA <- function(METcancer, 
                      METnormal, 
                      MAcancer, 
                      MAnormal = NULL, 
                      listOfGenes = NULL, 
                      filter = TRUE,
                      NoNormalMode = FALSE, 
                      OutputRoot = "") { 
    
    ### Step 0: Prepare data
    
    # Keep only those samples with both METcancer and MAcancer data
    OverlapSamples <- intersect(colnames(METcancer), colnames(MAcancer))
    MAcancer <- MAcancer[, OverlapSamples, drop = FALSE]
    METcancer <- METcancer[, OverlapSamples, drop = FALSE]
    # Same for normal samples if MAnormal is provided
    if (!is.null(MAnormal)) {
        OverlapSamples <- intersect(colnames(METnormal), colnames(MAnormal))
        if (length(OverlapSamples) == 0) {
            cat("Normal methylation and normal gene expression data come from different samples.\n")
            normalOverlap <- FALSE
        } else {
            MAnormal <- MAnormal[, OverlapSamples, drop = FALSE]
            METnormal <- METnormal[, OverlapSamples, drop = FALSE]
            normalOverlap <- TRUE
        }
    }
    
    # Keep only the genes provided by user
    if (!is.null(listOfGenes)) {
        listOfGenes <- intersect(listOfGenes, rownames(METcancer))
        METcancer <- METcancer[listOfGenes, , drop = F]
    }
    
    ### Step 1: modeling the gene expression using cancer methylation data (beta values scale)
    if (filter) {
        FunctionalGenesResults <- modelGeneExpression(METcancer, MAcancer, CovariateData = NULL)
        FunctionalGenes <- rownames(FunctionalGenesResults)
    } else {
        # No regression filter, but keep genes present in both METcancer and MAcancer data, as the ones returned by the regression filter
        METcancer.split.names  <- sapply(strsplit(rownames(METcancer),  '---'), function(x) x[1])
        genes.to.keep.MET <- METcancer.split.names %in% rownames(MAcancer)
        FunctionalGenes <- rownames(METcancer)[genes.to.keep.MET]
    }
 
    ### Step 2: modeling the methylation data as a mixture of beta/normal distributions 
    
    # Transform beta values into M values
    METcancer <- getMvalues(METcancer)
    METnormal <- getMvalues(METnormal)
    if (length(FunctionalGenes) > 0) {
        MixtureModelResults <- MixtureModel_Biv(METcancer, METnormal, MAcancer, MAnormal, FunctionalGenes, NoNormalMode)
        if (filter) {
            MixtureModelResults$FunctionalGenesResults <- FunctionalGenesResults[intersect(FunctionalGenes, MixtureModelResults$MethylationDrivers), ]
        }
    } else {
        cat("No transcriptionally predictive genes found.\n")
        return(NULL)
    }
    ### Step 3: write output to file
    if (OutputRoot != "") {
        saveRDS(MixtureModelResults, file = paste0(OutputRoot, "BIMEGA_Results.rds"))
    }
    return(MixtureModelResults)
}

#' The getMvalues function
#' 
#' Internal. Converts Beta values into M values
#' @param matr matrix of beta values.
#' @return matrix of M values
#' @keywords internal
#'
getMvalues <- function(matr) {   
    # Handle values out of range
    for (i in 1:nrow(matr)) {
        xmax = max(matr[i, ], na.rm = TRUE)
        xmin = min(matr[i, ], na.rm = TRUE)
        if (xmax >= 1) { # check this so we only run this where it's needed
            row = matr[i, ]
            max2 = max(row[row < 1], na.rm = TRUE)
            matr[i, ] = pmin(row, 1 - (1 - max2)/2)
        }
        if (xmin <= 0) {
            row = matr[i, ]
            min2 = min(row[row > 0], na.rm = TRUE)
            matr[i, ] = pmax(matr[i, ], min2/2)
        }
    }
    # Transform into M values
    matr <- log2(matr/(1 - matr))
    return (matr)
}

#' The modelGeneExpression function
#' 
#' Model gene expression as a function of gene expression with a polynomial robust regression model (robust base package). 
#' Up to grade 3 models are evaluated and compared. R-squared and estimated coefficients are returned for the best model only for the genes with significant models.
#' @param METcancer matrix with methylation data for cancer samples (genes in rows, samples in columns).
#' @param MAcancer matrix with gene expression data for cancer samples (genes in rows, samples in columns).
#' @param CovariateData vector (numeric or character) indicating a covariate to be included in the model to adjust for it. Not used in an standard run of BIMEGA.
#' It can be used if samples can from different tissue type, for example.
#' @return matrix with R-square, degree of the best polynomial model and estimated coefficients only for the significant genes.
#' @importFrom foreach %dopar%
#' @export
#' @examples
#' # load data sets
#' data(METcancer)
#' data(MAcancer)
#' 
#' # model gene expression
#' results <- modelGeneExpression(METcancer, MAcancer)
#' 
modelGeneExpression <- function(METcancer, MAcancer, CovariateData = NULL) {
    
    grade = 3
    
    # overlapping samples     
    OverlapSamples = intersect(colnames(METcancer), colnames(MAcancer))
    cat("Found", length(OverlapSamples), "samples with both methylation and expression data.\n")
    MAcancer = MAcancer[, OverlapSamples, drop = FALSE]
    METcancer = METcancer[, OverlapSamples, drop = FALSE]
    if (!is.null(CovariateData)) CovariateData = as.matrix(CovariateData[OverlapSamples, ])
    
    Genes = rownames(METcancer)  
    PvalueThreshold = 0.001  
    RsquareThreshold = 0.1
    
    cat("Modeling gene expression with methylation data...\n")
    
    i <- NULL # to avoid "no visible binding for global variable" in R CMD check
    Results = foreach::foreach(i = seq(nrow(METcancer)), .combine = 'rbind', .export = "fitRobustRegression_robustbase", .packages = "robustbase") %dopar% {
    # Results = matrix(0, nrow(METcancer), 8)
    # for (i in seq(nrow(METcancer))) {
        options(warn = -1)
        #cat(i, "\n")
        result = numeric(8)
        typeModel = 1 # 1 if it is robust, -1 if robust couldn't be fitted and it is OLS
        tmpGene = unlist(strsplit(Genes[i], '---'))[1]
        pos = which(rownames(MAcancer) == tmpGene)
        if (length(pos) > 0) {
            x = METcancer[Genes[i], ]
            y = MAcancer[pos, ]
            corrXY = cor(x, y)
            # Fit polynomial regression up to grade 3
            if (!is.null(CovariateData)) {
                res0 = robustbase::lmrob(y ~ CovariateData, control = robustbase::lmrob.control(max.it = 500, k.max = 500, refine.tol = 1e-05))
                res1 = robustbase::lmrob(y ~ x + CovariateData, control = robustbase::lmrob.control(max.it = 500, k.max = 500, refine.tol = 1e-05))
                res2 = robustbase::lmrob(y ~ x + I(x^2) + CovariateData, control = robustbase::lmrob.control(max.it = 500, k.max = 500, refine.tol = 1e-05))
                res3 = robustbase::lmrob(y ~ x + I(x^2) + I(x^3) + CovariateData, control = robustbase::lmrob.control(max.it = 500, k.max = 500, refine.tol = 1e-05))
            } else {       
                res0 = fitRobustRegression_robustbase(y ~ 1)
                res1 = fitRobustRegression_robustbase(y ~ x)
                res2 = fitRobustRegression_robustbase(y ~ x + I(x^2))
                res3 = fitRobustRegression_robustbase(y ~ x + I(x^2) + I(x^3))
                if (any(c(class(res0), class(res1), class(res2), class(res3)) == "lm")) {
                    # one of the robust models couldnt be fitted so I used OLS, I haven't seen this but just in case. Fit all as OLS lm
                    res0 = lm(y ~ 1); res1 = lm(y ~ x); res2 = lm(y ~ x + I(x^2)); res3 = lm(y ~ x + I(x^2) + I(x^3))
                    typeModel = -1
                    anovaRes = anova(res3, res2, res1, res0)$`Pr(>F)`[-1]
                    anovaRes[is.na(anovaRes)] = 1
                } else {
                    # the robust regression worked always fine
                    # there's a bug in the anova function in the robust package and to work correctly you have to do anova for each pair of models
                    anovaRes = try(c(anova(res3, res2)$`Pr(>chisq)`[2], anova(res2, res1)$`Pr(>chisq)`[2], anova(res1, res0)$`Pr(>chisq)`[2]), silent = T)
                    # Even though what I explained in the fitRobustRegression_robustbase function, still sometimes the anova throws an error saying that Models are not strictly nested
                    # I'll catch it here and replace it by a linear model
                    if (inherits(anovaRes, "try-error")) {
                        res0 = lm(y ~ 1); res1 = lm(y ~ x); res2 = lm(y ~ x + I(x^2)); res3 = lm(y ~ x + I(x^2) + I(x^3))
                        typeModel = -1
                        anovaRes = anova(res3, res2, res1, res0)$`Pr(>F)`[-1]
                        anovaRes[is.na(anovaRes)] = 1
                    }
                }
            }
            # Compare the nested models
            if (anovaRes[1] < PvalueThreshold) {
                result = c(summary(res3)$r.squared, corrXY, 3, typeModel, summary(res3)$coefficients[, 1])
            } else if (anovaRes[2] < PvalueThreshold) {
                result = c(summary(res2)$r.squared, corrXY, 2, typeModel, summary(res2)$coefficients[, 1])
            } else if (anovaRes[3] < PvalueThreshold) {
                result = c(summary(res1)$r.squared, corrXY, 1, typeModel, summary(res1)$coefficients[, 1])
            } else {
                result = rep(0, 8)
            }
            length(result) <- 8
        }
        # Results[i, ] = result
        options(warn = 0)
        result
    }
    colnames(Results) = c("Rsquare", "LinCorr", "Grade", "ModelType", "Intercept", "Coef1", "Coef2", "Coef3")
    Results <- Results[, -4] # remove ModelType column
    
    # Rsquare threshold
    rownames(Results) = Genes
    FunctionalGenes = Results[Results[, "Rsquare"] > RsquareThreshold, , drop = FALSE]
    cat("\nFound", nrow(FunctionalGenes), "transcriptionally predictive genes.\n")
    return(FunctionalGenes)
}

#' The fitRobustRegression_robustbase function
#' 
#' Internal. Helper function to fit the robust regression.
#' @param fitFormula formula for the model to be fitted.
#' @return fitted robust regression model object.
#' @keywords internal
#' 
fitRobustRegression_robustbase <- function(fitFormula) {
    mod = robustbase::lmrob(fitFormula, control = robustbase::lmrob.control(max.it = 500, k.max = 500, refine.tol = 1e-05))
    if (!mod$converged | any(is.na(mod$coefficients))) {
        mod = robustbase::lmrob(fitFormula, control = robustbase::lmrob.control(setting="KS2011", max.it = 500, k.max = 500, refine.tol = 1e-05))
        if (!mod$converged | any(is.na(mod$coefficients))) {
            mod = lm(fitFormula)
        }
    }
    mod
    # the any(is.na(mod$coefficients)) check is because when all the values of methylation are very close to 1, x^3 gives is the same as x^2, and it gets an NA
    # the coefficient is indetermined, this is also true for the lm model.
    # The robustbase anova gives an error in this case, so I just use the lm models, do the anova, and the respective pvalue is going to be NA
    # so in the main function a change that NA to 1.
}

#' The MixtureModel_Biv function
#' 
#' Internal. Prepares all the structures to store the results and calls in a foreach loop a function that fits the mixture model in each gene.
#' @param METcancer matrix with methylation data for cancer samples (genes in rows, samples in columns).
#' @param METnormal matrix with methylation data for normal samples (genes in rows, samples in columns).
#' @param MAcancer matrix with gene expression data for cancer samples (genes in rows, samples in columns).
#' @param MAnormal optional matrix with gene expression data for normal samples (genes in rows, samples in columns).
#' @param FunctionalGenes vector with genes names to be considered for the mixture models.
#' @param NoNormalMode logical, if TRUE no comparison to normal samples is performed. Defaults to FALSE.
#' @return MethylationStates matrix of DM values, with driver genes in the rows and samples in the columns.
#' @return NrComponents matrix with the number of components identified for each driver gene.
#' @return Models list with the mixture model fitted for each driver gene.
#' @return MethylationDrivers character vector with the genes found by BIMEGA as differentially methylated and transcriptionally predictive (driver genes).
#' @return MixtureStates a list with a matrix for each driver gene containing the DM values.
#' @return Classifications a vector indicating to which component each sample was assigned.
#' @importFrom foreach %dopar%
#' @keywords internal
#' 
MixtureModel_Biv <- function(METcancer, METnormal, MAcancer, MAnormal = NULL, FunctionalGenes, NoNormalMode = FALSE) {
    
    # overlap of genes between all data sets
    overlapMET = intersect(rownames(METcancer), rownames(METnormal))
    overlapMA = rownames(MAcancer)
    if (!is.null(MAnormal)) overlapMA = intersect(rownames(MAcancer), rownames(MAnormal))
    
    # Intersect genes in MET and MA data sets without loosing the "---Cluster" in MET data sets
    overlapMETsplit  <- sapply(strsplit(overlapMET,  '---'), function(x) x[1])
    genes.to.keep.MET = overlapMET[overlapMETsplit %in% overlapMA]
    genes.to.keep.MA = overlapMA[overlapMA %in% overlapMETsplit]
    METcancer = METcancer[genes.to.keep.MET, , drop = FALSE]
    METnormal = METnormal[genes.to.keep.MET, , drop = FALSE]
    MAcancer = MAcancer[genes.to.keep.MA, , drop = FALSE]
    if (!is.null(MAnormal)) MAnormal = MAnormal[genes.to.keep.MA, , drop = FALSE]
    
    # Keep only those genes that were identified as functional
    METcancer = METcancer[rownames(METcancer) %in% FunctionalGenes, , drop = FALSE]  
    METnormal = METnormal[rownames(METnormal) %in% FunctionalGenes, , drop = FALSE]
    FunctionalGenesMA  <- sapply(strsplit(FunctionalGenes,  '---'), function(x) x[1])
    MAcancer = MAcancer[rownames(MAcancer) %in% FunctionalGenesMA, , drop = FALSE]
    if (!is.null(MAnormal)) MAnormal = MAnormal[rownames(MAnormal) %in% FunctionalGenesMA, , drop = FALSE]
    # METcancer and METnormal will have more rows (genes) than MAcancer and MAnormal
    # if there are cases where there are more than one in MET for a gene in MA, 
    # like "name---Cluster1", "name---Cluster2" in MET and "name" in MA
    
    # dimnames = list(rownames(METcancer), colnames(METcancer))
    
    GeneNamesMET = rownames(METcancer)
    GeneNamesMA = sapply(strsplit(GeneNamesMET,  '---'), function(x) x[1])
    
    cat("Running Gaussian mixture model on", length(rownames(METcancer)), "genes and on", length(colnames(METcancer)), "samples.\n")
    options(warn = -1)

    i <- NULL # to avoid "no visible binding for global variable" in R CMD check
    res <- foreach::foreach(i = seq(rownames(METcancer)), .combine = "combineForEachOutput", .export = "ModelSingleGene_Biv", .packages = "mclust") %dopar% {
        geneNameMET = GeneNamesMET[i]
        geneNameMA = GeneNamesMA[i]
        if (!is.null(MAnormal)) MAdataNormalVector = MAnormal[geneNameMA, ] else MAdataNormalVector = NULL
        ModelSingleGene_Biv(geneNameMET, METcancer[geneNameMET, ], MAcancer[geneNameMA, ], METnormal[geneNameMET, ], MAdataNormalVector, NoNormalMode = NoNormalMode)
    }
    rownames(res$MethylationStates) <- rownames(res$Classifications) <- rownames(METcancer)
    colnames(res$MethylationStates) <- colnames(res$Classifications) <- colnames(METcancer)
    options(warn = 0)
    
    # Removing the genes without any differential methylation. 
    if (!NoNormalMode) {
        # If NoNormalMode == T, no comparison to normal is made, and we don't remove genes with methylation states equal to 0
        # (for example, in pancancer analysis, running MethylMix only with normal samples, we use NoNormalMode = TRUE, and genes with only one state equal to 0 are kept.)
        NonZeroPositions = rowSums(res$MethylationStates) != 0
        res$NrComponents = res$NrComponents[NonZeroPositions]
        res$MixtureStates = res$MixtureStates[NonZeroPositions]
        res$Models = res$Models[NonZeroPositions]
        res$MethylationStates = res$MethylationStates[NonZeroPositions, , drop=FALSE]
        res$Classifications = res$Classifications[NonZeroPositions, , drop=FALSE]
    }
    
    # Adding names
    res$MethylationDrivers = rownames(res$MethylationStates)
    names(res$AllNrComponents) = names(res$MethylationStates)
    names(res$MixtureStates) = rownames(res$MethylationStates)
    names(res$Models) = rownames(res$MethylationStates)
    return(res)
}

#' The ModelSingleGene_Biv function
#' 
#' Internal. For a given gene, this function fits the mixture model, selects the number of components and defines the respective methylation states. 
#' @param GeneName character string with the name of the gene to model
#' @param METdataVector vector with methylation data for cancer samples.
#' @param MAdataVector vector with gene expression data for cancer samples
#' @param METdataNormalVector vector with methylation data for normal samples.
#' @param MAdataNormalVector vector with gene expression data for normal samples (optional).
#' @param test how to do the comparison between cancer samples in one state and the normal samples. Default is "ttest", the other option is "wilcoxon".
#' @param NoNormalMode logical, if TRUE no comparison to normal samples is performed. Defaults to FALSE.
#' @param maxComp maximum number of mixture components admitted in the model (3 by default).
#' @param PvalueThreshold threshold to consider results significant.
#' @param METDiffThreshold threshold in beta value scale from which two methylation means are considered different.
#' @param minSamplesPerGroup minimum number of samples required to belong to a new mixture component in order to accept it. Defaults to 5\% of all cancer samples.
#' @return NrComponents number of components identified.
#' @return Models an object of class 'Mclust' with the output model
#' @return MethylationStates vector with DM values for each sample.
#' @return MixtureStates vector with DMvalues for each component.
#' @return Classifications a vector indicating to which component each sample was assigned.
#' @details test, maxComp, PvalueThreshold, METDiffThreshold, minSamplesPerGroup are arguments for this function but are fixed in their default values for the user
#' because they are not available in the main BIMEGA function, to keep it simple. It would be easy to make them available to the user if we want to.
#' @keywords internal
#' 
ModelSingleGene_Biv <- function(GeneName, METdataVector, MAdataVector, METdataNormalVector, MAdataNormalVector = NULL,
                                NoNormalMode = FALSE, test = c("ttest", "wilcoxon"), maxComp = 3, 
                                PvalueThreshold = 0.01, METDiffThreshold = 0.10, minSamplesPerGroup = -1) {
    
    test <- match.arg(test)
    
    ## CHOOSE NUMBER OF COMPONENTS / STATES
    
    # 1 component model
    mods = vector("list", maxComp + 1)
    mods[[1]] = mclust::Mclust(data.frame(METdataVector, MAdataVector), G = 1, modelNames = "VVV")
    bic = numeric(maxComp + 1)
    bic[1] = - mods[[1]]$bic # mclust provides BIC with opposite sign                     
    
    # 2- to maxComp components model
    for (comp in 2:maxComp) {
        res = mclust::Mclust(data.frame(METdataVector, MAdataVector), G = comp, modelNames = "VVV")
        if (is.null(res)) mods[[comp]] = NA  else mods[[comp]] = res
        if (all(is.na(mods[[comp]]))) {
            # If it wasn't able to fit the model, try again with a conjugate prior on the variances
            res = mclust::Mclust(data.frame(METdataVector, MAdataVector), G = comp, modelNames = "VVV", prior = mclust::priorControl(scale = diag(1, comp)))
            if (is.null(res)) mods[[comp]] = NA  else mods[[comp]] = res
            # If it didn't converge again, quit trying this comp and keep comp-1 as the number of components.
            if (all(is.na(mods[[comp]]))) {
                NrComponents = comp - 1
                break
            }                
        }
        bic[comp] = - mods[[comp]]$bic
        
        modelMETmeans = summary(mods[[comp]])$mean["METdataVector", ]
        modelMETmeansBeta = sort((2^modelMETmeans) / (2^modelMETmeans + 1))
        
        DiffMETmeans = ifelse(all(abs(diff(modelMETmeansBeta)) > METDiffThreshold), T, F)
        
        # Check if smallest group has at least minSamplesPerGroup observations:
        if (minSamplesPerGroup < 2) {
            minSamplesPerGroup = max(5, 0.05 * length(METdataVector))
        }
        minOK = ifelse(min(table(mods[[comp]]$classification)) >= minSamplesPerGroup, TRUE, FALSE)
        
        # In original MethylMix, we try adding another component if the following 2 conditions are satisfied:
        #   A: Adding one component reduces BIC
        #   B: All absolute differences between methylation means in model with one extra component are above the MeanDifferenceThreshold
        # Here, A is the same, and for B I check differences in both methyltion and GE as described above
        # But I also check C = smallest group has at least minSamplesPerGroup observations:
        # So, continue with another component if A & B & C, else not continue, which is the same as saying not continue if !A OR !B OR !C        
        if ( (bic[comp] >= bic[comp - 1]) | !DiffMETmeans | !minOK) {
            NrComponents = comp - 1
            break
        } else {
            # It improved, try one more (unless comp already is maxComp, in that case we end here)
            NrComponents = comp
        }
    }                
    
    mod = mods[[NrComponents]]
    
    # Results of the classification
    classification = mod$classification
    if (NrComponents == 1) names(classification) = names(METdataVector)
    
    ## COMPARE EACH COMPONENT/STATE TO NORMAL
    MethylationState = matrix(0, 1, length(METdataVector)) # to hold the DMvalue in beta scale for each sample
    MixtureStates = matrix(0, NrComponents, 1) # to hold the DM-value for each component

    for (comp in 1:NrComponents) {
        
        METdataVector_comp = METdataVector[classification == comp]
        MAdataVector_comp = MAdataVector[classification == comp]
        
        # ttest or wilcoxon to test MET individually
        if (test == "ttest") {
            METres = t.test(METdataVector_comp, METdataNormalVector)
        } else if (test == "wilcoxon") {
            METres = wilcox.test(METdataVector_comp, METdataNormalVector)
        }
        
        pvalue = METres$p.value < PvalueThreshold
        
        # Methylation means in M-value and Beta-value
        METmeans = c(mean(METdataVector_comp), mean(METdataNormalVector))
        METmeansBeta = (2^METmeans) / (2^METmeans + 1)
        
        # Difference in methylation
        DiffMETmeans = ifelse(abs(METmeansBeta[1] - METmeansBeta[2]) > METDiffThreshold, T, F)
        
        # Difference in methylation? Use methylation to define both MixtureStates and METStates
        if ((pvalue & DiffMETmeans) | NoNormalMode) {
            DMvalueBeta = METmeansBeta[1] - METmeansBeta[2]
            MethylationState[1, classification == comp] = DMvalueBeta
            MixtureStates[comp, 1] = DMvalueBeta
        }
    }  
    
    message = ifelse(NrComponents == 1, " component is best.\n", " components are best.\n")
    cat(c(GeneName, ": ", NrComponents, message))
    
    return(list(NrComponents = NrComponents,
                Models = list(mod),
                MethylationStates = MethylationState,
                MixtureStates = list(MixtureStates), 
                Classifications = classification))
}

#' The combineForEachOutput function
#' 
#' Internal. Function to combine results from the foreach loop.
#' @param out1 result from one foreach loop.
#' @param out2 result from another foreach loop.
#' @return List with the combined results.
#' @keywords internal
#'
combineForEachOutput <- function(out1, out2) {
    NrComponents <- c(out1$NrComponents, out2$NrComponents)
    MethylationStates <- rbind(out1$MethylationStates, out2$MethylationStates)
    MixtureStates <- append(out1$MixtureStates, out2$MixtureStates)
    Models <- append(out1$Models, out2$Models)
    Classifications <- rbind(out1$Classifications, out2$Classifications)
    return(list(NrComponents = NrComponents,
                MethylationStates = MethylationStates,
                MixtureStates = MixtureStates,
                Models = Models,
                Classifications = Classifications))
}

#' The BIMEGA_Plot function.
#' 
#' Produces plots to represent BIMEGA's output.
#' @param GeneName Name of the gene for which to create a BIMEGA plot.
#' @param MixtureModelResults List returned by BIMEGA function.
#' @param METcancer	Matrix with the methylation data of cancer tissue with genes in rows and samples in columns.
#' @param MAcancer Gene expression data for cancer tissue with genes in rows and samples in columns.
#' @param METnormal	Matrix with the normal methylation data of the same genes as in METcancer (optional). Again genes in rows and samples in columns.
#' @param MAnormal Gene expression data for normal tissue with genes in rows and samples in columns (optional).
#' @param title A title for the plot.
#' @return BIMEGA plot, a scatterplot between DNA methylation and gene expression, showing the different mixture components identified.
#' @importFrom dplyr %>%
#' @export
#' @examples 
#' # load the data sets needed for BIMEGA
#' data(METcancer)
#' data(METnormal)
#' data(MAcancer)
#' data(MAnormal)
#' 
#' # run BIMEGA on a small set of example data
#' BIMEGAresults <- BIMEGA(METcancer, METnormal, MAcancer, MAnormal)
#' 
#' # Produce plots of differentially methylated genes
#' for (gene in BIMEGAresults$MethylationDrivers) {
#'      g <- BIMEGA_Plot(gene, BIMEGAresults, METcancer, MAcancer, METnormal, MAnormal)
#'      plot(g)
#' }
#' 
BIMEGA_Plot <- function(GeneName, MixtureModelResults, METcancer, MAcancer, METnormal = NULL, MAnormal = NULL, title = NULL) {
    
    met <- ge <- group <- NULL # to avoid "no visible binding for global variable" in R CMD check
    
    GeneNameMA <- unlist(strsplit(GeneName, '---'))[1]
    
    if (GeneName %in% MixtureModelResults$MethylationDrivers && GeneName %in% rownames(METcancer) && GeneNameMA %in% rownames(MAcancer)) {
        # Prepare data for plot
        grade <- MixtureModelResults$FunctionalGenesResults[GeneName, "Grade"]
        OverlapSamples <- intersect(colnames(METcancer), colnames(MAcancer))
        data <- data.frame(met = METcancer[GeneName, OverlapSamples],
                           ge = MAcancer[GeneNameMA, OverlapSamples],
                           group = factor(MixtureModelResults$Classification[GeneName, OverlapSamples], levels = 0:length(unique(MixtureModelResults$Classification[GeneName, ]))))
        data$met[data$met < 0] <- 0
        data$met[data$met > 1] <- 1
        cols <- c("black", RColorBrewer::brewer.pal(8, "Dark2")[1:length(unique(data$group))])
        names(cols) <- c("Normal", 1:length(unique(data$group)))
        # Create plot
        g <- ggplot2::ggplot(data, ggplot2::aes(x = met, y = ge, color = group)) + 
            ggplot2::geom_point(shape = 19, alpha = 0.4, size = 4) +
            ggplot2::scale_color_manual(values = cols, name="Mixture\ncomponent") + #, labels=c("Normal", 1:length(unique(data$group)))) +
            ggplot2::scale_y_continuous(name = "Gene expression") +
            ggplot2::scale_x_continuous(name = "Methylation (Beta values)", limits = 0:1) + 
            # ggplot2::stat_ellipse(type = "norm", size = 1.3) +
            ggplot2::theme(axis.text = ggplot2::element_text(size=16),
                            axis.title = ggplot2::element_text(size=16),
                            legend.title = ggplot2::element_text(size=14),
                            legend.text = ggplot2::element_text(size=14),
                            plot.title = ggplot2::element_text(size=18, face="bold")) +
            ggplot2::stat_smooth(method = "lm", se = F, formula =y ~ poly(x, grade, raw = TRUE), colour = "red")
        
        # Add symbols for centroids
        centroids <- data %>% dplyr::group_by(group) %>% dplyr::summarise(met = mean(met), ge = mean(ge))
        g <- g + 
            ggplot2::geom_point(data = centroids, ggplot2::aes(fill = group), shape = 21, size = 5, col = "white", stroke = 2) + 
            ggplot2::scale_fill_manual(values = cols, guide = "none")
        
        if (!is.null(title)) {
            if (title != "") {
                g <- g + ggplot2::ggtitle(title) # Add provided title
            } # else add no title
        } else {
            g <- g + ggplot2::ggtitle(paste0("Mixture components for ", GeneName)) # add default title
        }
        
        # If METnormal and MAnormal available, add ellipse
        if (!is.null(METnormal) && !is.null(MAnormal)) {
            OverlapSamples <- intersect(colnames(METnormal), colnames(MAnormal))
            # only add ellipse for samples with both MET and MA data
            if (length(OverlapSamples) > 1) {
                dataEllipse <- data.frame(met = METnormal[GeneName, OverlapSamples],
                                          ge =  MAnormal[GeneNameMA, OverlapSamples])
                dataEllipse <- data.frame(ellipse::ellipse(cov(dataEllipse), centre = colMeans(dataEllipse, na.rm = T), level = 0.95))
                dataEllipse$group <- factor(rep("Normal", nrow(dataEllipse), levels = c("Normal", 1:length(table(data$group)))))
                g <- g + ggplot2::geom_path(data = dataEllipse, ggplot2::aes(x = met, y = ge), size = 1.5)
                centroidNormal <- data.frame(met = mean(dataEllipse$met), ge = mean(dataEllipse$ge))
            } else {
                centroidNormal <- data.frame(met = mean(METnormal[GeneName, ]), ge = mean(MAnormal[GeneNameMA, ]))
            }
            g <- g + ggplot2::geom_point(data = centroidNormal, shape = 21, size = 3, col = "black", fill = "black", stroke = 3)
        }
        
        # Add line with confidence interval for methylation with normal data
        if (!is.null(METnormal)) {
            metnormal <- METnormal[GeneName, , drop=FALSE]
            if (length(metnormal) > 1) {
                tmpTtest1 <- t.test(metnormal)
                g <- g + ggplot2::geom_segment(ggplot2::aes(x = tmpTtest1$conf.int[1], y = min(data$ge), xend = tmpTtest1$conf.int[2], yend = min(data$ge)), color = "black", size = 1.5)
            }
        }
        # Add line with confidence interval for gene expression with normal data
        if (!is.null(MAnormal)) {
            manormal <- MAnormal[GeneNameMA, ]
            if (length(manormal) > 1) {
                tmpTtest2 <- t.test(manormal)
                g <- g + ggplot2::geom_segment(ggplot2::aes(x = min(data$met), y = tmpTtest2$conf.int[1], xend = min(data$met), yend = tmpTtest2$conf.int[2]), color = "black", size = 1.5)
            }
        }
        return(g)
    } else {
        cat(paste(GeneName, "not found.\n"))
        return(NULL)
    }
}

