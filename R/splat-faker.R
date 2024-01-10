#' Splat faker
#'
#' Fake count data from a fictional single-cell RNA-seq experiment using
#' the Splat method.
#'
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' @param method which simulation method to use. Options are "single" which
#'        produces a single population, "groups" which produces distinct groups
#'        (eg. cell types), or "paths" which selects cells from continuous
#'        trajectories (eg. differentiation processes).
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' Parameters can be set in a variety of ways. If no parameters are provided
#' the default parameters are used. Any parameters in \code{params} can be
#' overridden by supplying additional arguments through a call to
#' \code{\link{setParams}}. This design allows the user flexibility in
#' how they supply parameters and allows small adjustments without creating a
#' new \code{SplatParams} object. See examples for a demonstration of how this
#' can be used.
#'
#' The simulation involves the following steps:
#' \enumerate{
#'     \item Set up simulation object
#'     \item Simulate library sizes
#'     \item Simulate gene means
#'     \item Simulate groups/paths
#'     \item Simulate BCV adjusted cell means
#'     \item Simulate true counts
#'     \item Simulate dropout
#'     \item Create final dataset
#' }
#'
#' The final output is a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object that
#' contains the simulated counts but also the values for various intermediate
#' steps. These are stored in the \code{\link{colData}} (for cell specific
#' information), \code{\link{rowData}} (for gene specific information) or
#' \code{\link{assays}} (for gene by cell matrices) slots. This additional
#' information includes:
#' \describe{
#'     \item{\code{colData}}{
#'         \describe{
#'             \item{Cell}{Unique cell identifier.}
#'             \item{Group}{The group or path the cell belongs to.}
#'             \item{ExpLibSize}{The expected library size for that cell.}
#'             \item{Step (paths only)}{how far along the path each cell is.}
#'         }
#'     }
#'     \item{\code{rowData}}{
#'         \describe{
#'             \item{Gene}{Unique gene identifier.}
#'             \item{BaseGeneMean}{The base expression level for that gene.}
#'             \item{OutlierFactor}{Expression outlier factor for that gene.
#'             Values of 1 indicate the gene is not an expression outlier.}
#'             \item{GeneMean}{Expression level after applying outlier factors.}
#'             \item{BatchFac[Batch]}{The batch effects factor for each gene for
#'             a particular batch.}
#'             \item{DEFac[Group]}{The differential expression factor for each
#'             gene in a particular group. Values of 1 indicate the gene is not
#'             differentially expressed.}
#'             \item{SigmaFac[Path]}{Factor applied to genes that have
#'             non-linear changes in expression along a path.}
#'         }
#'     }
#'     \item{\code{assays}}{
#'         \describe{
#'             \item{BatchCellMeans}{The mean expression of genes in each cell
#'             after adding batch effects.}
#'             \item{BaseCellMeans}{The mean expression of genes in each cell
#'             after any differential expression and adjusted for expected
#'             library size.}
#'             \item{BCV}{The Biological Coefficient of Variation for each gene
#'             in each cell.}
#'             \item{CellMeans}{The mean expression level of genes in each cell
#'             adjusted for BCV.}
#'             \item{TrueCounts}{The simulated counts before dropout.}
#'             \item{Dropout}{Logical matrix showing which values have been
#'             dropped in which cells.}
#'         }
#'     }
#' }
#'
#' Values that have been added by Splatter are named using \code{UpperCamelCase}
#' in order to differentiate them from the values added by analysis packages
#' which typically use \code{underscore_naming}.
#'
#' @return SingleCellExperiment object containing the simulated counts and
#' intermediate values.
#'
#' @references
#' Zappia L, Phipson B, Oshlack A. Splatter: simulation of single-cell RNA
#' sequencing data. Genome Biology (2017).
#'
#' Paper: \url{10.1186/s13059-017-1305-0}
#'
#' Code: \url{https://github.com/Oshlack/splatter}
#'
#' @seealso
#' \code{\link{splatSimLibSizes}}, \code{\link{splatSimGeneMeans}},
#' \code{\link{splatSimBatchEffects}}, \code{\link{splatSimBatchCellMeans}},
#' \code{\link{splatSimDE}}, \code{\link{splatSimCellMeans}},
#' \code{\link{splatSimBCVMeans}}, \code{\link{splatSimTrueCounts}},
#' \code{\link{splatSimDropout}}
#'
#' @examples
#' # Simulation with default parameters
#' sim <- splatSimulate()
#' \dontrun{
#' # Simulation with different number of genes
#' sim <- splatSimulate(nGenes = 1000)
#' # Simulation with custom parameters
#' params <- newSplatParams(nGenes = 100, mean.rate = 0.5)
#' sim <- splatSimulate(params)
#' # Simulation with adjusted custom parameters
#' sim <- splatSimulate(params, mean.rate = 0.6, out.prob = 0.2)
#' # Simulate groups
#' sim <- splatSimulate(method = "groups")
#' # Simulate paths
#' sim <- splatSimulate(method = "paths")
#' }
#' @importFrom SummarizedExperiment rowData colData colData<- assays
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods validObject
#' @export
splatFake <- function(gene_means, params = newSplatParams(),
                          method = c("single", "groups", "paths"),
                          cell_name_prefix = "Cell",
                          sparsify = TRUE, verbose = TRUE, ...) {

    checkmate::assertClass(params, "SplatParams")
    checkmate::assertVector(gene_means, strict=TRUE, any.missing=FALSE, min.len=1, names="strict")

    method <- match.arg(method)

    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)

    # Set random seed
    seed <- getParam(params, "seed")
    withr::with_seed(seed, {

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- len(gene_means)
    nBatches <- getParam(params, "nBatches")
    batch.cells <- getParam(params, "batchCells")
    nGroups <- getParam(params, "nGroups")
    group.prob <- getParam(params, "group.prob")

    if (nGroups == 1 && method == "groups") {
        warning("nGroups is 1, switching to single mode")
        method <- "single"
    }

    # Set up name vectors
    if (verbose) {message("Creating simulation object...")}
    cell.names <- paste0(cell_name_prefix, seq_len(nCells))
    gene.names <- names(gene_means)
    batch.names <- paste0("Batch", seq_len(nBatches))
    if (method == "groups") {
        group.names <- paste0("Group", seq_len(nGroups))
    } else if (method == "paths") {
        group.names <- paste0("Path", seq_len(nGroups))
    }

    # Create SingleCellExperiment to store simulation
    cells <-  data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    features <- data.frame(Gene = gene.names)
    rownames(features) <- gene.names
    sim <- SingleCellExperiment(rowData = features, colData = cells,
                                metadata = list(Params = params))

    # Make batches vector which is the index of param$batchCells repeated
    # params$batchCells[index] times
    batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                      b = batch.cells)
    batches <- unlist(batches)
    colData(sim)$Batch <- batch.names[batches]

    if (method != "single") {
        groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                         replace = TRUE)
        colData(sim)$Group <- factor(group.names[groups], levels = group.names)
    }

    if (verbose) {message("Simulating library sizes...")}
    sim <- splatSimLibSizes(sim, params)
    if (verbose) {message("Simulating gene means...")}
    sim <- splatSetGeneMeans(sim, params, gene_means)
    if (nBatches > 1) {
        if (verbose) {message("Simulating batch effects...")}
        sim <- splatSimBatchEffects(sim, params)
    }
    sim <- splatSimBatchCellMeans(sim, params)
    if (method == "single") {
        sim <- splatSimSingleCellMeans(sim, params)
    } else if (method == "groups") {
        if (verbose) {message("Simulating group DE...")}
        sim <- splatSimGroupDE(sim, params)
        if (verbose) {message("Simulating cell means...")}
        sim <- splatSimGroupCellMeans(sim, params)
    } else {
        if (verbose) {message("Simulating path endpoints...")}
        sim <- splatSimPathDE(sim, params)
        if (verbose) {message("Simulating path steps...")}
        sim <- splatSimPathCellMeans(sim, params)
    }
    if (verbose) {message("Simulating BCV...")}
    sim <- splatSimBCVMeans(sim, params)
    if (verbose) {message("Simulating counts...")}
    sim <- splatSimTrueCounts(sim, params)
    if (verbose) {message("Simulating dropout (if needed)...")}
    sim <- splatSimDropout(sim, params)

    if (sparsify) {
        if (verbose) {message("Sparsifying assays...")}
        assays(sim) <- sparsifyMatrices(assays(sim), auto = TRUE,
                                        verbose = verbose)
    }
    })

    if (verbose) {message("Done!")}
    return(sim)

}

#' Set gene means
#'
#' Set gene means from a given vector. Also simulates outlier expression factors.
#' Genes with an outlier factor not equal to 1 are replaced
#' with the median mean expression multiplied by the outlier factor.
#'
#' @param sim SingleCellExperiment to add gene means to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated gene means.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom stats rgamma median
splatSetGeneMeans <- function(sim, params, gene_means) {

    # Note: splatPopSimGeneMeans in splatPop-simulate.R mirrors this function.
    # If changes are made to the "add expression outliers" method here, please
    # make the same changes in splatPopSimGeneMeans.

    nGenes <- len(gene_means)
    out.prob <- getParam(params, "out.prob")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")

    # Simulate base gene means, just assign
    base.means.gene <- gene_means

    # Add expression outliers
    outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc,
                                    out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]

    rowData(sim)$BaseGeneMean <- base.means.gene
    rowData(sim)$OutlierFactor <- outlier.facs
    rowData(sim)$GeneMean <- means.gene

    return(sim)
}
