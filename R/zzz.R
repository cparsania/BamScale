#' @keywords internal
#' @useDynLib BamScale, .registration = TRUE
#' @import methods
#' @importFrom BiocParallel bplapply bpnworkers bpworkers bpparam
#' @importFrom BiocGenerics path unlist
#' @importFrom Rcpp evalCpp
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
    # The compatible seq/qual path builds DNAStringSet/PhredQuality in C via the
    # C-callable constructors exported by IRanges and XVector, resolved at runtime
    # with R_GetCCallable. Those symbols are only registered once the providing
    # package's shared object is loaded. Force-load the namespaces here so the
    # callables are available in every process -- in particular fresh SnowParam
    # SOCK workers, which start clean and would otherwise fail with
    # "function '_new_IRanges' not provided by package 'IRanges'".
    # IRanges/XVector register the C-callable constructors; Biostrings defines the
    # DNAStringSet/PhredQuality S4 classes the constructors instantiate. All must
    # be loaded in the worker before the builder runs.
    #
    # GenomicAlignments must be loaded for the same reason: the GAlignments /
    # GAlignmentPairs output path builds objects via methods::new("GAlignments"),
    # which fails with "'GAlignments' is not a defined class" unless the class is
    # registered. BamScale reaches GenomicAlignments only via `::`, so loading
    # BamScale alone does not register the class in the current (or a worker)
    # process.
    for (pkg in c("S4Vectors", "IRanges", "XVector", "Biostrings", "GenomicAlignments")) {
        requireNamespace(pkg, quietly = TRUE)
    }
    invisible()
}
