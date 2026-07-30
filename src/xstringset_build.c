/*
 * xstringset_build.c — construct Biostrings XStringSet objects (DNAStringSet,
 * PhredQuality) directly from a contiguous byte buffer, using the
 * R_GetCCallable-exported XVector/IRanges C constructors (the same interface
 * Rsamtools itself uses). The buffer becomes the XStringSet's shared pool tag
 * with zero per-element copy.
 *
 * Compiled as C so the un-`extern "C"` Bioconductor interface headers and their
 * stub translation units link cleanly. bamscale.cpp calls the two `bs_*`
 * entries below through an `extern "C"` declaration.
 *
 * The `_*_stubs.c` includes define lazy R_GetCCallable resolvers (no link-time
 * dependency on the XVector/IRanges shared objects) and MUST appear in exactly
 * one translation unit — this one.
 */
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "IRanges_interface.h"
#include "XVector_interface.h"

#include "_IRanges_stubs.c"
/* Both stub files define these helper macros identically; undo before the
 * second include to avoid a benign -Wmacro-redefined warning. */
#undef DEFINE_CCALLABLE_STUB
#undef DEFINE_NOVALUE_CCALLABLE_STUB
#include "_XVector_stubs.c"

/*
 * Multi-buffer builder.
 *   tags               : VECSXP of RAWSXP pool buffers, each length < INT_MAX
 *   start, width, group: INTSXP length n; start is 1-based within its tag,
 *                        group is the 1-based index of the tag each range lives in.
 * Returns the finished S4 object (e.g. DNAStringSet / PhredQuality).
 */
SEXP bs_new_xrawlist(const char *classname, const char *element_type,
                     SEXP tags, SEXP start, SEXP width, SEXP group)
{
    SEXP ranges = PROTECT(new_IRanges("IRanges", start, width, R_NilValue));
    SEXP ans    = PROTECT(new_XRawList_from_tags(classname, element_type,
                                                 tags, ranges, group));
    UNPROTECT(2);
    return ans;
}

/*
 * Single-buffer builder (the < 2 GB case → the entire benchmark). Uses
 * new_XRawList_from_tag (no group vector), which yields the exact single-pool,
 * group-all-1 layout Rsamtools::scanBam produces — so the result is
 * byte-for-byte identical() to scanBam.
 *   tag         : one RAWSXP pool buffer (length < INT_MAX)
 *   start, width: INTSXP length n; start is 1-based within the tag.
 */
SEXP bs_new_xrawlist_single(const char *classname, const char *element_type,
                            SEXP tag, SEXP start, SEXP width)
{
    SEXP ranges = PROTECT(new_IRanges("IRanges", start, width, R_NilValue));
    SEXP ans    = PROTECT(new_XRawList_from_tag(classname, element_type,
                                                tag, ranges));
    UNPROTECT(2);
    return ans;
}
