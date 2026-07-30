#include <Rcpp.h>
using namespace Rcpp;

// Required to print ompBAM logging to the R console.
#define cout Rcpp::Rcout

// [[Rcpp::depends(ompBAM)]]
#include <ompBAM.hpp>

#include <algorithm>
#include <cstring>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Defined in xstringset_build.c: construct Biostrings XStringSet objects
// (DNAStringSet / PhredQuality) directly from a shared byte buffer via the
// R_GetCCallable-exported XVector/IRanges constructors.
extern "C" {
SEXP bs_new_xrawlist_single(const char* classname, const char* element_type,
                            SEXP tag, SEXP start, SEXP width);
SEXP bs_new_xrawlist(const char* classname, const char* element_type,
                     SEXP tags, SEXP start, SEXP width, SEXP group);
}

namespace {

static const std::string TAG_NA_SENTINEL = "__BAMSCALE_TAG_NA__";

unsigned int clamp_threads(const int requested) {
#ifdef _OPENMP
  if (requested <= 1) return 1U;
  const unsigned int max_threads = static_cast<unsigned int>(omp_get_max_threads());
  if (static_cast<unsigned int>(requested) > max_threads) return max_threads;
  return static_cast<unsigned int>(requested);
#else
  (void)requested;
  return 1U;
#endif
}

struct PackedStringColumn;

inline void append_decoded_seq_ascii(PackedStringColumn& out, const uint8_t* packed_seq, const uint32_t n_bases);
inline void append_decoded_seq_encoded(PackedStringColumn& out, const uint8_t* packed_seq, const uint32_t n_bases);
inline void append_encoded_qual_ascii(PackedStringColumn& out, const uint8_t* qual_ptr, const uint32_t n_bases);
inline void append_encoded_qual_native(PackedStringColumn& out, const uint8_t* qual_ptr, const uint32_t n_bases);

template <typename T>
std::string join_numeric_vector(const std::vector<T>& x) {
  std::ostringstream oss;
  for (size_t i = 0; i < x.size(); ++i) {
    if (i > 0) oss << ',';
    oss << static_cast<double>(x[i]);
  }
  return oss.str();
}

std::string extract_tag_as_string(pbam1_t& read, const std::string& tag) {
  const char type = read.Tag_Type(tag);
  switch (type) {
    case 'A': {
      const char v = read.tagVal_A(tag);
      if (v == '\0') return TAG_NA_SENTINEL;
      return std::string(1, v);
    }
    case 'c':
      return std::to_string(static_cast<int>(read.tagVal_c(tag)));
    case 'C':
      return std::to_string(static_cast<unsigned int>(read.tagVal_C(tag)));
    case 's':
      return std::to_string(static_cast<int>(read.tagVal_s(tag)));
    case 'S':
      return std::to_string(static_cast<unsigned int>(read.tagVal_S(tag)));
    case 'i':
      return std::to_string(read.tagVal_i(tag));
    case 'I':
      return std::to_string(read.tagVal_I(tag));
    case 'f': {
      std::ostringstream oss;
      oss << read.tagVal_f(tag);
      return oss.str();
    }
    case 'Z': {
      std::string out;
      const int ret = read.tagVal_Z(tag, out);
      if (ret < 0) return TAG_NA_SENTINEL;
      if (!out.empty() && out.back() == '\0') out.pop_back();
      return out;
    }
    case 'B': {
      const char subtype = read.Tag_Subtype(tag);
      switch (subtype) {
        case 'c': {
          std::vector<int8_t> out;
          if (read.tagVal_B(tag, out) < 0) return TAG_NA_SENTINEL;
          return join_numeric_vector(out);
        }
        case 'C': {
          std::vector<uint8_t> out;
          if (read.tagVal_B(tag, out) < 0) return TAG_NA_SENTINEL;
          return join_numeric_vector(out);
        }
        case 's': {
          std::vector<int16_t> out;
          if (read.tagVal_B(tag, out) < 0) return TAG_NA_SENTINEL;
          return join_numeric_vector(out);
        }
        case 'S': {
          std::vector<uint16_t> out;
          if (read.tagVal_B(tag, out) < 0) return TAG_NA_SENTINEL;
          return join_numeric_vector(out);
        }
        case 'i': {
          std::vector<int32_t> out;
          if (read.tagVal_B(tag, out) < 0) return TAG_NA_SENTINEL;
          return join_numeric_vector(out);
        }
        case 'I': {
          std::vector<uint32_t> out;
          if (read.tagVal_B(tag, out) < 0) return TAG_NA_SENTINEL;
          return join_numeric_vector(out);
        }
        case 'f': {
          std::vector<float> out;
          if (read.tagVal_B(tag, out) < 0) return TAG_NA_SENTINEL;
          return join_numeric_vector(out);
        }
        default:
          return TAG_NA_SENTINEL;
      }
    }
    default:
      return TAG_NA_SENTINEL;
  }
}

List raw_list_from_bytes(const std::vector<std::vector<uint8_t> >& buffers) {
  const R_xlen_t n = static_cast<R_xlen_t>(buffers.size());
  List out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const std::vector<uint8_t>& src = buffers[static_cast<size_t>(i)];
    RawVector dst(static_cast<R_xlen_t>(src.size()));
    if (!src.empty()) {
      std::copy(src.begin(), src.end(), dst.begin());
    }
    out[i] = dst;
  }
  return out;
}

struct PackedStringColumn {
  std::vector<size_t> offsets;
  std::vector<size_t> lengths;
  std::vector<char> bytes;

  void clear() {
    offsets.clear();
    lengths.clear();
    bytes.clear();
  }

  void reserve(const size_t n_strings, const size_t n_bytes) {
    if (n_strings > 0U) {
      offsets.reserve(n_strings);
      lengths.reserve(n_strings);
    }
    if (n_bytes > 0U) {
      bytes.reserve(n_bytes);
    }
  }

  size_t size() const {
    return offsets.size();
  }

  size_t bytes_size() const {
    return bytes.size();
  }

  void append(const char* ptr, const size_t len) {
    const size_t offset = bytes.size();
    offsets.push_back(offset);
    lengths.push_back(len);
    if (ptr != NULL && len > 0U) {
      bytes.resize(offset + len);
      std::memcpy(bytes.data() + offset, ptr, len);
    }
  }

  void append(const std::string& value) {
    append(value.data(), value.size());
  }

  void append_from(const PackedStringColumn& other) {
    const size_t base = bytes.size();
    offsets.reserve(offsets.size() + other.offsets.size());
    lengths.reserve(lengths.size() + other.lengths.size());
    bytes.reserve(bytes.size() + other.bytes.size());

    for (size_t i = 0; i < other.offsets.size(); ++i) {
      offsets.push_back(base + other.offsets[i]);
      lengths.push_back(other.lengths[i]);
    }
    bytes.insert(bytes.end(), other.bytes.begin(), other.bytes.end());
  }
};

inline void append_decoded_seq_ascii(PackedStringColumn& out, const uint8_t* packed_seq, const uint32_t n_bases) {
  static const char nt16_map[16] = {
    '=', 'A', 'C', 'M',
    'G', 'R', 'S', 'V',
    'T', 'W', 'Y', 'H',
    'K', 'D', 'B', 'N'
  };

  const size_t offset = out.bytes.size();
  const size_t len = static_cast<size_t>(n_bases);
  out.offsets.push_back(offset);
  out.lengths.push_back(len);

  if (len == 0U || packed_seq == NULL) return;

  out.bytes.resize(offset + len);
  char* dst = out.bytes.data() + offset;

  const uint32_t n_pairs = n_bases >> 1;
  uint32_t j = 0U;
  for (uint32_t i = 0; i < n_pairs; ++i) {
    const uint8_t byte = packed_seq[i];
    dst[j++] = nt16_map[static_cast<uint8_t>((byte >> 4) & 0x0FU)];
    dst[j++] = nt16_map[static_cast<uint8_t>(byte & 0x0FU)];
  }

  if ((n_bases & 1U) != 0U) {
    const uint8_t byte = packed_seq[n_pairs];
    dst[j] = nt16_map[static_cast<uint8_t>((byte >> 4) & 0x0FU)];
  }
}

// Fill `out` with DNAString *internal* encoding: BAM nt16 codes 1-15 map to
// themselves and code 0 ('=') maps to 16, exactly matching the byte layout a
// DNAStringSet's shared pool holds (verified against Rsamtools::scanBam). This
// skips the ASCII intermediate entirely — the nibble becomes the pool byte.
inline void append_decoded_seq_encoded(PackedStringColumn& out, const uint8_t* packed_seq, const uint32_t n_bases) {
  static const char enc[16] = {
    16, 1, 2, 3,
    4, 5, 6, 7,
    8, 9, 10, 11,
    12, 13, 14, 15
  };

  const size_t offset = out.bytes.size();
  const size_t len = static_cast<size_t>(n_bases);
  out.offsets.push_back(offset);
  out.lengths.push_back(len);

  if (len == 0U || packed_seq == NULL) return;

  out.bytes.resize(offset + len);
  char* dst = out.bytes.data() + offset;

  const uint32_t n_pairs = n_bases >> 1;
  uint32_t j = 0U;
  for (uint32_t i = 0; i < n_pairs; ++i) {
    const uint8_t byte = packed_seq[i];
    dst[j++] = enc[static_cast<uint8_t>((byte >> 4) & 0x0FU)];
    dst[j++] = enc[static_cast<uint8_t>(byte & 0x0FU)];
  }

  if ((n_bases & 1U) != 0U) {
    const uint8_t byte = packed_seq[n_pairs];
    dst[j] = enc[static_cast<uint8_t>((byte >> 4) & 0x0FU)];
  }
}

inline void append_encoded_qual_ascii(PackedStringColumn& out, const uint8_t* qual_ptr, const uint32_t n_bases) {
  const size_t offset = out.bytes.size();
  const size_t len = static_cast<size_t>(n_bases);

  if (len == 0U || qual_ptr == NULL) {
    out.offsets.push_back(offset);
    out.lengths.push_back(0U);
    return;
  }

  bool missing = true;
  out.bytes.resize(offset + len);
  char* dst = out.bytes.data() + offset;
  for (uint32_t i = 0; i < n_bases; ++i) {
    const uint8_t q = qual_ptr[i];
    if (q != 255U) missing = false;
    const uint8_t phred = (q == 255U) ? 0U : q;
    dst[static_cast<size_t>(i)] = static_cast<char>(phred + 33U);
  }

  out.offsets.push_back(offset);
  if (missing) {
    out.lengths.push_back(1U);
    out.bytes.resize(offset + 1U);
    out.bytes[offset] = '*';
  } else {
    out.lengths.push_back(len);
  }
}

// Fill `out` with PhredQuality bytes exactly as Rsamtools::scanBam does: full
// read width, each byte (q + 33) & 0xFF. A missing quality string (ompBAM fills
// it with 0xFF to l_seq) therefore becomes width-many 0x20 (space) bytes, NOT a
// width-1 "*". This matches scanBam byte-for-byte.
inline void append_encoded_qual_native(PackedStringColumn& out, const uint8_t* qual_ptr, const uint32_t n_bases) {
  const size_t offset = out.bytes.size();
  const size_t len = static_cast<size_t>(n_bases);
  out.offsets.push_back(offset);
  out.lengths.push_back(len);

  if (len == 0U || qual_ptr == NULL) return;

  out.bytes.resize(offset + len);
  char* dst = out.bytes.data() + offset;
  for (uint32_t i = 0; i < n_bases; ++i) {
    dst[static_cast<size_t>(i)] =
        static_cast<char>((static_cast<unsigned int>(qual_ptr[i]) + 33U) & 0xFFU);
  }
}

CharacterVector packed_strings_to_character(const PackedStringColumn& packed) {
  const R_xlen_t n = static_cast<R_xlen_t>(packed.size());
  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP empty = Rf_mkCharLen("", 0);

  for (R_xlen_t i = 0; i < n; ++i) {
    const size_t idx = static_cast<size_t>(i);
    const size_t offset = packed.offsets[idx];
    const size_t len = packed.lengths[idx];
    if (len == 0U) {
      SET_STRING_ELT(out, i, empty);
    } else {
      SET_STRING_ELT(
        out,
        i,
        Rf_mkCharLen(packed.bytes.data() + offset, static_cast<int>(len))
      );
    }
  }

  UNPROTECT(1);
  return CharacterVector(out);
}

// Split threshold for the shared pool. A DNAStringSet/PhredQuality pool buffer,
// like R's SharedRaw and the IRanges start/width, is indexed by int, so a single
// buffer can hold up to INT_MAX bytes. Rsamtools::scanBam uses exactly this
// limit, so matching it means any result scanBam can produce (< INT_MAX bytes)
// takes our single-buffer path and is byte-for-byte identical() to scanBam;
// only a > INT_MAX result (which scanBam cannot return either) splits into
// multiple buffers. `total < BS_TAG_CAP` keeps every start+1 within int.
static const size_t BS_TAG_CAP = 2147483647U;  // INT_MAX

// Build a DNAStringSet / PhredQuality from a merged PackedStringColumn, sharing
// its `bytes` buffer as the XStringSet pool (one memcpy, no per-element CHARSXP).
// `classname`/`element_type` select the S4 type: ("DNAStringSet","DNAString")
// or ("PhredQuality","BString"). Single buffer when total < BS_TAG_CAP (the
// whole benchmark, identical() to scanBam); multi-buffer only above that.
Rcpp::RObject xstringset_from_packed(const PackedStringColumn& col,
                                     const char* classname,
                                     const char* element_type) {
  const R_xlen_t n = static_cast<R_xlen_t>(col.size());
  const size_t total = col.bytes.size();

  if (n == 0) {
    SEXP tag = PROTECT(Rf_allocVector(RAWSXP, 0));
    SEXP st = PROTECT(Rf_allocVector(INTSXP, 0));
    SEXP wd = PROTECT(Rf_allocVector(INTSXP, 0));
    SEXP ans = PROTECT(bs_new_xrawlist_single(classname, element_type, tag, st, wd));
    Rcpp::RObject result(ans);
    UNPROTECT(4);
    return result;
  }

  SEXP st = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP wd = PROTECT(Rf_allocVector(INTSXP, n));
  int* stp = INTEGER(st);
  int* wdp = INTEGER(wd);

  if (total < BS_TAG_CAP) {
    // Single shared buffer — identical() to Rsamtools::scanBam.
    SEXP tag = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(total)));
    if (total > 0U) std::memcpy(RAW(tag), col.bytes.data(), total);  // the only copy
    for (R_xlen_t i = 0; i < n; ++i) {
      const size_t idx = static_cast<size_t>(i);
      stp[i] = static_cast<int>(col.offsets[idx]) + 1;  // 1-based within the tag
      wdp[i] = static_cast<int>(col.lengths[idx]);
    }
    SEXP ans = PROTECT(bs_new_xrawlist_single(classname, element_type, tag, st, wd));
    Rcpp::RObject result(ans);
    UNPROTECT(4);  // st, wd, tag, ans
    return result;
  }

  // Multi-buffer: total decoded bytes exceed BS_TAG_CAP (> ~2 GB). Split at read
  // boundaries into pool buffers each < BS_TAG_CAP. Content-correct; no single-
  // buffer identical() reference exists (scanBam cannot return > 2 GB either).
  SEXP grp = PROTECT(Rf_allocVector(INTSXP, n));
  int* gp = INTEGER(grp);
  std::vector<size_t> tag_first_read;  // first read index of each pool buffer
  tag_first_read.push_back(0U);
  size_t cur = 0U;       // bytes accumulated in the current tag
  int gcur = 1;          // 1-based current tag index
  size_t tag_base = 0U;  // byte offset (in col.bytes) where the current tag starts
  for (R_xlen_t i = 0; i < n; ++i) {
    const size_t idx = static_cast<size_t>(i);
    const size_t w = col.lengths[idx];
    if (cur > 0U && cur + w > BS_TAG_CAP) {  // start a new tag at read i
      ++gcur;
      tag_base = col.offsets[idx];
      cur = 0U;
      tag_first_read.push_back(idx);
    }
    stp[i] = static_cast<int>(col.offsets[idx] - tag_base) + 1;  // within-tag 1-based
    wdp[i] = static_cast<int>(w);
    gp[i] = gcur;
    cur += w;
  }
  const int ntags = gcur;
  SEXP tags = PROTECT(Rf_allocVector(VECSXP, ntags));
  for (int t = 0; t < ntags; ++t) {
    const size_t r0 = tag_first_read[static_cast<size_t>(t)];
    const size_t b0 = col.offsets[r0];
    const size_t b1 = (static_cast<size_t>(t) + 1U < tag_first_read.size())
                          ? col.offsets[tag_first_read[static_cast<size_t>(t) + 1U]]
                          : total;
    const size_t sz = b1 - b0;
    SEXP tag = Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(sz));
    SET_VECTOR_ELT(tags, t, tag);
    if (sz > 0U) std::memcpy(RAW(tag), col.bytes.data() + b0, sz);
  }
  SEXP ans = PROTECT(bs_new_xrawlist(classname, element_type, tags, st, wd, grp));
  Rcpp::RObject result(ans);
  UNPROTECT(5);  // st, wd, grp, tags, ans
  return result;
}

// ---------------------------------------------------------------------------
// Typed tag accumulation.
//
// Rsamtools::scanBam returns numeric SAM tags with their native R type (e.g.
// NM:i as an integer vector, a float tag as a numeric vector), reserving
// character output for A/Z/H/B tags. BamScale mirrors that here: each tag column
// is accumulated as the widest type it needs, so a homogeneous numeric tag never
// pays for per-value string formatting or per-value CHARSXP allocation. Columns
// only fall back to character when a genuinely non-numeric value appears.
// ---------------------------------------------------------------------------

enum TagKind : uint8_t {
  TAG_EMPTY = 0,  // no concrete value observed yet (only absent/NA)
  TAG_INT   = 1,  // integer-valued numeric tag
  TAG_DBL   = 2,  // real-valued numeric tag (float seen, or int out of R range)
  TAG_STR   = 3   // character tag (A/Z/H/B), or numeric promoted after a string
};

// One probed value from a single read.
struct TagValue {
  uint8_t kind;       // 0 = absent/NA, 1 = numeric, 2 = string
  double num;         // valid when kind == 1
  bool is_integer;    // numeric came from an integer BAM type
  bool int_safe;      // numeric fits in an R integer
  std::string str;    // valid when kind == 2
};

inline std::string format_tag_num(const double v, const bool as_integer) {
  if (as_integer) {
    return std::to_string(static_cast<long long>(v));
  }
  std::ostringstream oss;
  oss << static_cast<float>(v);
  return oss.str();
}

TagValue probe_tag(pbam1_t& read, const std::string& tag) {
  TagValue tv;
  tv.kind = 0;
  tv.num = NA_REAL;
  tv.is_integer = true;
  tv.int_safe = true;

  const char type = read.Tag_Type(tag);
  switch (type) {
    case 'c': tv.kind = 1; tv.num = static_cast<double>(static_cast<int>(read.tagVal_c(tag))); break;
    case 'C': tv.kind = 1; tv.num = static_cast<double>(static_cast<unsigned int>(read.tagVal_C(tag))); break;
    case 's': tv.kind = 1; tv.num = static_cast<double>(static_cast<int>(read.tagVal_s(tag))); break;
    case 'S': tv.kind = 1; tv.num = static_cast<double>(static_cast<unsigned int>(read.tagVal_S(tag))); break;
    case 'i': tv.kind = 1; tv.num = static_cast<double>(read.tagVal_i(tag)); break;
    case 'I': {
      const uint32_t v = read.tagVal_I(tag);
      tv.kind = 1;
      tv.num = static_cast<double>(v);
      if (v > 2147483647u) tv.int_safe = false;
      break;
    }
    case 'f':
      tv.kind = 1;
      tv.is_integer = false;
      tv.int_safe = false;
      tv.num = static_cast<double>(read.tagVal_f(tag));
      break;
    case 'A':
    case 'Z':
    case 'H':
    case 'B':
      tv.kind = 2;
      tv.str = extract_tag_as_string(read, tag);
      break;
    default:
      tv.kind = 0;  // tag absent / unknown type -> NA
      break;
  }
  return tv;
}

struct TagColumn {
  uint8_t kind = TAG_EMPTY;
  size_t n = 0U;          // records represented (num/str size once kind is known)
  bool int_safe = true;   // all integer values fit in an R integer
  std::vector<double> num;
  PackedStringColumn str;

  void clear() {
    kind = TAG_EMPTY;
    n = 0U;
    int_safe = true;
    num.clear();
    str.clear();
  }

  void reserve(const size_t hint) {
    if (hint == 0U) return;
    if (kind == TAG_STR) {
      str.reserve(hint, 0U);
    } else {
      num.reserve(hint);
    }
  }

  // Convert already-accumulated entries to string form (numeric -> formatted,
  // deferred absents -> NA sentinel). Only hit on the rare mixed-type column.
  void promote_to_str() {
    if (kind == TAG_STR) return;
    PackedStringColumn s;
    s.reserve(n, 0U);
    if (kind == TAG_EMPTY) {
      for (size_t i = 0; i < n; ++i) s.append(TAG_NA_SENTINEL);
    } else {
      const bool as_integer = (kind == TAG_INT);
      for (size_t i = 0; i < num.size(); ++i) {
        if (R_IsNA(num[i])) {
          s.append(TAG_NA_SENTINEL);
        } else {
          s.append(format_tag_num(num[i], as_integer));
        }
      }
    }
    str = std::move(s);
    num.clear();
    num.shrink_to_fit();
    kind = TAG_STR;
  }

  // Bring the column up to a numeric kind, materializing deferred absents.
  void ensure_numeric(const uint8_t target) {
    if (kind == TAG_EMPTY) {
      num.assign(n, NA_REAL);
      kind = target;
    } else {
      kind = target;  // INT -> DBL widening; stored values are already double
    }
  }

  void append(const TagValue& tv) {
    if (tv.kind == 1) {           // numeric value
      if (kind == TAG_STR) {
        str.append(format_tag_num(tv.num, tv.is_integer));
      } else {
        const uint8_t want = tv.is_integer ? TAG_INT : TAG_DBL;
        if (kind == TAG_EMPTY) {
          ensure_numeric(want);
        } else if (want == TAG_DBL && kind == TAG_INT) {
          kind = TAG_DBL;         // widen INT -> DBL, never narrow
        }
        if (!tv.int_safe) int_safe = false;
        num.push_back(tv.num);
      }
    } else if (tv.kind == 2) {    // string value
      if (kind != TAG_STR) promote_to_str();
      str.append(tv.str);
    } else {                      // absent / NA
      if (kind == TAG_STR) {
        str.append(TAG_NA_SENTINEL);
      } else if (kind != TAG_EMPTY) {
        num.push_back(NA_REAL);
      }
      // TAG_EMPTY: defer; materialized when the kind is decided.
    }
    ++n;
  }

  // Concatenate another thread's column onto this one, widening as needed.
  void merge(TagColumn& other) {
    const uint8_t target = std::max(kind, other.kind);
    if (target == TAG_EMPTY) {
      n += other.n;
      return;
    }
    if (target == TAG_STR) {
      promote_to_str();
      other.promote_to_str();
      str.append_from(other.str);
    } else {
      // target is TAG_INT or TAG_DBL
      if (kind == TAG_EMPTY) ensure_numeric(target); else kind = target;
      if (other.kind == TAG_EMPTY) other.ensure_numeric(target);
      if (!other.int_safe) int_safe = false;
      num.insert(num.end(), other.num.begin(), other.num.end());
    }
    n += other.n;
  }

  RObject finalize() const {
    if (kind == TAG_STR) {
      return packed_strings_to_character(str);
    }
    if (kind == TAG_INT && int_safe) {
      IntegerVector v(static_cast<R_xlen_t>(n));
      for (size_t i = 0; i < n; ++i) {
        v[static_cast<R_xlen_t>(i)] = R_IsNA(num[i]) ? NA_INTEGER : static_cast<int>(num[i]);
      }
      return v;
    }
    if (kind == TAG_INT || kind == TAG_DBL) {
      return NumericVector(num.begin(), num.end());
    }
    // TAG_EMPTY: tag never present -> all-NA character column (prior behaviour).
    CharacterVector v(static_cast<R_xlen_t>(n));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(n); ++i) v[i] = NA_STRING;
    return v;
  }
};

// Build an R factor directly from BAM reference ids. Levels are the header
// chromosome names, so the integer ref id (+1) is already the factor code.
// Unmapped / unset references (id < 0) become <NA>, matching Rsamtools::scanBam
// (which returns rname/mrnm as factors with NA for absent references).
IntegerVector factor_from_ids(
    const std::vector<int>& ref_ids,
    const std::vector<std::string>& chr_names
) {
  const R_xlen_t n = static_cast<R_xlen_t>(ref_ids.size());
  const int n_levels = static_cast<int>(chr_names.size());
  IntegerVector out(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const int ref_id = ref_ids[static_cast<size_t>(i)];
    out[i] = (ref_id >= 0 && ref_id < n_levels) ? (ref_id + 1) : NA_INTEGER;
  }

  out.attr("levels") = wrap(chr_names);
  out.attr("class") = "factor";
  return out;
}

// Build an R factor from strand codes (0 = "+", 1 = "-", 2 = "*"), matching the
// strand factor levels used by Rsamtools::scanBam / GenomicAlignments.
IntegerVector strand_factor_from_codes(const std::vector<int>& strand_codes) {
  const R_xlen_t n = static_cast<R_xlen_t>(strand_codes.size());
  IntegerVector out(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    switch (strand_codes[static_cast<size_t>(i)]) {
      case 1:  out[i] = 2; break;  // "-"
      case 2:  out[i] = 3; break;  // "*"
      default: out[i] = 1; break;  // "+"
    }
  }

  out.attr("levels") = CharacterVector::create("+", "-", "*");
  out.attr("class") = "factor";
  return out;
}

struct QueryInterval {
  int start;
  int end;
  std::string label;
};

using IntervalMap = std::unordered_map<std::string, std::vector<QueryInterval> >;

IntervalMap build_interval_map(
    const CharacterVector& seqnames,
    const IntegerVector& starts,
    const IntegerVector& ends,
    const CharacterVector& labels
) {
  IntervalMap out;
  if (seqnames.size() == 0) return out;

  for (R_xlen_t i = 0; i < seqnames.size(); ++i) {
    if (CharacterVector::is_na(seqnames[i]) || starts[i] == NA_INTEGER || ends[i] == NA_INTEGER) {
      continue;
    }

    QueryInterval interval;
    interval.start = std::max(1, static_cast<int>(starts[i]));
    interval.end = std::max(interval.start, static_cast<int>(ends[i]));

    if (labels.size() == seqnames.size() && !CharacterVector::is_na(labels[i])) {
      interval.label = as<std::string>(labels[i]);
    } else {
      std::ostringstream oss;
      oss << as<std::string>(seqnames[i]) << ':' << interval.start << '-' << interval.end;
      interval.label = oss.str();
    }

    out[as<std::string>(seqnames[i])].push_back(interval);
  }

  return out;
}

bool match_interval(
    const IntervalMap& intervals,
    const std::string& seqname,
    const int pos,
    std::string& label
) {
  if (intervals.empty()) return true;

  const auto it = intervals.find(seqname);
  if (it == intervals.end()) return false;

  const std::vector<QueryInterval>& vec = it->second;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (pos >= vec[i].start && pos <= vec[i].end) {
      label = vec[i].label;
      return true;
    }
  }
  return false;
}

}  // namespace

// [[Rcpp::export]]
CharacterVector decode_compact_seq_cpp(const List& seq_raw, const IntegerVector& qwidth) {
  static const char nt16_map[16] = {
    '=', 'A', 'C', 'M',
    'G', 'R', 'S', 'V',
    'T', 'W', 'Y', 'H',
    'K', 'D', 'B', 'N'
  };

  const R_xlen_t n = seq_raw.size();
  if (qwidth.size() != n) {
    stop("`seq` and `qwidth` must have the same length");
  }

  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP empty = Rf_mkCharLen("", 0);

  for (R_xlen_t i = 0; i < n; ++i) {
    const int width = qwidth[i];

    if (width == NA_INTEGER) {
      SET_STRING_ELT(out, i, NA_STRING);
      continue;
    }
    if (width < 0) {
      UNPROTECT(1);
      stop("`qwidth` must be non-negative");
    }
    if (width == 0) {
      SET_STRING_ELT(out, i, empty);
      continue;
    }

    SEXP elt = seq_raw[i];
    if (TYPEOF(elt) != RAWSXP) {
      UNPROTECT(1);
      stop("Compact `seq` must be a list of raw vectors");
    }

    RawVector src(elt);
    const int needed = (width + 1) / 2;
    if (src.size() < needed) {
      UNPROTECT(1);
      stop("Compact `seq` entry is shorter than required by `qwidth`");
    }

    std::string decoded;
    decoded.resize(static_cast<size_t>(width));
    int j = 0;
    for (int k = 0; k < needed && j < width; ++k) {
      const uint8_t byte = src[static_cast<R_xlen_t>(k)];
      decoded[static_cast<size_t>(j++)] = nt16_map[(byte >> 4) & 0x0F];
      if (j < width) {
        decoded[static_cast<size_t>(j++)] = nt16_map[byte & 0x0F];
      }
    }

    SET_STRING_ELT(out, i, Rf_mkCharLen(decoded.data(), width));
  }

  UNPROTECT(1);
  return CharacterVector(out);
}

// [[Rcpp::export]]
CharacterVector decode_compact_qual_cpp(const List& qual_raw) {
  const R_xlen_t n = qual_raw.size();
  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP empty = Rf_mkCharLen("", 0);
  SEXP star = Rf_mkCharLen("*", 1);

  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP elt = qual_raw[i];
    if (TYPEOF(elt) != RAWSXP) {
      UNPROTECT(1);
      stop("Compact `qual` must be a list of raw vectors");
    }

    RawVector src(elt);
    const R_xlen_t len = src.size();

    if (len == 0) {
      SET_STRING_ELT(out, i, empty);
      continue;
    }

    bool missing = true;
    std::string decoded;
    decoded.resize(static_cast<size_t>(len));
    for (R_xlen_t j = 0; j < len; ++j) {
      const uint8_t q = src[j];
      if (q != 255U) {
        missing = false;
      }
      const uint8_t phred = (q == 255U) ? 0U : q;
      decoded[static_cast<size_t>(j)] = static_cast<char>(phred + 33U);
    }

    if (missing) {
      SET_STRING_ELT(out, i, star);
    } else {
      SET_STRING_ELT(out, i, Rf_mkCharLen(decoded.data(), static_cast<int>(len)));
    }
  }

  UNPROTECT(1);
  return CharacterVector(out);
}

// [[Rcpp::export]]
List read_bam_cpp(
    const std::string& bam_file,
    const int n_threads,
    const int min_mapq,
    const bool include_unmapped,
    const bool include_seq,
    const bool include_qual,
    const bool compact_seqqual,
    const int flag_require_set,
    const int flag_require_unset,
    const CharacterVector& tag_names_r,
    const CharacterVector& which_seqnames,
    const IntegerVector& which_starts,
    const IntegerVector& which_ends,
    const CharacterVector& which_labels,
    const bool with_which_label,
    const int field_mask
) {
  const unsigned int threads = clamp_threads(n_threads);
  const int min_mapq_clamped = std::max(0, min_mapq);
  const uint32_t flag_set = static_cast<uint32_t>(std::max(0, flag_require_set));
  const uint32_t flag_unset = static_cast<uint32_t>(std::max(0, flag_require_unset));

  const int mask = std::max(0, field_mask);
  const bool need_qname = (mask & (1 << 0)) != 0;
  const bool need_flag = (mask & (1 << 1)) != 0;
  const bool need_rname = (mask & (1 << 2)) != 0;
  const bool need_strand = (mask & (1 << 3)) != 0;
  const bool need_pos = (mask & (1 << 4)) != 0;
  const bool need_qwidth = (mask & (1 << 5)) != 0;
  const bool need_mapq = (mask & (1 << 6)) != 0;
  const bool need_cigar = (mask & (1 << 7)) != 0;
  const bool need_mrnm = (mask & (1 << 8)) != 0;
  const bool need_mpos = (mask & (1 << 9)) != 0;
  const bool need_isize = (mask & (1 << 10)) != 0;
  const bool need_seq = include_seq && ((mask & (1 << 11)) != 0);
  const bool need_qual = include_qual && ((mask & (1 << 12)) != 0);
  const bool need_seq_compact = need_seq && compact_seqqual;
  const bool need_qual_compact = need_qual && compact_seqqual;

  std::vector<std::string> tag_names;
  tag_names.reserve(static_cast<size_t>(tag_names_r.size()));
  for (R_xlen_t i = 0; i < tag_names_r.size(); ++i) {
    if (CharacterVector::is_na(tag_names_r[i])) continue;
    tag_names.push_back(as<std::string>(tag_names_r[i]));
  }

  const IntervalMap intervals = build_interval_map(
    which_seqnames,
    which_starts,
    which_ends,
    which_labels
  );
  const bool need_interval_match = !intervals.empty();

  pbam_in inbam;
  if (inbam.openFile(bam_file, threads) != 0) {
    stop("Failed to open BAM file");
  }

  std::vector<std::string> chr_names;
  std::vector<uint32_t> chr_lens;
  const int chrom_count = inbam.obtainChrs(chr_names, chr_lens);
  if (chrom_count <= 0) {
    stop("Failed to read BAM header");
  }

  PackedStringColumn qname_out;
  std::vector<int> flag_out;
  std::vector<int> rname_out;
  std::vector<int> pos_out;
  std::vector<int> strand_out;
  std::vector<int> qwidth_out;
  std::vector<int> mapq_out;
  PackedStringColumn cigar_out;
  std::vector<int> mrnm_out;
  std::vector<int> mpos_out;
  std::vector<int> isize_out;
  PackedStringColumn seq_out;
  PackedStringColumn qual_out;
  std::vector< std::vector<uint8_t> > seq_raw_out;
  std::vector< std::vector<uint8_t> > qual_raw_out;
  std::vector<std::string> which_label_out;
  std::vector<TagColumn> tag_out(tag_names.size());
  size_t n_records_total = 0U;

  struct alignas(64) ThreadChunk {
    PackedStringColumn qname;
    std::vector<int> flag;
    std::vector<int> rname;
    std::vector<int> pos;
    std::vector<int> strand;
    std::vector<int> qwidth;
    std::vector<int> mapq;
    PackedStringColumn cigar;
    std::vector<int> mrnm;
    std::vector<int> mpos;
    std::vector<int> isize;
    PackedStringColumn seq;
    PackedStringColumn qual;
    std::vector< std::vector<uint8_t> > seq_raw;
    std::vector< std::vector<uint8_t> > qual_raw;
    std::vector<std::string> which_label;
    std::vector<TagColumn> tag;
    size_t n_records = 0U;
  };

  std::vector<ThreadChunk> chunk_data(threads);
  for (unsigned int tid = 0; tid < threads; ++tid) {
    chunk_data[tid].tag.resize(tag_names.size());
  }

  const auto reserve_chunk = [&](ThreadChunk& local, const size_t hint) {
    if (hint == 0U) return;

    if (need_qname) local.qname.reserve(hint, 0U);
    if (need_flag) local.flag.reserve(hint);
    if (need_rname) local.rname.reserve(hint);
    if (need_pos) local.pos.reserve(hint);
    if (need_strand) local.strand.reserve(hint);
    if (need_qwidth) local.qwidth.reserve(hint);
    if (need_mapq) local.mapq.reserve(hint);
    if (need_cigar) local.cigar.reserve(hint, 0U);
    if (need_mrnm) local.mrnm.reserve(hint);
    if (need_mpos) local.mpos.reserve(hint);
    if (need_isize) local.isize.reserve(hint);
    if (need_seq) {
      if (need_seq_compact) {
        local.seq_raw.reserve(hint);
      } else {
        local.seq.reserve(hint, 0U);
      }
    }
    if (need_qual) {
      if (need_qual_compact) {
        local.qual_raw.reserve(hint);
      } else {
        local.qual.reserve(hint, 0U);
      }
    }
    if (with_which_label) local.which_label.reserve(hint);
    for (size_t j = 0; j < local.tag.size(); ++j) {
      local.tag[j].reserve(hint);
    }
  };

  while (true) {
    const int state = inbam.fillReads();
    if (state == 1) break;
    if (state == -1 || inbam.GetErrorState() == -1) {
      stop("BAM decompression failed while reading alignments");
    }

    for (unsigned int tid = 0; tid < threads; ++tid) {
      ThreadChunk& local = chunk_data[tid];
      local.n_records = 0U;
      local.qname.clear();
      local.flag.clear();
      local.rname.clear();
      local.pos.clear();
      local.strand.clear();
      local.qwidth.clear();
      local.mapq.clear();
      local.cigar.clear();
      local.mrnm.clear();
      local.mpos.clear();
      local.isize.clear();
      local.seq.clear();
      local.qual.clear();
      local.seq_raw.clear();
      local.qual_raw.clear();
      local.which_label.clear();
      for (size_t j = 0; j < local.tag.size(); ++j) {
        local.tag[j].clear();
      }

      const size_t bytes_hint = inbam.remainingThreadReadsBuffer(tid);
      const size_t record_hint = std::max(static_cast<size_t>(1U), bytes_hint / 128U);
      reserve_chunk(local, record_hint);
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
    for (unsigned int tid = 0; tid < threads; ++tid) {
      ThreadChunk& local = chunk_data[tid];
      std::string cigar_scratch;
      while (true) {
        pbam1_t read(inbam.supplyRead(tid));
        if (!read.validate()) break;

        const uint32_t flag = read.flag();
        if ((flag & flag_set) != flag_set) continue;
        if ((flag & flag_unset) != 0U) continue;

        int this_mapq = 0;
        if (min_mapq_clamped > 0 || need_mapq) {
          this_mapq = static_cast<int>(read.mapq());
          if (this_mapq < min_mapq_clamped) continue;
        }

        const int32_t ref_id = read.refID();
        const bool unmapped = (ref_id < 0) || ((flag & 0x4U) != 0U);
        if (!include_unmapped && unmapped) continue;

        int this_pos = NA_INTEGER;
        int this_ref_id = -1;
        if (!unmapped && ref_id >= 0 && ref_id < chrom_count && (need_rname || need_pos || need_interval_match)) {
          this_ref_id = static_cast<int>(ref_id);
          this_pos = static_cast<int>(read.pos()) + 1;
        }

        std::string this_which_label("*");
        if (need_interval_match) {
          if (this_ref_id < 0 || this_pos == NA_INTEGER) continue;
          const std::string& this_rname_ref = chr_names[static_cast<size_t>(this_ref_id)];
          if (!match_interval(intervals, this_rname_ref, this_pos, this_which_label)) {
            continue;
          }
        }

        local.n_records++;

        const uint32_t read_len = read.l_seq();

        if (need_qname) {
          const char* qname_ptr = read.read_name();
          const uint8_t qname_len = read.l_read_name();
          local.qname.append(
            qname_ptr,
            qname_len > 0 ? static_cast<size_t>(qname_len - 1) : 0U
          );
        }
        if (need_flag) {
          local.flag.push_back(static_cast<int>(flag));
        }
        if (need_rname) {
          local.rname.push_back(this_ref_id);
        }
        if (need_pos) {
          local.pos.push_back(this_pos);
        }
        if (need_strand) {
          if (unmapped) {
            local.strand.push_back(2);
          } else if ((flag & 0x10U) != 0U) {
            local.strand.push_back(1);
          } else {
            local.strand.push_back(0);
          }
        }
        if (need_qwidth) {
          local.qwidth.push_back(static_cast<int>(read_len));
        }
        if (need_mapq) {
          local.mapq.push_back(this_mapq);
        }
        if (need_cigar) {
          cigar_scratch.clear();
          read.cigar(cigar_scratch);
          local.cigar.append(cigar_scratch);
        }
        if (need_mrnm) {
          const int32_t next_ref_id = read.next_refID();
          local.mrnm.push_back(
            (next_ref_id >= 0 && next_ref_id < chrom_count) ? static_cast<int>(next_ref_id) : -1
          );
        }
        if (need_mpos) {
          const int32_t next_pos = read.next_pos();
          if (next_pos >= 0) {
            local.mpos.push_back(static_cast<int>(next_pos) + 1);
          } else {
            local.mpos.push_back(NA_INTEGER);
          }
        }
        if (need_isize) {
          local.isize.push_back(static_cast<int>(read.tlen()));
        }
        if (need_seq) {
          if (need_seq_compact) {
            local.seq_raw.emplace_back();
            const uint8_t* seq_ptr = read.seq();
            const size_t n_bytes = static_cast<size_t>((read_len + 1U) / 2U);
            if (seq_ptr != NULL && n_bytes > 0U) {
              local.seq_raw.back().assign(seq_ptr, seq_ptr + n_bytes);
            }
          } else {
            append_decoded_seq_encoded(local.seq, read.seq(), read_len);
          }
        }
        if (need_qual) {
          if (need_qual_compact) {
            local.qual_raw.emplace_back();
            const uint8_t* qual_ptr = reinterpret_cast<const uint8_t*>(read.qual());
            if (qual_ptr != NULL && read_len > 0U) {
              local.qual_raw.back().assign(
                qual_ptr,
                qual_ptr + static_cast<size_t>(read_len)
              );
            }
          } else {
            append_encoded_qual_native(
              local.qual,
              reinterpret_cast<const uint8_t*>(read.qual()),
              read_len
            );
          }
        }
        if (with_which_label) {
          local.which_label.push_back(this_which_label);
        }

        for (size_t j = 0; j < tag_names.size(); ++j) {
          local.tag[j].append(probe_tag(read, tag_names[j]));
        }
      }
    }

    size_t batch_records = 0U;
    for (unsigned int tid = 0; tid < threads; ++tid) {
      batch_records += chunk_data[tid].n_records;
    }

    if (batch_records > 0U) {
      const size_t target = n_records_total + batch_records;
      if (need_qname) {
        size_t batch_qname_bytes = 0U;
        for (unsigned int tid = 0; tid < threads; ++tid) {
          batch_qname_bytes += chunk_data[tid].qname.bytes_size();
        }
        qname_out.reserve(target, qname_out.bytes_size() + batch_qname_bytes);
      }
      if (need_flag) flag_out.reserve(target);
      if (need_rname) rname_out.reserve(target);
      if (need_pos) pos_out.reserve(target);
      if (need_strand) strand_out.reserve(target);
      if (need_qwidth) qwidth_out.reserve(target);
      if (need_mapq) mapq_out.reserve(target);
      if (need_cigar) {
        size_t batch_cigar_bytes = 0U;
        for (unsigned int tid = 0; tid < threads; ++tid) {
          batch_cigar_bytes += chunk_data[tid].cigar.bytes_size();
        }
        cigar_out.reserve(target, cigar_out.bytes_size() + batch_cigar_bytes);
      }
      if (need_mrnm) mrnm_out.reserve(target);
      if (need_mpos) mpos_out.reserve(target);
      if (need_isize) isize_out.reserve(target);
      if (need_seq) {
        if (need_seq_compact) {
          seq_raw_out.reserve(target);
        } else {
          size_t batch_seq_bytes = 0U;
          for (unsigned int tid = 0; tid < threads; ++tid) {
            batch_seq_bytes += chunk_data[tid].seq.bytes_size();
          }
          seq_out.reserve(target, seq_out.bytes_size() + batch_seq_bytes);
        }
      }
      if (need_qual) {
        if (need_qual_compact) {
          qual_raw_out.reserve(target);
        } else {
          size_t batch_qual_bytes = 0U;
          for (unsigned int tid = 0; tid < threads; ++tid) {
            batch_qual_bytes += chunk_data[tid].qual.bytes_size();
          }
          qual_out.reserve(target, qual_out.bytes_size() + batch_qual_bytes);
        }
      }
      if (with_which_label) which_label_out.reserve(target);
      for (size_t j = 0; j < tag_names.size(); ++j) {
        tag_out[j].reserve(target);
      }
    }

    for (unsigned int tid = 0; tid < threads; ++tid) {
      ThreadChunk& local = chunk_data[tid];
      n_records_total += local.n_records;

      if (need_qname) {
        qname_out.append_from(local.qname);
      }
      if (need_flag) {
        flag_out.insert(flag_out.end(), local.flag.begin(), local.flag.end());
      }
      if (need_rname) {
        rname_out.insert(
          rname_out.end(),
          std::make_move_iterator(local.rname.begin()),
          std::make_move_iterator(local.rname.end())
        );
      }
      if (need_pos) {
        pos_out.insert(pos_out.end(), local.pos.begin(), local.pos.end());
      }
      if (need_strand) {
        strand_out.insert(
          strand_out.end(),
          std::make_move_iterator(local.strand.begin()),
          std::make_move_iterator(local.strand.end())
        );
      }
      if (need_qwidth) {
        qwidth_out.insert(qwidth_out.end(), local.qwidth.begin(), local.qwidth.end());
      }
      if (need_mapq) {
        mapq_out.insert(mapq_out.end(), local.mapq.begin(), local.mapq.end());
      }
      if (need_cigar) {
        cigar_out.append_from(local.cigar);
      }
      if (need_mrnm) {
        mrnm_out.insert(
          mrnm_out.end(),
          std::make_move_iterator(local.mrnm.begin()),
          std::make_move_iterator(local.mrnm.end())
        );
      }
      if (need_mpos) {
        mpos_out.insert(mpos_out.end(), local.mpos.begin(), local.mpos.end());
      }
      if (need_isize) {
        isize_out.insert(isize_out.end(), local.isize.begin(), local.isize.end());
      }
      if (need_seq) {
        if (need_seq_compact) {
          seq_raw_out.insert(
            seq_raw_out.end(),
            std::make_move_iterator(local.seq_raw.begin()),
            std::make_move_iterator(local.seq_raw.end())
          );
        } else {
          seq_out.append_from(local.seq);
        }
      }
      if (need_qual) {
        if (need_qual_compact) {
          qual_raw_out.insert(
            qual_raw_out.end(),
            std::make_move_iterator(local.qual_raw.begin()),
            std::make_move_iterator(local.qual_raw.end())
          );
        } else {
          qual_out.append_from(local.qual);
        }
      }
      if (with_which_label) {
        which_label_out.insert(
          which_label_out.end(),
          std::make_move_iterator(local.which_label.begin()),
          std::make_move_iterator(local.which_label.end())
        );
      }

      for (size_t j = 0; j < tag_names.size(); ++j) {
        tag_out[j].merge(local.tag[j]);
      }
    }
  }

  inbam.closeFile();

  const int n = static_cast<int>(n_records_total);
  List out;
  if (need_qname) out["qname"] = packed_strings_to_character(qname_out);
  if (need_flag) out["flag"] = wrap(flag_out);
  if (need_rname) out["rname"] = factor_from_ids(rname_out, chr_names);
  if (need_pos) out["pos"] = wrap(pos_out);
  if (need_strand) out["strand"] = strand_factor_from_codes(strand_out);
  if (need_qwidth) out["qwidth"] = wrap(qwidth_out);
  if (need_mapq) out["mapq"] = wrap(mapq_out);
  if (need_cigar) out["cigar"] = packed_strings_to_character(cigar_out);
  if (need_mrnm) out["mrnm"] = factor_from_ids(mrnm_out, chr_names);
  if (need_mpos) out["mpos"] = wrap(mpos_out);
  if (need_isize) out["isize"] = wrap(isize_out);

  if (need_seq) {
    if (need_seq_compact) {
      out["seq"] = raw_list_from_bytes(seq_raw_out);
    } else {
      out["seq"] = xstringset_from_packed(seq_out, "DNAStringSet", "DNAString");
    }
  }
  if (need_qual) {
    if (need_qual_compact) {
      out["qual"] = raw_list_from_bytes(qual_raw_out);
    } else {
      out["qual"] = xstringset_from_packed(qual_out, "PhredQuality", "BString");
    }
  }
  if (with_which_label) {
    out["which_label"] = wrap(which_label_out);
  }

  for (size_t j = 0; j < tag_names.size(); ++j) {
    out[tag_names[j]] = tag_out[j].finalize();
  }

  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -n);
  out.attr("seqnames_header") = wrap(chr_names);
  out.attr("seqlengths_header") = wrap(chr_lens);
  out.attr("tag_na_sentinel") = TAG_NA_SENTINEL;
  out.attr("bam_file") = bam_file;

  return out;
}

// [[Rcpp::export]]
DataFrame count_bam_cpp(
    const std::string& bam_file,
    const int n_threads,
    const int min_mapq,
    const bool include_unmapped,
    const int flag_require_set,
    const int flag_require_unset,
    const CharacterVector& which_seqnames,
    const IntegerVector& which_starts,
    const IntegerVector& which_ends
) {
  const unsigned int threads = clamp_threads(n_threads);
  const int min_mapq_clamped = std::max(0, min_mapq);
  const uint32_t flag_set = static_cast<uint32_t>(std::max(0, flag_require_set));
  const uint32_t flag_unset = static_cast<uint32_t>(std::max(0, flag_require_unset));

  CharacterVector labels(which_seqnames.size());
  for (R_xlen_t i = 0; i < labels.size(); ++i) labels[i] = NA_STRING;
  const IntervalMap intervals = build_interval_map(
    which_seqnames,
    which_starts,
    which_ends,
    labels
  );

  pbam_in inbam;
  if (inbam.openFile(bam_file, threads) != 0) {
    stop("Failed to open BAM file");
  }

  std::vector<std::string> chr_names;
  std::vector<uint32_t> chr_lens;
  const int chrom_count = inbam.obtainChrs(chr_names, chr_lens);
  if (chrom_count <= 0) {
    stop("Failed to read BAM header");
  }

  std::vector<unsigned long long> total(static_cast<size_t>(chrom_count), 0ULL);
  unsigned long long total_unmapped = 0ULL;

  const size_t n_chr = static_cast<size_t>(chrom_count);
  const size_t n_cells = static_cast<size_t>(threads) * n_chr;
  std::vector<unsigned long long> local_counts(n_cells, 0ULL);
  std::vector<unsigned long long> local_unmapped(threads, 0ULL);

  while (true) {
    const int state = inbam.fillReads();
    if (state == 1) break;
    if (state == -1 || inbam.GetErrorState() == -1) {
      stop("BAM decompression failed while counting alignments");
    }

    std::fill(local_counts.begin(), local_counts.end(), 0ULL);
    std::fill(local_unmapped.begin(), local_unmapped.end(), 0ULL);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
    for (unsigned int tid = 0; tid < threads; ++tid) {
      unsigned long long* local = local_counts.data() + (static_cast<size_t>(tid) * n_chr);
      unsigned long long local_unmapped_count = 0ULL;

      while (true) {
        pbam1_t read(inbam.supplyRead(tid));
        if (!read.validate()) break;

        const uint32_t flag = read.flag();
        if ((flag & flag_set) != flag_set) continue;
        if ((flag & flag_unset) != 0U) continue;

        const int this_mapq = static_cast<int>(read.mapq());
        if (this_mapq < min_mapq_clamped) continue;

        const int32_t ref_id = read.refID();
        const bool unmapped = (ref_id < 0) || ((flag & 0x4U) != 0U);

        if (!intervals.empty()) {
          if (unmapped || ref_id < 0 || ref_id >= chrom_count) continue;
          const std::string& seq = chr_names[static_cast<size_t>(ref_id)];
          const int pos = static_cast<int>(read.pos()) + 1;
          std::string label;
          if (!match_interval(intervals, seq, pos, label)) continue;
        }

        if (unmapped) {
          if (include_unmapped) local_unmapped_count++;
        } else if (ref_id >= 0 && ref_id < chrom_count) {
          local[static_cast<size_t>(ref_id)]++;
        }
      }

      local_unmapped[tid] = local_unmapped_count;
    }

    for (unsigned int tid = 0; tid < threads; ++tid) {
      const unsigned long long* local = local_counts.data() + (static_cast<size_t>(tid) * n_chr);
      for (size_t i = 0; i < n_chr; ++i) {
        total[i] += local[i];
      }
      total_unmapped += local_unmapped[tid];
    }
  }

  inbam.closeFile();

  std::vector<std::string> seqname(chr_names.begin(), chr_names.end());
  std::vector<int> seqlength;
  seqlength.reserve(chr_lens.size() + (include_unmapped ? 1 : 0));
  for (size_t i = 0; i < chr_lens.size(); ++i) {
    seqlength.push_back(static_cast<int>(chr_lens[i]));
  }

  NumericVector count(static_cast<R_xlen_t>(total.size() + (include_unmapped ? 1 : 0)));
  for (size_t i = 0; i < total.size(); ++i) {
    count[static_cast<R_xlen_t>(i)] = static_cast<double>(total[i]);
  }

  if (include_unmapped) {
    seqname.push_back("*");
    seqlength.push_back(NA_INTEGER);
    count[count.size() - 1] = static_cast<double>(total_unmapped);
  }

  return DataFrame::create(
    _["seqname"] = seqname,
    _["seqlength"] = seqlength,
    _["count"] = count,
    _["stringsAsFactors"] = false
  );
}
