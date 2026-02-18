#include <Rcpp.h>
using namespace Rcpp;

// Required to print ompBAM logging to the R console.
#define cout Rcpp::Rcout

// [[Rcpp::depends(ompBAM)]]
#include <ompBAM.hpp>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

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

inline void decode_seq_ascii(const uint8_t* packed_seq, const uint32_t n_bases, std::string& out) {
  static const char nt16_map[16] = {
    '=', 'A', 'C', 'M',
    'G', 'R', 'S', 'V',
    'T', 'W', 'Y', 'H',
    'K', 'D', 'B', 'N'
  };

  out.clear();
  if (n_bases == 0U || packed_seq == NULL) return;

  out.resize(static_cast<size_t>(n_bases));
  for (uint32_t i = 0; i < n_bases; ++i) {
    const uint8_t byte = packed_seq[i >> 1];
    const uint8_t code = (i & 1U) ? static_cast<uint8_t>(byte & 0x0FU)
                                  : static_cast<uint8_t>((byte >> 4) & 0x0FU);
    out[static_cast<size_t>(i)] = nt16_map[code];
  }
}

inline void encode_qual_ascii(const uint8_t* qual_ptr, const uint32_t n_bases, std::string& out) {
  out.clear();
  if (n_bases == 0U || qual_ptr == NULL) return;

  bool missing = true;
  out.resize(static_cast<size_t>(n_bases));
  for (uint32_t i = 0; i < n_bases; ++i) {
    const uint8_t q = qual_ptr[i];
    if (q != 255U) missing = false;
    const uint8_t phred = (q == 255U) ? 0U : q;
    out[static_cast<size_t>(i)] = static_cast<char>(phred + 33U);
  }

  if (missing) out.assign("*");
}

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
List read_bam_cpp(
    const std::string& bam_file,
    const int n_threads,
    const int min_mapq,
    const bool include_unmapped,
    const bool include_seq,
    const bool include_qual,
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

  std::vector<std::string> qname_out;
  std::vector<int> flag_out;
  std::vector<std::string> rname_out;
  std::vector<int> pos_out;
  std::vector<std::string> strand_out;
  std::vector<int> qwidth_out;
  std::vector<int> mapq_out;
  std::vector<std::string> cigar_out;
  std::vector<std::string> mrnm_out;
  std::vector<int> mpos_out;
  std::vector<int> isize_out;
  std::vector<std::string> seq_out;
  std::vector<std::string> qual_out;
  std::vector<std::string> which_label_out;
  std::vector< std::vector<std::string> > tag_out(tag_names.size());
  size_t n_records_total = 0U;

  struct ThreadChunk {
    std::vector<std::string> qname;
    std::vector<int> flag;
    std::vector<std::string> rname;
    std::vector<int> pos;
    std::vector<std::string> strand;
    std::vector<int> qwidth;
    std::vector<int> mapq;
    std::vector<std::string> cigar;
    std::vector<std::string> mrnm;
    std::vector<int> mpos;
    std::vector<int> isize;
    std::vector<std::string> seq;
    std::vector<std::string> qual;
    std::vector<std::string> which_label;
    std::vector< std::vector<std::string> > tag;
    size_t n_records = 0U;
  };

  std::vector<ThreadChunk> chunk_data(threads);
  for (unsigned int tid = 0; tid < threads; ++tid) {
    chunk_data[tid].tag.resize(tag_names.size());
  }

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
      local.which_label.clear();
      for (size_t j = 0; j < local.tag.size(); ++j) {
        local.tag[j].clear();
      }
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
    for (unsigned int tid = 0; tid < threads; ++tid) {
      ThreadChunk& local = chunk_data[tid];
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

        std::string this_rname("*");
        int this_pos = NA_INTEGER;
        if (!unmapped && ref_id >= 0 && ref_id < chrom_count && (need_rname || need_pos || need_interval_match)) {
          this_rname = chr_names.at(static_cast<size_t>(ref_id));
          this_pos = static_cast<int>(read.pos()) + 1;
        }

        std::string this_which_label("*");
        if (need_interval_match) {
          if (unmapped || this_pos == NA_INTEGER) continue;
          if (!match_interval(intervals, this_rname, this_pos, this_which_label)) {
            continue;
          }
        }

        local.n_records++;

        const uint32_t read_len = read.l_seq();

        if (need_qname) {
          const char* qname_ptr = read.read_name();
          const uint8_t qname_len = read.l_read_name();
          local.qname.emplace_back(
            qname_ptr,
            qname_ptr + (qname_len > 0 ? static_cast<size_t>(qname_len - 1) : 0U)
          );
        }
        if (need_flag) {
          local.flag.push_back(static_cast<int>(flag));
        }
        if (need_rname) {
          local.rname.push_back(this_rname);
        }
        if (need_pos) {
          local.pos.push_back(this_pos);
        }
        if (need_strand) {
          if (unmapped) {
            local.strand.push_back("*");
          } else if ((flag & 0x10U) != 0U) {
            local.strand.push_back("-");
          } else {
            local.strand.push_back("+");
          }
        }
        if (need_qwidth) {
          local.qwidth.push_back(static_cast<int>(read_len));
        }
        if (need_mapq) {
          local.mapq.push_back(this_mapq);
        }
        if (need_cigar) {
          local.cigar.emplace_back();
          read.cigar(local.cigar.back());
        }
        if (need_mrnm) {
          const int32_t next_ref_id = read.next_refID();
          if (next_ref_id >= 0 && next_ref_id < chrom_count) {
            local.mrnm.push_back(chr_names.at(static_cast<size_t>(next_ref_id)));
          } else {
            local.mrnm.push_back("*");
          }
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
          local.seq.emplace_back();
          decode_seq_ascii(read.seq(), read_len, local.seq.back());
        }
        if (need_qual) {
          local.qual.emplace_back();
          encode_qual_ascii(
            reinterpret_cast<const uint8_t*>(read.qual()),
            read_len,
            local.qual.back()
          );
        }
        if (with_which_label) {
          local.which_label.push_back(this_which_label);
        }

        for (size_t j = 0; j < tag_names.size(); ++j) {
          local.tag[j].push_back(extract_tag_as_string(read, tag_names[j]));
        }
      }
    }

    for (unsigned int tid = 0; tid < threads; ++tid) {
      ThreadChunk& local = chunk_data[tid];
      n_records_total += local.n_records;

      if (need_qname) {
        qname_out.insert(
          qname_out.end(),
          std::make_move_iterator(local.qname.begin()),
          std::make_move_iterator(local.qname.end())
        );
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
        cigar_out.insert(
          cigar_out.end(),
          std::make_move_iterator(local.cigar.begin()),
          std::make_move_iterator(local.cigar.end())
        );
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
        seq_out.insert(
          seq_out.end(),
          std::make_move_iterator(local.seq.begin()),
          std::make_move_iterator(local.seq.end())
        );
      }
      if (need_qual) {
        qual_out.insert(
          qual_out.end(),
          std::make_move_iterator(local.qual.begin()),
          std::make_move_iterator(local.qual.end())
        );
      }
      if (with_which_label) {
        which_label_out.insert(
          which_label_out.end(),
          std::make_move_iterator(local.which_label.begin()),
          std::make_move_iterator(local.which_label.end())
        );
      }

      for (size_t j = 0; j < tag_names.size(); ++j) {
        tag_out[j].insert(
          tag_out[j].end(),
          std::make_move_iterator(local.tag[j].begin()),
          std::make_move_iterator(local.tag[j].end())
        );
      }
    }
  }

  inbam.closeFile();

  const int n = static_cast<int>(n_records_total);
  List out;
  if (need_qname) out["qname"] = wrap(qname_out);
  if (need_flag) out["flag"] = wrap(flag_out);
  if (need_rname) out["rname"] = wrap(rname_out);
  if (need_pos) out["pos"] = wrap(pos_out);
  if (need_strand) out["strand"] = wrap(strand_out);
  if (need_qwidth) out["qwidth"] = wrap(qwidth_out);
  if (need_mapq) out["mapq"] = wrap(mapq_out);
  if (need_cigar) out["cigar"] = wrap(cigar_out);
  if (need_mrnm) out["mrnm"] = wrap(mrnm_out);
  if (need_mpos) out["mpos"] = wrap(mpos_out);
  if (need_isize) out["isize"] = wrap(isize_out);

  if (need_seq) {
    out["seq"] = wrap(seq_out);
  }
  if (need_qual) {
    out["qual"] = wrap(qual_out);
  }
  if (with_which_label) {
    out["which_label"] = wrap(which_label_out);
  }

  for (size_t j = 0; j < tag_names.size(); ++j) {
    out[tag_names[j]] = wrap(tag_out[j]);
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

  while (true) {
    const int state = inbam.fillReads();
    if (state == 1) break;
    if (state == -1 || inbam.GetErrorState() == -1) {
      stop("BAM decompression failed while counting alignments");
    }

    std::vector< std::vector<unsigned long long> > local_counts(
      threads,
      std::vector<unsigned long long>(static_cast<size_t>(chrom_count), 0ULL)
    );
    std::vector<unsigned long long> local_unmapped(threads, 0ULL);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
    for (unsigned int tid = 0; tid < threads; ++tid) {
      std::vector<unsigned long long>& local = local_counts[tid];
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
          const std::string seq = chr_names.at(static_cast<size_t>(ref_id));
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
      const std::vector<unsigned long long>& local = local_counts[tid];
      for (size_t i = 0; i < local.size(); ++i) {
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
