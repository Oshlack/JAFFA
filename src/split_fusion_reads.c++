// split_fusion_reads_mask.cpp  (C++11 compatible)
// Produce TWO full-length masked FASTA files so minimap2 keeps original query coordinates.
//
// Usage:
//   split_fusion_reads_mask <breakpoint_table> <reads_fasta> <buffer_bp> <left_out_fa> <right_out_fa>
//
// Breakpoint table format (whitespace/tab-separated):
//   col1: readID
//   col2: bp_left   (e.g. 293)  approximate breakpoint region start
//   col3: bp_right  (e.g. 296)  approximate breakpoint region end
//   col4+: ignored (gene etc.)
//   col5: read_length (optional; used only for sanity checks)
//
// Masking strategy (0-based, end-exclusive coordinates):
//   left_out_fa  : "only-left-maps"  -> mask RIGHT side [right_start, L) with 'N'
//   right_out_fa : "only-right-maps" -> mask LEFT  side [0, left_end)     with 'N'
//
// Where:
//   left_end    = clamp(bp_right + buffer, 0, L)
//   right_start = clamp(bp_left  - buffer, 0, L)
//
// Notes:
// - Reads not present in the breakpoint table are skipped.
// - If a breakpoint ID is not found in FASTA, it is reported at end.
// - FASTA header ID is taken as first token before whitespace.

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

struct BPInfo {
  long bp_left = -1;
  long bp_right = -1;
  long declared_len = -1; // optional
};

static inline std::string trim(const std::string& s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) b++;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) e--;
  return s.substr(b, e - b);
}

static inline long clamp_long(long v, long lo, long hi) {
  if (v < lo) return lo;
  if (v > hi) return hi;
  return v;
}

static bool parse_breakpoint_line(const std::string& line, std::string& read_id, BPInfo& info) {
  std::istringstream iss(line);
  std::string col1, col2, col3, col4, col5;
  if (!(iss >> col1 >> col2 >> col3)) return false;

  read_id = col1;

  try {
    info.bp_left = std::stol(col2);
    info.bp_right = std::stol(col3);
  } catch (...) {
    return false;
  }

  // Optional: col5 is declared length
  if (iss >> col4) {
    if (iss >> col5) {
      try { info.declared_len = std::stol(col5); } catch (...) { info.declared_len = -1; }
    }
  }

  if (info.bp_left > info.bp_right) std::swap(info.bp_left, info.bp_right);
  return true;
}

static std::unordered_map<std::string, BPInfo> read_breakpoints(const std::string& path) {
  std::unordered_map<std::string, BPInfo> bps;
  std::ifstream in(path.c_str());
  if (!in) throw std::runtime_error("Could not open breakpoint table: " + path);

  std::string line;
  size_t lineno = 0;
  while (std::getline(in, line)) {
    lineno++;
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;

    std::string read_id;
    BPInfo info;
    if (!parse_breakpoint_line(line, read_id, info)) {
      std::cerr << "Warning: could not parse breakpoint line " << lineno << " (skipping)\n";
      continue;
    }
    bps[read_id] = info;
  }
  return bps;
}

struct FastaRec {
  std::string header; // without '>'
  std::string seq;    // concatenated sequence
};

static bool read_next_fasta(std::istream& in, FastaRec& rec) {
  rec = FastaRec();
  std::string line;

  // Find next header
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '>') {
      rec.header = line.substr(1);
      break;
    }
  }
  if (rec.header.empty()) return false;

  // Read sequence lines until next header or EOF
  std::streampos pos;
  while (true) {
    pos = in.tellg();
    if (!std::getline(in, line)) break;
    if (!line.empty() && line[0] == '>') {
      in.seekg(pos); // rewind so next call sees this header
      break;
    }
    line = trim(line);
    if (!line.empty()) rec.seq += line;
  }
  return true;
}

static inline std::string header_id_token(const std::string& header) {
  size_t i = 0;
  while (i < header.size() && !std::isspace(static_cast<unsigned char>(header[i]))) i++;
  return header.substr(0, i);
}

static void write_fasta(std::ostream& out, const std::string& header, const std::string& seq, size_t wrap = 80) {
  out << ">" << header << "\n";
  for (size_t i = 0; i < seq.size(); i += wrap) {
    out << seq.substr(i, std::min(wrap, seq.size() - i)) << "\n";
  }
}

int main(int argc, char** argv) {
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0]
              << " <breakpoint_table> <reads_fasta> <buffer_bp> <left_out_fa> <right_out_fa>\n";
    return 1;
  }

  const std::string bp_path = argv[1];
  const std::string fa_path = argv[2];

  long buffer = 0;
  try { buffer = std::stol(argv[3]); } catch (...) { buffer = -1; }
  if (buffer < 0) {
    std::cerr << "Error: buffer_bp must be >= 0\n";
    return 1;
  }

  const std::string left_out = argv[4];
  const std::string right_out = argv[5];

  std::unordered_map<std::string, BPInfo> bps;
  try {
    bps = read_breakpoints(bp_path);
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  std::ifstream fin(fa_path.c_str());
  if (!fin) {
    std::cerr << "Error: could not open fasta: " << fa_path << "\n";
    return 1;
  }

  std::ofstream fout_left(left_out.c_str());
  std::ofstream fout_right(right_out.c_str());
  if (!fout_left) {
    std::cerr << "Error: could not write left output: " << left_out << "\n";
    return 1;
  }
  if (!fout_right) {
    std::cerr << "Error: could not write right output: " << right_out << "\n";
    return 1;
  }

  size_t n_total = 0, n_matched = 0, n_written = 0;
  FastaRec rec;

  while (read_next_fasta(fin, rec)) {
    n_total++;

    const std::string id = header_id_token(rec.header);
    std::unordered_map<std::string, BPInfo>::const_iterator it = bps.find(id);
    if (it == bps.end()) continue; // skip reads without breakpoints
    n_matched++;

    const BPInfo& info = it->second;
    const long L = static_cast<long>(rec.seq.size());

    if (info.declared_len > 0 && info.declared_len != L) {
      std::cerr << "Warning: length mismatch for " << id
                << " table_len=" << info.declared_len
                << " fasta_len=" << L << "\n";
    }

    // Define masked regions (0-based, end-exclusive)
    const long left_end = clamp_long(info.bp_left - buffer, 0L, L);      // [0, left_end)
    const long right_start = clamp_long(info.bp_right + buffer, 0L, L);    // [right_start, L)

    // Create full-length masked sequences:
    // - only_left : mask right side, so only left prefix maps
    // - only_right: mask left side, so only right suffix maps
    std::string only_left = rec.seq;
    if (right_start < L) {
      std::fill(only_left.begin() + static_cast<size_t>(right_start),
                only_left.end(),
                'N');
    }

    std::string only_right = rec.seq;
    if (left_end > 0) {
      std::fill(only_right.begin(),
                only_right.begin() + static_cast<size_t>(left_end),
                'N');
    }

    // Output headers (keep provenance)
    std::ostringstream hL, hR;
    hL << id << " side=left_only"
       << " bp_left=" << info.bp_left
       << " bp_right=" << info.bp_right
       << " buffer=" << buffer
       << " masked=" << right_start << "-" << L
       << " orig_len=" << L;

    hR << id << " side=right_only"
       << " bp_left=" << info.bp_left
       << " bp_right=" << info.bp_right
       << " buffer=" << buffer
       << " masked=0-" << left_end
       << " orig_len=" << L;

    write_fasta(fout_left, hL.str(), only_left);
    write_fasta(fout_right, hR.str(), only_right);
    n_written++;
  }

  std::cerr << "FASTA records read: " << n_total << "\n";
  std::cerr << "Reads with breakpoints: " << n_matched << "\n";
  std::cerr << "Reads written: " << n_written << "\n";

  if (n_matched < bps.size()) {
    std::cerr << "Note: " << (bps.size() - n_matched)
              << " breakpoint IDs were not found in the FASTA.\n";
  }

  return 0;
}
