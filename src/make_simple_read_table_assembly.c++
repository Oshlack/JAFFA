/**********************************************************************
 * compute_spanning.cpp
 *
 * Port of R script computing:
 *   - spanning_pairs
 *   - spanning_reads
 *
 * C++11, STL-only, standalone.
 * NOTE/WARNING: This file was auto-ported from R to C++ using AI (chatGPT)
 *********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>

using namespace std;

/**********************************************************************
 * Utility functions
 *********************************************************************/

static inline vector<string> split(const string& s, char delim) {
    vector<string> out;
    string tmp;
    for (char c : s) {
        if (c == delim) {
            out.push_back(tmp);
            tmp.clear();
        } else tmp.push_back(c);
    }
    out.push_back(tmp);
    return out;
}

static inline long long to_ll(const string& s) {
    return s.empty() ? 0 : stoll(s);
}

static inline double to_double(const string& s) {
    return s.empty() ? 0.0 : stod(s);
}

/**********************************************************************
 * Structures
 *********************************************************************/

struct Fusion {
    string transcript;
    long long break_min;
    long long break_max;
    string fusion_genes;
    long long contig_length;

    long long spanning_pairs = 0;
    long long spanning_reads = 0;
};

struct Alignment {
    string transcript;
    long long start;   // V2
    long long end;     // V3
    long long read_length; // from read_length_list
};

/**********************************************************************
 * MAIN
 *********************************************************************/

int main(int argc, char* argv[]) {

    if (argc != 7) {
        cerr << "Usage: compute_spanning sample candidates bam_data output read_length_list hang\n";
        return 1;
    }

    string sample              = argv[1];
    string candidates_file     = argv[2];
    string bam_file            = argv[3];
    string output_file         = argv[4];
    string read_length_file    = argv[5];
    long long hang             = stoll(argv[6]);

    cerr << "[INFO] Loading fusion candidates\n";

    /******************************************************************
     * Load fusion candidates
     ******************************************************************/
    vector<Fusion> fusions;
    {
        ifstream in(candidates_file);
        string line;
        while (getline(in, line)) {
            if (line.empty()) continue;
            auto fields = split(line, '\t');
            if (fields.size() < 5) continue;

            Fusion f;
            f.transcript     = fields[0];
            f.break_min      = to_ll(fields[1]);
            f.break_max      = to_ll(fields[2]);
            f.fusion_genes   = fields[3];
            f.contig_length  = to_ll(fields[4]);
            fusions.push_back(f);
        }
    }

    cerr << "[INFO] Loaded " << fusions.size() << " fusions\n";

    /******************************************************************
     * Load alignments
     * Equivalent to:
     * alignments_table = cbind(bam_data_file, read_length_list)
     ******************************************************************/
    cerr << "[INFO] Loading alignment data\n";

    vector<Alignment> alignments;
    {
        ifstream bam_in(bam_file);
        ifstream len_in(read_length_file);

        string line_bam, line_len;

        while (getline(bam_in, line_bam) &&
               getline(len_in, line_len)) {

            auto f_bam = split(line_bam, '\t');
            auto f_len = split(line_len, '\t');

            if (f_bam.size() < 3) continue;

            Alignment a;
            a.transcript = f_bam[0];
            a.start      = to_ll(f_bam[1]);
            a.end        = to_ll(f_bam[2]);
            a.read_length = to_ll(f_len[0]);

            alignments.push_back(a);
        }
    }

    cerr << "[INFO] Loaded " << alignments.size() << " alignments\n";

    /******************************************************************
     * Build transcript index (equivalent to split(..., V1))
     ******************************************************************/
    unordered_map<string, vector<int>> sal;
    for (size_t i = 0; i < alignments.size(); ++i) {
        sal[alignments[i].transcript].push_back((int)i);
    }

    cerr << "[INFO] Computing spanning metrics\n";

    /******************************************************************
     * Compute spanning_pairs and spanning_reads
     ******************************************************************/
    for (auto& fusion : fusions) {

        auto it = sal.find(fusion.transcript);
        if (it == sal.end()) {
            fusion.spanning_pairs = 0;
            fusion.spanning_reads = 0;
            continue;
        }

        const vector<int>& idxs = it->second;

        // Compute mean read length
        double mean_length = 0.0;
        for (int i : idxs)
            mean_length += alignments[i].read_length;

        if (!idxs.empty())
            mean_length /= idxs.size();

        if (mean_length == 0) {
            fusion.spanning_pairs = 0;
            fusion.spanning_reads = 0;
            continue;
        }

        long long pos = fusion.break_min;
        long long contig_length = fusion.contig_length;

        // Short contig filter (exact R logic)
        if ((pos < mean_length * 2) ||
            ((contig_length - pos) < mean_length * 2)) {

            fusion.spanning_pairs = 0;
        } else {

            long long count_pairs = 0;
            for (int i : idxs) {
                const Alignment& a = alignments[i];

                bool this_below = a.start < (pos - a.read_length);
                bool pair_above = a.end > pos;

                if (this_below && pair_above)
                    count_pairs++;
            }

            fusion.spanning_pairs = count_pairs;
        }

        // Spanning reads (R logic exactly)
        long long count_reads = 0;
        for (int i : idxs) {
            const Alignment& a = alignments[i];

            bool below_fusion = a.start > (pos + hang - a.read_length);
            bool above_fusion = a.start < (pos - hang);

            if (below_fusion && above_fusion)
                count_reads++;
        }

        fusion.spanning_reads = count_reads;
    }

    /******************************************************************
     * Write output
     ******************************************************************/
    cerr << "[INFO] Writing output\n";

    ofstream out(output_file);
    out << "transcript\tbreak_min\tbreak_max\tfusion_genes\tspanning_pairs\tspanning_reads" << endl;
    for (const auto& f : fusions) {
        out << f.transcript << "\t"
            << f.break_min << "\t"
            << f.break_max << "\t"
            << f.fusion_genes << "\t"
            << f.spanning_pairs << "\t"
            << f.spanning_reads << "\n";
    }

    cerr << "[INFO] Done\n";

    return 0;
}
