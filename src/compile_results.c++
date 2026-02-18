/**********************************************************************
 * compile_summaries.cpp
 *
 * Port of R script that:
 *  - Collects JAFFAL summary files
 *  - Filters missing/empty files
 *  - Adds sample column
 *  - Renames/reorders columns
 *  - Sorts by:
 *      1) gap (ascending)
 *      2) supporting reads (descending)
 *      3) classification priority
 *  - Writes combined CSV
 *
 * C++11, STL-only
 * NOTE/WARNING: This code was auto-ported from R to C++ using AI (ChatGPT)
 *********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <sys/stat.h>

using namespace std;

/**********************************************************************
 * Utilities
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

static inline bool file_exists_and_nonempty(const string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return st.st_size > 0;
}

static inline double to_double(const string& s) {
    if (s.empty() || s == "-") return 0.0;
    return stod(s);
}

/**********************************************************************
 * Row container
 *********************************************************************/

struct Row {
    unordered_map<string,string> data;
};

/**********************************************************************
 * MAIN
 *********************************************************************/

int main(int argc, char* argv[]) {

    if (argc < 3) {
        cerr << "Usage: compile_summaries out_name summary1 summary2 ...\n";
        return 1;
    }

    string out_name = argv[1];

    vector<string> summary_files;
    vector<string> sample_names;

    // Extract sample name = second-to-last directory
    for (int i = 2; i < argc; ++i) {
        string full_path = argv[i];
        auto parts = split(full_path, '/');
        if (parts.size() >= 2) {
            string sample = parts[parts.size()-2];
            summary_files.push_back(full_path);
            sample_names.push_back(sample);
        }
    }

    cout << "Compiling the results from:\n";
    for (const auto& s : sample_names)
        cout << s << " ";
    cout << "\n";

    /******************************************************************
     * Filter missing/empty files
     ******************************************************************/
    vector<string> valid_files;
    vector<string> valid_samples;

    for (size_t i = 0; i < summary_files.size(); ++i) {
        if (file_exists_and_nonempty(summary_files[i])) {
            valid_files.push_back(summary_files[i]);
            valid_samples.push_back(sample_names[i]);
        } else {
            cout << "No fusions found for sample: "
                 << sample_names[i] << "\n";
        }
    }

    if (valid_files.empty()) {
        cerr << "No valid summary files found.\n";
        return 1;
    }

    /******************************************************************
     * Read and combine files
     ******************************************************************/
    vector<Row> full_list;

    for (size_t i = 0; i < valid_files.size(); ++i) {

        ifstream in(valid_files[i]);
        string header_line;
        getline(in, header_line);

        auto headers = split(header_line, '\t');

        string line;
        while (getline(in, line)) {
            if (line.empty()) continue;

            auto fields = split(line, '\t');
            if (fields.size() != headers.size()) continue;

            Row r;
            for (size_t j = 0; j < headers.size(); ++j) {
                string col = headers[j];

                // Replace "_" with space
                replace(col.begin(), col.end(), '_', ' ');

                if (col == "transcript")
                    col = "contig";
                if (col == "gap")
                    col = "gap (kb)";

                r.data[col] = fields[j];
            }

            r.data["sample"] = valid_samples[i];
            full_list.push_back(r);
        }
    }

    /******************************************************************
     * Column order (matches R)
     ******************************************************************/
    vector<string> column_order = {
        "sample","fusion genes","chrom1","base1","strand1",
        "chrom2","base2","strand2",
        "gap (kb)","spanning pairs","spanning reads",
        "inframe","aligns","rearrangement",
        "contig","contig break",
        "classification","geneCounts1","geneCounts2",
        "known mitelman","known cosmic",
        "cosmic tier","gtex samples"
    };

    /******************************************************************
     * Sorting
     ******************************************************************/

    // First sort by gap ascending
    stable_sort(full_list.begin(), full_list.end(),
        [](const Row& a, const Row& b) {

	const string& sa = a.data.at("gap (kb)");
        const string& sb = b.data.at("gap (kb)");

	if(sa == "inf") return false;
	if(sb == "inf") return true;
	// Otherwise numeric compare
	return stod(sa) < stod(sb);
	}
    );

    // Then sort by supporting reads descending
    stable_sort(full_list.begin(), full_list.end(),
        [](const Row& a, const Row& b) {

            double sp1 = to_double(a.data.at("spanning pairs"))
                       + to_double(a.data.at("spanning reads"));
            double sp2 = to_double(b.data.at("spanning pairs"))
                       + to_double(b.data.at("spanning reads"));

            return sp1 > sp2;
        }
    );

    // Classification priority (R rbind order)
    unordered_map<string,int> class_rank = {
        {"HighConfidence", 0},
        {"MediumConfidence", 1},
        {"LowConfidence", 2},
        {"PotentialTransSplicing", 3},
        {"PotentialReadThrough", 4}
    };

    stable_sort(full_list.begin(), full_list.end(),
        [&](const Row& a, const Row& b) {

            int ra = class_rank.count(a.data.at("classification"))
                     ? class_rank[a.data.at("classification")] : 100;
            int rb = class_rank.count(b.data.at("classification"))
                     ? class_rank[b.data.at("classification")] : 100;

            return ra < rb;
        }
    );

    /******************************************************************
     * Write CSV
     ******************************************************************/
    string output_file = out_name + ".csv";
    ofstream out(output_file);

    // Header
    for (size_t i = 0; i < column_order.size(); ++i) {
        out << column_order[i];
        if (i + 1 < column_order.size()) out << ",";
    }
    out << "\n";

    // Rows
    for (const auto& r : full_list) {
        for (size_t i = 0; i < column_order.size(); ++i) {
            auto it = r.data.find(column_order[i]);
            if (it != r.data.end())
                out << it->second;
            if (i + 1 < column_order.size()) out << ",";
        }
        out << "\n";
    }

    cout << "Done writing output " << output_file << "\n";

    return 0;
}
