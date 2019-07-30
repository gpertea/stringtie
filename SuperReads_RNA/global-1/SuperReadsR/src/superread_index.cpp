#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>

#include "timer.h"
#include "superread_index.hpp"

//Minimum weight each read is given in coverage computation
#define BASE_WEIGHT 1.0

//Delimiters for parsing files
const std::string SR_DELIM = "_\n"; 
const std::string READ_DELIM = " \t\n";

Read::Read(std::string read_str) {

    char *cstr = new char[read_str.length()+1];
    strcpy(cstr, read_str.c_str());

    //Make sure line contains anything    
    if (read_str.length() == 0) {
        is_valid = false;
        return;
    }

    //Determine how much space to reserve for k-unitigs
    int delim_ct = 0;
    for (unsigned int i = 0; i < read_str.length(); i++) {
        for (unsigned int j = 0; j < READ_DELIM.size(); j++) {
            if (cstr[i] == READ_DELIM[j]) {
                delim_ct++;
                break;
            }
        }
    }

    ku_ids.reserve((delim_ct-1)/3);
    ku_strands.reserve(ku_ids.capacity());
    
    char *read_name = strtok(cstr, READ_DELIM.c_str());

    //Make sure read name is long enough
    if (strlen(read_name) < PREFIX_LEN+1) {
        is_valid = false;
        return;
    }

    read_id = atoi(&read_name[PREFIX_LEN]);

    length = atoi(strtok(NULL, READ_DELIM.c_str()));
    
    char *ku_id;

    is_valid = false;

    while ( (ku_id = strtok(NULL, READ_DELIM.c_str())) != NULL ) {
        ku_ids.push_back(atoi(ku_id));

        //Store last ku position
        last_ku_pos = atoi(strtok(NULL, READ_DELIM.c_str()));

        ku_strands.push_back(*strtok(NULL, READ_DELIM.c_str()));

        //Read must contain at least on k-unitig
        //Read was read by the read reader
        is_valid = true;
    }

    delete[] cstr;
}

//Copy constructor
Read::Read(const Read &r) {
    read_id = r.read_id;
    ku_ids = r.ku_ids;
    ku_strands = r.ku_strands;
    last_ku_pos = r.last_ku_pos;
}

bool Read::is_mate(Read r) {
    if (read_id % 2 == 0)
        return read_id+1 == r.read_id;
    return read_id-1 == r.read_id;
}

IndexQuery::IndexQuery(int buffer_size) {
    reads.reserve(buffer_size);
    results = new std::vector< std::vector<int>* >();    
    matched_reads = new std::vector<Read*>();
    finished = true;
}

IndexQuery::~IndexQuery() {
    delete results;
}

bool IndexQuery::add_read(Read *read) {
    if (is_full())
        return false;

    reads.push_back(read);

    return true;
}

bool IndexQuery::is_full() {
    return reads.size() >= reads.capacity()-1;
}

bool IndexQuery::is_empty() {
    return reads.size() == 0;
}

void IndexQuery::clear() {
    for (unsigned int r = 0; r < results->size(); r++)
        delete results->at(r);                                    
    results->clear();
    reads.clear();
}


ReadIO::ReadIO(std::string read_ku_fname, 
               char *sr_matches_fname) {
    read_ku_file.open(read_ku_fname.c_str());

    //Open file to output matches if filename provided
    if (sr_matches_fname != NULL)
        sr_matches_file.open(sr_matches_fname);
}

ReadIO::~ReadIO() {
}

bool ReadIO::fill_buffer(IndexQuery *q) {
    std::string read_ku_line;

    Read *read;
        
    //Add more reads until full or EOF
    while (!(q->is_full()) && getline(read_ku_file, read_ku_line)) {
        read = new Read(read_ku_line);

        if (!read->is_valid) {
            continue;
        }

        if (!q->add_read(read)) {
            break;
        }
    }
   
    //No reads added - must have reached EOF
    if (q->is_empty())
        return false;

    q->finished = false;

    return true;
}

void ReadIO::update_matched(IndexQuery *q) {
    
    for (unsigned int i = 0; i < q->results->size(); i++) {
        Read *read = q->matched_reads->at(i);
        
        if (sr_matches_file.is_open()) {
            sr_matches_file << read->read_id;

            std::vector<int> &sr_list = *(q->results->at(i));
            for (unsigned int j = 0; j < sr_list.size(); j++) {
                sr_matches_file << "\t" << SuperreadIndex::sr_id_lens[sr_list[j]].first;
            }

            sr_matches_file << "\n";
        }
    }

}

bool ReadIO::copy_fastq(std::ifstream &in, std::ofstream &out) {
    std::string s;

    if (!getline(in, s))
        return false;
    
    out << s << "\n";

    for (unsigned int i = 0; i < 3; i++) {
        getline(in, s);
        out << s << "\n";
    }
    
    return true;
}

bool ReadIO::skip_fastq(std::ifstream &in) {
    bool not_eof;
    std::string s;

    for (unsigned int i = 0; i < 4; i++) 
        not_eof = (bool)getline(in, s);
    
    return not_eof;
}

//Declare static variables
std::vector< std::vector<SuperreadIndex::KUnitig> > SuperreadIndex::ku_list;
std::map<int, int> SuperreadIndex::sr_indices;
std::vector< std::pair<int, int> > SuperreadIndex::sr_id_lens;
std::map< std::vector<int>, int > SuperreadIndex::sr_frag_counts;
std::vector<double> SuperreadIndex::sr_covs;
std::string SuperreadIndex::sr_sam_fname;

bool SuperreadIndex::load_files(std::string sr_name_path, std::string sr_sam_path) {
    
    //Save fasta filename for later output
    sr_sam_fname = sr_sam_path;

    //Only want to index SRs present in sr_seq_file
    //SR KUs read from sr_name_file
    std::ifstream sr_name_file(sr_name_path.c_str()),
                  sr_sam_file(sr_sam_path.c_str());

    
    std::string sr_name_line, sr_sam_line, ku_fasta_line;

    if (!sr_name_file.is_open()) {
        std::cerr << "Error: unable to open sr file '" << sr_name_path << "'\nExiting\n";
        return false;
    } else if (!sr_sam_file.is_open()) {
        std::cerr << "Error: unable to open sam file '" << sr_sam_path << "'\nExiting\n";
        return false;
    }


    std::cerr << "Parsing SAM file...\n";

    //Iterate thru sr_seq file to find SRs to store
    int sr_id, flag;
    while (getline(sr_sam_file, sr_sam_line)) {
        if (sr_sam_line[0] == '@') {
            continue;
        }

        int i = sr_sam_line.find('\t'), j = sr_sam_line.find('\t', i+1);

        sr_id = atoi(sr_sam_line.substr(0, i).c_str());
        i++;
        flag = atoi(sr_sam_line.substr(i, j-i).c_str());

        if (!(flag & 0x4) & !(flag & 0x900)) {
            for (int k = 0; k < 8; k++) {
                i = j+1;
                j = sr_sam_line.find('\t', i);
            }
            sr_id_lens.push_back(std::pair<int, int>(sr_id, j - i));
        }
    }

    std::sort(sr_id_lens.begin(), sr_id_lens.end());

    std::cerr << "Removing duplicates...\n";
    unsigned int i = 1, j = 1;
    for (; j < sr_id_lens.size(); j++) {
        if (sr_id_lens[j].first != sr_id_lens[j-1].first) {
            sr_id_lens[i++] = sr_id_lens[j];
        }
    }

    while (i < sr_id_lens.size()) {
        sr_id_lens.pop_back();
    }

    std::cerr << "Parsing SR file...\n";
    unsigned int sr_i = 0;
    sr_id = 0;
    while (getline(sr_name_file, sr_name_line) && sr_i < sr_id_lens.size()) {
        if (sr_id == sr_id_lens[sr_i].first) {
            sr_indices[sr_id] = sr_i;
            add_superread_kus(sr_i, sr_name_line);
            sr_i++;
        } //else std::cout << sr_i << " " << sr_id << " "  << sr_id_lens[sr_i].first << "\n";
        sr_id++;
    }

    sr_covs = std::vector<double>(sr_id_lens.size(), 0.0);

    sr_name_file.close();
    sr_sam_file.close();

    return true;
}
 
void SuperreadIndex::destroy() {
    for (unsigned int i = 0; i < ku_list.size(); i++) {
        for (unsigned int j = 0; j < ku_list[i].size(); j++) {
            delete ku_list[i][j];
        }
        ku_list[i].clear();
    }
    ku_list.clear();
}


void SuperreadIndex::update_frag_counts(IndexQuery *q) {
    for (unsigned int i = 0; i < q->results->size(); i++) {
        std::vector<int> &sr_list = *(q->results->at(i));
        
        //Get the read length
        int read_len;

        read_len = q->matched_reads->at(i)->length;

        //Add read length to fragment counts
        if (sr_frag_counts.count(sr_list) == 0) {
            sr_frag_counts[sr_list] = read_len;
        } else {
            sr_frag_counts[sr_list] += read_len;
        }


        //Update initial super-read coverages
        if (sr_list.size() == 1) {
            sr_covs[sr_list[0]] += (double) read_len / sr_id_lens[sr_list[0]].second;
        }
    }
}

double SuperreadIndex::iter_sr_covs() {
    int count;
    double weight_sum;
    std::vector<double> new_sr_covs(sr_covs.size(), 0.0);
    
    //Iterate thru each fragment overlap
    std::map<std::vector<int>, int>::iterator frag_iter = sr_frag_counts.begin();
    for (; frag_iter != sr_frag_counts.end(); ++frag_iter) {
        const std::vector<int> &sr_list = frag_iter->first;
        count = frag_iter->second;
        
        //Calc denominator
        weight_sum = 0;
        for (unsigned int i = 0; i < sr_list.size(); i++)
            weight_sum += BASE_WEIGHT + sr_covs[sr_list[i]];
        
        //Add to total read length
        for (unsigned int i = 0; i < sr_list.size(); i++) {
            new_sr_covs[sr_list[i]] += 
                    count * (BASE_WEIGHT + sr_covs[sr_list[i]]) / weight_sum;

        }
    }
    
    double diff, max_diff = 0;
    for (unsigned int i = 0; i < sr_covs.size(); i++) {
        
        //Get final coverage
        new_sr_covs[i] /= sr_id_lens[i].second;
        
        //Get absolute difference between estimates
        diff = sr_covs[i]>new_sr_covs[i] ? 
                  sr_covs[i]-new_sr_covs[i] 
                : new_sr_covs[i]-sr_covs[i];
        
        //Update max difference
        if (diff > max_diff)
            max_diff = diff;
    }
    
    //Set new coverage
    sr_covs.swap(new_sr_covs);

    return max_diff;
}

void SuperreadIndex::write_sam(std::string fname, std::string delim) {
    std::ifstream sr_sam_in(sr_sam_fname.c_str());
    std::ofstream sr_sam_out(fname.c_str());

    std::string sam_line;

    while (getline(sr_sam_in, sam_line)) {

        if (sam_line[0] == '@') {
            sr_sam_out << sam_line << std::endl;
            continue;
        }

        int sr_id = atoi(sam_line.substr(0, sam_line.find('\t')).c_str());

        std::map<int, int>::iterator i = sr_indices.find(sr_id);


        if (i != sr_indices.end() && sr_covs[i->second] > 0) {
            sr_sam_out << sam_line << "\t" << delim << sr_covs[i->second] << std::endl;
        }
    }
}

//Private functions
void SuperreadIndex::add_superread_kus(int sr_index, std::string sr_str) {
    char *cstr = new char[sr_str.length()+1];
    strcpy(cstr, sr_str.c_str());

    KUnitig ku;

    char *ku_str;   //Stores each std::string representing a KU (ex 12345F)
    unsigned int ku_id,      //Stores int of KU (ex 12345)
                 ku_len,     //Stores length of KU std::string
                 sr_idx = 0;

    ku_str = strtok(cstr, SR_DELIM.c_str());

    do {

        ku_len = strlen(ku_str);

        //Create new K-Unitig record
        ku = new struct kunitig;
        ku->sr = sr_index;
        ku->strand = ku_str[ku_len-1]; //Strand is last char
        ku->sr_idx = sr_idx;

        ku_str[ku_len] = '\0'; //Remove strand char
        ku_id = atoi(ku_str);

        //Make room for KU (if needed) and store it
        if (ku_id > ku_list.size())
            ku_list.resize(ku_id+1);
        ku_list[ku_id].push_back(ku);

        
        sr_idx++;

    } while ( (ku_str = strtok(NULL, SR_DELIM.c_str())) != NULL );

    delete[] cstr;
}


int SuperreadIndex::match_reads(IndexQuery *q) { 
    
    #ifdef DEBUG_TIMER
    Timer timer;
    #endif

    std::vector<int> *read_srs;

    Read *read;
    for (unsigned int r = 0; r < q->reads.size(); r += 1) {
        read = q->reads[r];
        read_srs = new std::vector<int>(match_read(read));
        
        if (read_srs->size() > 0) {
            q->results->push_back(read_srs);
            q->matched_reads->push_back(read);
        } else {
            delete read_srs;
            delete read;
        }

    }

    q->reads.clear();
    
    #ifdef DEBUG_TIMER
    q->thread_time = timer.get();
    #endif
     
    q->finished = true;

    return 0;
}

std::vector<int> SuperreadIndex::match_read(Read *read) {
    std::vector<KUnitig> prev_kus, next_kus;

    prev_kus = ku_list[read->ku_ids[0]]; //Stores previously matched SR KUs
    next_kus.reserve(prev_kus.size());   //Stores next matching of SR KUs
    
    //Loops thru every read KU
    for (unsigned int ku = 1; ku < read->ku_ids.size(); ku++) {

        //Get SRs containing current KU
        std::vector<KUnitig> &sr_kus = ku_list[read->ku_ids[ku]];

        //Match up SR KUs to previously matched KUs
        unsigned int i = 0, j = 0, k;
        while (i < prev_kus.size() && j < sr_kus.size()) {

            //Previous SR doesn't contain current KU
            if (prev_kus[i]->sr < sr_kus[j]->sr) {
                i++;

            //Current KU not present in previous SR
            } else if (prev_kus[i]->sr > sr_kus[j]->sr) {
                j++;

            //KU present in previous and current SRs - make sure positions match
            } else {
                k = j;
                while (k < sr_kus.size() && sr_kus[k]->sr == prev_kus[i]->sr) {
                    int pos_shift = 1;
                    if (sr_kus[k]->strand == read->ku_strands[ku]) 
                        pos_shift = -1;

                    if (prev_kus[i]->sr_idx == sr_kus[k]->sr_idx + pos_shift) {
                        next_kus.push_back(sr_kus[k]);
                        break;
                    }
                    k++;
                }
                if (k == j)
                    j++;
                i++;
            }
        }

        if (next_kus.size() == 0) {
            prev_kus.clear();
            break;
        }

        prev_kus.swap(next_kus);
        next_kus.clear();
    }
    
    //Fill result vector
    std::vector<int> result;
    for (unsigned int i = 0; i < prev_kus.size(); i++) {
        result.push_back(prev_kus[i]->sr);
    }

    return result;
}



