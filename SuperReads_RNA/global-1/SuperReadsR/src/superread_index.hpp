#ifndef SUPERREAD_INDEX_H
#define SUPERREAD_INDEX_H

#include <string.h>
#include <vector>
#include <map>
#include <set>

//Length of read dataset prefix
#define PREFIX_LEN 2

//Stores k-unitigs of a single read
//Should should be read from read/k-unitig ("newTestOutput...") file
class Read {
    public:

    Read() {}
    Read(const Read &r); 
    
    //Main constructor - parses line from "newTestOutput..." file
    Read(std::string read_str);

    //Returns true if instance and r are mates based on read id
    //Note: does not check prefixes
    bool is_mate(Read r);

    //Returns true if reads are part of same dataset (have same prefix)
    bool same_prefix(Read r);
    
    int read_id,
        length; //Read length in basepairs

    //Set to true if the read was successfully parsed
    bool is_valid;

    //Stores read prefix - note not a proper string (no null character)
    char prefix[PREFIX_LEN];

    //Stores k-unitig information
    std::vector<int> ku_ids;
    std::vector<char> ku_strands;

    int last_ku_pos;
};

//Stores a set of reads to be passed to index thread
//Also stores results for thread to return
//Should only contain reads from a single dataset (sharing a prefix)
class IndexQuery {
    public:
    //Creates a buffer to store specified number of paired or unpaired reads
    IndexQuery(int buffer_size);

    ~IndexQuery();

    //Adds a single or paired-end read
    //Will only treat reads as paired if "is_paired" is true
    bool add_read(Read *read);

    bool is_full();
    bool is_empty();

    //Empties read buffer and all results
    void clear();
	
    //Stores reads loaded from a single dataset
    std::vector<Read*> reads;

    //Stores sets of super-reads matched to each read
    std::vector< std::vector<int>* > *results;

    //Stores read_ids corrasponding to each element in "results"
    std::vector<Read*> *matched_reads;

    bool finished, //Set to true once thread is finished
         is_paired; //Will treat reads as paired if true
    
    //Stores which dataset buffer has been loaded from
    //Set and read by ReadIO
    int prefix_id;

    #ifdef DEBUG_TIMER
    float thread_time;
    #endif
};

//Reads read/k-unitig file and writes read fastq files
//Keeps track of which dataset each read comes from by prefix
class ReadIO {
    public:

    //Sets read/k-unitig file and the fastq(s) for each dataset
    ReadIO(std::string read_ku_fname, 
           char *match_out_fname);
    ~ReadIO();

    //Loads reads into buffer and sets prefix and paired-end info
    bool fill_buffer(IndexQuery *q);

    //Gets successfully matched reads from a finished query
    void update_matched(IndexQuery *q);

    //Writes all paired and unpaired reads from all datasets to fastq files
    void write_unmatched(std::string pe1_fname, 
                         std::string pe2_fname, 
                         std::string unpaired_fname);
    
    private:

    //Translates read prefix into dataset ID
    int get_prefix_id(Read r);

    //Copies a fastq entry from one file to another
    bool copy_fastq(std::ifstream &in, std::ofstream &out);

    //Reads past a fastq entry
    bool skip_fastq(std::ifstream &in);
    
    //Buffer for read/k-unitig file
    std::ifstream read_ku_file;

    //Buffer to output super-read/read matches
    std::ofstream sr_matches_file;

    //Stores fastqs for each dataset
    std::vector<std::string> fastq_fnames;

    //Stores matched reads for each dataset
    std::vector< std::vector<int> > matched_reads;

    //Stores which datasets contain paired or unpaired reads
    std::vector<bool> is_paired;

    //Stores all dataset prefixes
    //Locations corrasponds to prefix_ids
    std::string prefixes;
};

//Stores map from k-unitigs to super-read locations
//and functions to match reads to super-reads
//All methods are static to make passing to threads easier
class SuperreadIndex {
    public:

    //Creates a new index given a list of superread names (K-unitig sequences) 
    //and a fasta of reduced SRs (only SRs present in fasta will be included).
    //Will reserve space for 'max_ku' K-unitigs if included
    //Returns true if successful, false otherwise
    static bool load_files(std::string sr_name_path, 
                           std::string sr_sam_path);

    //Frees all memory used by index
    static void destroy();

    //Matches reads to superreads and stores results in 'q'
    //Designed to be called as a thread
    //Always returns 0 (required for 'thread_pool')
    static int match_reads(IndexQuery *q);
    
    //Increments fragment (reads or mates) counts for matched superread groups
    static void update_frag_counts(IndexQuery *q);
    
    //Perform iteration of estimation/maximization to update superread sr_covs
    //Weights initially set to number of reads/fragments uniquely assigned
    //Returns maximum absolute difference of weight estimates for each superread
    static double iter_sr_covs();

    //Writes superreads to fastq with superread sr_covs
    //Should be called after iter_sr_covs()
    static void write_sam(std::string fname, std::string delim=":");
    
    //Stores superread IDs corrasponding to other SR lists
    static std::map<int, int> sr_indices;
    static std::vector< std::pair<int, int> > sr_id_lens;

    private:

    //Stores a record of a KU in a superread
    typedef struct kunitig {
        int sr,     //Which superread it is in
            sr_idx; //Index in superread, starting at 0 from the leftmost KU
        char strand;//Forward (F) or reverse (R) strand
    } *KUnitig;
   
    static std::string sr_sam_fname;
    
    //Stores locations of all K-unitigs
    //Index in list is equivilent to KU number
    static std::vector< std::vector<KUnitig> > ku_list;

    //Stores lengths of superreads
    static std::vector<int> sr_lens;

    //Stores counts of reads/fragments assigned to sets of superreads
    //Used to calculate sr_covs
    static std::map< std::vector<int>, int > sr_frag_counts;

    //Stores current suprread weight estimates
    //Before iter_sr_covs() is called, stores number of reads/fragments uniquely
    //assigned to each superread (initial estimate)
    static std::vector<double> sr_covs;

    //Parses superread std::string and indexes locations K-unitigs
    //'sr_index' should corraspond to element in 'sr_ids' and 'sr_lens'
    static void add_superread_kus(int sr_index, std::string sr_str);


    //Returns all superreads that contain given read
    static std::vector<int> match_read(Read *read);

    static double get_trimmed_len(int sr);
};


#endif


