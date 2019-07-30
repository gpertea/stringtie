#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <time.h>
#include <unistd.h>

#include "thread_pool.hpp"
#include "superread_index.hpp"
#include "timer.h"


//Command line arguments
const char OPTIONS[] = "w:s:r:o:p:b:m:k:";
const char WORK_DIR_OPT       = 'w',
           SR_SAM_OPT         = 's',
           READ_KU_OPT        = 'r',
           COV_DELIM_OPT      = 'd',
           OUT_FNAME_OPT     = 'o',
           NUM_THREADS_OPT    = 'p',
           BUFFER_SIZE_OPT    = 'b',
           MATCHES_OUT_OPT    = 'm';

//Argument defaults
const std::string DEF_WORK_DIR = "./work1/",
             DEF_READ_KU_FNAME = "newTestOutput.nucmerLinesOnly",
             DEF_COV_DELIM     = "YK:f:",
             SR_FASTA_FNAME    = "superReadSequences.fasta",
             SR_NAME_FNAME     = "superReadNames.txt"; 

int BUFFER_SIZE = 1000;
int NUM_THREADS = 8;

bool build_index(std::string sr_name_fname, std::string sr_seq_fname, int max_ku);

int main(int argc, char **argv) {
    
    //Keeps track of total time for various tasks
    //DEBUG_TIMER shoud be set in "superread_index.hpp"
    #ifdef DEBUG_TIMER
    double index_time = 0,
           input_time = 0,
           submit_time = 0,
           output_time = 0,
           total_thread_time = 0;
    Timer timer;
    int thread_count = 0;
    #endif

    //Stores option arguments
    std::string work_dir = DEF_WORK_DIR, 
           sr_cov_delim = DEF_COV_DELIM,
           read_ku_path,
           sr_sam_path, 
           sr_name_path,
           out_fname;

    char *sr_matches_fname = NULL;

    //Parse arguments
    char o;
    while ( (o = getopt(argc, argv, OPTIONS)) != -1 ) {
        switch (o) {
            case WORK_DIR_OPT:
                work_dir = std::string(optarg);
                if (work_dir[work_dir.size()-1] != '/')
                    work_dir += "/";
                break;
            case SR_SAM_OPT:
                sr_sam_path = std::string(optarg);
                break;
            case READ_KU_OPT:
                read_ku_path = std::string(optarg);
                break;
            case COV_DELIM_OPT:
                sr_cov_delim = std::string(optarg);
                break;
            case OUT_FNAME_OPT:
                out_fname = std::string(optarg);
                break;
            case MATCHES_OUT_OPT:
                sr_matches_fname = optarg;
                break;
            case NUM_THREADS_OPT:
                NUM_THREADS = atoi(optarg);
                break;
            case BUFFER_SIZE_OPT:
                BUFFER_SIZE = atoi(optarg);
                break;
        }
    }

    //Set default filenames
    //sr_fasta_path = work_dir + SR_FASTA_FNAME;
    sr_name_path = work_dir + SR_NAME_FNAME;

    if (read_ku_path.empty()) {
        read_ku_path = work_dir + DEF_READ_KU_FNAME;
    }


    //Initialze read io
    ReadIO short_read_io(read_ku_path, sr_matches_fname);

    //Build index    
    std::cerr << "Building super-read index...\n";
    if (!SuperreadIndex::load_files(sr_name_path, sr_sam_path))
        return 1;

    #ifdef DEBUG_TIMER
    index_time = timer.lap();
    #endif

    //Init thread pool
    int thread_ret; //Unused thread return value
    thread_pool<IndexQuery*, int> threads(NUM_THREADS, 
                                  SuperreadIndex::match_reads);

    //Init thread input buffers
    std::vector<IndexQuery*> queries(NUM_THREADS);
    for (int q = 0; q < NUM_THREADS; q++)
        queries[q] = new IndexQuery(BUFFER_SIZE);
    
    //Query buffer - where input is actually read
    IndexQuery *query_buf = new IndexQuery(BUFFER_SIZE);

    //true if buffer is ready to be filled
    //false if buffer points to thread parameter (not yet swapped)
    bool buffer_swapped = true;

    //True when all reads have been read and all threads are finished
    bool all_reads_assigned = false;

    std::cerr << "Index built, assigning reads...\n";

    //Read K Unitig file line-by-line
    while(!all_reads_assigned) {
    
        #ifdef DEBUG_TIMER
	    timer.reset();
	    #endif

        all_reads_assigned = !short_read_io.fill_buffer(query_buf);
        
        #ifdef DEBUG_TIMER
        input_time += timer.lap();
        #endif

    	//Submit job
        if (!all_reads_assigned) {
            threads.submit_job(&query_buf, &thread_ret);
            buffer_swapped = false; 
        }
        
        #ifdef DEBUG_TIMER
        submit_time += timer.lap();
        #endif

        //Check status of running threads
        for (int q = 0; q < NUM_THREADS; q++) {
            if (queries[q]->finished) {

                #ifdef DEBUG_TIMER
        		total_thread_time += queries[q]->thread_time;
                #endif

                //Parse results if present
                if (queries[q]->results->size() > 0) {

                    #ifdef DEBUG_TIMER
                    thread_count++;
                    #endif

                    //Store super-read results
                    SuperreadIndex::update_frag_counts(queries[q]);
                    
                    //Store matched reads
                    short_read_io.update_matched(queries[q]);
                    
                    //Empty buffer
                    queries[q]->clear();
                }

                //Store buffer results if needed
                if (!buffer_swapped) {
                    std::swap(query_buf, queries[q]);
                    query_buf->finished = false;
                    buffer_swapped = true;
                }
            
            //Threads still running
            } else {
                all_reads_assigned = false;
            }
        }
        
        #ifdef DEBUG_TIMER
        output_time += timer.get();
        #endif
    }
    
    std::cerr << "All reads assigned, calculating coverage..." << "\n";
    
    //for (int i = 0; i < 100; i++) 
    //    if (SuperreadIndex::iter_sr_covs() < 0.01)
    //        break;

    int iters = 0;
    while (SuperreadIndex::iter_sr_covs() > 0.01) {
        iters++;
    }
    std::cerr << "Completed after " << iters << " iterations\n";
    
    std::cerr << "Weights calculated, writing super-read SAM\n";
    SuperreadIndex::write_sam(out_fname, sr_cov_delim);

    std::cerr << "Cleaning up" << "\n";
    
    threads.release_workers();
    SuperreadIndex::destroy();
    delete query_buf;

    for (int q = 0; q < NUM_THREADS; q++) 
        delete queries[q];

    #ifdef DEBUG_TIMER
    std::cerr << "\nMain thread\n";
    std::cerr << "Index time:  " << index_time  << "s\n"; 
    std::cerr << "Parse time:  " << input_time  << "s\n";
    std::cerr << "Wait time:   " << submit_time << "s\n";
    std::cerr << "Output time: " << output_time << "s\n";
    
    std::cerr << "\nIndex threads\n";
    std::cerr << "Total thread time: " << total_thread_time <<  "\n";
    std::cerr << "Threads run:       " << thread_count << "\n";
    std::cerr << "Time per thread:   " << (total_thread_time / thread_count) << "\n";
    #endif

    return 0;
}

