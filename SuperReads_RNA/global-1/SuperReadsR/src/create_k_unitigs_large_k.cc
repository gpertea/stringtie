/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>

#include <iostream>
#include <memory>
#include <algorithm>
#include <set>

#include <src/create_k_unitigs_common.hpp>
#include <gzip_stream.hpp>
#include <jflib/multiplexed_io.hpp>
#include <src/create_k_unitigs_large_k_cmdline.hpp>

// GLOBAL: command line switches
cmdline_parse args;

class mer_counter {
  mer_dna_bloom_counter counter_;
public:
  mer_counter(double fp, size_t size) : counter_(fp, size) { }
  mer_dna_bloom_counter::element_proxy operator[](const mer_dna& m) {
    return mer_dna_bloom_counter::element_proxy(counter_, m);
  }
  mer_dna_bloom_counter::const_element_proxy operator[](const mer_dna& m) const {
    return mer_dna_bloom_counter::const_element_proxy(counter_, m);
  }
  void done() { }
};

typedef populate_mer_set<mer_counter, sequence_parser, read_mers> mer_populate;


typedef create_k_unitig<mer_counter, mer_bloom, mer_bloom,
                        sequence_parser, read_mers,
                        cmdline_parse> unitiger_type;

std::ostream* open_output() {
  if(!args.output_given)
    return new std::ostream(std::cout.rdbuf());
  if(args.gzip_flag)
    return new gzipstream(args.output_arg);
  return new std::ofstream(args.output_arg);
}

int main(int argc, char *argv[])
{
  args.parse(argc, argv);

  mer_dna::k(args.mer_arg);
  if(!args.min_len_given)
    args.min_len_arg = args.mer_arg + 1;

  // Populate Bloom filter with k-mers
  std::unique_ptr<mer_counter> kmers;
  std::unique_ptr<jellyfish::mapped_file> dbf;

  {
    // if(args.load_given) {
    //   dbf.reset(new jellyfish::mapped_file(args.load_arg));
    //   uint64_t* base = (uint64_t*)dbf->base();
    //   kmers.reset(new mer_dna_bloom_counter(base[0], base[1], (unsigned char*)(base + 2)));
    // } else {
      kmers.reset(new mer_counter(args.false_positive_arg, args.nb_mers_arg));
      stream_manager manager(args.input_arg.begin(), args.input_arg.end());\
      sequence_parser parser(mer_dna::k(), manager.nb_streams(), 3 * args.threads_arg, 4096, manager);
      mer_populate populate(args.mer_arg, *kmers, parser);
      populate.exec_join(args.threads_arg);
      //    }
  }

  {
    std::unique_ptr<std::ostream> output_ostream(open_output());
    mer_bloom used(args.false_positive_arg, args.nb_mers_arg);
    mer_bloom ends(args.false_positive_arg, args.nb_mers_arg);
    stream_manager manager(args.input_arg.begin(), args.input_arg.end());
    sequence_parser parser(mer_dna::k(), manager.nb_streams(), 3 * args.threads_arg, 4096, manager);
    unitiger_type unitiger(*kmers, used, ends, args.threads_arg, parser,
                           *output_ostream, args);
    unitiger.exec_join(args.threads_arg);
  }

  return EXIT_SUCCESS;
}
