#include <iostream>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/mer_dna.hpp>
#include <src/bloom_counter2.hpp>
#include <src/bloom_query_cmdline.hpp>

bloom_query_cmdline args;

struct mer_dna_hash {
  void operator()(const jellyfish::mer_dna& m, uint64_t *hashes) const {
    MurmurHash3_T_128(m, (m.k() / 4) + (m.k() % 4 != 0), 0, hashes);
  }
};
typedef bloom_counter2<jellyfish::mer_dna, mer_dna_hash> mer_bloom_counter2;

int main(int argc, char *argv[])
{
  args.parse(argc, argv);

  jellyfish::mer_dna::k(args.mer_arg);

  jellyfish::mapped_file dbf(args.input_arg);
  uint64_t* base = (uint64_t*)dbf.base();
  mer_bloom_counter2 mers(base[0], base[1], (unsigned char*)(base + 2));

  std::string word;
  jellyfish::mer_dna m;
  jellyfish::mer_dna mq;
  while(true) {
    std::cin >> word;
    if(!std::cin) 
      break;
    m = word;
    mq = m;
    mq.reverse_complement();
    std::cout << word << " " << mers[m < mq ? m : mq] << "\n";
  }

  return 0;
}
