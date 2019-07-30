#include <jellyfish/err.hpp>

#include <src/mer_database.hpp>

namespace err = jellyfish::err;

int main(int argc, char *argv[])
{
  if(argc < 3)
    err::die(err::msg() << "Usage: " << argv[0] << " db mer ...");

  database_query mer_database(argv[1]);
  mer_dna::k(mer_database.header().key_len() / 2);
  std::cout << mer_dna::k() << "\n";
  mer_dna m;

  for(int i = 2; i < argc; ++i) {
    m = argv[i];
    m.canonicalize();
    auto v = mer_database[m];
    std::cout << argv[i] << ":" << m << " val:" << v.first << " qual:" << v.second << "\n";
  }

  return 0;
}
