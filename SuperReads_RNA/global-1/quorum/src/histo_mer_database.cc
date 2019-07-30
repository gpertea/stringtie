#include <algorithm>

#include <jellyfish/err.hpp>
#include <src/mer_database.hpp>

namespace err = jellyfish::err;

int main(int argc, char* argv[])
{
  if(argc != 2)
    err::die(err::msg() << "Usage: " << argv[0] << " db");

  static const size_t hlen = 1001;
  uint64_t histos[hlen][2];
  memset(histos, '\0', sizeof(histos));
  const database_query mer_database(argv[1]);
  for(auto it = mer_database.begin(); it != mer_database.end(); ++it) {
    auto& vals = it->second;
    const uint64_t val = std::min(vals.first, hlen - 1);
    ++histos[val][vals.second];
  }

  for(size_t i = 0; i < hlen; ++i) {
    if(histos[i][0] || histos[i][1])
      std::cout << i << " " << histos[i][0] << " " << histos[i][1] << "\n";
  }

  return 0;
}
