#include <iostream>
#include <fstream>
#include <string>

#include <src/split_mate_pairs_cmdline.hpp>

void print_sequence(std::string& line, std::ostream& os, const std::string& file) {
  os << line << '\n';
  std::getline(std::cin, line);
  os << line << '\n';
  if(!os.good())
    cmdline::error() << "Error while writing to file '" << file << "'";
}

int main(int argc, char *argv[]) {
  cmdline args(argc, argv);

  std::string file1(args.prefix_arg), file2(args.prefix_arg);
  file1 += "_1.fa";
  file2 += "_2.fa";
  std::ofstream out1(file1), out2(file2);
  if(!out1.good())
    cmdline::error() << "Failed to open output file '" << file1 << "'";
  if(!out2.good())
    cmdline::error() << "Failed to open output file '" << file2 << "'";

  std::string line;
  bool first = true;
  for( ; std::getline(std::cin, line); first = !first)
    print_sequence(line, first ? out1 : out2, first ? file1 : file2);

  return 0;
}
