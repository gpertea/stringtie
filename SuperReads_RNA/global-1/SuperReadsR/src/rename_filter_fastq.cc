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

#include <iostream>
#include <fstream>

#include <jellyfish/err.hpp>

#include <src/rename_filter_fastq_cmdline.hpp>

namespace err = jellyfish::err;

static args_t args;


const char* const error_header  = "Missing read header '@'";
const char* const error_qheader = "Missing qual value header '+'";
const char* const error_eof     = "Unexpected end of file";

void output_N_read(std::ostream& os, const size_t read_number) {
  os << "@" << args.library_arg << read_number << "\nN\n+\n" << args.base_quality_arg << "\n";
}

const char* parse_write(std::istream& is, std::ostream& os, const size_t read_number) {
  static std::string header;
  static std::string sequence;
  static std::string qheader;
  static std::string quals;

  if(!std::getline(is, header))
    return error_eof;
  if(header[0] != '@') return error_header;
  std::getline(is, sequence);
  std::getline(is, qheader);
  if(qheader[0] != '+') return error_qheader;
  std::getline(is, quals);
  if(sequence.size() == quals.size() && sequence.find_first_not_of("ACGTNacgtn") == std::string::npos) {
    os << "@" << args.library_arg << read_number << "\n"
       << sequence << "\n+\n" << quals << "\n";
  } else {
    output_N_read(os, read_number);
  }
  return 0;
}

int main(int argc, char *argv[])
{
  args.parse(argc, argv);
  std::string& path1 = args.path1_arg;
  std::string& path2 = args.path2_arg;

  bool unmated = (path1 == path2) || path2.empty();
  std::ifstream infile1(path1);
  if(!infile1.good())
    err::die(err::msg() << "Error while opening sequence file '" << path1 << "'");
  std::ifstream infile2;
  if(!unmated) {
    infile2.open(path2);
    if(!infile2.good())
      err::die(err::msg() << "Error while opening sequence file '" << path2 << "'");
  }

  std::ofstream outfile;
  if(args.output_given) {
    outfile.open(args.output_arg);
    if(!outfile.good())
      err::die(err::msg() << "Error while opening output file '" << args.output_arg << "'");
  }
  std::ostream& out = args.output_given ? outfile : std::cout;

  std::string line;
  size_t read_number = 0;
  while(infile1.good()) {
    const char* res = parse_write(infile1, out, read_number++);
    if(res == error_eof) break;
    if(res)
      err::die(err::msg() << "Invalid fastq format (" << res << ") in file '" << path1 << "' around position " << infile1.tellg());
    if(!unmated) {
      res = parse_write(infile2, out, read_number++);
      if(res)
        err::die(err::msg() << "Invalid fastq format (" << res << ") in file '" << path2 << "' around position " << infile2.tellg());
    } else {
      output_N_read(out, read_number++);
    }
  }

  return 0;
}
