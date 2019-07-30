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


#include <gtest/gtest.h>
#include <unittests/misc.hpp>
#include <src/read_parser.hpp>
#include <tmpstream.hpp>

static const int nb_sequences = 100000;
static const int nb_bases_dev = 100;
static const int nb_threads   = 5;

class ReadParserRandomSeq : public ::testing::Test {
protected:
  static void SetUpTestCase() {
    static const char bases[4] = { 'A', 'C', 'G', 'T' };

    fasta_total_bases = 0;
    for(int i = 0; i < nb_sequences; ++i) {
      random_fasta << ">" << i << "\n";
      int nb_bases = 50 + (random() % nb_bases_dev);
      for(int j = 0; j < nb_bases; ++j) {
        long r = random();
        random_fasta << bases[random() % 4];
        if(j != nb_bases - 1 && (r >> 2) % 100 == 0)
          random_fasta << "\n";
      }
      random_fasta << "\n";
      fasta_total_bases += nb_bases;
    }

    fastq_total_bases = 0;
    for(int i = 0; i < nb_sequences; ++i) {
      random_fastq << "@" << i << "\n";
      int nb_bases = 50 + (random() % nb_bases_dev);
      for(int j = 0; j < nb_bases; ++j) {
        long r = random();
        random_fastq << bases[random() % 4];
        if((r >> 2) % 100 == 0)
          random_fastq << "\n";
      }
      random_fastq << "\n+\n";
      for(int j = 0; j < nb_bases; ++j) {
        long r = random();
        random_fastq << (char)('!' + (r % 40));
        if(j != nb_bases - 1 && (r >> 6) % 100 == 0)
          random_fastq << "\n";
      }
      random_fastq << "\n";
      fastq_total_bases += nb_bases;
    }
  }

  std::istream fasta_stream;
  std::istream fastq_stream;
public:
  ReadParserRandomSeq() :
    fasta_stream(random_fasta.rdbuf()),
    fastq_stream(random_fastq.rdbuf())
  { }

  virtual void SetUp() {
    fasta_stream.seekg(0);
    fasta_stream.clear();
    fastq_stream.seekg(0);
    fastq_stream.clear();
  }


  static tmpstream random_fasta;
  static tmpstream random_fastq;
  static int fasta_total_bases;
  static int fastq_total_bases;
};

tmpstream ReadParserRandomSeq::random_fasta;
tmpstream ReadParserRandomSeq::random_fastq;
int ReadParserRandomSeq::fasta_total_bases = 0;
int ReadParserRandomSeq::fastq_total_bases = 0;

struct sum_line_lengths_data {
  read_parser              parser;
  jflib::atomic_field<int> sum;
  sum_line_lengths_data(std::istream& is) : parser(is), sum(0) { }
};
void* sum_line_lengths(void* d) {
  auto data = (sum_line_lengths_data*)d;
  read_parser::stream read_stream(data->parser);
  int bases = 0;

  for( ; read_stream; ++read_stream)
    bases += strlen(read_stream->sequence);

  data->sum += bases;
  return 0;
}

TEST_F(ReadParserRandomSeq, FastaLengths) {
  sum_line_lengths_data data(fasta_stream);

  EXPECT_TRUE(data.parser.good());
  EXPECT_FALSE(data.parser.eof());
  EXPECT_FALSE(data.parser.fail());
  
  pdo(nb_threads, sum_line_lengths, (void*)&data);

  EXPECT_EQ(fasta_total_bases, jflib::a_load(data.sum));
  EXPECT_TRUE(data.parser.eof());
  EXPECT_FALSE(data.parser.fail());
  EXPECT_FALSE(data.parser.good());
}

TEST_F(ReadParserRandomSeq, FastqLengths) {
  sum_line_lengths_data data(fastq_stream);

  EXPECT_TRUE(data.parser.good());
  EXPECT_FALSE(data.parser.eof());
  EXPECT_FALSE(data.parser.fail());
  
  pdo(nb_threads, sum_line_lengths, (void*)&data);

  EXPECT_EQ(fastq_total_bases, jflib::a_load(data.sum));
  EXPECT_TRUE(data.parser.eof());
  EXPECT_FALSE(data.parser.fail());
  EXPECT_FALSE(data.parser.good());
}

// TODO: test invalid input
// TEST(ReadParser, "Invalid input") {
  
// }
