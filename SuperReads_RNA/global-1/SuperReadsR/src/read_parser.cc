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


#include <src/read_parser.hpp>

read_parser::~read_parser() {
  if(filebufs_.empty()) {
    if(close_input_)
      delete input_.rdbuf();
  } else {
    for(auto it = filebufs_.begin(); it != filebufs_.end(); ++it)
      delete *it;
  }
}

void read_parser::parser_loop() {
  if(filebufs_.empty())
    parse_input_stream();
  else {
    for(auto it = filebufs_.begin(); it != filebufs_.end(); ++it) {
      input_.rdbuf(*it);
      if(!parse_input_stream())
        break;
    }
  }
}


bool read_parser::parse_input_stream() {
  switch(input_.peek()) {
  case '>': return fasta_reader_loop();
  case '@': return fastq_reader_loop();
  default:
    super::report_error("Invalid file format");
    return false;
  }
}

bool read_parser::fasta_reader_loop() {
  int nextc = input_.peek();
  while(nextc != EOF) {
    elt e(elt_init());
    if(e.is_empty())
      break;

    size_type& i = e->nb_filled;
    for(i = 0 ; nextc != EOF && i < group_size(); ++i) {
      read& r = e->elements[i];
      getline(input_, r.header);
      r.sequence.clear();
      nextc = input_.peek();
      while(nextc != EOF && nextc != '>') {
        getline_append(input_, r.sequence);
        r.sequence.chomp();
        nextc = input_.peek();
      }
    }
  }
  return true;
}

bool read_parser::fastq_reader_loop() {
  charb unused_line;
  int nextc = input_.peek();
  while(nextc != EOF) {
    elt e(elt_init());
    if(e.is_empty())
      break;

    size_type& i = e->nb_filled;
    for(i= 0; nextc != EOF && i < group_size(); ++i) {
      read& r = e->elements[i];
      getline(input_, r.header);
      if(r.header[0] != '@') {
        report_error("Found bad sequence header");
        nextc = EOF;
        return false;
      }
      r.sequence.clear();
      r.quals.clear();
      nextc = input_.peek();
      while(nextc != EOF && nextc != '+') {
        getline_append(input_, r.sequence);
        r.sequence.chomp();
        nextc = input_.peek();
      }
      getline(input_, unused_line); // read quals header: ignored
      while(nextc != EOF && r.quals.size() < r.sequence.size()) {
        getline_append(input_, r.quals);
        r.quals.chomp();
        nextc = input_.peek();
      }
      if(r.quals.size() != r.sequence.size()) { // Invalid input file
        report_error("Number of qual values != number of bases");
        nextc = EOF;
        return false;
      }
    }
  }

  return true;
}

