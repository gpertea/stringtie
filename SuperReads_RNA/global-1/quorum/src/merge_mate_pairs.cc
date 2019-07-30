#include <vector>
#include <iterator>
#include <iostream>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <src/merge_mate_pairs_cmdline.hpp>

template<typename ITERATOR>
class skip_iterator :
  public std::iterator<std::input_iterator_tag, typename std::iterator_traits<ITERATOR>::value_type>
{
  typedef typename std::iterator<std::input_iterator_tag, typename std::iterator_traits<ITERATOR>::value_type> base;
  ITERATOR m_base;
public:
  typedef typename base::value_type value_type;
  typedef typename base::reference reference;
  typedef typename base::pointer pointer;

  skip_iterator(const ITERATOR& it) : m_base(it) { }
  skip_iterator& operator=(const ITERATOR& it) {
    m_base = it;
  }
  bool operator==(const skip_iterator& rhs) const {
    if(m_base == rhs.m_base) return true;
    ITERATOR next(m_base);
    std::advance(next, 1);
    if(m_base == next) return true;
    next = rhs.m_base;
    std::advance(next, 1);
    return m_base == next;
  }
  bool operator!=(const skip_iterator& rhs) const {
    return !(*this == rhs);
  }
  reference operator*() { return *m_base; }
  pointer operator->() { return &*m_base; }
  skip_iterator& operator++() {
    std::advance(m_base, 2);
    return *this;
  }
  skip_iterator operator++(int) {
    skip_iterator res(*this);
    ++*this;
    return res;
  }
};
typedef skip_iterator<cmdline::file_arg_it>              file_iterator_type;
typedef jellyfish::stream_manager<file_iterator_type>    stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

template<typename SEQ>
void print_seq(const SEQ& s) {
  std::cout << '@' << s.header << '\n'
            << s.seq
            << "\n+\n";
  if(!s.qual.empty())
    std::cout << s.qual << '\n';
  else
    std::cout << std::string(s.seq.size(), '*') << '\n';
}

int main(int argc, char *argv[]) {
  cmdline args(argc, argv);

  if(args.file_arg.size() % 2 != 0)
    cmdline::error() << "Must give a even number files";

  file_iterator_type even_file(args.file_arg.begin());
  file_iterator_type odd_file(args.file_arg.begin() + 1);
  const file_iterator_type end_file(args.file_arg.end());

  stream_manager even_stream(even_file, end_file);
  stream_manager odd_stream(odd_file, end_file);
  sequence_parser even_sequence(4, 100, 1, even_stream);
  sequence_parser odd_sequence(4, 100, 1, odd_stream);

  while(true) {
    sequence_parser::job even_j(even_sequence);
    sequence_parser::job odd_j(odd_sequence);
    if(even_j.is_empty() != odd_j.is_empty())
      throw std::runtime_error("Input files are not paired reads.");
    if(even_j.is_empty()) break;
    if(even_j->nb_filled != odd_j->nb_filled)
      throw std::runtime_error("Input files are not paired reads.");
    for(size_t i = 0; i < even_j->nb_filled; ++i) {
      print_seq(even_j->data[i]);
      print_seq(odd_j->data[i]);
    }
  }

  return 0;
}
