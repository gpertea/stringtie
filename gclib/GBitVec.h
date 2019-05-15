#ifndef __GBITVEC_H__
#define __GBITVEC_H__
#include "GBase.h"
//this code is lifted from LLVM (llvm.org, BitVector.h)

/// bitCount_32 - this function counts the number of set bits in a value.
/// Ex. CountPopulation(0xF000F000) = 8
/// Returns 0 if the word is zero.
inline uint bitCount_32(uint32_t Value) {
#if __GNUC__ >= 4
    return __builtin_popcount(Value);
#else
    uint32_t v = Value - ((Value >> 1) & 0x55555555);
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
    return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
#endif
  }

/// bitCount_64 - this function counts the number of set bits in a value,
/// (64 bit edition.)
inline uint bitCount_64(uint64_t Value) {
#if __GNUC__ >= 4
    return __builtin_popcountll(Value);
#else
    uint64_t v = Value - ((Value >> 1) & 0x5555555555555555ULL);
    v = (v & 0x3333333333333333ULL) + ((v >> 2) & 0x3333333333333333ULL);
    v = (v + (v >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return uint((uint64_t)(v * 0x0101010101010101ULL) >> 56);
#endif
  }

/// CountTrailingZeros_32 - this function performs the platform optimal form of
/// counting the number of zeros from the least significant bit to the first one
/// bit.  Ex. CountTrailingZeros_32(0xFF00FF00) == 8.
/// Returns 32 if the word is zero.
inline unsigned bitCountTrailingZeros_32(uint32_t Value) {
#if __GNUC__ >= 4
  return Value ? __builtin_ctz(Value) : 32;
#else
  static const unsigned Mod37BitPosition[] = {
    32, 0, 1, 26, 2, 23, 27, 0, 3, 16, 24, 30, 28, 11, 0, 13,
    4, 7, 17, 0, 25, 22, 31, 15, 29, 10, 12, 6, 0, 21, 14, 9,
    5, 20, 8, 19, 18
  };
  return Mod37BitPosition[(-Value & Value) % 37];
#endif
}

// CountTrailingZeros_64 - This function performs the platform optimal form
/// of counting the number of zeros from the least significant bit to the first
/// one bit (64 bit edition.)
/// Returns 64 if the word is zero.
inline unsigned bitCountTrailingZeros_64(uint64_t Value) {
#if __GNUC__ >= 4
  return Value ? __builtin_ctzll(Value) : 64;
#else
  static const unsigned Mod67Position[] = {
    64, 0, 1, 39, 2, 15, 40, 23, 3, 12, 16, 59, 41, 19, 24, 54,
    4, 64, 13, 10, 17, 62, 60, 28, 42, 30, 20, 51, 25, 44, 55,
    47, 5, 32, 65, 38, 14, 22, 11, 58, 18, 53, 63, 9, 61, 27,
    29, 50, 43, 46, 31, 37, 21, 57, 52, 8, 26, 49, 45, 36, 56,
    7, 48, 35, 6, 34, 33, 0
  };
  return Mod67Position[(-Value & Value) % 67];
#endif
}

class GBitVec {
  typedef unsigned long BitWord;

  enum { BITWORD_SIZE = (uint)sizeof(BitWord) * CHAR_BIT };

  BitWord  *fBits;        // Actual bits.
  uint Size;         // Size of GBitVec in bits.
  uint Capacity;     // Size of allocated memory in BitWord.

public:
  // Encapsulation of a single bit.
  class GBitRef {
    friend class GBitVec;
    BitWord *WordRef;
    uint BitPos;
    GBitRef();  // Undefined
  public:
    GBitRef(GBitVec &b, uint Idx) {
      WordRef = &b.fBits[Idx / BITWORD_SIZE];
      BitPos = Idx % BITWORD_SIZE;
    }

    ~GBitRef() {}

    GBitRef &operator=(GBitRef t) {
      *this = bool(t);
      return *this;
    }

    GBitRef& operator=(bool t) {
      if (t)
        *WordRef |= 1L << BitPos;
      else
        *WordRef &= ~(1L << BitPos);
      return *this;
    }

    operator bool() const {
      return ((*WordRef) & (1L << BitPos)) ? true : false;
    }
  };


  /// GBitVec default ctor - Creates an empty GBitVec.
  GBitVec() : Size(0), Capacity(0) {
    fBits = 0;
  }

  /// GBitVec ctor - Creates a GBitVec of specified number of bits. All
  /// bits are initialized to the specified value.
  explicit GBitVec(uint s, bool value = false) : Size(s) {
    Capacity = NumBitWords(s);
    //fBits = (BitWord *)std::malloc(Capacity * sizeof(BitWord));
    GMALLOC(fBits, Capacity * sizeof(BitWord));
    init_words(fBits, Capacity, value);
    if (value)
      clear_unused_bits();
  }
  unsigned long getMemorySize() const {
	   unsigned long r = ((unsigned long) Capacity) * sizeof(BitWord);
	   return r;
  }

  /// GBitVec copy ctor.
  GBitVec(const GBitVec &RHS) : Size(RHS.size()) {
    if (Size == 0) {
      fBits = 0;
      Capacity = 0;
      return;
    }

    Capacity = NumBitWords(RHS.size());
    GMALLOC(fBits, Capacity * sizeof(BitWord));
    memcpy(fBits, RHS.fBits, Capacity * sizeof(BitWord));
  }

  ~GBitVec() {
    GFREE(fBits);
  }

  /// empty - Tests whether there are no bits in this GBitVec.
  bool empty() const { return Size == 0; }

  /// size - Returns the number of bits in this GBitVec.
  uint size() const { return Size; }


  void bitSizeError() {
    GError("Error at GBitVec: unsupported BitWord size (%d)!\n", 
        sizeof(BitWord));
    }
  /// count - Returns the number of bits which are set.
  uint count() {
    uint NumBits = 0;
    for (uint i = 0; i < NumBitWords(size()); ++i)
      if (sizeof(BitWord) == 4)
        NumBits += bitCount_32((uint32_t)fBits[i]);
      else if (sizeof(BitWord) == 8)
        NumBits += bitCount_64(fBits[i]);
      else
        bitSizeError();
    return NumBits;
  }

  /// any - Returns true if any bit is set.
  bool any() {
    for (uint i = 0; i < NumBitWords(size()); ++i)
      if (fBits[i] != 0)
        return true;
    return false;
  }

  /// all - Returns true if all bits are set.
  bool all() {
    // TODO: Optimize this.
    return count() == size();
  }

  /// none - Returns true if none of the bits are set.
  bool none() {
    return !any();
  }

  /// find_first - Returns the index of the first set bit, -1 if none
  /// of the bits are set.
  int find_first() {
    for (uint i = 0; i < NumBitWords(size()); ++i)
      if (fBits[i] != 0) {
        if (sizeof(BitWord) == 4)
          return i * BITWORD_SIZE + bitCountTrailingZeros_32((uint32_t)fBits[i]);
        else if (sizeof(BitWord) == 8)
          return i * BITWORD_SIZE + bitCountTrailingZeros_64(fBits[i]);
        else
          bitSizeError();
      }
    return -1;
  }

  /// find_next - Returns the index of the next set bit following the
  /// "Prev" bit. Returns -1 if the next set bit is not found.
  int find_next(uint Prev) {
    ++Prev;
    if (Prev >= Size)
      return -1;

    uint WordPos = Prev / BITWORD_SIZE;
    uint BitPos = Prev % BITWORD_SIZE;
    BitWord Copy = fBits[WordPos];
    // Mask off previous bits.
    Copy &= ~0UL << BitPos;

    if (Copy != 0) {
      if (sizeof(BitWord) == 4)
        return WordPos * BITWORD_SIZE + bitCountTrailingZeros_32((uint32_t)Copy);
      else if (sizeof(BitWord) == 8)
        return WordPos * BITWORD_SIZE + bitCountTrailingZeros_64(Copy);
      else
        bitSizeError();
    }

    // Check subsequent words.
    for (uint i = WordPos+1; i < NumBitWords(size()); ++i)
      if (fBits[i] != 0) {
        if (sizeof(BitWord) == 4)
          return i * BITWORD_SIZE + bitCountTrailingZeros_32((uint32_t)fBits[i]);
        else if (sizeof(BitWord) == 8)
          return i * BITWORD_SIZE + bitCountTrailingZeros_64(fBits[i]);
        else
          bitSizeError();
      }
    return -1;
  }

  /// clear - Clear all bits; does NOT release memory
  void clear() {
    Size = 0;
  }

  /// resize - Grow or shrink the GBitVec.
  void resize(uint N, bool value = false) {
    if (N > Capacity * BITWORD_SIZE) {
      uint OldCapacity = Capacity;
      grow(N);
      init_words(&fBits[OldCapacity], (Capacity-OldCapacity), value);
    }

    // Set any old unused bits that are now included in the GBitVec. This
    // may set bits that are not included in the new vector, but we will clear
    // them back out below.
    if (N > Size)
      set_unused_bits(value);

    // Update the size, and clear out any bits that are now unused
    uint OldSize = Size;
    Size = N;
    if (value || N < OldSize)
      clear_unused_bits();
  }

  void reserve(uint N) {
    if (N > Capacity * BITWORD_SIZE)
      grow(N);
  }

  // Set, reset, flip
  GBitVec &set() {
    init_words(fBits, Capacity, true);
    clear_unused_bits();
    return *this;
  }

  GBitVec &set(uint Idx) {
    fBits[Idx / BITWORD_SIZE] |= 1L << (Idx % BITWORD_SIZE);
    return *this;
  }

  GBitVec &reset() {
    init_words(fBits, Capacity, false);
    return *this;
  }

  GBitVec &reset(uint Idx) {
    fBits[Idx / BITWORD_SIZE] &= ~(1L << (Idx % BITWORD_SIZE));
    return *this;
  }

  GBitVec &flip() {
    for (uint i = 0; i < NumBitWords(size()); ++i)
      fBits[i] = ~fBits[i];
    clear_unused_bits();
    return *this;
  }

  GBitVec &flip(uint Idx) {
    fBits[Idx / BITWORD_SIZE] ^= 1L << (Idx % BITWORD_SIZE);
    return *this;
  }

  // No argument flip.
  GBitVec operator~() const {
    return GBitVec(*this).flip();
  }

  inline static void indexCheck(uint vIdx, uint vSize) {
    if (vIdx >= vSize)
      GError("Error at GBitVec: index %d out of bounds (size %d)\n", 
        (int)vIdx, vSize);
   }

  // Indexing.
  GBitRef operator[](uint Idx) {
    //assert (Idx < Size && "Out-of-bounds Bit access.");
    indexCheck(Idx, Size);
    return GBitRef(*this, Idx);
  }

  bool operator[](uint Idx) const {
    indexCheck(Idx, Size);
    BitWord Mask = 1L << (Idx % BITWORD_SIZE);
    return (fBits[Idx / BITWORD_SIZE] & Mask) != 0;
  }

  bool test(uint Idx) const {
    return (*this)[Idx];
  }

  // Comparison operators.
  bool operator==(const GBitVec &RHS) const {
    uint ThisWords = NumBitWords(size());
    uint RHSWords  = NumBitWords(RHS.size());
    uint i;
    uint imax=GMIN(ThisWords, RHSWords);
    for (i = 0; i != imax; ++i)
      if (fBits[i] != RHS.fBits[i])
        return false;

    // Verify that any extra words are all zeros.
    if (i != ThisWords) {
      for (; i != ThisWords; ++i)
        if (fBits[i])
          return false;
    } else if (i != RHSWords) {
      for (; i != RHSWords; ++i)
        if (RHS.fBits[i])
          return false;
    }
    return true;
  }

  bool operator!=(const GBitVec &RHS) const {
    return !(*this == RHS);
  }

  // Intersection, union, disjoint union.
  GBitVec &operator&=(const GBitVec &RHS) {
    uint ThisWords = NumBitWords(size());
    uint RHSWords  = NumBitWords(RHS.size());
    uint i;
    uint imax=GMIN(ThisWords, RHSWords);
    for (i = 0; i != imax; ++i)
      fBits[i] &= RHS.fBits[i];

    // Any bits that are just in this GBitVec become zero, because they aren't
    // in the RHS bit vector.  Any words only in RHS are ignored because they
    // are already zero in the LHS.
    for (; i != ThisWords; ++i)
      fBits[i] = 0;

    return *this;
  }

  GBitVec &operator|=(const GBitVec &RHS) {
    if (size() < RHS.size())
      resize(RHS.size());
    for (size_t i = 0, e = NumBitWords(RHS.size()); i != e; ++i)
      fBits[i] |= RHS.fBits[i];
    return *this;
  }

  GBitVec &operator^=(const GBitVec &RHS) {
    if (size() < RHS.size())
      resize(RHS.size());
    for (size_t i = 0, e = NumBitWords(RHS.size()); i != e; ++i)
      fBits[i] ^= RHS.fBits[i];
    return *this;
  }

  // Assignment operator.
  const GBitVec &operator=(const GBitVec &RHS) {
    if (this == &RHS) return *this;

    Size = RHS.size();
    uint RHSWords = NumBitWords(Size);
    if (Size <= Capacity * BITWORD_SIZE) {
      if (Size)
        memcpy(fBits, RHS.fBits, RHSWords * sizeof(BitWord));
      clear_unused_bits();
      return *this;
    }

    // Grow the GBitVec to have enough elements.
    Capacity = RHSWords;
    //BitWord *NewBits = (BitWord *)std::malloc(Capacity * sizeof(BitWord));
    BitWord *NewBits = NULL;
    GMALLOC(NewBits, Capacity * sizeof(BitWord));
    memcpy(NewBits, RHS.fBits, Capacity * sizeof(BitWord));

    // Destroy the old bits.
    GFREE(fBits);
    fBits = NewBits;

    return *this;
  }

  void swap(GBitVec &RHS) {
    Gswap(fBits, RHS.fBits);
    Gswap(Size, RHS.Size);
    Gswap(Capacity, RHS.Capacity);
  }

private:
  uint NumBitWords(uint S) const {
    return (S + BITWORD_SIZE-1) / BITWORD_SIZE;
  }

  // Set the unused bits in the high words.
  void set_unused_bits(bool value = true) {
    //  Set high words first.
    uint UsedWords = NumBitWords(Size);
    if (Capacity > UsedWords)
      init_words(&fBits[UsedWords], (Capacity-UsedWords), value);

    //  Then set any stray high bits of the last used word.
    uint ExtraBits = Size % BITWORD_SIZE;
    
    if (ExtraBits) {
      BitWord ExtraBitMask = ~0UL << ExtraBits;
      if (value)
        fBits[UsedWords-1] |= ExtraBitMask;
      else
        fBits[UsedWords-1] &= ~ExtraBitMask;
    }
  }

  // Clear the unused bits in the high words.
  void clear_unused_bits() {
    set_unused_bits(false);
  }

  void grow(uint NewSize) {
    Capacity = GMAX(NumBitWords(NewSize), Capacity * 2);
    //fBits = (BitWord *)std::realloc(fBits, Capacity * sizeof(BitWord));
    GREALLOC(fBits, Capacity * sizeof(BitWord));
    clear_unused_bits();
  }

  void init_words(BitWord *B, uint NumWords, bool value) {
    memset(B, 0 - (int)value, NumWords*sizeof(BitWord));
  }
};


inline GBitVec operator&(const GBitVec &LHS, const GBitVec &RHS) {
  GBitVec Result(LHS);
  Result &= RHS;
  return Result;
}

inline GBitVec operator|(const GBitVec &LHS, const GBitVec &RHS) {
  GBitVec Result(LHS);
  Result |= RHS;
  return Result;
}

inline GBitVec operator^(const GBitVec &LHS, const GBitVec &RHS) {
  GBitVec Result(LHS);
  Result ^= RHS;
  return Result;
}

inline void Gswap(GBitVec &LHS, GBitVec &RHS) {
  LHS.swap(RHS);
  }

#endif
