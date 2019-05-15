#ifndef _GHASHT_HH
#define _GHASHT_HH
#include "GBase.h"
//----------------------------------------------
//  Hash table templates based on Jeff Preshing's code
// ---------------------------------------------
//  Maps 32-bit integers to user data
//  Uses open addressing with linear probing.
//  In the m_cells array, key = 0 is reserved to indicate an unused cell.
//  Actual value for key 0 (if any) is stored in m_zeroCell.
//  The hash table automatically doubles in size when it becomes 75% full.
//  The hash table never shrinks in size, even after Clear(),
//  unless you explicitly call Compact().
//----------------------------------------------
inline uint32_t upper_power_of_two(uint32_t v) {
    v--;
    v |= v >> 1; v |= v >> 2;
    v |= v >> 4; v |= v >> 8;
    v |= v >> 16; v++;
    return v;
}

inline uint64_t upper_power_of_two(uint64_t v) {
  v--;
  v |= v >> 1; v |= v >> 2;
  v |= v >> 4; v |= v >> 8;
  v |= v >> 16; v |= v >> 32;
  v++;
  return v;
}

template <typename OBJ, typename CELL> class GHashT {
public:
//protected:
    CELL* m_cells;
    uint32 m_arraySize;
    uint32 m_population;
    bool m_zeroUsed;
    CELL m_zeroCell;
    void Resize(uint32 desiredSize);
    //for iteration over elements
//public:
    void init(uint32 initialSize = 32);
    GHashT(uint32 initialSize = 32):m_cells(NULL),m_arraySize(0),
       m_population(0),m_zeroUsed(false), m_zeroCell(), m_cur(NULL) { init(initialSize); }
    ~GHashT() {  delete[] m_cells; }
    uint32 Count() { return m_population; }
    // Basic operations
    CELL* Lookup(uint32 key);
    CELL* Insert(uint32 key);//Important: set the value to Insert()->value
    void Delete(CELL* cell);
    void Clear(uint32 initSize = 32) {
      delete[] m_cells;
      init(initSize);
    }
    void Compact() {
      Resize(upper_power_of_two((m_population * 4 + 3) / 3));
    }
 //----------------------------------------------
 //  Iteration
 //----------------------------------------------
 //protected:
    CELL* m_cur;
 //public:
    void startIterate() {
       m_cur = &m_zeroCell;
       if (!m_zeroUsed) NextCell();
    }
    CELL* NextCell();
};

template <class OBJ> class GIntHash { //OBJ requires a copy operator=
protected:
	struct Cell {
		uint32 key;
		OBJ value;
		Cell():key(0) { }
	};
	GHashT<OBJ, Cell> ghash;
public:
  GIntHash():ghash() {}
  OBJ* Add(uint32 key, OBJ val) {
  	Cell* c=ghash.Insert(key);
  	c->value = val;
  	return &(c->value);
  }
  OBJ* set(uint32 key, OBJ val) {
  	Cell* c=ghash.Insert(key);
  	c->value = val;
  	return &(c->value);
  }
  uint32 Count() { return ghash.Count(); }
  OBJ Replace(uint32 key, OBJ val) {
     //just like set() but returns a copy of the *old* value, if any
	 Cell* c=ghash.Insert(key);
	 OBJ oldv=c->value;
	 c->value=val;
	 return oldv;
  }
  void Clear() { ghash.Clear(); }
  void Compact() { ghash.Compact(); }
  void startIterate() { ghash.startIterate(); }
  void Delete(uint32 key) {
      Cell* cell = ghash.Lookup(key);
      if (cell) Delete(cell);
  }
  OBJ* Find(uint32 key) {
  	Cell* cell = ghash.Lookup(key);
  	return (cell ? & cell->value : NULL);
  }
  OBJ* get(uint32 key) {
  	Cell* cell = ghash.Lookup(key);
  	return (cell ? & cell->value : NULL);
  }
  OBJ* operator[](const uint32 ky) {
  	Cell* cell = ghash.Lookup(ky);
  	return (cell ? & cell->value : NULL);
  }
  OBJ* Next(uint32& nextky) {
      Cell* cell=ghash.NextCell();
      if (cell) {
        nextky=cell->key;
        return & (cell->value);
      }
      else {
        nextky=0;
        return NULL;
      }
  }
  uint32 NextKey() {
     Cell* cell=ghash.NextCell();
     if (cell) return cell->key;
       else return 0;
  }
  OBJ* NextValue() {
     Cell* cell=ghash.NextCell();
     if (cell) return & (cell->value);
       else return NULL;
  }
};

template <class OBJ> class GIntHashP {
protected:
	struct Cell {
		uint32 key;
		OBJ* value;
		Cell():key(0),value(NULL) { }
	};
	GHashT<OBJ, Cell> ghash;
	bool doFreeItems;
public:
  GIntHashP(bool freeItems=true):ghash(),doFreeItems(freeItems) {}
  ~GIntHashP() { Clear(); }
  OBJ* Add(uint32 key, OBJ* val) {
  	Cell* c=ghash.Insert(key);
  	c->value = val;
  	return c->value;
  }
  OBJ* set(uint32 key, OBJ* val) {
  	Cell* c=ghash.Insert(key);
  	c->value = val;
  	return c->value;
  }
  uint32 Count() { return ghash.Count(); }
  OBJ* Replace(uint32 key, OBJ* val) {
     //just like set() but returns a copy of the *old* value, if any
	 Cell* c=ghash.Insert(key);
	 OBJ* oldv=c->value;
	 c->value=val;
	 return oldv;
  }

  void startIterate() { ghash.startIterate(); }
  void Compact() { ghash.Compact(); }

  void Clear() {
	if (doFreeItems) {
		if (ghash.m_zeroUsed) delete ghash.m_zeroCell.value;
		ghash.startIterate();
		while (Cell* cell=ghash.NextCell()) {
			delete cell->value;
		}
	}
	ghash.Clear();
  }
  void Delete(uint32 key) {
      Cell* cell = ghash.Lookup(key);
      if (cell) {
    	  if (doFreeItems) {
    		  delete cell->value;
    	  }
    	  Delete(cell);
      }
  }
  OBJ* Find(uint32 key) {
  	Cell* cell = ghash.Lookup(key);
  	return (cell ? cell->value : NULL);
  }
  OBJ* get(uint32 key) {
  	Cell* cell = ghash.Lookup(key);
  	return (cell ? cell->value : NULL);
  }
  OBJ* operator[](const uint32 ky) {
  	Cell* cell = ghash.Lookup(ky);
  	return (cell ? cell->value : NULL);
  }
  OBJ* Next(uint32& nextky) {
      Cell* cell=ghash.NextCell();
      if (cell) {
        nextky=cell->key;
        return cell->value;
      }
      else {
        nextky=0;
        return NULL;
      }
  }
  uint32 NextKey() {
     Cell* cell=ghash.NextCell();
     if (cell) return cell->key;
       else return 0;
  }
  OBJ* NextValue() {
     Cell* cell=ghash.NextCell();
     if (cell) return cell->value;
       else return NULL;
  }
};


// from code.google.com/p/smhasher/wiki/MurmurHash3
inline uint32_t integerHash(uint32_t h)
{
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}

// from code.google.com/p/smhasher/wiki/MurmurHash3
inline uint64_t integerHash(uint64_t k)
{
	k ^= k >> 33;
	k *= 0xff51afd7ed558ccd;
	k ^= k >> 33;
	k *= 0xc4ceb9fe1a85ec53;
	k ^= k >> 33;
	return k;
}

#define GIHASH_FIRST_CELL(hash) (m_cells + ((hash) & (m_arraySize - 1)))
#define GIHASH_CIRCULAR_NEXT(c) ((c) + 1 != m_cells + m_arraySize ? (c) + 1 : m_cells)
#define GIHASH_CIRCULAR_OFFSET(a, b) ((b) >= (a) ? (b) - (a) : m_arraySize + (b) - (a))


//----------------------------------------------
//  constructor
//----------------------------------------------
template <typename OBJ, typename CELL> void GHashT<OBJ,CELL>::init(uint32 initialSize) {
    // Initialize regular cells
    m_arraySize = initialSize;
    GASSERT((m_arraySize & (m_arraySize - 1)) == 0);   // Must be a power of 2
    m_cells = new CELL[m_arraySize];
    memset(m_cells, 0, sizeof(CELL) * m_arraySize);
    m_population = 0;

    // Initialize zero cell
    m_zeroUsed = 0;
    m_zeroCell.key = 0;
    //m_zeroCell.value = 0;
}

//----------------------------------------------
//  Lookup key
//----------------------------------------------
template <typename OBJ, typename CELL> CELL* GHashT<OBJ, CELL>::Lookup(uint32 key) {
  if (key) {
    // Check regular cells
    for (CELL* cell = GIHASH_FIRST_CELL(integerHash(key));;
                cell = GIHASH_CIRCULAR_NEXT(cell)) {
      if (cell->key == key) return cell;
      if (!cell->key) return NULL;
    }
  }
  else {
    // Check zero cell
    if (m_zeroUsed) return &m_zeroCell;
    return NULL;
  }
};

//-----------------------------------------------------------------------
// Adding a key pair to the hash table, returns CELL
//IMPORTANT: Caller is responsible of setting the value into CELL->value
//-----------------------------------------------------------------------
template <typename OBJ, typename CELL> CELL* GHashT<OBJ, CELL>::Insert(uint32 key) {
    if (key) {
        // Check regular cells
        for (;;) {
            for (CELL* cell = GIHASH_FIRST_CELL(integerHash(key));;
                    cell = GIHASH_CIRCULAR_NEXT(cell)) {
                if (cell->key == key) {
                  // Found
                  //cell->value=val;
                  return cell;
                }
                if (cell->key == 0) {
                    // Insert here
                    if ((m_population + 1) * 4 >= m_arraySize * 3) {
                        // Time to resize
                        Resize(m_arraySize * 2);
                        break;
                    }
                    ++m_population;
                    cell->key = key;
                    //cell->value = val;
                    return cell;
                }
            }
        }
    }
    else   {
      // Check zero cell
      if (!m_zeroUsed) {
        // Insert here
        m_zeroUsed = true;
        if (++m_population * 4 >= m_arraySize * 3) {
          // Even though we didn't use a regular slot, let's keep the sizing rules consistent
          Resize(m_arraySize * 2);
        }
      }
      //m_zeroCell.value=val;
      return &m_zeroCell;
    }
}

//----------------------------------------------
//  Delete a key-value pair in the hash table
//----------------------------------------------
template <typename OBJ, typename CELL> void GHashT<OBJ,CELL>::Delete(CELL* cell) {
    if (cell != &m_zeroCell) {
        // Delete from regular cells
        GASSERT(cell >= m_cells && cell - m_cells < m_arraySize);
        GASSERT(cell->key);
        // Remove this cell by shuffling neighboring cells so there are no gaps in anyone's probe chain
        for (CELL* neighbor = GIHASH_CIRCULAR_NEXT(cell);;
                  neighbor = GIHASH_CIRCULAR_NEXT(neighbor)) {
            if (!neighbor->key) {
                // There's nobody to swap with. Go ahead and clear this cell, then return
                cell->key = 0;
                //cell->value = 0;
                m_population--;
                if (m_population<m_arraySize/2) Compact();
                return;
            }
            CELL* ideal = GIHASH_FIRST_CELL(integerHash(neighbor->key));
            if (GIHASH_CIRCULAR_OFFSET(ideal, cell) <
                     GIHASH_CIRCULAR_OFFSET(ideal, neighbor)) {
                // Swap with neighbor, then make neighbor the new cell to remove.
                *cell = *neighbor;
                cell = neighbor;
            }
        }
    }
    else {
        // Delete zero cell
        GASSERT(m_zeroUsed);
        m_zeroUsed = false;
        //cell->value = 0;
        m_population--;
        if (m_population<m_arraySize/2) Compact();
        return;
    }
}

//----------------------------------------------
//  Resize hash table
//----------------------------------------------
template <typename OBJ, typename CELL> void GHashT<OBJ,CELL>::Resize(uint32 desiredSize) {
    GASSERT((desiredSize & (desiredSize - 1)) == 0);   // Must be a power of 2
    GASSERT(m_population * 4  <= desiredSize * 3);
    // Get start/end pointers of old array
    CELL* oldCells = m_cells;
    CELL* end = m_cells + m_arraySize;
    // Allocate new array
    m_arraySize = desiredSize;
    m_cells = new CELL[m_arraySize];
    memset(m_cells, 0, sizeof(CELL) * m_arraySize);
    // Iterate through old array
    for (CELL* c = oldCells; c != end; c++)  {
        if (c->key) {
            // Insert this element into new array
            for (CELL* cell = GIHASH_FIRST_CELL(integerHash(c->key));;
            		cell = GIHASH_CIRCULAR_NEXT(cell))  {
                if (!cell->key) {
                    // Insert here
                    *cell = *c;
                    break;
                }
            }
        }
    }
    delete[] oldCells; // Delete old array
}

//--------------------------------------------------
//  return next cell (requires startIterate() first)
//--------------------------------------------------
template <typename OBJ, typename CELL> CELL* GHashT<OBJ,CELL>::NextCell() {
    // Already finished?
    if (!m_cur) return m_cur;
    // Iterate past zero cell
    if (m_cur == &m_zeroCell) m_cur = & (m_cells[-1]);
    // Iterate through the regular cells
    CELL* end = m_cells + m_arraySize;
    while (++m_cur != end) {
      if (m_cur->key) return m_cur;
    }
    // Finished
    return m_cur = NULL;
}


#endif
