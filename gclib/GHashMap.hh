/********************************************************************************
 *                  Hash map class templates
 *********************************************************************************/

#ifndef GHashMap_HH
#define GHashMap_HH
#include "GBase.h"
#include "khashl.hh"
#include <type_traits>
#include <typeinfo>

#define XXH_INLINE_ALL 1
#include "xxhash.h"
#include "wyhash.h"

template <typename K> struct GHashKey_xxHash32 { //K generic (class, primitive, pointer except const char* )
  //template <typename T=K> inline typename std::enable_if< std::is_trivial<T>::value, uint32_t>::type
 uint32_t operator()(const K& s) const { //only works for trivial types!
      static_assert(std::is_trivial<K>::value, "Error: cannot use this for non-trivial types!\n");
      return XXH32((const void *) &s, sizeof(K), 0);
    }
};

template <> struct GHashKey_xxHash32<const char*> {
   inline uint32_t operator()(const char* s) const {
      return XXH32(s, strlen(s), 0);
   }
};

template <typename K> struct GHashKey_xxHash { //K generic (class, primitive, pointer except const char* )
  //template <typename T=K> inline typename std::enable_if< std::is_trivial<T>::value, uint32_t>::type
 uint64_t operator()(const K& s) const { //only works for trivial types!
      static_assert(std::is_trivial<K>::value, "Error: cannot use this for non-trivial types!\n");
      return XXH64((const void *) &s, sizeof(K), 0);
    }
};

template <> struct GHashKey_xxHash<const char*> {
   inline uint32_t operator()(const char* s) const {
      return XXH64(s, strlen(s), 0);
   }
};


template <typename K> struct GHashKey_wyHash { //K generic (class, primitive, pointer except const char* )
  //template <typename T=K> inline typename std::enable_if< std::is_trivial<T>::value, uint32_t>::type
 uint64_t operator()(const K& s) const { //only works for trivial types!
      static_assert(std::is_trivial<K>::value, "Error: cannot use this for non-trivial types!\n");
      return wyhash((const void *) &s, sizeof(K), 0, _wyp);
    }
};

template <> struct GHashKey_wyHash<const char*> {
   inline uint32_t operator()(const char* s) const {
      return wyhash(s, strlen(s), 0, _wyp);
   }
};

template <typename K> struct GHashKey_Eq { //K is a type having the == operator defined
    inline bool operator()(const K& x, const K& y) const {
      return (x == y); //requires == operator to be defined for K
    }
};

template <> struct GHashKey_Eq<const char*> {
    inline bool operator()(const char* x, const char* y) const {
      return (strcmp(x, y) == 0);
    }
};

// GHashSet<KType> never makes a deep copy of a char* key, it only stores the pointer
//  - for pointer keys like char*, key allocation must be managed separately (and should always survive the GHashSet)
template <typename K=const char*, class Hash=GHashKey_wyHash<K>, class Eq=GHashKey_Eq<K>, typename khInt_t=uint64_t >
  class GHashSet: public std::conditional< is_char_ptr<K>::value,
    klib::KHashSetCached< K, Hash,  Eq, khInt_t >,
	klib::KHashSet< K, Hash,  Eq, khInt_t > >::type  {
protected:
	khInt_t i_iter=0;
public:
	inline khInt_t Add(const K ky) { // return -1 if the key already exists
		int absent=-1;
		khInt_t i=this->put(ky, &absent);
		if (absent==1) //key was actually added
		   return i;
		return -1;
	}

	inline khInt_t Remove(K ky) { //return index being removed, or -1 if no such key exists
		khInt_t i=this->get(ky);
		  if (i!=this->end()) {
			  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	inline void Clear() {
		this->clear(); //does not shrink
	}

	inline void Reset() {
		this->clear();
		GFREE(this->used); GFREE(this->keys);
		this->bits=0; this->count=0;
	}

    ~GHashSet() {
    	this->Reset();
    }

	inline bool operator[](K ky) { //RH only (read-only), cannot assign (use Add instead)
		return (this->get(ky)!=this->end());
	}

	inline bool hasKey(K ky) {
		return (this->get(ky)!=this->end());
	}

	int Find(K ky) {//return internal slot location if found,
	                // or -1 if not found
	  khInt_t r=this->get(ky);
  	  if (r==this->end()) return -1;
  	  return (int)r;
	}

	void startIterate() { //iterator-like initialization
	  i_iter=0;
	}

	K* Next() {
		//returns a pointer to next valid key in the table (NULL if no more)
		if (this->count==0) return NULL;
		uint32_t nb=this->n_buckets();
		while (i_iter<nb && !this->_used(i_iter)) i_iter++;
		if (i_iter==nb) return NULL;
		K* k=&(this->key(i_iter-1));
		++i_iter;
		return k;
	}

	inline uint32_t Count() { return this->count; }

};

// GStrSet always allocates a new copy of each added string;
//  if you don't want that, just use GHashSet<const char*> instead and manage the key allocation separately
template <class Hash=GHashKey_wyHash<const char*>, class Eq=GHashKey_Eq<const char*>, typename khInt_t=uint64_t>
  class GStrSet: public GHashSet<const char*, Hash, Eq, khInt_t> {
  protected:
	const char* lastKey=NULL;
  public:
	inline int Add(const char* ky) { // return -1 if the key already exists
		int absent=-1;
		khInt_t i=this->put(ky, &absent);
		if (absent==1) {//key was actually added
		   const char* s=Gstrdup(ky);
		   this->key(i)=s; //store a copy of the key string
		   lastKey=s;
		   return i;
		}
		//key was already there
		return -1;
	}

	inline const char* getLastKey() { return lastKey; }

	int Remove(const char* ky) { //return index being removed, or -1 if no such key exists
		  khInt_t i=this->get(ky);
		  if (i!=this->end()) {
			  const char* s=this->key(i);
			  if (s==lastKey) lastKey=NULL;
			  GFREE(s); //free string copy
			  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	inline void Clear() {
		khInt_t nb=this->n_buckets();
		for (khInt_t i = 0; i != nb; ++i) {
			if (!this->_used(i)) continue;
			//deallocate string copy
			GFREE(this->key(i));
		}
		lastKey=NULL;
		this->clear(); //does not shrink !
	}

	inline void Reset() {
		this->Clear();
		GFREE(this->used); GFREE(this->keys);
		lastKey=NULL;
		this->bits=0; this->count=0;
	}

    ~GStrSet() {
    	this->Reset();
    }

};

// Generic hash map where keys and values can be of any type
// Note: keys are always copied (shared) as simple value, there is no deep copy/allocation for pointers
//     so pointer keys must me managed separately
// Note: pointer values are automatically deallocated on container destruction by default,
//         use GHashMap(false) to disable that when V is a pointer
template <class K, class V, class Hash=GHashKey_wyHash<K>, class Eq=GHashKey_Eq<K>, typename khInt_t=uint64_t>
  class GHashMap:public std::conditional< is_char_ptr<K>::value,
    klib::KHashMapCached< K, V, Hash,  Eq, khInt_t>,
    klib::KHashMap< K, V, Hash,  Eq, khInt_t> >::type  {
protected:
	khInt_t i_iter=0;
	bool freeItems=false;
public:
	//---- these should be reimplemented for GHash
	inline int Add(const K ky, const V val) { // if a key does not exist allocate a copy of the key
		// return -1 if the key already exists
		int absent=-1;
		khInt_t i=this->put(ky, &absent);
		if (absent==1) { //key was actually added
			this->value(i)=val; //value is always copied
			return i;
		}
		return -1;
	}
	template <typename T=V> inline
		typename std::enable_if< std::is_pointer<T>::value, int>::type
			Remove(K ky) { //return index being removed
		khInt_t i=this->get(ky);
		  if (i!=this->end()) {
			  if (freeItems) delete this->value(i);
			  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	template <typename T=V> inline
		typename std::enable_if< !std::is_pointer<T>::value, int>::type
			Remove(K ky) { //return index being removed
		  khInt_t i=this->get(ky);
		  if (i!=this->end()) {
			  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	template <typename T=V> inline
		typename std::enable_if< std::is_pointer<T>::value, void>::type
		Clear() {
		if (!freeItems) {
			this->clear(); //does not shrink !
			return;
		}
		khInt_t nb=this->n_buckets();
		for (khInt_t i = 0; i != nb; ++i) {
			if (!this->_used(i)) continue;
			if (freeItems) delete this->value(i);
		}
		this->clear();
	}

	template <typename T=V> inline
		typename std::enable_if< !std::is_pointer<T>::value, void>::type
		Clear() { this->clear(); }

	inline void Reset() {
		this->Clear();
		GFREE(this->used); GFREE(this->keys);
		this->bits=0; this->count=0;
	}

    ~GHashMap() {
    	this->Reset();
    }

  // -- these can be shared with GHash:
	GHashMap(bool doFree=std::is_pointer<V>::value):freeItems(doFree) {
		static_assert(std::is_trivial<K>::value,
				"Error: cannot use this for non-trivial types!\n");
		if (!std::is_pointer<V>::value) doFree=false;
	};
	//return pointer to stored value if found, NULL otherwise
	// if the stored value is a pointer, it's going to be a pointer to that
	template <typename T=V> inline
			typename std::enable_if< std::is_pointer<T>::value, T>::type
	        Find(const K ky) {
	  khInt_t r=this->get(ky);
	  if (r==this->end()) return NULL;
	  return this->value(r);
	}

	template <typename T=V> inline
			typename std::enable_if< !std::is_pointer<T>::value, T*>::type
	        Find(const K ky) {
	  khInt_t r=this->get(ky);
	  if (r==this->end()) return NULL;
	  return &(this->value(r));
	}

	//-- operator[] should be defined just like Find?
	template <typename T=V> inline
			typename std::enable_if< std::is_pointer<T>::value, T>::type
	        operator[](const K ky) {
	  khInt_t r=this->get(ky);
	  if (r==this->end()) return NULL;
	  return this->value(r);
	}

	template <typename T=V> inline
			typename std::enable_if< !std::is_pointer<T>::value, T*>::type
	        operator[](const K ky) {
	  khInt_t r=this->get(ky);
	  if (r==this->end()) return NULL;
	  return &(this->value(r));
	}

	inline bool hasKey(K ky) {
		return (this->get(ky)!=this->end());
	}

	inline void startIterate() { //iterator-like initialization
	  i_iter=0;
	}

	template <typename T=K> inline
				typename std::enable_if< !std::is_pointer<T>::value, T*>::type
	  Next (V& val) {
		//returns a pointer to next key entry in the table (NULL if no more)
		if (this->count==0) return NULL;
		khInt_t nb=this->n_buckets();
		while (i_iter<nb && !this->_used(i_iter)) i_iter++;
		if (i_iter==nb) return NULL;
		val=this->value(i_iter);
		K* k=&(this->key(i_iter));
		++i_iter;
		return k;
	}

	template <typename T=K> inline
				typename std::enable_if< std::is_pointer<T>::value, T>::type
	  Next (V& val) {
		//returns a pointer to next key entry in the table (NULL if no more)
		if (this->count==0) return NULL;
		khInt_t nb=this->n_buckets();
		while (i_iter<nb && !this->_used(i_iter)) i_iter++;
		if (i_iter==nb) return NULL;
		val=this->value(i_iter);
		K k = this->key(i_iter);
		++i_iter;
		return k;
	}

	template <typename T=V> inline
				typename std::enable_if< !std::is_pointer<T>::value, T*>::type
	  NextData () {
		//returns a pointer to next key entry in the table (NULL if no more)
		if (this->count==0) return NULL;
		khInt_t nb=this->n_buckets();
		while (i_iter<nb && !this->_used(i_iter)) i_iter++;
		if (i_iter==nb) return NULL;
		T* val=&(this->value(i_iter));
		++i_iter;
		return val;
	}

	template <typename T=V> inline
				typename std::enable_if< std::is_pointer<T>::value, T>::type
	  NextData () {
		//returns a pointer to next key entry in the table (NULL if no more)
		if (this->count==0) return NULL;
		khInt_t nb=this->n_buckets();
		while (i_iter<nb && !this->_used(i_iter)) i_iter++;
		if (i_iter==nb) return NULL;
		T val=this->value(i_iter);
		++i_iter;
		return val;
	}

	inline uint32_t Count() { return this->count; }

};

// GHash<VType>(doFree=true) -- basic string hashmap
// Note: this hash map always makes a copy of the string key which can be costly
//     use GHashMap<const char*, VTYPE> for a faster alternative
template <class V, class Hash=GHashKey_wyHash<const char*>, class Eq=GHashKey_Eq<const char*>, typename khInt_t=uint64_t >
  class GHash:public GHashMap<const char*, V, Hash, Eq, khInt_t>  {
protected:
  const char* lastKey=NULL;
public:
	GHash(bool doFree=true) {
		this->freeItems=doFree;
	};
	//---- these should be now reimplemented
	inline int Add(const char* ky, const V val) { // if a key does not exist allocate a copy of the key
		// return -1 if the key already exists
		int absent=-1;
		khInt_t i=this->put(ky, &absent);
		if (absent==1) { //key was actually added
			const char* s=Gstrdup(ky);
			this->key(i)=s; //store a copy of the key string
			lastKey=s;
			this->value(i)=val; //value is always copied
			return i;
		}
		return -1;
	}

	inline const char* getLastKey() { return lastKey; }

	template <typename T=V> inline
		typename std::enable_if< std::is_pointer<T>::value, int>::type
	Remove(const char* ky) { //return index being removed
		  khInt_t i=this->get(ky);
		  if (i!=this->end()) {
			  const char* s=this->key(i);
			  if (s==lastKey) lastKey=NULL;
			  GFREE(s); //free string copy
			  if (this->freeItems) delete this->value(i);
	    	  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	template <typename T=V> inline
		typename std::enable_if< !std::is_pointer<T>::value, int>::type
	Remove(const char* ky) { //return index being removed
		  khInt_t i=this->get(ky);
		  if (i!=this->end()) {
			  const char* s=this->key(i);
			  if (s==lastKey) lastKey=NULL;
			  GFREE(s); //free string copy
	    	  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	template <typename T=V> inline
		typename std::enable_if< std::is_pointer<T>::value, void>::type
		Clear() {
		khInt_t nb=this->n_buckets();
		for (khInt_t i = 0; i != nb; ++i) {
			if (!this->_used(i)) continue;
			if (this->freeItems) delete this->value(i);
			GFREE(this->key(i));
		}
		lastKey=NULL;
		this->clear();
	}

	template <typename T=V> inline
		typename std::enable_if< !std::is_pointer<T>::value, void>::type
		Clear() {
		khInt_t nb=this->n_buckets();
		for (khInt_t i = 0; i != nb; ++i) {
			if (! this->_used(i) ) continue;
			GFREE(this->key(i));
		}
		lastKey=NULL;
		this->clear();
	}

	inline void Reset() {
		this->Clear();
		GFREE(this->used); GFREE(this->keys);
		this->bits=0; this->count=0;
	}

    ~GHash() {
    	this->Reset();
    }
};

template<typename T>
 using GIntHash = GHashMap<int, T, GHashKey_xxHash32<int>, GHashKey_Eq<int>, uint32_t>;

#endif
