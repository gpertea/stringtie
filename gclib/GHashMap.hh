/********************************************************************************
 *                  Hash map class templates
 *********************************************************************************/

#ifndef GHashMap_HH
#define GHashMap_HH
#include "GBase.h"
#include "city.h"
#include "khashl.hh"
#include <type_traits>
#include <typeinfo>
//#include <type_traits> //for std::is_trivial

/*
struct cstr_eq {
    inline bool operator()(const char* x, const char* y) const {
      return (strcmp(x, y) == 0);
    }
};

struct cstr_hash {
   inline uint32_t operator()(const char* s) const {
      return CityHash32(s, std::strlen(s));
   }
};
*/

#define GSTR_HASH(s) strhash(s)
//#define GSTR_HASH(s) djb_hash(s)
//#define GSTR_HASH(s) fnv1a_hash(s)
//#define GSTR_HASH(s) murmur3(s)

//template<typename T>
// using T_Ptr = typename std::conditional< std::is_pointer<T>::value, T, T* >::type;

/* defined it in GBase.h
template< typename T >
 struct is_char_ptr
  : std::enable_if< std::is_same< T, char* >::value ||
                    std::is_same< T, const char* >::value >
  {};
*/

template <typename K> struct GHashKey_Hash { //K generic (class, primitive, pointer except const char* )
  //template <typename T=K> inline typename std::enable_if< std::is_trivial<T>::value, uint32_t>::type
 uint32_t operator()(const K& s) const { //only works for trivial types!
      static_assert(std::is_trivial<K>::value, "Error: cannot use this for non-trivial types!\n");
      return CityHash32((const char *) &s, sizeof(K));
    }
};

template <> struct GHashKey_Hash<const char*> {
   inline uint32_t operator()(const char* s) const {
      return CityHash32(s, strlen(s));
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

//GHashSet is never making a deep copy of the char* key, it only stores the pointer
template <typename K=const char*, class Hash=GHashKey_Hash<K>, class Eq=GHashKey_Eq<K> >
  class GHashSet: public std::conditional< is_char_ptr<K>::value,
    klib::KHashSetCached< K, Hash,  Eq >,
	klib::KHashSet< K, Hash,  Eq > >::type  {
protected:
	uint32_t i_iter=0;
public:
	inline int Add(const K ky) { // return -1 if the key already exists
		int absent=-1;
		uint32_t i=this->put(ky, &absent);
		if (absent==1) //key was actually added
		   return i;
		return -1;
	}

	int Remove(K ky) { //return index being removed, or -1 if no such key exists
		  uint32_t i=this->get(ky);
		  if (i!=this->end()) {
			  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	inline void Clear() {
		this->clear(); //does not shrink !
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
  	  uint32_t r=this->get(ky);
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
		while (i_iter<nb && !this->occupied(i_iter)) i_iter++;
		if (i_iter==nb) return NULL;
		K* k=&(this->key(i_iter-1));
		++i_iter;
		return k;
	}

	inline uint32_t Count() { return this->count; }

};

//GStrSet always allocates a copy of each added string;
// if you don't want that (keys are shared), just use GHashSet<const char*> instead
template <class Hash=GHashKey_Hash<const char*>, class Eq=GHashKey_Eq<const char*> >
  class GStrSet: public GHashSet<const char*, Hash, Eq> {
  public:
	inline int Add(const char* ky) { // return -1 if the key already exists
		int absent=-1;
		uint32_t i=this->put(ky, &absent);
		if (absent==1) {//key was actually added
		   const char* s=Gstrdup(ky);
		   this->key(i)=s; //store a copy of the key string
		   return i;
		}
		//key was already there
		return -1;
	}

	int Remove(const char* ky) { //return index being removed, or -1 if no such key exists
		  uint32_t i=this->get(ky);
		  if (i!=this->end()) {
			  GFREE(this->key(i)); //free string copy
			  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	inline void Clear() {
		int nb=this->n_buckets();
		for (int i = 0; i != nb; ++i) {
			if (!this->__kh_used(this->used, i)) continue;
			//deallocate string copy
			GFREE(this->key(i));
		}
		this->clear(); //does not shrink !
	}

	inline void Reset() {
		this->Clear();
		GFREE(this->used); GFREE(this->keys);
		this->bits=0; this->count=0;
	}

    ~GStrSet() {
    	this->Reset();
    }

};

//generic hash map where keys and values can be of any type
template <class K, class V, class Hash=GHashKey_Hash<K>, class Eq=GHashKey_Eq<K> >
  class GHashMap:public std::conditional< is_char_ptr<K>::value,
    klib::KHashMapCached< K, V, Hash,  Eq >,
    klib::KHashMap< K, V, Hash,  Eq > >::type  {
protected:
	uint32_t i_iter=0;
	bool freeItems=false;
public:
	//---- these should be reimplemented for GHash
	inline int Add(const K ky, const V val) { // if a key does not exist allocate a copy of the key
		// return -1 if the key already exists
		int absent=-1;
		uint32_t i=this->put(ky, &absent);
		if (absent==1) { //key was actually added
			this->value(i)=val; //value is always copied
			return i;
		}
		return -1;
	}
	template <typename T=V> inline
		typename std::enable_if< std::is_pointer<T>::value, int>::type
			Remove(K ky) { //return index being removed
		  uint32_t i=this->get(ky);
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
		  uint32_t i=this->get(ky);
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
		int nb=this->n_buckets();
		for (int i = 0; i != nb; ++i) {
			if (!this->__kh_used(this->used, i)) continue;
			if (freeItems) delete this->value(i);
		}
		this->clear();
	}

	template <typename T=V> inline
		typename std::enable_if< !std::is_pointer<T>::value, void>::type
		Clear() {
		if (!freeItems) {
			this->clear(); //does not shrink !
			return;
		}
		int nb=this->n_buckets();
		for (int i = 0; i != nb; ++i) {
			if (!this->__kh_used(this->used, i)) continue;
		}
		this->clear();
	}

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
	  uint32_t r=this->get(ky);
	  if (r==this->end()) return NULL;
	  return this->value(r);
	}

	template <typename T=V> inline
			typename std::enable_if< !std::is_pointer<T>::value, T*>::type
	        Find(const K ky) {
	  uint32_t r=this->get(ky);
	  if (r==this->end()) return NULL;
	  return &(this->value(r));
	}

	//-- operator[] should be defined just like Find?
	template <typename T=V> inline
			typename std::enable_if< std::is_pointer<T>::value, T>::type
	        operator[](const K ky) {
	  uint32_t r=this->get(ky);
	  if (r==this->end()) return NULL;
	  return this->value(r);
	}

	template <typename T=V> inline
			typename std::enable_if< !std::is_pointer<T>::value, T*>::type
	        operator[](const K ky) {
	  uint32_t r=this->get(ky);
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
		uint32_t nb=this->n_buckets();
		while (i_iter<nb && !this->occupied(i_iter)) i_iter++;
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
		uint32_t nb=this->n_buckets();
		while (i_iter<nb && !this->occupied(i_iter)) i_iter++;
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
		uint32_t nb=this->n_buckets();
		while (i_iter<nb && !this->occupied(i_iter)) i_iter++;
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
		uint32_t nb=this->n_buckets();
		while (i_iter<nb && !this->occupied(i_iter)) i_iter++;
		if (i_iter==nb) return NULL;
		T val=this->value(i_iter);
		++i_iter;
		return val;
	}



	inline uint32_t Count() { return this->count; }

};

template <class V, class Hash=GHashKey_Hash<const char*>, class Eq=GHashKey_Eq<const char*> >
  class GHash:public GHashMap<const char*, V, Hash, Eq>  {
protected:

public:
	GHash(bool doFree=true) {
		this->freeItems=doFree;
	};
	//---- these should be now reimplemented
	inline int Add(const char* ky, const V val) { // if a key does not exist allocate a copy of the key
		// return -1 if the key already exists
		int absent=-1;
		uint32_t i=this->put(ky, &absent);
		if (absent==1) { //key was actually added
			const char* s=Gstrdup(ky);
			this->key(i)=s; //store a copy of the key string
			this->value(i)=val; //value is always copied
			return i;
		}
		return -1;
	}
	template <typename T=V> inline
		typename std::enable_if< std::is_pointer<T>::value, int>::type
	Remove(const char* ky) { //return index being removed
		  uint32_t i=this->get(ky);
		  if (i!=this->end()) {
			  GFREE(this->key(i)); //free string copy
			  if (this->freeItems) delete this->value(i);
	    	  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	template <typename T=V> inline
		typename std::enable_if< !std::is_pointer<T>::value, int>::type
	Remove(const char* ky) { //return index being removed
		  uint32_t i=this->get(ky);
		  if (i!=this->end()) {
			  GFREE(this->key(i)); //free string copy
	    	  this->del(i);
	    	  return i;
		  }
		  return -1;
	}

	template <typename T=V> inline
		typename std::enable_if< std::is_pointer<T>::value, void>::type
		Clear() {
		int nb=this->n_buckets();
		for (int i = 0; i != nb; ++i) {
			if (!this->__kh_used(this->used, i)) continue;
			if (this->freeItems) delete this->value(i);
			GFREE(this->key(i));
		}
		this->clear();
	}

	template <typename T=V> inline
		typename std::enable_if< !std::is_pointer<T>::value, void>::type
		Clear() {
		int nb=this->n_buckets();
		for (int i = 0; i != nb; ++i) {
			if (!this->__kh_used(this->used, i)) continue;
			GFREE(this->key(i));
		}
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

#endif
