//---------------------------------------------------------------------------
/*
Sortable collection of pointers to objects
*/

#ifndef _GVec_HH
#define _GVec_HH

#include "GBase.h"
#include <algorithm>

#define GVEC_INDEX_ERR "GVec error: invalid index: %d\n"
 #if defined(NDEBUG) || defined(NODEBUG) || defined(_NDEBUG) || defined(NO_DEBUG)
 #define TEST_INDEX(x)
#else
 #define TEST_INDEX(x) \
 if (x<0 || x>=fCount) GError(GVEC_INDEX_ERR, x)
#endif

#define GVEC_CAPACITY_ERR "GVec error: invalid capacity: %d\n"
#define GVEC_COUNT_ERR "GVec error: invalid count: %d\n"

#define MAXLISTSIZE INT_MAX-1

#define FREEDATA (fFreeProc!=NULL)

template <class OBJ> int DefLTCompareProc(pointer p1, pointer p2) {
 OBJ& o1 = *((OBJ*) p1);
 OBJ& o2 = *((OBJ*) p2);
 if (o1 < o2) return -1;
   else return ((o2 < o1) ? 1 : 0 );
}

//basic template for array of objects;
//so it doesn't require comparison operators to be defined
template <class OBJ> class GVec {
  protected:
    OBJ* fArray;
    int fCount;
    int fCapacity;
    void qSort(int L, int R, GCompareProc* cmpFunc);
  public:
    GVec(int init_capacity=2);
    GVec(int init_count, const OBJ init_val);
    GVec(int init_count, OBJ* init_val, bool delete_initval=true); //convenience constructor for complex vectors
    GVec(const GVec<OBJ>& array); //copy constructor
    GVec(GVec<OBJ>&& array); //move constructor
    GVec<OBJ>& operator=(const GVec<OBJ>& array); //copy operator
    GVec<OBJ>& operator=(GVec<OBJ>&& array); //move operator
    virtual ~GVec();
    void Insert(int idx, OBJ item) { Insert(idx, &item); }
    void Insert(int idx, OBJ* item);
    void idxInsert(int idx, OBJ& item) { Insert(idx, &item); }
    void Grow();
    void Grow(int newCap);
    void Grow(int idx, OBJ& item); //grow and add/insert item copy
    void Reverse(); //WARNING: will break the sort order if SORTED!
    int Add(OBJ* item); // simply append to the end of fArray, reallocating as needed
    int Add(OBJ& item) { return Add(&item); }
    int cAdd(OBJ item) { return Add(&item); } //all these will CREATE a new OBJ and COPY to it
    //                   // using OBJ copy operator=
    // -- stack/queue usage:
    //int Push(OBJ& item) { return Add(&item); }
    int Push(OBJ& item) { return Add(&item); }
    int cPush(OBJ item) { return Add(&item); }
    OBJ Pop();// Stack use; removes and returns a copy of the last item
    OBJ Shift(); //Queue use: removes and returns a copy of the first item
    void Shift(int idx); //Queue use: removes first idx elements from array

    void Add(GVec<OBJ>& list); //append copies of all items from another list

    OBJ& Get(int idx) {
          TEST_INDEX(idx);
          return fArray[idx];
          }
    inline OBJ& operator[](int i) {
          TEST_INDEX(i);
          return fArray[i];
          }
    OBJ& Last() {
         TEST_INDEX(fCount-1);
         return fArray[fCount-1];
         }
    OBJ& First() {
         TEST_INDEX(0);
         return fArray[0];
         }
    void Clear();
    void Delete(int index);
    void Replace(int idx, OBJ& item); //Put, use operator= to copy
    void Exchange(int idx1, int idx2);
    void Swap(int idx1, int idx2)  { Exchange(idx1, idx2); }
    int  Capacity() { return fCapacity; }
    //this will reject identical items in sorted lists only!
    inline void setCapacity(int NewCapacity) {
      if (NewCapacity < fCount || NewCapacity > MAXLISTSIZE)
          GError(GVEC_CAPACITY_ERR, NewCapacity);
       //error: NewCapacity MUST be > fCount
       //if you want to shrink it use Resize() or setCount()
      if (NewCapacity!=fCapacity) {
        if (NewCapacity==0) {
           delete[] fArray;
           fArray=nullptr;
        }
        else {
            OBJ* oldArray=fArray;
    		fArray=new OBJ[NewCapacity];
    		if (oldArray) {
    		  std::move(&oldArray[0], &oldArray[fCount], & fArray[0]);
              delete[] oldArray;
    		}
    		//this relies on in-class initializers, default constructors etc.
        }
        fCapacity=NewCapacity;
      }
   }
    /*
    void insertGrow(int NewCapacity, int idx, OBJ& item) {
        OBJ* newList=new OBJ[NewCapacity];
        // operator= required!
        std::move(& fArray[0], & fArray[idx], & newList[0]);
        newList[idx]=item;
        //copy data after idx
        std::move(& fArray[idx], & fArray[fCount], & newList[idx+1]);
        delete[] fArray;
        fArray=newList;
        fCapacity=NewCapacity;
        fCount++;
    }
    */

    template <typename T=OBJ>
         typename std::enable_if< std::is_trivial<T>::value, void>::type
      topOff(int NewCount) {
        memset(fArray+fCount, 0, (NewCount-fCount)*sizeof(OBJ));
    }

    template <typename T=OBJ>
             typename std::enable_if< !std::is_trivial<T>::value, void>::type
      topOff(int NewCount) {
    	OBJ v;
    	if (NewCount>fCount)
    	 	for (int i=fCount;i<NewCount;i++)
    	   	  fArray[i]=v;
    }

    int  Count() { return fCount; }
    void setCount(int NewCount);         // will trim or expand the array as needed
    void setCount(int NewCount, OBJ v);
    void Resize(int NewCount) { setCount(NewCount); }
    void Resize(int NewCount, OBJ v) { setCount(NewCount, v); }

    //void Move(int curidx, int newidx);
    bool isEmpty() { return fCount==0; }
    bool notEmpty() { return fCount>0; }

    void Sort(GCompareProc* cmpFunc);
    void Sort();
};

//---- template for dynamic array of object pointers
//---- it's faster than GVec<OBJ*> and has item deallocation awareness
template <class OBJ> class GPVec {
  protected:
    OBJ* *fList; //pointer to an array of pointers to objects
    int fCount; //total number of entries in list
    int fCapacity; //current allocated size
    GFreeProc* fFreeProc; //useful for deleting objects
    //---
    void Expand();
    void Grow();
    void Grow(int idx, OBJ* item);
    void qSort(int L, int R, GCompareProc* cmpFunc);
  public:
    static void DefaultFreeProc(pointer item) {
      delete (OBJ*)item;
    }
    virtual ~GPVec();
    GPVec(int init_capacity=2, bool free_elements=true); //also the default constructor
    GPVec(bool free_elements);
    GPVec(const GPVec<OBJ>& list); //copy constructor
    GPVec(GPVec<OBJ>&& list); //move construstor
    GPVec(GPVec<OBJ>* list); //similar to a copy constructor
    GPVec<OBJ>& operator=(const GPVec<OBJ>& list);
    GPVec<OBJ>& operator=(GPVec<OBJ>&& list);//move assignment operator
    inline OBJ* Get(int i) {
    	TEST_INDEX(i);
        return fList[i];
    }
     //OBJ* operator[](int i) { return this->Get(i); }
    inline OBJ*& operator[](int i) {
    	 TEST_INDEX(i); return fList[i];
    }
    void Reverse(); //reverse pointer array; WARNING: will break(reverse) the sort order if sorted!
    void freeItem(int idx); //calls fFreeProc (or DefaultFreeProc) on fList[idx] and sets NULL there, doesn't pack!
                      //it will free even if fFreeProc is NULL!
    void setFreeItem(GFreeProc *freeProc) { fFreeProc=freeProc; }
    void setFreeItem(bool doFree) {
       if (doFree) fFreeProc=DefaultFreeProc;
             else  fFreeProc=NULL;
       }
    // -- stack usage:
    int Push(OBJ* item) { return Add(item); }
    OBJ* Pop();// Stack use; removes and returns last item,but does NOT FREE it
    OBJ* Shift(); //Queue use: removes and returns first item, but does NOT FREE it
    void deallocate_item(OBJ*& item); //forcefully call fFreeProc or delete on item
    void Clear();
    void Exchange(int idx1, int idx2);
    void Swap(int idx1, int idx2)  { Exchange(idx1, idx2); }
    OBJ* First() { return (fCount>0)?fList[0]:nullptr; }
    OBJ* Last()  { return (fCount>0)?fList[fCount-1]:nullptr;}
    bool isEmpty() { return fCount==0; }
    bool notEmpty() { return fCount>0; }
    int Capacity() { return fCapacity; }
    int Count()   { return fCount; }
    void setCapacity(int NewCapacity);
    void setCount(int NewCount); //the same as setCapacity() but the new item range is filled with NULLs
    int Add(OBJ* item); //simply append the pointer copy
    void Add(GPVec<OBJ>& list); //add all pointers from another list
    void addNew(GPVec<OBJ>& list); //add new OBJ copies of items in another list
    void Insert(int idx, OBJ* item);
    void Move(int curidx, int newidx);
    void Put(int idx, OBJ* item);
    void Pack();
    void Delete(int index); //also frees the item if fFreeProc!=NULL, and shifts the successor items
    void Forget(int idx); //simply places a NULL at fList[idx], nothing else
    int RemovePtr(pointer item); //always use linear search to find the pointer! calls Delete() if found
    int IndexOf(pointer item); //a linear search for pointer address!
    void Sort(GCompareProc* cmpFunc);
    void Sort();
 };

//-------------------- TEMPLATE IMPLEMENTATION-------------------------------

template <class OBJ> GVec<OBJ>::GVec(int init_capacity) {
  fCount=0;
  fCapacity=0;
  fArray=nullptr;
  setCapacity(init_capacity);
  //if (set_count) fCount = init_capacity;
}


template <class OBJ> GVec<OBJ>::GVec(int init_count, const OBJ init_val) {
  fCount=0;
  fCapacity=0;
  fArray=nullptr;
  setCapacity(init_count);
  fCount = init_count;
  for (int i=0;i<fCount;i++)
    fArray[i]=init_val;
}

template <class OBJ> GVec<OBJ>::GVec(int init_count, OBJ* init_val, bool delete_initval) {
	  fCount=0;
	  fCapacity=0;
	  fArray=nullptr;
	  setCapacity(init_count);
	  fCount = init_count;
	  for (int i=0;i<fCount;i++)
	    fArray[i]=*init_val;
	  if (delete_initval) { delete init_val; }
}


template <class OBJ> GVec<OBJ>::GVec(const GVec<OBJ>& array) { //copy constructor
 fCount=array.fCount;
 fCapacity=array.fCapacity;
 fArray=nullptr;
 if (fCapacity>0) {
     fArray=new OBJ[this->fCapacity];
     std::copy(& array.fArray[0], & array.fArray[fCount], & fArray[0]);
   }
 }

template <class OBJ> GVec<OBJ>::GVec(GVec<OBJ>&& array) { //move constructor
  fCount=array.fCount;
  fCapacity=array.fCapacity;
  fArray=array.fArray;
  array.fArray=nullptr;
  array.fCount=0;
  array.fCapacity=0;
 }

template <class OBJ> GVec<OBJ>& GVec<OBJ>::operator=(const GVec<OBJ>& array) {
 if (&array==this) return *this;
 Clear();
 fCapacity=array.fCapacity;
 fCount=array.fCount;
 if (fCapacity>0) {
    fArray=new OBJ[this->fCapacity];
    // uses OBJ operator=
    std::copy(& array.fArray[0], & array.fArray[fCount], & fArray[0]);
 }
 return *this;
}

template <class OBJ> GVec<OBJ>& GVec<OBJ>::operator=(GVec<OBJ>&& array) {
 if (&array==this) return *this;
 Clear();
 fCapacity=array.fCapacity;
 fCount=array.fCount;
 fArray=array.fArray;
 array.fArray=nullptr;
 array.fCapacity=0;
 array.fCount=0;
 return *this;
}

template <class OBJ> GVec<OBJ>::~GVec() {
 this->Clear();
}

template <class OBJ> void GVec<OBJ>::Clear() {
  fCount=0;
  delete[] fArray;
  fArray=nullptr;
  fCapacity=0;
}

template <class OBJ> void GVec<OBJ>::Grow() {
 int delta = (fCapacity>8) ? (fCapacity>>2) : 2 ;
 setCapacity(fCapacity + delta);
}

template <class OBJ> void GVec<OBJ>::Grow(int newCap) {
 if (newCap<fCapacity) return; //never shrink with this method
 int delta=newCap-fCapacity+2;
 setCapacity(fCapacity + delta);
}

template <class OBJ> void GVec<OBJ>::Reverse() {
  int l=0;
  int r=fCount-1;
  OBJ c;
  while (l<r) {
     c=fArray[l];fArray[l]=fArray[r];
     fArray[r]=c;
     l++;r--;
  }
}

template <class OBJ> void GVec<OBJ>::Grow(int idx, OBJ& item) {
 int delta = (fCapacity>8) ? (fCapacity>>2) : 1 ;
 int NewCapacity=fCapacity+delta;
 if (NewCapacity <= fCount || NewCapacity >= MAXLISTSIZE)
    GError(GVEC_CAPACITY_ERR, NewCapacity);
    //error: capacity not within range
 if (idx==fCount) { //append item
         setCapacity(NewCapacity);
         fArray[idx]=item;
         fCount++;
         return;
 }
 //expands and inserts item at idx at the same time
 //insertGrow(NewCapacity, idx, item);
 OBJ* newList=new OBJ[NewCapacity];
 // move operator= required!
 std::move(& fArray[0], & fArray[idx], & newList[0]);
 newList[idx]=item; //operator=
 //copy data after idx
 std::move(& fArray[idx], & fArray[fCount], & newList[idx+1]);
 delete[] fArray;
 fArray=newList;
 fCapacity=NewCapacity;
 fCount++;
}
template <class OBJ> int GVec<OBJ>::Add(OBJ* item) {
 if (item==NULL) return -1;
 if (fCount==fCapacity) Grow();
 fArray[fCount] = *item; //OBJ::operator= must copy OBJ properly!
 fCount++;
 return fCount-1;
}


template <class OBJ> void GVec<OBJ>::Add(GVec<OBJ>& list) {
  if (list.Count()==0) return;
  //simply copy
  setCapacity(fCapacity+list.fCount);
  if (std::is_trivial<OBJ>::value) {
    memcpy( &fArray[fCount], list.fArray, list.fCount*sizeof(OBJ));
    }
   else {
    for (int i=0;i<list.fCount;i++)
          fArray[fCount+i]=list.fArray[i];
    }
  fCount+=list.fCount;
}

//Stack usage:
template <class OBJ> OBJ GVec<OBJ>::Pop() {
 if (fCount<=0) GError("Error: invalid GVec::Pop() operation!\n");
 fCount--;
 //OBJ o(fArray[fCount]); //copy constructor
 //o=fList[fCount];
 //fArray[fCount]=NULL;
 return fArray[fCount]; //copy of the last element (copy constructor called)
}

//Queue usage:
template <class OBJ> OBJ GVec<OBJ>::Shift() {
 if (fCount<=0) GError("Error: invalid GVec::Shift() operation!\n");
 fCount--;
 OBJ o(fArray[0]); //copy constructor
 if (fCount>0)
   memmove(&fArray[0], &fArray[1], (fCount)*sizeof(OBJ));
 //fList[fCount]=NULL; //not that it matters..
 return o;
}

template <class OBJ> void GVec<OBJ>::Shift(int idx) {
 if (idx<=0 || fCount-idx<=0) GError("Error: invalid GVec::Shift() operation!\n");
 fCount-=idx;
 if (fCount>0)
   memmove(&fArray[0], &fArray[idx], (fCount)*sizeof(OBJ));
}

template <class OBJ> void GVec<OBJ>::Insert(int idx, OBJ* item) {
 //idx must be the new position this new item must have
 //so the allowed range is [0..fCount]
 //the old idx item all the above will be shifted to idx+1
 if (idx<0 || idx>fCount) GError(GVEC_INDEX_ERR, idx);
 if (fCount==fCapacity) { //need to resize the array
    //expand and also copy/move data and insert the new item
    //more efficient than copying everything and then moving elements
    Grow(idx, *item);
    return;
 }
 //move data around to make room for the new item
 if (idx<fCount) {
   //copy after-idx items (shift up)
   //this->shiftUp(idx);
	std::move_backward(& fArray[idx],  & fArray[fCount], & fArray[fCount+1]);
 }
 fArray[idx]=*item;
 fCount++;
}


/*template <class OBJ> void GVec<OBJ>::Move(int curidx, int newidx) { //swap
 if (curidx!=newidx || newidx>=fCount)
     GError(GVEC_INDEX_ERR, newidx);
 OBJ tmp=fArray[curidx]; //copy constructor here
 fArray[curidx]=fArray[newidx];
 fArray[newidx]=tmp;
}*/


template <class OBJ> void GVec<OBJ>::Replace(int idx, OBJ& item) {
 TEST_INDEX(idx);
 fArray[idx]=item;
}

template <class OBJ> void GVec<OBJ>::Exchange(int idx1, int idx2) {
 TEST_INDEX(idx1);
 TEST_INDEX(idx2);
 OBJ item=fArray[idx1];
 fArray[idx1]=fArray[idx2];
 fArray[idx2]=item;
}


template <class OBJ> void GVec<OBJ>::Delete(int idx) {
 TEST_INDEX(idx);
 std::move(& fArray[idx+1], & fArray[fCount], & fArray[idx]);
 fCount--;
 /*
 if (std::is_trivial<OBJ>::value) {
   if (index<fCount)
    //move higher elements if any (shift down)
      memmove(&fArray[index], &fArray[index+1], (fCount-index)*sizeof(OBJ));
   }
 else {
   while (index<fCount) {
      fArray[index]=fArray[index+1];
      index++;
      }
  }
  */
}

template <class OBJ> void GVec<OBJ>::setCount(int NewCount) {
	if (NewCount<0 || NewCount > MAXLISTSIZE)
	   GError(GVEC_COUNT_ERR, NewCount);
	//if (NewCount > fCapacity) setCapacity(NewCount);
	Grow(NewCount);
	if (NewCount>fCount) this->topOff(NewCount);
    fCount = NewCount; //new items should be populated by the default object constructor(!)
}
/*
template <class OBJ> void GVec<OBJ>::setCount(int NewCount, OBJ* v) {
	if (NewCount<0 || NewCount > MAXLISTSIZE)
	  GError(GVEC_COUNT_ERR, NewCount);
	while (NewCount > fCapacity) Grow();
	if (NewCount>fCount) {
		for (int i=fCount;i<NewCount;i++)
		  fArray[i]=*v;
	}
	fCount = NewCount;
}
*/
template <class OBJ> void GVec<OBJ>::setCount(int NewCount, OBJ v) {
	if (NewCount<0 || NewCount > MAXLISTSIZE)
	   GError(GVEC_COUNT_ERR, NewCount);
	Grow(NewCount);
	if (NewCount>fCount) {
		for (int i=fCount;i<NewCount;i++)
		  fArray[i]=v;
	}
	fCount = NewCount;
}

template <class OBJ> void GVec<OBJ>::qSort(int l, int r, GCompareProc* cmpFunc) {
 int i, j;
 OBJ p,t;
 do {
    i = l; j = r;
    p = this->fArray[l+ ((r - l) >> 1)];
    do {
      while (cmpFunc(&(this->fArray[i]), &p) < 0) i++;
      while (cmpFunc(&(this->fArray[j]), &p) > 0) j--;
      if (i <= j) {
        t = this->fArray[i];
        this->fArray[i] = this->fArray[j];
        this->fArray[j] = t;
        i++; j--;
        }
      } while (i <= j);
    if (l < j) qSort(l, j, cmpFunc);
    l = i;
    } while (i < r);
}

template <class OBJ> void GVec<OBJ>::Sort(GCompareProc* cmpFunc) {
 if (cmpFunc==NULL) {
   GMessage("Warning: NULL compare function given, useless Sort() call.\n");
   return;
 }
 if (this->fArray!=nullptr && this->fCount>0)
     qSort(0, this->fCount-1, cmpFunc);
}

template <class OBJ> void GVec<OBJ>::Sort() {
  GCompareProc* cmpFunc = DefLTCompareProc<OBJ>;
  Sort(cmpFunc);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//*=> GPVec implementation

template <class OBJ> GPVec<OBJ>::GPVec(const GPVec& list) { //copy constructor
 fCount=list.fCount;
 fCapacity=list.fCapacity;
 fList=nullptr;
 fFreeProc=list.fFreeProc;
 fCount=list.fCount;
 if (fCapacity>0) {
    //GMALLOC(fList, fCapacity*sizeof(OBJ*));
	fList=new OBJ*[fCapacity];
	std::copy(&list.fList[0], &list.fList[fCount], &fList[0]);
    //memcpy(fList, list.fList, fCount*sizeof(OBJ*));
 }
}

template <class OBJ> GPVec<OBJ>::GPVec(GPVec&& list) { //copy constructor
 fCount=list.fCount;
 fCapacity=list.fCapacity;
 fList=list.fList;
 fFreeProc=list.fFreeProc;
 list.fCount=0;
 list.fCapacity=0;
 list.fList=nullptr;
}

template <class OBJ> GPVec<OBJ>::GPVec(GPVec* plist) { //another copy constructor
  fCount=0;
  fCapacity=plist->fCapacity;
  fList=nullptr;
  fFreeProc=plist->fFreeProc;
  fCount=plist->fCount;
  if (fCapacity>0) {
    //GMALLOC(fList, fCapacity*sizeof(OBJ*));
    fList=new OBJ*[fCapacity];
    std::copy(& (plist->fList[0]), & (plist->fList[fCount]), &fList[0]);
    //memcpy(fList, plist->fList, fCount*sizeof(OBJ*));
  }
}

template <class OBJ> GPVec<OBJ>& GPVec<OBJ>::operator=(const GPVec& list) {
 if (&list!=this) {
     Clear();
     fFreeProc=list.fFreeProc;
     //Attention: only the *POINTERS* are copied,
     // the actual objects are NOT duplicated
     fCount=list.fCount;
     fCapacity=list.fCapacity;
     if (fCapacity>0) {
        //GMALLOC(fList, fCapacity*sizeof(OBJ*));
        //memcpy(fList, list.fList, fCount*sizeof(OBJ*));
   		fList=new OBJ*[fCapacity];
   		std::copy(&list.fList[0], &list.fList[fCount], &fList[0]);
        }
     }
 return *this;
}

template <class OBJ> GPVec<OBJ>& GPVec<OBJ>::operator=(GPVec&& list) {
 if (&list!=this) {
     Clear();
     fFreeProc=list.fFreeProc;
     //Attention: only the *POINTERS* are copied,
     // the actual objects are NOT duplicated
     fCount=list.fCount;
     fCapacity=list.fCapacity;
     fList=list.fList;
     list.fList=nullptr;
     list.fCapacity=0;
     list.fCount=0;
     }
 return *this;
}


template <class OBJ> void GPVec<OBJ>::Add(GPVec<OBJ>& list) {
  if (list.Count()==0) return;
  //simply copy the pointers! -- the objects will be shared
  setCapacity(fCapacity+list.fCount);
  memcpy( & (fList[fCount]), list.fList, list.fCount*sizeof(OBJ*));
  fCount+=list.fCount;
}

template <class OBJ> void GPVec<OBJ>::addNew(GPVec<OBJ>& list) {
  if (list.Count()==0) return;
  //requires OBJ copy constructor
  setCapacity(fCapacity+list.fCount);
  //memcpy( & (fList[fCount]), list.fList, list.fCount*sizeof(OBJ*));
  for (int i=fCount;i<list.fCount;i++)
	  fList[i]=new OBJ(*list.fList[i-fCount]);
  fCount+=list.fCount;
}


template <class OBJ> void GPVec<OBJ>::Reverse() {
  int l=0;
  int r=fCount-1;
  OBJ* c;
  while (l<r) {
     c=fList[l];fList[l]=fList[r];
     fList[r]=c;
     l++;r--;
     }
}

template <class OBJ> GPVec<OBJ>::GPVec(int init_capacity, bool free_elements) {
  fCount=0;
  fCapacity=0;
  fList=nullptr;
  fFreeProc=(free_elements) ? DefaultFreeProc : NULL;
  if (init_capacity>0)
    setCapacity(init_capacity);
}

template <class OBJ> GPVec<OBJ>::GPVec(bool free_elements) {
  fCount=0;
  fCapacity=0;
  fList=nullptr;
  fFreeProc=(free_elements) ? DefaultFreeProc : NULL;
}

template <class OBJ> GPVec<OBJ>::~GPVec() {
 this->Clear();//this will free the items if fFreeProc is defined
}

template <class OBJ> void GPVec<OBJ>::setCapacity(int NewCapacity) {
  if (NewCapacity < fCount || NewCapacity > MAXLISTSIZE)
    GError(GVEC_CAPACITY_ERR, NewCapacity);
    //error: capacity not within range
  if (NewCapacity!=fCapacity) {
   if (NewCapacity==0) {
      //GFREE(fList);
	   delete[] fList;
       fList=nullptr;
      }
    else {
      //GREALLOC(fList, NewCapacity*sizeof(OBJ*));
        OBJ** oldList=fList;
  	    fList=new OBJ*[NewCapacity];
        if (oldList) {
		   std::move(&oldList[0], &oldList[fCount], & fList[0]);
		   delete[] oldList;
        }
     }
   fCapacity=NewCapacity;
   }
}

template <class OBJ> void GPVec<OBJ>::deallocate_item(OBJ* &item) {
 if (item==NULL) return;
 if (FREEDATA) {
   (*fFreeProc)(item);
   item=NULL;
   }
 else {
  delete item;
  item=NULL;
  }
}

template <class OBJ> void GPVec<OBJ>::Clear() {
 if (FREEDATA) {
   for (int i=0; i<fCount; i++) {
     (*fFreeProc)(fList[i]);
     }
   }
 //GFREE(fList);
 delete[] fList;
 fList=nullptr;
 fCount=0;
 fCapacity=0;
}

template <class OBJ> void GPVec<OBJ>::Exchange(int idx1, int idx2) {
 TEST_INDEX(idx1);
 TEST_INDEX(idx2);
 OBJ* item=fList[idx1];
 fList[idx1]=fList[idx2];
 fList[idx2]=item;
}

template <class OBJ> void GPVec<OBJ>::Expand() {
 if (fCount==fCapacity) Grow();
 //return this;
}

template <class OBJ> void GPVec<OBJ>::Grow() {
 	int delta = (fCapacity>8) ? (fCapacity>>2) : 1;
	setCapacity(fCapacity + delta);
}

template <class OBJ> void GPVec<OBJ>::Grow(int idx, OBJ* item) {
 int delta = (fCapacity>8) ? (fCapacity>>2) : 1 ;
 int NewCapacity=fCapacity+delta;
 if (NewCapacity <= fCount || NewCapacity > MAXLISTSIZE)
    GError(GVEC_CAPACITY_ERR, NewCapacity);
 if (idx==fCount) {
    //GREALLOC(fList, NewCapacity*sizeof(OBJ*));
	setCapacity(NewCapacity);
    fList[idx]=item;
    fCount++;
    return;
 }
 OBJ** newList=new OBJ*[NewCapacity];
 std::move(& fList[0], & fList[idx], & newList[0]);
 newList[idx]=item;
 //copy data after idx
 std::move(& fList[idx], & fList[fCount], & newList[idx+1]);
 delete[] fList;
 fList=newList;
 fCapacity=NewCapacity;
 fCount++;
}
 /*
 GMALLOC(newList, NewCapacity*sizeof(OBJ*));
 //copy data before idx
 memcpy(&newList[0],&fList[0], idx*sizeof(OBJ*));
 newList[idx]=newitem;
 //copy data after idx
 memmove(&newList[idx+1],&fList[idx], (fCount-idx)*sizeof(OBJ*));
 memset(&newList[fCount+1], 0, (NewCapacity-fCount-1)*sizeof(OBJ*));
 //data copied:
 GFREE(fList);
 fList=newList;
 fCount++;
 fCapacity=NewCapacity;
 */

template <class OBJ> int GPVec<OBJ>::IndexOf(pointer item) {
 for (int i=0;i<fCount;i++) {
     if (item==(pointer)fList[i]) return i;
     }
 return -1;
 }

template <class OBJ> int GPVec<OBJ>::Add(OBJ* item) {
 int result;
 if (item==NULL) return -1;
 result = fCount;
 if (result==fCapacity) this->Grow();
 fList[result]=item;
 fCount++;
 return fCount-1;
}

template <class OBJ> void GPVec<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 if (idx<0 || idx>fCount) GError(GVEC_INDEX_ERR, idx);
 if (fCount==fCapacity) {
   Grow(idx, item);
   return;
   }
 if (idx<fCount)
      memmove(&fList[idx+1], &fList[idx], (fCount-idx)*sizeof(OBJ*));
 fList[idx]=item;
 fCount++;
}

template <class OBJ> void GPVec<OBJ>::Move(int curidx, int newidx) { //s
 //BE_UNSORTED; //cannot do that in a sorted list!
 if (curidx!=newidx || newidx>=fCount)
     GError(GVEC_INDEX_ERR, newidx);
 OBJ* p;
 p=Get(curidx);
 //this is a delete:
 fCount--;
 if (curidx<fCount)
    memmove(&fList[curidx], &fList[curidx+1], (fCount-curidx)*sizeof(OBJ*));
 //-this was instead of delete
 Insert(newidx, p);
}

template <class OBJ> void GPVec<OBJ>::Put(int idx, OBJ* item) {
 //WARNING: this will never free the replaced item!
 TEST_INDEX(idx);
 fList[idx]=item;
}

template <class OBJ> void GPVec<OBJ>::Forget(int idx) {
 TEST_INDEX(idx);
 fList[idx]=NULL; //user should free that somewhere else
}

template <class OBJ> void GPVec<OBJ>::freeItem(int idx) {
  TEST_INDEX(idx);
  if (fFreeProc!=NULL) {
      (*fFreeProc)(fList[idx]);
      }
    else this->DefaultFreeProc(fList[idx]);
  fList[idx]=NULL;
}

template <class OBJ> void GPVec<OBJ>::Delete(int index) {
 TEST_INDEX(index);
 if (fFreeProc!=NULL && fList[index]!=NULL) {
   (*fFreeProc)(fList[index]); //freeItem
   }
 fList[index]=NULL;
 fCount--;
 if (index<fCount) //move higher elements if any
   memmove(&fList[index], &fList[index+1], (fCount-index)*sizeof(OBJ*));
}

//Stack usage:
template <class OBJ> OBJ* GPVec<OBJ>::Pop() {
 if (fCount<=0) return NULL;
 fCount--;
 OBJ* o=fList[fCount];
 fList[fCount]=NULL;
 return o;
}

//Queue usage:
template <class OBJ> OBJ* GPVec<OBJ>::Shift() {
 if (fCount<=0) return NULL;
 fCount--;
 OBJ* o=fList[0];
 if (fCount>0)
   memmove(&fList[0], &fList[1], (fCount)*sizeof(OBJ*));
 fList[fCount]=NULL; //not that it matters..
 return o;
}

//linear search for the pointer address
template <class OBJ> int GPVec<OBJ>::RemovePtr(pointer item) {
if (item==NULL) return -1;
for (int i=0;i<fCount;i++)
   if ((pointer)fList[i] == item) {
       Delete(i);
       return i;
       }
return -1; //not found
}

template <class OBJ> void GPVec<OBJ>::Pack()  {
 for (int i=fCount-1; i>=0; i--)
    if (fList[i]==NULL) Delete(i); //shift rest of fList content accordingly
}

template <class OBJ> void GPVec<OBJ>::setCount(int NewCount) {
  if (NewCount<0 || NewCount > MAXLISTSIZE)
     GError(GVEC_COUNT_ERR, NewCount);
  if (NewCount > fCapacity) setCapacity(NewCount);
  if (NewCount > fCount) //pad with NULL pointers!
    memset(& fList[fCount], 0, (NewCount - fCount) * sizeof(OBJ*));
  fCount = NewCount;
}

template <class OBJ> void GPVec<OBJ>::qSort(int L, int R, GCompareProc* cmpFunc) {
 int I, J;
 OBJ* P;
 OBJ* T;
 do {
    I = L;
    J = R;
    P = this->fList[L + ((R-L)>>1)];
    do {
      while (cmpFunc(this->fList[I], P) < 0) I++;
      while (cmpFunc(this->fList[J], P) > 0) J--;
      if (I <= J) {
        T = this->fList[I];
        this->fList[I] = this->fList[J];
        this->fList[J] = T;
        I++;
        J--;
        }
      }
    while (I <= J);
    if (L < J) qSort(L, J, cmpFunc);
    L = I;
    }
 while (I < R);
}

template <class OBJ> void GPVec<OBJ>::Sort(GCompareProc* cmpFunc) {
 if (cmpFunc==NULL) {
    GMessage("Warning: NULL compare function given, useless Sort() call.\n");
    return;
    }
 if (this->fList!=NULL && this->fCount>0)
     qSort(0, this->fCount-1, cmpFunc);
}


template <class OBJ> void GPVec<OBJ>::Sort() {
  GCompareProc* cmpFunc = DefLTCompareProc<OBJ>;
  Sort(cmpFunc);
}


//---------------------------------------------------------------------------
#endif
