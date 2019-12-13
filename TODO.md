* BundleData::__refseq__ and BundleData::__gseq__ should be bundled together as a shared pointer to a reference seq data structure (e.g. similar to `GRefData`) , _reference counted_ , as in a multi-threaded environment there could be many/multiple bundles sharing that info and it is wasteful to duplicate it for every bundle

* we should get rid of that terrible hashing hack (appending a dummy suffix, `id+=`) in `rlink.cpp` processRead(), by using a better hash function 

* implement a better SPMC threading model with less locking -- perhaps one with a queue per consumer/worker? However, multiple small bundles should be packaged together (say 100 tiny bundles at once) and passed on to a worker thread all at once.

