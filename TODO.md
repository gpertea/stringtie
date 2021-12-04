* implement a better SPMC threading model with less locking. Multiple small bundles should be packaged together (say 50 tiny bundles at once), with a dynamic limit based on the total number of junctions across the bundles, and passed on to a worker thread all at once.

