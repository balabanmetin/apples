* `apples/PoolQuerkWorker.py`: the main worker each thread runs
  * the function `runquery` is the entry point 
  * One worker per thread
  * Each worker gets its private copy of all the data structure including reference 
* `apples/reference.py`: the reference tree, with all the preprocessing results included
  * `build_appledb.py`: builds an apple reference without any queries; it can be pickled.  
* `FM.py`, `BE.py`, etc.: the algorithm implementations. 
  * Entry point is `dp_frag`.
  * `placement_per_edge`: computes the error for each edge
  * `placement`: finds the best edge
* `subtree.py`: computes the subtree relevant per query for a given alignment. 
  * `validate_edges` is the main function that decides what to ignore. This works based on `obs-dist` calculated before. 
  * similar function in `PoolQuerkWorker.py` for when distances are given as input. 
  * `unroll_changes`: removes changes before going to the next query
