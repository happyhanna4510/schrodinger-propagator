Update Rk4Evolver buffers in include/src to add preallocated sum_ scratch.
Refactor step() in evolve_rk4.cpp to reuse buffers, noalias operations, and shared profiling logic.
Run existing profiling command to capture updated parity metrics and record them in RESULTS.txt.
