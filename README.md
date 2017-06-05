# Density-of-States
A research tool to investigate different algorithms to determine the density of states of a thermodynamic system.

## License
 * the code under gengraph folder is licensed under GPL2 (see gengraph/README.txt)
 * the remaining code is licensed under a MIT license

## Prerequisites
 * cmake >= 2.6
 * GSL (GNU Scientific Library)
 * Boost
 * open-mpi (optional)
 * doxygen (optional)

## Compile
    mkdir debug && cd debug
    cmake -DCMAKE_BUILD_TYPE=DEBUG ..
    
    mkdir release && cd release
    cmake -DCMAKE_BUILD_TYPE=RELEASE ..
