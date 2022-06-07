# GraphMatching
A graph matching library which provides stable and popular matching algorithms for the HR problem.


## Dependencies
clang compiler, cmake, ninja (if desired)


## Installation
    $ mkdir build; cd build
    $ cmake -G "Ninja" ../
    $ ninja

You can replace Ninja with "Unix Makefiles" in the above command.
Then type make to build.

This should build executables named *graphmatching* and *tester* inside the build directory.


## Usage
The executable takes a set of parameters to compute the desired matching:

	-s -- compute a stable matching
	-p -- compute a maximum cardinality popular matching

To run the code in test mode:

	-t -- test mode

To provide an input graph, and the output matching and/or signature filename:

	-i -- /path/to/graphfile
	-o -- /path/to/store/the/matching
	-g -- /path/to/store/signature

Also, for the -s and -p parameters, you could specify the resident/hospital
proposing algorithm (by default it runs the resident proposing algorithm).

	-A -- run the resident proposing algorithm
	-B -- run the hospital proposing algorithm

### Computing
For e.g., to compute a stable matching with the hospitals proposing (assuming inside the build directory):

	$ ./graphmatching -B -s -i ../resources/hrlq_m6.txt -o ../resources/hrlq_m6_stable.txt

### Testing
For e.g., to test the popular matching algorithm with the residents proposing (assuming inside the build directory):

	$ ./tester -A -p -o ../resources/tester_output.txt

## License
MIT
