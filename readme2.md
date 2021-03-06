# All the implemented changes are listed here.

1. Added "string" header file in Tdefs.h file

2. Replaced all "throw" cases with error_message function to output errors in
   ostream.

3. Implemented a method "compatibility_checker" 
   in the class GraphReader, which checks for the compatibility of input

4. Implemented signature for output
   in Utils.cc, which prints the number of allotments with the given signature
   - Changed some minor parts in print_signature definition and added new argument for choosing stream

5. Created test suite

6. Included parameter -g for choosing between type of outputs
   in main.cc, which lets us choose between resident ranks type or signature type

7. Implemented dual signature for one-to-one popularity matching.

8. Added -t parameter for testing if matching is popular 
   - Runs testcases in the suite on p. Logs on the output file whether the testcases ran successfully, and gives dual signature for popular matchings.

9. Added new executable for tester
   - Run using `./tester -A -s -o <output_file>` for testing stable matching and `./tester -A -p -o <output_file>` for testing popular matching.

10. Removed extra files other than what's needed for compilation and comment out code other than popular and stable matching.

11. While printing signature it prints features; cardinality of A, the cardinality of B, max size of preference list in A size and also prints out only max-pref-list-size of A ranks.

12. Modified signature so that we are able to print signature in one file, matching in another file in the same command

13. Prints the usage if u run just ./graphmatching or give incorrect arguments.
