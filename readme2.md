1) Added "string" header file in Tdefs.h file

2) Replaced all "throw" cases with error_message function to output errors in
   ostream.

3) Implemented a method "compatibility_checker" 
   in the class GraphReader, which checks for the compatibility of input

4) Implemented signature for output
   in Utils.cc, which prints the number of allotments with the given signature
   Changed some minor parts in print_signature definition and added new argument for choosing stream

5) Created test suite, for real this time

6) Included parameter for choosing between type of outputs
   in main.cc, which lets us choose between resident ranks type or signature type

7) Implemented dual signature for one-to-one popularity matching

./graphmatching -A -p -i ../resources/one2one.txt -o ../resources/outone2one1.txt

8) Added -t parameter for testing if matching is popular

   Runs testcases in the suite on p. Logs on the output file whether the testcases ran successfully, and gives dual signature for popular matchings.