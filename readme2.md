1) Added "string" header file in Tdefs.h file

2) Replaced all "throw" cases with error_message function to output errors in
   ostream.

3) Implemented a method "compatibility_checker" 
   in the class GraphReader, which checks for the compatibility of input

4) Implemented signature for output
   in Utils.cc, which prints the number of allotments with the given signature
   Changed some minor parts in print_signature definition and added new argument for choosing stream

5) Included parameter for choosing between type of outputs
   in main.cc, which lets us choose between resident ranks type or signature type

./graphmatching -A -p -i ../resources/our_input.txt -o ../resources/outout2.txt