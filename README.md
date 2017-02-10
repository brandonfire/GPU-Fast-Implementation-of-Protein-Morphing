The initial version of GPU implement MORPH-PRO
This package including three programs:
input, intermed and distangle


To complie the program simple type make:
NOTICE: the make file is original from cuda 7.0. if your version cannot complie it by make. you may change the directory in make file or complie the program by gcc and nvcc for each program

INPUT is for processing the original pdb file the alpha Carbon atoms only file.

$ ./input {start.pdb} {pstart.pdb}
$ ./input {end.pdb} {pend.pdb}
the first parameter is the input file. the second is the file name you want to store new file.

Intermed is for processing the start and end pose to generate the intermediate pose

$ ./intermed {pstart.pdb} {pend.pdb} {interpose.pdb}

first argument is for start pose pdb file. second argument is for end pose pdb file. 
the third file is the output file name you want to store you results. Right now I set the total pose generate as 20. 
You could change it in input.h line 7

distangle is to use GPU to process the distance and angles.

$ ./distangle {p.pdb}

The argument is the pdb file you want to calculate. 

