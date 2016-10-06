# IBMGenerator
IBM Synthetic Data Generator for Itemsets and Sequences

Type make, which will create the executable file 'gen'

type ./gen -help for general help

For itemsets, type ./gen lit -help
An example run can be: 
./gen lit -ntrans 100 -tlen 10 -nitems 1 -npats 1000 -patlen 4 -fname T10I4D100K -ascii

This will generate a datafile named "T10I4D100K.data" which is in the IBM format
