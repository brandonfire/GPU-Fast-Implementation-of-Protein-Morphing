#ifndef INPUT_H_
#define INPUT_H_
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define PROTEINLENGTH 125 //protein has 125 AA in our sample input file 
#define FLAME 20.0 //K = 20 for the out put
// RAM block size used for reading in data
#define RAMusage 95
typedef struct calpha{
	double x_pos;
	double y_pos;
	double z_pos;
} Proteinbone;
void Readpdb(char *file_name, char *output);

#endif
