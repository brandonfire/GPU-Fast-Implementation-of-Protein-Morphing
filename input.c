#include "input.h"
void Readpdb(char *file_name, char *output)
{
//printf ("Opening PDB file \"%s\".\n", file_name);
	
	// array must be (num. characters + 1) length for end-of-line byte!
	// strangely enough, declaring character arrays without defining them creates errors
	
	//char molType[7], atmNumTmp[6],  atomName[5], altLoc[2], resName[4];
	//char resNumTmp[5], Xtmp[9], Ytmp[9], Ztmp[9], chainID[2], occupTmp[7];
	//char betaTmp[7], comments[15];
	
	// declare and fill all output character arrays to eliminate errors
	char molType[7] = { ' ', ' ', ' ', ' ', ' ', ' ' };
	char atmNumTmp[6] = { ' ', ' ', ' ', ' ', ' ' };
	char atomName[5] = { ' ', ' ', ' ', ' ' };
	char altLoc[2] = { ' ' };
	char resName[4] = { ' ', ' ', ' ' };
	char resNumTmp[5] = { ' ', ' ', ' ', ' ' };
	char Xtmp[9] = { ' ', ' ',' ', ' ', ' ', ' ', ' ', ' ' };
	char Ytmp[9] = { ' ', ' ',' ', ' ', ' ', ' ', ' ', ' ' };
	char Ztmp[9] = { ' ', ' ',' ', ' ', ' ', ' ', ' ', ' ' };
	char chainID[2] = { ' ' };
	char occupTmp[7] = { ' ', ' ', ' ', ' ', ' ', ' ' };
	char betaTmp[7] = { ' ', ' ', ' ', ' ', ' ', ' ' };
	//I do not read comments
	//char comments[15] = { ' ', ' ',' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ' };
	
	// logical and input variables
	char fileBuffer[RAMusage]; // input character array from file
	int atomNum, resNum; // output integers
	double X, Y, Z, occup, beta; // output floats
	int i, j, endOfLine, endOfComment, size; // logical integers
	
	// open file, check file valid, check line of text exists
	FILE *fileID;
	FILE *pdbfiletowrite;
	
	fileID = fopen( file_name, "r" );
	pdbfiletowrite = fopen( output, "w" );
	if ( fileID != NULL ) {
		
		while ( fgets( fileBuffer, sizeof (fileBuffer), fileID ) != NULL ) {
			
			// check if text contains valid data
			j = 0; for ( i =  0; i <  6; i++ ) { molType[j]   = fileBuffer[i]; j++; }
			for ( i = 0; fileBuffer[i] != '\0' ; i++ ) {} endOfLine = i - 1;// find end of line
			
			if ( memcmp( molType, "ATOM" , 4 ) == 0 && endOfLine > 66  ) {
				
				// split line into various data types
				j = 0; for ( i =  6; i < 11; i++ ) { atmNumTmp[j] = fileBuffer[i]; j++; }
				j = 0; for ( i = 12; i < 16; i++ ) { atomName[j]  = fileBuffer[i]; j++; }
				//printf("%s\n",atomName);
				if(memcmp(atomName," CA ",4)!=0){continue;}
				j = 0;       i = 16;                 altLoc[j]    = fileBuffer[i];
				j = 0; for ( i = 17; i < 20; i++ ) { resName[j]   = fileBuffer[i]; j++; }
				j = 0;       i = 21;                 chainID[j]   = fileBuffer[i];
				j = 0; for ( i = 22; i < 26; i++ ) { resNumTmp[j] = fileBuffer[i]; j++; }
				j = 0; for ( i = 30; i < 38; i++ ) { Xtmp[j]      = fileBuffer[i]; j++; }
				j = 0; for ( i = 38; i < 46; i++ ) { Ytmp[j]      = fileBuffer[i]; j++; }
				j = 0; for ( i = 46; i < 54; i++ ) { Ztmp[j]      = fileBuffer[i]; j++; }
				j = 0; for ( i = 54; i < 60; i++ ) { occupTmp[j]  = fileBuffer[i]; j++; }
				j = 0; for ( i = 60; i < 66; i++ ) { betaTmp[j]   = fileBuffer[i]; j++; }
				
				// deal with unknown length of comment section
				endOfComment = endOfLine; 
				//if ( endOfLine > ( sizeof( comments ) + 65 ) ) 
				//	{ endOfComment = ( sizeof( comments ) + 65 ); }
				//j = 0; for ( i = 66; i < endOfComment; i++ ) { comments[j]  = fileBuffer[i]; j++; }
				
				// convert text to numbers where appropriate
				atomNum      = atof( atmNumTmp );
				resNum       = atof( resNumTmp );
				X     = atof( Xtmp );
				Y     = atof( Ytmp );
				Z     = atof( Ztmp );
				occup = atof( occupTmp );
				beta  = atof( betaTmp );
				
				// print out data of each line displaying interpreted format
				fprintf(pdbfiletowrite,"%-5s%5d %4s%1s%3s %1s%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", 
					   molType, atomNum, atomName, altLoc, resName, chainID, 
					   resNum, X, Y, Z, occup, beta);//, comments);%s
				
			} //else { printf( "%s", fileBuffer ); } // don't interpret non-recognized lines
		

		}
	printf("File processing complete. please check %s for the output file. \n", output);
	} else { printf("File could not be opened\n"); } // if file opening failed
	
	fclose( fileID );
	fclose( pdbfiletowrite );
	

}




int main ( int argc, char *argv[] ) {
	if (argc != 3){
	printf("ERROR: please include only two arguments: the input pdb filename and output pdbfilename.\n");
	return 1;
	}
	Readpdb(argv[1], argv[2]);
	return 0;
}
