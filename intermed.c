#include "intermed.h"
//#include "input.h"

struct calpha* start;
struct calpha* end;
//struct calpha* intermed
void pdbtoarray(char *file_name,Proteinbone *atom)//function that will read pdb file and store to array
{
	//for pure proccessed pdb I only read XYZ to my array
	char Xtmp[9] = { ' ', ' ',' ', ' ', ' ', ' ', ' ', ' ' };
	char Ytmp[9] = { ' ', ' ',' ', ' ', ' ', ' ', ' ', ' ' };
	char Ztmp[9] = { ' ', ' ',' ', ' ', ' ', ' ', ' ', ' ' };
	char fileBuffer[RAMusage]; // input character array from file
	double X, Y, Z;// output floats
	int line, i, j, endOfLine;//, endOfComment, size; // logical integers
	FILE *fileID;
	fileID = fopen( file_name, "r" );
	if ( fileID != NULL ) 
	{	
		line = 0;
		while ( fgets( fileBuffer, sizeof (fileBuffer), fileID ) != NULL )
		{
			for ( i = 0; fileBuffer[i] != '\0' ; i++ ) {} endOfLine = i - 1;// find end of line
			if (endOfLine > 54  )
			{
				j = 0; for ( i = 30; i < 38; i++ ) { Xtmp[j]      = fileBuffer[i]; j++; }
				j = 0; for ( i = 38; i < 46; i++ ) { Ytmp[j]      = fileBuffer[i]; j++; }
				j = 0; for ( i = 46; i < 54; i++ ) { Ztmp[j]      = fileBuffer[i]; j++; }
				X     = atof( Xtmp );
				Y     = atof( Ytmp );
				Z     = atof( Ztmp );
				atom[line].x_pos = X;
				atom[line].y_pos = Y;
				atom[line].z_pos = Z;
				//printf("%8.3lf%8.3lf%8.3lf\n", X, Y, Z);
				line++;
			}
		}
	}else { printf("File could not be opened\n"); } // if file opening failed
	fclose( fileID );
}

int main ( int argc, char *argv[])
{
	if(argc != 4){
		printf("Error: only 3 arguments allowed, please input start phase, end phase and the final result pdb file\n");
		return 1;
	}
	char molType[7] = { 'A', 'T', 'O', 'M', ' ', ' ' };
	char atmNumTmp[6] = { ' ', ' ', ' ', ' ', ' ' };
	char atomName[5] = { ' ', 'C', 'A', ' ' };
	char altLoc[2] = { ' ' };
	char resName[4] = { ' ', ' ', ' ' };
	char resNumTmp[5] = { ' ', ' ', ' ', ' ' };
	char chainID[2] = { 'A' };
	char occupTmp[7] = { ' ', ' ', '1', '.', '0', '0' };
	char betaTmp[7] = { ' ', ' ', ' ', ' ', ' ', ' ' };
	int atomNum, resNum; // output integers
	double X,Y,Z, occup, beta; 
	atomNum      = atof( atmNumTmp );
	resNum       = atof( resNumTmp );
	occup = atof( occupTmp );
	beta  = atof( betaTmp );

	start = (Proteinbone *)malloc(sizeof(Proteinbone)*PROTEINLENGTH);
	end = (Proteinbone *)malloc(sizeof(Proteinbone)*PROTEINLENGTH);
	pdbtoarray(argv[1],start);
	pdbtoarray(argv[2],end);
	int i;	
	/*printf(" The start protein backbone:\n");
	for(i=0;i<PROTEINLENGTH;i++){
		printf("%8.3lf%8.3lf%8.3lf\n", start[i].x_pos, start[i].y_pos, start[i].z_pos);
	}
	printf(" The end protein backbone:\n");
	for(i=0;i<PROTEINLENGTH;i++){
		printf("%8.3lf%8.3lf%8.3lf\n", end[i].x_pos, end[i].y_pos, end[i].z_pos);
	}*/

	int m;
	//const double FLAME = 200.0;
	double alpha;
	FILE *pdbfiletowrite;
	pdbfiletowrite = fopen( argv[3], "w" );
	for(m=0;m<FLAME;m++){
		fprintf(pdbfiletowrite,"MODEL %d\n",m);
		alpha = (FLAME-m)/FLAME;
		for(i=0;i<PROTEINLENGTH;i++){
			X = alpha*start[i].x_pos  +  (1-alpha)*end[i].x_pos;
			Y = alpha*start[i].y_pos  +  (1-alpha)*end[i].y_pos;
			Z = alpha*start[i].z_pos  +  (1-alpha)*end[i].z_pos;
			fprintf(pdbfiletowrite,"%-5s%5d %4s%1s%3s %1s%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", 
					   molType, atomNum, atomName, altLoc, resName, chainID, 
					   resNum, X, Y, Z, occup, beta);
		}
		fprintf(pdbfiletowrite,"END\n");
		fprintf(pdbfiletowrite,"ENDMDL\n");
	}
	printf("Output the results at %s. \n", argv[3]);
	fclose( pdbfiletowrite );
	return 0;
}
