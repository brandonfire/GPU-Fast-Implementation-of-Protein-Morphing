
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
//#include <cuda.h>
#include <cuda_runtime.h>
#include "input.h"
//#include "intermed.h"
#define DISTANCERSIZE 15625
typedef struct dist_recorder{ //structure to record distance of each position we can use index to indicate poisiton
	//int pos1;//position 1
	//int pos2;//position 2
	double d;//this is the distance
} Distancer;

typedef struct ang_recorder{ //structure to record angler of each position
	//int pos;//position use index as pos
	double d;//this is the angle
} Angler;

Proteinbone * protein;
Distancer * distrecorder;
Angler * anglecorder;

/*device function to caculate distance*/
__device__
double d_p2p_distance(double x1, double x2, double y1, double y2, double z1, double z2) {
	return sqrt((x1 - x2)*(x1-x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

/* chengbin: kenrel function to caculate histogram*/


__global__
void D_dist_ang_function(double * x,double *y, double * z, Distancer * distrecorder, Angler * anglerecorder, int proteinlen){
	/*shared tiling input*/
	__shared__ double ix1[34];
	__shared__ double iy1[34];
	__shared__ double iz1[34];
	__shared__ double ix2[34];
	__shared__ double iy2[34];
	__shared__ double iz2[34];

	/*shared ouput
	extern __shared__ unsigned long long p_hist[];*/
	int i, j, ti,k;
	//input[threadIdx.x]=atomlist[];
	int gd =gridDim.x;
	int bd = blockDim.x;
	int bdx = blockIdx.x;
	ti = threadIdx.x;
	i = bdx * bd + ti;
	int didx;
	
	//for(j=ti;j<n_buckets;j+=bd)p_hist[j]=0;//iniatilize the ouput histogram
	//copy the anchor tile data to ix1,iy1,iz1 according to i
	if(i<proteinlen){	
	ix1[ti] = x[i];
	iy1[ti] = y[i];
	iz1[ti] = z[i];	
	}
	if(i=proteinlen-1){
		if(bdx!=0){
			ix1[32] = x[(bdx+1) * bd];
			iy1[32] = y[(bdx+1) * bd];
			iz1[32] = z[(bdx+1) * bd];
			ix1[33] = x[(bdx-1) * bd+31];
			iy1[33] = y[(bdx-1) * bd+31];
			iz1[33] = z[(bdx-1) * bd+31];		
		}else{
			ix1[32] = x[(bdx+1) * bd];
			iy1[32] = y[(bdx+1) * bd];
			iz1[32] = z[(bdx+1) * bd];
		}	
	}
	//double COS_70;
	//double COS_120;
	//COS_70    = cos(70. * 3.15159 / 180.);
    	//COS_120   = -0.5; 
	//angle part
	int i_1;
	if(ti == 0){
	i_1 = 33;} else {i_1=ti-1;}
	double aa = d_p2p_distance(ix1[ti],ix1[i_1],iy1[ti],iy1[i_1],iz1[ti],iz1[i_1]);
   	double bb = d_p2p_distance(ix1[ti],ix1[ti+1],iy1[ti],iy1[ti+1],iz1[ti],iz1[ti+1]);
   	double cc = d_p2p_distance(ix1[i_1],ix1[ti+1],iy1[i_1],iy1[ti+1],iz1[i_1],iz1[ti+1]);

   	double dd = 2*sqrt(aa*bb);

   	double cos_v1v2v3 = (aa + bb - cc)/dd;
	anglerecorder[i].d = cos_v1v2v3;
	

	//ix2[ti] = x[i];
	//iy2[ti] = y[i];
	//iz2[ti] = z[i];
	__syncthreads();
	//distance part
	double dist;
	
	int lastblock = gd-1;
	int lastblocklength = proteinlen - bd*(gd-1);
	if(bdx<lastblock)
	{
		
		for(j=ti+1; j<bd;j++)
		{
			dist = d_p2p_distance(ix1[ti],ix1[j],iy1[ti],iy1[j],iz1[ti],iz1[j]);
			//d_pos = (int) (dist / PDH_w);
			//atomicAdd(&(p_hist[d_pos]),1);
			didx = i*125+(bdx * bd + j);
			distrecorder[didx].d = dist;
		}
		

		
	__syncthreads();
	} else 
	{
		
		if(i<proteinlen)
		{

			
			
			for(j=ti+1; j<lastblocklength;j++)
			{
				dist = d_p2p_distance(ix1[ti],ix1[j],iy1[ti],iy1[j],iz1[ti],iz1[j]);
				//d_pos = (int) (dist / PDH_w);
				//atomicAdd(&(p_hist[d_pos]),1);
				didx = i*125+(bdx * bd + j);
				distrecorder[didx].d = dist;
			}
		}
		__syncthreads();
	}
	__syncthreads();
	
	//calcute the points between blocks.
	int cycle = ceil(gd/2.0);//becareful the last block	
	for(k=1;k<cycle;k++)//caculate points between blocks
		{
			j = (bdx+k)%gd;
			if(j<lastblock) // j is not the last block
			{
				ix2[ti] = x[j* bd + ti];
				iy2[ti] = y[j* bd + ti];
				iz2[ti] = z[j* bd + ti];
				__syncthreads();
				if(i<proteinlen)
				{
					for(int m = 0; m<bd; m++)
					{
						dist = d_p2p_distance(ix1[ti],ix2[m],iy1[ti],iy2[m],iz1[ti],iz2[m]);
						//d_pos = (int) (dist / PDH_w);
						//atomicAdd(&(p_hist[d_pos]),1);
						didx = i*125+ j* bd + m;
						distrecorder[didx].d = dist;
					}
				}
				__syncthreads();
			} else //J is the last block
			{
				
				if(ti<lastblocklength)
				{
					ix2[ti] = x[j* bd + ti];
					iy2[ti] = y[j* bd + ti];
					iz2[ti] = z[j* bd + ti];
				}
				__syncthreads();
				if(i<proteinlen)
				{
					for(int m = 0; m<lastblocklength; m++)
					{
						dist = d_p2p_distance(ix1[ti],ix2[m],iy1[ti],iy2[m],iz1[ti],iz2[m]);
						//d_pos = (int) (dist / PDH_w);
						//atomicAdd(&(p_hist[d_pos]),1);
						didx = i*125+ j* bd + m;
						distrecorder[didx].d = dist;
					}
				}
				__syncthreads();
			}
	
		}//last half cycle for gridDim.x%2==0
		if(gd%2==0)
		{
			
			if(bdx<gd/2)
			{
				j = (bdx+cycle)%gd;
				if(j<lastblock) // j is not the last block
				{
					ix2[ti] = x[j* bd + ti];
					iy2[ti] = y[j* bd + ti];
					iz2[ti] = z[j* bd + ti];
					__syncthreads();
					if(i<proteinlen)
					{
						for(int m = 0; m<bd; m++)
						{
							dist = d_p2p_distance(ix1[ti],ix2[m],iy1[ti],iy2[m],iz1[ti],iz2[m]);
							//d_pos = (int) (dist / PDH_w);
							//atomicAdd(&(p_hist[d_pos]),1);
							didx = i*125+ j* bd + m;
							distrecorder[didx].d = dist;
						}
					}
					__syncthreads();
				} else //J is the last block
				{
					
					if(ti<lastblocklength)
					{
						ix2[ti] = x[j* bd + ti];
						iy2[ti] = y[j* bd + ti];
						iz2[ti] = z[j* bd + ti];
						
					}
					__syncthreads();
					if(i<proteinlen)
					{
						for(int m = 0; m<lastblocklength; m++)
						{
							dist = d_p2p_distance(ix1[ti],ix2[m],iy1[ti],iy2[m],iz1[ti],iz2[m]);
							//d_pos = (int) (dist / PDH_w);
							//atomicAdd(&(p_hist[d_pos]),1);
							didx = i*125+ j* bd + m;
							distrecorder[didx].d = dist;
						}
					
					}
					__syncthreads();
				}



			}

		}

	__syncthreads();
	
	


}

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

/*
__global__
void D_initialize(bucket * h, int n_buckets){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i<n_buckets) h[i].d_cnt=0;

}*/
int main(int argc, char **argv)
{
	if(argc != 2) {
        printf("ERROR please input one argument: %s {the input pdb file}\n",argv[0]);
        exit(1);
    	}
    	protein = (Proteinbone *)malloc(sizeof(Proteinbone)*PROTEINLENGTH);
    	pdbtoarray(argv[1],protein);
    	double * h_x, * h_y, * h_z;//seperate host input array
    	double * d_x,* d_y,* d_z;//seperate device input array	
    	h_x = (double *)malloc(sizeof(double)*PROTEINLENGTH);
	h_y = (double *)malloc(sizeof(double)*PROTEINLENGTH);
	h_z = (double *)malloc(sizeof(double)*PROTEINLENGTH);
	cudaMalloc((void**)&d_x, sizeof(double)*PROTEINLENGTH);
	cudaMalloc((void**)&d_y, sizeof(double)*PROTEINLENGTH);	
	cudaMalloc((void**)&d_z, sizeof(double)*PROTEINLENGTH);
	int i;
	for(i = 0;  i < PROTEINLENGTH; i++) {
		h_x[i] = protein[i].x_pos;
		h_y[i] = protein[i].y_pos;
		h_z[i] = protein[i].z_pos;
	}
	distrecorder = (Distancer *)malloc(sizeof(Distancer)*DISTANCERSIZE);
	anglecorder = (Angler *)malloc(sizeof(Angler)*PROTEINLENGTH);
	Distancer * d_distrecorder;
	Angler * d_anglerecorder;
	cudaMalloc((void**)&d_distrecorder, sizeof(Distancer)*DISTANCERSIZE);
	cudaMalloc((void**)&d_anglerecorder, sizeof(Angler)*PROTEINLENGTH);	
	cudaMemcpy(d_x,h_x,sizeof(double)*PROTEINLENGTH, cudaMemcpyHostToDevice);
	cudaMemcpy(d_y,h_y,sizeof(double)*PROTEINLENGTH, cudaMemcpyHostToDevice);
	cudaMemcpy(d_z,h_z,sizeof(double)*PROTEINLENGTH, cudaMemcpyHostToDevice);
	/*chengbin: defince grid and block parameter*/
	dim3 dimGrid((int)ceil(PROTEINLENGTH/(float)32),1,1);
	dim3 dimBlock(32,1,1);
	D_dist_ang_function<<<dimGrid,dimBlock>>>(d_x, d_y, d_z, d_distrecorder,d_anglerecorder, PROTEINLENGTH);
	cudaMemcpy(distrecorder,d_distrecorder,sizeof(Distancer)*DISTANCERSIZE, cudaMemcpyDeviceToHost);
	cudaMemcpy(anglecorder,d_anglerecorder,sizeof(Angler)*PROTEINLENGTH, cudaMemcpyDeviceToHost);
	for(i = 0; i < DISTANCERSIZE; i ++ ){

	printf("%f\n",distrecorder[i].d);
	}
	for(i = 0; i < PROTEINLENGTH; i ++ ){

	printf("angle at %d: %f\n",i,anglecorder[i].d);
	}


	cudaFree(d_distrecorder);
	cudaFree(d_anglerecorder);
	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_z);
	free(protein);
	free(anglecorder);
	free(h_x);
	free(h_y);
	free(h_z);	
	free(distrecorder);
	return 0;

}

