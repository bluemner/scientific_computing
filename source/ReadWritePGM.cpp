/********************************************************
***IMPORTANT NOTE***
If you have problems with the provided sample code,
part of the reason might be due to the function "fopen".
Please try changing "r/w" to "rb/wb" or the other way
when you use this function.

Note: this file is not part of the betacore project(s) work and is not bound by the same licence

LICENSE: University of Wisconsin - Milwaukee Proprietary
*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>  
#include <memory.h>
#include "../headers/linear_system.hpp"
#include "../headers/matrix.hpp"
#define max(x, y) ((x>y) ? (x):(y))
#define min(x, y) ((x<y) ? (x):(y))

typedef double T;

int xdim;
int ydim;
int maxraw;
unsigned char *image;

void ReadPGM(FILE*);
void WritePGM(FILE*);


int main(int argc, char **argv)
{
	int i, j;
	FILE *fp;

	if (argc != 3){
		printf("Usage: MyProgram <input_ppm> <output_ppm> \n");
		printf("       <input_ppm>: PGM file \n");
		printf("       <output_ppm>: PGM file \n");
		exit(0);
	}
 
	/* begin reading PGM.... */
	printf("begin reading PGM.... \n");
	if ((fp=fopen(argv[1], "rb"))==NULL){
		printf("read error...\n");
		exit(0);
	}
	ReadPGM(fp);
	 T A[9][6]=
	 {// a, b, c, d, e, f
		 { 0, 0, 0, 0, 0, 1}, // p(0,0) => p(i,j)
		 { 0, 1, 0, 0, 1, 1}, // p(0,1) => p(i, j+i)
		 { 0, 1, 0, 0,-1, 1}, // p(0,-1)=> p(i, j-i)
		 { 1, 0, 0, 1, 0, 1}, // p(1,0) => p(i+1, j)
		 { 1, 1, 1, 1, 1, 1}, // p(1,1) => p(i+1, j+i)
		 { 1, 1,-1,-1, 1, 1}, // p(-1,1) => p(i-1, j+i)
		 { 1, 1,-1, 1,-1, 1}, // p(1,-1) => p(i+1, j-1)
		 { 1, 0, 0,-1, 0, 1}, // p(-1,0) => p(i-1, j)
		 { 1, 1, 1,-1,-1, 1}  // p(-1,-1) => p(i-1, j-1)
	 };
	T x[6]={1,1,1,1,1,1};
	T b[9]={0,0,0,0,0,0};


	betacore::Matrix<T> mA(A);
	betacore::Matrix<T> mAT(A);
	betacore::Matrix<T> mAP(1,1);
	mAT.transpose();
	mAT.print();
	mAP = mAT * mA;
 T lower[6][6]=	{{1,0,0,0,0,0},
{0.66666666666,1,0,0,0,0},
{0,0,1,0,0,0},
{0,0,0,1,0,0},
{0,0,0,0,1,0},
{1,0.6,0,0,0,1}};
//betacore::Matrix<T> mAP(lower);
	mAP.print();
	//
	// your application here
	for (j=0; j<ydim; j++){
		// x: 0, 1, ..., i, ..., xdim;
		for (i=0; i < xdim; i++) {

		 /*
	

	 
		 */ 
		//3x3 
		// 
		//   -j   0  +j
		//-i 0. -1. -2.
		//    |   |   | 
		//  0 3. -4. -5.
		//    |   |   |
		//+i 6. -7. -8.
		//      image [i][j] = image[j*xdim+i]
		int _i = (i-1<0)? i: i-1 ;
		int i_ = (i+1>=xdim)?i: i+1;

		int _j = (j-1<0)? j: j-1 ;
		int j_ = (j+1>=xdim)?j: j+1;
		//0 p(0,0) => p(i,j)
		//1 p(0,1) => p(i, j+i)
		//2 p(0,-1)=> p(i, j-i)
		//3 p(1,0) => p(i+1, j)
		//4 p(1,1) => p(i+1, j+i)
		//5 p(-1,1) => p(i-1, j+i)
		//6 p(1,-1) => p(i+1, j-1)
		//7 p(-1,0) => p(i-1, j)
		//8 p(-1,-1) => p(i-1, j-1)

		b[0]= image[j*xdim+i]; // image[_j*xdim+_i];
		b[1]= image[j_*xdim+i] ;//  image[j*xdim+_i] ;
		b[2]= image[_j*xdim+i];// image[j_*xdim+_i] ;

		b[3]=image[j*xdim+i_];//  image[_j*xdim+i] ;
		b[4]=image[j_*xdim+i_];// 
		b[5]=image[j_*xdim+_i];//  

		b[6]=image[_j*xdim+i_];//  image[_j*xdim+i_] ;
		b[7]=image[j*xdim+_i];//  image[j*xdim+i_] ;
		b[8]=image[_j*xdim+_i];//  image[j_*xdim+i_] ;
		
		betacore::Matrix<T> mB(b);
		betacore::Matrix<T> mX(x);
		
		mB=  mAT * mB;	
		// At * A  x = At * b
		betacore::Linear_System<T> ls(mAP, mX, mB);
		ls.forward_substitution();
		ls.get_x(mX);	
		image[j*xdim+i] = (char) mX.at(0,0);
		}
	// printf("\n\n\n");
	}

	// As an example, let's just make an inversion of the input image.
	// for (j=0; j<ydim; j++){
	//   for (i=0; i<xdim; i++) {
	//     image[j*xdim+i] = 255 - image[j*xdim+i];
	//   }
	// }
	
	/* Begin writing PGM.... */
	printf("Begin writing PGM.... \n");
	if ((fp=fopen(argv[2], "wb")) == NULL){
		 printf("write pgm error....\n");
		 exit(0);
	 }
	WritePGM(fp);

	free(image);

	return (1);
}



void ReadPGM(FILE* fp)
{
		int c;
		int i,j;
		int val;
		unsigned char *line;
		char buf[1024];

		while ((c=fgetc(fp)) == '#')
				fgets(buf, 1024, fp);
		 ungetc(c, fp);
		 if (fscanf(fp, "P%d\n", &c) != 1) {
			 printf ("read error ....");
			 exit(0);
		 }
		 if (c != 5 && c != 2) {
			 printf ("read error ....");
			 exit(0);
		 }

		 if (c==5) {
			 while ((c=fgetc(fp)) == '#')
				 fgets(buf, 1024, fp);
			 ungetc(c, fp);
			 if (fscanf(fp, "%d%d%d",&xdim, &ydim, &maxraw) != 3) {
				 printf("failed to read width/height/max\n");
				 exit(0);
			 }
			 printf("Width=%d, Height=%d \nMaximum=%d\n",xdim,ydim,maxraw);

			 image = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
			 getc(fp);

			 line = (unsigned char *)malloc(sizeof(unsigned char)*xdim);
			 for (j=0; j<ydim; j++) {
					fread(line, 1, xdim, fp);
					for (i=0; i<xdim; i++) {
						image[j*xdim+i] = line[i];
				 }
			 }
			 free(line);

		 }

		 else if (c==2) {
			 while ((c=fgetc(fp)) == '#')
				 fgets(buf, 1024, fp);
			 ungetc(c, fp);
			 if (fscanf(fp, "%d%d%d", &xdim, &ydim, &maxraw) != 3) {
				 printf("failed to read width/height/max\n");
				 exit(0);
			 }
			 printf("Width=%d, Height=%d \nMaximum=%d,\n",xdim,ydim,maxraw);

			 image = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
			 getc(fp);

			 for (j=0; j<ydim; j++)
				 for (i=0; i<xdim; i++) {
						fscanf(fp, "%d", &val);
						image[j*xdim+i] = val;
				 }

		 }

		 fclose(fp);
}


void WritePGM(FILE* fp)
{
	int i,j;
	

	fprintf(fp, "P5\n%d %d\n%d\n", xdim, ydim, 255);
	for (j=0; j<ydim; j++)
		for (i=0; i<xdim; i++) {
			fputc(image[j*xdim+i], fp);
		}

	fclose(fp);
	
}
