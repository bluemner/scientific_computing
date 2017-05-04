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
#include <string>
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
void getB(unsigned char * image,int i,int j,int filter_size, betacore::Matrix<T> &mB){
	int side =  (int) filter_size / (int) 2;
	unsigned int size = filter_size * filter_size;
	int x = -side;
	int y = -side;
	for(unsigned int s=0; s< size; s++){
		//check array bounds 
		

		int _j = j+y;
		int _i = i+x;

		bool __T__Y = _j >-1 && _i>-1 && _j < ydim && _i < xdim;

		unsigned char value = ( __T__Y  )? image[ _j*xdim +_i] :0; //image[j*xdim+i];
		if( __T__Y){
			int jjjj =0;
		}
		mB.set(s,0, value);
		if( s != 0 && s % filter_size ==0 ){
			++x;
		}
		y = y+1 > side? -side: y+1;
	}
}
void getFilter(int filter_size, betacore::Matrix<T> &mA){
	unsigned int size = filter_size * filter_size;
	betacore::Matrix<T> temp(size, 6);
	int side =  (int) filter_size / (int) 2;
	int x = -side;
	int y = -side;
	for( unsigned int i=0; i< size; i++){		
		// + +  + dx + ey +f, 
		temp.set( i , 0 , x*x); //ax2
		temp.set( i , 1 , y*y); //by2 
		temp.set( i , 2 , x*y); //cxy
		temp.set( i , 3 , x); //dx 
		temp.set( i , 4 , y); //ey
		temp.set( i , 5 ,1 ); //f
		
		if( i != 0 && i % filter_size ==0 ){
			++x;
		}
		y = y+1 > side? -side: y+1;
	}
	mA=temp;
}

int main(int argc, char **argv)
{
	int i, j;
	FILE *fp;

	if (argc != 4){
		printf("Usage: MyProgram <input_ppm> <output_ppm> \n");
		printf("       <input_ppm>: PGM file \n");
		printf("       <output_ppm>: PGM file \n");
		printf("       filter size\n");
		exit(0);
	}
 
	/* begin reading PGM.... */
	printf("begin reading PGM.... \n");
	if ((fp=fopen(argv[1], "rb"))==NULL){
		printf("read error...\n");
		exit(0);
	}
	ReadPGM(fp);
	int filter_size = std::stoi(argv[3]);
	if(filter_size<3){
		std::cout <<"Don't be a dick, filter size > 2 "<<std::endl;
		return (-1);
	}
	//  T A[9][6]=
	//  {// a, b, c, d, e, f
	// 	 { 0, 0, 0, 0, 0, 1}, // p(0,0) => p(i,j)
	// 	 { 0, 1, 0, 0, 1, 1}, // p(0,1) => p(i, j+i)
	// 	 { 0, 1, 0, 0,-1, 1}, // p(0,-1)=> p(i, j-i)
	// 	 { 1, 0, 0, 1, 0, 1}, // p(1,0) => p(i+1, j)
	// 	 { 1, 1, 1, 1, 1, 1}, // p(1,1) => p(i+1, j+i)
	// 	 { 1, 1,-1,-1, 1, 1}, // p(-1,1) => p(i-1, j+i)
	// 	 { 1, 1,-1, 1,-1, 1}, // p(1,-1) => p(i+1, j-1)
	// 	 { 1, 0, 0,-1, 0, 1}, // p(-1,0) => p(i-1, j)
	// 	 { 1, 1, 1,-1,-1, 1}  // p(-1,-1) => p(i-1, j-1)
	//  };
	T x[6]={1,1,1,1,1,1};
	//T b[9]={0,0,0,0,0,0};


	betacore::Matrix<T> mA(1);
	getFilter(filter_size,mA);
	mA.print();
	betacore::Matrix<T> mAT(1);
	mAT = mA;
	betacore::Matrix<T> mAP(1,1);
	mAT.transpose();
	mAT.print();
	mAP = mAT * mA;
	mAP.print();
	//
	// your application here
	for (j=0; j<ydim; j++){
		// x: 0, 1, ..., i, ..., xdim;
		for (i=0; i < xdim; i++) {

		//3x3 
		// 
		//   -j   0  +j
		//-i 0. -1. -2.
		//    |   |   | 
		//  0 3. -4. -5.
		//    |   |   |
		//+i 6. -7. -8.
		//      image [i][j] = image[j*xdim+i]
	
		
		betacore::Matrix<T> mB(filter_size*filter_size,1);
		getB(image,i,j,filter_size,mB);
		betacore::Matrix<T> mX(filter_size*filter_size,1);
		//mB.print();
		mB=  mAT * mB;

		// At * A  x = At * b
		betacore::Linear_System<T> ls(mAP, mX, mB);
		//ls.print();
		ls.solveCrout();
		//ls.forward_substitution();
		//ls.print();
		ls.get_x(mX);
		//ls.print();
		//mX.print();
		unsigned int __mid__ = filter_size * filter_size / 2;
		unsigned char __x___ = (unsigned char) mX.at( 5 ,0);
		//std::cout<< "before\t" << image[j*xdim+i] << "\tafter\t" << __x___ << std::endl;
		image[j*xdim+i] = __x___;
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
