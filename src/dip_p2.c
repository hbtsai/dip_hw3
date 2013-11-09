#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define Size 64

int write_pgm_image(char* filename, int x_dim, int y_dim, unsigned char* image)
{
	unsigned char* y = image;
	FILE* filehandle = NULL;
	filehandle = fopen(filename, "wb");
	if (filehandle) 
	{
		fprintf(filehandle, "P5\n\n%d %d 255\n", x_dim, y_dim);
		fwrite(y, 1, x_dim * y_dim, filehandle);
		fclose(filehandle);
		return 0;
	} 
	else
	  return 1;
}

void write_pgm_image_double(char* filename, int x_dim, int y_dim, double* image)
{
	unsigned char image_c[Size*Size]={};
	int i=0;
	double max=0, min=0;

	for(i=0; i<x_dim*y_dim; i++)
	{
		if( image[i] > max)
			max = image[i] ;
		if(image[i]  < min)
			min = image[i] ;
	}

	fprintf(stderr, "max=%f min=%f\n", max, min);

	int dist=(int)(max-min)+5;
	int check_val = 0;

	for(i=0; i<Size*Size; i++)
	{
		check_val = 254* ((image[i]+abs(min)+1)/dist);

		if(check_val<0 || check_val>255)
			fprintf(stderr, "warning!! value %d incorrect\n", check_val);

		image_c[i]=  (unsigned char) check_val ;
	}

	write_pgm_image(filename, x_dim, y_dim, image_c);
}


void readfile(char* fname, unsigned char* image_r)
{
	FILE *file = NULL;
	if (!(file=fopen(fname,"rb")))
	{
		fprintf(stderr, "Cannot open file!\n");
		exit(1);
	}
	fread(image_r, sizeof(unsigned char), Size*Size, file);
	fclose(file);
}


// 1/36
int laws1[9]={ 1, 2, 1, 
						 2, 4, 2,
						 1, 2, 1};

// 1/12
int laws2[9]={ 1, 0, -1,
						 2, 0, -2,
						 1, 0, -1};

// 1/12
int laws3[9]={ -1, 2, -1,
						 -2, 4, -2,
						 -1, 2, -1};

// 1/12
int laws4[9]={ -1, -2, -1, 
						 0 , 0, 0, 
						  1, 2, 1};

// 1/4
int laws5[9]={ 1, 0, -1, 
						 0, 0, 0, 
						 -1, 0, 1};

// 1/4
int laws6[9]={ -1, 2, -1, 
						  0, 0, 0, 
						  1, -2, 1};

// 1/12
int laws7[9]={ -1, -2, -1, 
						  2, 4, 2, 
						 -1, -2, -1};

// 1/4
int laws8[9]={ -1, 0, 1, 
						  2, 0, -2, 
						  -1, 0, 1};

// 1/4
int laws9[9]={ 1, -2, 1, 
						 -2, 4, -2, 
						 1, -2, 1};


void convolve(int width, int height, int* kernel, unsigned char* image, double* image_r, double divider)
{
	int i=0, j=0;
	double pix0=0;
	double pix1=0;
	double pix2=0;
	double pix3=0;
	double pix4=0;
	double pix5=0;
	double pix6=0;
	double pix7=0;
	double pix8=0;
	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{

			if(i==0 || j==0 || i+1==height || j+1==width)
				continue;

			pix0=(double)image[i*width+(j-1)]/255;
			pix0*=kernel[0];
			pix1=(double)image[(i-1)*width+j]/255;
			pix1*=kernel[1];
			pix2=(double)image[(i-1)*width+(j+1)]/255;
			pix2*=kernel[2];
			pix3=(double)image[i*width+(j-1)]/255;
			pix3*=kernel[3];
			pix4=(double)image[i*width+j]/255;
			pix4*=kernel[4];
			pix5=(double)image[i*width+(j+1)]/255;
			pix5*=kernel[5];
			pix6=(double)image[(i+1)*width+(j-1)]/255;
			pix6*=kernel[6];
			pix7=(double)image[(i+1)*width+j]/255;
			pix7*=kernel[7];
			pix8=(double)image[(i+1)*width+(j+1)]/255;
			pix8*=kernel[8];
			image_r[i*width+j]=(pix0 + pix1 + pix2 + pix3 + pix4 + pix5 + pix6 + pix7 + pix8 )*255/divider;
		}
	
	}
	return ;
}

void create_features(int sample, unsigned char *image, double *sample_laws)
{
	char filename[1024]={};
	convolve(Size, Size, laws1, image, sample_laws, 36);
	convolve(Size, Size, laws2, image, sample_laws+Size*Size, 12);
	convolve(Size, Size, laws3, image, sample_laws+Size*Size*2, 12);
	convolve(Size, Size, laws4, image, sample_laws+Size*Size*3, 12);
	convolve(Size, Size, laws5, image, sample_laws+Size*Size*4, 4);
	convolve(Size, Size, laws6, image, sample_laws+Size*Size*5, 4);
	convolve(Size, Size, laws7, image, sample_laws+Size*Size*6, 12);
	convolve(Size, Size, laws8, image, sample_laws+Size*Size*7, 4);
	convolve(Size, Size, laws9, image, sample_laws+Size*Size*8, 4);
	int i=0; 
	for(i=0; i<9; i++)
	{
		sprintf(filename, "sample%02d_laws%d.pgm", sample, i+1);
		write_pgm_image_double(filename, Size, Size, &sample_laws[i*Size*Size]);
	}
	
}

void calculate_energy(double* sample, long long * energy)
{
	int idx=0; /* 0-8 */
	int i=0; /* 0-4095 */
	double val = 0;
	for(idx = 0; idx < 9; idx++)
	{
		for(i=0; i<Size*Size; i++)
			val += sample[idx*Size*Size+i]*sample[idx*Size*Size+i];

		energy[idx]=(long long)val;
	}
}

int main(int argc, char** argv)
{
	unsigned char image01[Size*Size] = {};
	unsigned char image02[Size*Size] = {};
	unsigned char image03[Size*Size] = {};
	unsigned char image04[Size*Size] = {};
	unsigned char image05[Size*Size] = {};
	unsigned char image06[Size*Size] = {};
	unsigned char image07[Size*Size] = {};
	unsigned char image08[Size*Size] = {};
	unsigned char image09[Size*Size] = {};
	unsigned char image10[Size*Size] = {};
	unsigned char image11[Size*Size] = {};
	unsigned char image12[Size*Size] = {};
	unsigned char image13[Size*Size] = {};
	unsigned char image14[Size*Size] = {};
	unsigned char image15[Size*Size] = {};
	unsigned char image16[Size*Size] = {};


	readfile("sample01.raw", image01);
	readfile("sample02.raw", image02);
	readfile("sample03.raw", image03);
	readfile("sample04.raw", image04);
	readfile("sample05.raw", image05);
	readfile("sample06.raw", image06);
	readfile("sample07.raw", image07);
	readfile("sample08.raw", image08);
	readfile("sample09.raw", image09);
	readfile("sample10.raw", image10);
	readfile("sample11.raw", image11);
	readfile("sample12.raw", image12);
	readfile("sample13.raw", image13);
	readfile("sample14.raw", image14);
	readfile("sample15.raw", image15);
	readfile("sample16.raw", image16);

	write_pgm_image("sample01.pgm", Size, Size, image01);
	write_pgm_image("sample02.pgm", Size, Size, image02);
	write_pgm_image("sample03.pgm", Size, Size, image03);
	write_pgm_image("sample04.pgm", Size, Size, image04);
	write_pgm_image("sample05.pgm", Size, Size, image05);
	write_pgm_image("sample06.pgm", Size, Size, image06);
	write_pgm_image("sample07.pgm", Size, Size, image07);
	write_pgm_image("sample08.pgm", Size, Size, image08);
	write_pgm_image("sample09.pgm", Size, Size, image09);
	write_pgm_image("sample10.pgm", Size, Size, image10);
	write_pgm_image("sample11.pgm", Size, Size, image11);
	write_pgm_image("sample12.pgm", Size, Size, image12);
	write_pgm_image("sample13.pgm", Size, Size, image13);
	write_pgm_image("sample14.pgm", Size, Size, image14);
	write_pgm_image("sample15.pgm", Size, Size, image15);
	write_pgm_image("sample16.pgm", Size, Size, image16);

	/* image01 */
	double* sample01_laws=calloc(9*Size*Size, sizeof(double));
	create_features(1, image01, sample01_laws);

	/* image02 */
	double* sample02_laws=calloc(9*Size*Size, sizeof(double));
	create_features(2, image02, sample02_laws);

	/* image03 */
	double* sample03_laws=calloc(9*Size*Size, sizeof(double));
	create_features(3, image03, sample03_laws);

	/* image04 */
	double* sample04_laws=calloc(9*Size*Size, sizeof(double));
	create_features(4, image04, sample04_laws);

	/* image05 */
	double* sample05_laws=calloc(9*Size*Size, sizeof(double));
	create_features(5, image05, sample05_laws);

	/* image06 */
	double* sample06_laws=calloc(9*Size*Size, sizeof(double));
	create_features(6, image06, sample06_laws);

	/* image07 */
	double* sample07_laws=calloc(9*Size*Size, sizeof(double));
	create_features(7, image07, sample07_laws);

	/* image08 */
	double* sample08_laws=calloc(9*Size*Size, sizeof(double));
	create_features(8, image08, sample08_laws);

	/* image09 */
	double* sample09_laws=calloc(9*Size*Size, sizeof(double));
	create_features(9, image09, sample09_laws);

	/* image10 */
	double* sample10_laws=calloc(9*Size*Size, sizeof(double));
	create_features(10, image10, sample10_laws);

	/* image11 */
	double* sample11_laws=calloc(9*Size*Size, sizeof(double));
	create_features(11, image11, sample11_laws);
	
	/* image12 */
	double* sample12_laws=calloc(9*Size*Size, sizeof(double));
	create_features(12, image12, sample12_laws);

	/* image13 */
	double* sample13_laws=calloc(9*Size*Size, sizeof(double));
	create_features(13, image13, sample13_laws);

	/* image14 */
	double* sample14_laws=calloc(9*Size*Size, sizeof(double));
	create_features(14, image14, sample14_laws);
	
	/* image15 */
	double* sample15_laws=calloc(9*Size*Size, sizeof(double));
	create_features(15, image15, sample15_laws);

	/* image16 */
	double* sample16_laws=calloc(9*Size*Size, sizeof(double));
	create_features(16, image16, sample16_laws);

	long long energy01[9]={};
	long long energy02[9]={};
	long long energy03[9]={};
	long long energy04[9]={};
	long long energy05[9]={};
	long long energy06[9]={};
	long long energy07[9]={};
	long long energy08[9]={};
	long long energy09[9]={};
	long long energy10[9]={};
	long long energy11[9]={};
	long long energy12[9]={};
	long long energy13[9]={};
	long long energy14[9]={};
	long long energy15[9]={};
	long long energy16[9]={};

	calculate_energy(sample01_laws, energy01);
	calculate_energy(sample02_laws, energy02);
	calculate_energy(sample03_laws, energy03);
	calculate_energy(sample04_laws, energy04);
	calculate_energy(sample05_laws, energy05);
	calculate_energy(sample06_laws, energy06);
	calculate_energy(sample07_laws, energy07);
	calculate_energy(sample08_laws, energy08);
	calculate_energy(sample09_laws, energy09);
	calculate_energy(sample10_laws, energy10);
	calculate_energy(sample11_laws, energy11);
	calculate_energy(sample12_laws, energy12);
	calculate_energy(sample13_laws, energy13);
	calculate_energy(sample14_laws, energy14);
	calculate_energy(sample15_laws, energy15);
	calculate_energy(sample16_laws, energy16);
	
	int law=0;
	for(law=0; law<9; law++)
	{
		fprintf(stderr, "following is law's feature %d\n\n", law+1);
		fprintf(stderr, "%lld\n", energy01[law]);
		fprintf(stderr, "%lld\n", energy02[law]);
		fprintf(stderr, "%lld\n", energy03[law]);
		fprintf(stderr, "%lld\n", energy04[law]);
		fprintf(stderr, "%lld\n", energy05[law]);
		fprintf(stderr, "%lld\n", energy06[law]);
		fprintf(stderr, "%lld\n", energy07[law]);
		fprintf(stderr, "%lld\n", energy08[law]);
		fprintf(stderr, "%lld\n", energy09[law]);
		fprintf(stderr, "%lld\n", energy10[law]);
		fprintf(stderr, "%lld\n", energy11[law]);
		fprintf(stderr, "%lld\n", energy12[law]);
		fprintf(stderr, "%lld\n", energy13[law]);
		fprintf(stderr, "%lld\n", energy14[law]);
		fprintf(stderr, "%lld\n", energy15[law]);
		fprintf(stderr, "%lld\n", energy16[law]);
		fprintf(stderr, "\n");
	}

	free(sample01_laws);
	free(sample02_laws);
	free(sample03_laws);
	free(sample04_laws);
	free(sample05_laws);
	free(sample06_laws);
	free(sample07_laws);
	free(sample08_laws);
	free(sample09_laws);
	free(sample10_laws);
	free(sample11_laws);
	free(sample12_laws);
	free(sample13_laws);
	free(sample14_laws);
	free(sample15_laws);
	free(sample16_laws);

	sample01_laws=NULL;
	sample02_laws=NULL;
	sample03_laws=NULL;
	sample04_laws=NULL;
	sample05_laws=NULL;
	sample06_laws=NULL;
	sample07_laws=NULL;
	sample08_laws=NULL;
	sample09_laws=NULL;
	sample10_laws=NULL;
	sample11_laws=NULL;
	sample12_laws=NULL;
	sample13_laws=NULL;
	sample14_laws=NULL;
	sample15_laws=NULL;
	sample16_laws=NULL;
	exit(0);
	return 0;
}


