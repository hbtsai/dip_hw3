#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define Size 512

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
	{
	  return 1;
	}
}

int paint_histogram(int width, int height, unsigned char* image, char* filename)
{

	unsigned char *histoimg = calloc(1, 256*256);
	double histo[256]={};
	int i=0, j=0;
	for(i=0; i<256; i++)
		histo[i]=0;

	for(i=0; i<width*height; i++)
		histo[image[i]]++;

	for(i=0; i<256; i++)
	{
		histo[i]/=65536;
		histo[i]*=255*40;
	}

	for(i=0; i<256; i++)
	{
		for(j=0; j<histo[i]; j++)
		{
			if(j>255)
				break;
			histoimg[i*256+j]=238;
		}

	}

	write_pgm_image(filename, 256, 256, histoimg);
	free(histoimg);
	histoimg=NULL;
	return 0;
}


int main(int argc, char** argv)
{
	FILE *file = NULL;
	// image data array
	unsigned char Imagedata[Size*Size] = {};

	char fname[1024]={};
	if(argv[1] != NULL && strlen(argv[1])>0)
		strcpy(fname, argv[1]);
	else
	{
		fprintf(stderr, "please specify filename of raw input image.\n");
		exit(-1);
	}

	if (!(file=fopen(fname,"rb")))
	{
		fprintf(stderr, "Cannot open file!\n");
		exit(1);
	}
	fread(Imagedata, sizeof(unsigned char), Size*Size, file);
	fclose(file);


	int i=0;
	for(i=0; i<Size*Size; i++)
	{
		if(Imagedata[i]>128)
			Imagedata[i]=255;
		else
			Imagedata[i]=0;
	}
	/* save the original image for comparision */
	write_pgm_image("sample3.pgm", Size, Size, Imagedata);



	exit(0);
	return 0;
}




	
