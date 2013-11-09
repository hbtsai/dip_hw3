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


unsigned char erode_filter[9]={ 1, 1, 1,
								1, 1, 1,
	     					    1, 1, 1 };

/* 3x3 filter, make sure image is not at boundary */
int compare(int x, int y, int width, unsigned char* image, unsigned char* filter, int type)
{
	int f_width=3;
	int ret = 0;
	int xc=image[x*width+(y)];
	int x0=image[x*width+(y+1)];
	int x1=image[(x-1)*width+(y+1)];
	int x2=image[(x-1)*width+y];
	int x3=image[(x-1)*width+(y-1)];
	int x4=image[x*width+(y-1)];
	int x5=image[(x+1)*width+(y-1)];
	int x6=image[(x+1)*width+y];
	int x7=image[(x+1)*width+(y+1)];

	int mc=255*erode_filter[ 1*f_width+   (1)];
	int m0=255*erode_filter[ 1*f_width+   (1+1)];
	int m1=255*erode_filter[(1-1)*f_width+(1+1)];
	int m2=255*erode_filter[(1-1)*f_width+ 1];
	int m3=255*erode_filter[(1-1)*f_width+(1-1)];
	int m4=255*erode_filter[ 1*f_width+   (1-1)];
	int m5=255*erode_filter[(1+1)*f_width+(1-1)];
	int m6=255*erode_filter[(1+1)*f_width+ 1];
	int m7=255*erode_filter[(1+1)*f_width+(1+1)];

	if(type)
	{
		if(xc==mc && x0==m0 && 
		   x1==m1 && x2==m2 &&  
		   x3==m3 && x4==m4 && 
		   x5==m5 && x6==m6 && 
		   x7==m7)
			ret = 1;
		else
			ret = 0;

		if(!x0 && !x1 && !x2 && !x3 && !x4 && !x5 && !x6 && !x7 && xc )
			ret =1;

		if(xc && x0 && x1 && x2 && !x3 && !x4 && !x5 && !x6 && !x7)
			ret = 1;

		if(xc && x7 && !x6 && !x2 && !x0 && !x3 && !x5 && !x1 && !x4)
			ret = 1;

	}
	else
	{
		if(xc!=mc && x0!=m0 && 
		   x1!=m1 && x2!=m2 &&  
		   x3!=m3 && x4!=m4 && 
		   x5!=m5 && x6!=m6 && 
		   x7!=m7)
			ret = 0;
		else
			ret = 1;
	}

	return ret;
}


/* type=1: erosion, type=0: dilation*/
void erosion(int width, int height, unsigned char* image, unsigned char* image_r, int type)
{
	int i=0, j=0;
	int pix=0;
	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			if((i-1)<0 || (j-1)<0 || (i+1)>=width || (j+1)>=height)
				continue;

			pix = compare (i, j, width, image, erode_filter, type);

			if(pix)
				image_r[i*width+j]=255;
			else
				image_r[i*width+j]=0;

		}
	}
}

// analysis 3x3 patch
void dfs(int x, int y, int width, unsigned char* patch, int color)
{
	int pc=0,p0=0,p1=0,p2=0,p3=0,p4=0,p5=0,p6=0,p7=0;

	pc=patch[x*width+(y)];
	p0=patch[x*width+(y+1)];
	p1=patch[(x-1)*width+(y+1)];
	p2=patch[(x-1)*width+y];
	p3=patch[(x-1)*width+(y-1)];
	p4=patch[x*width+(y-1)];
	p5=patch[(x+1)*width+(y-1)];
	p6=patch[(x+1)*width+y];
	p7=patch[(x+1)*width+(y+1)];

	// set center color first
	if(pc==255)
	{
		patch[x*width+y] = color;
		if(p0==255)
			dfs(x, y+1, width, patch, color);
		if(p1==255)
			dfs(x-1, y+1, width, patch, color);
		if(p2==255)
			dfs(x-1, y, width, patch, color);
		if(p3==255)
			dfs(x-1, y-1, width, patch, color);
		if(p4==255)
			dfs(x, y+1, width, patch, color);
		if(p5==255)
			dfs(x+1, y-1, width, patch, color);
		if(p6==255)
			dfs(x+1, y, width, patch, color);
		if(p7==255)
			dfs(x+1, y+1, width, patch, color);
	}

}

void conn_label(int width, int height, unsigned char* image)
{
	int x=0, y=0;
	int gradient=30, count=0;
	char filename[256]={};
	for(x=0; x<height; x++)
	{
		for(y=0; y<width; y++)
		{
			if((x-1)<0 || (y-1)<0 || (x+1)>=(height) || (y+1) >= (width))
				continue;

			if(image[x*width+y] == 255)
			{
				dfs(x, y, width, image, gradient);
				count++;
				gradient+=25;
				sprintf(filename, "conn_label_%d.pgm", gradient);
				write_pgm_image(filename, Size, Size, image);
			}
		}
	}

	fprintf(stderr, "count=%d\n", count);
	paint_histogram( width, height, image, "conn_label_histogram.pgm");
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
	write_pgm_image("sample1.pgm", Size, Size, Imagedata);


	unsigned char tmp[Size*Size] = {};
	unsigned char ero[Size*Size] = {};
	unsigned char ero2[Size*Size] = {};
	erosion(Size, Size, Imagedata, ero, 1);
	erosion(Size, Size, ero, ero2, 1);
	int iter=30;
	while(iter>0)
	{
		memcpy(tmp, ero2, sizeof(ero2));
		memset(ero2, 0, sizeof(ero2));
		erosion(Size, Size, tmp, ero2, 1);
		iter--;
	}
	write_pgm_image("erosion.pgm", Size, Size, ero);
	write_pgm_image("erosion2.pgm", Size, Size, ero2);

	int z=0, circle=0;
	for(z=0; z<Size*Size; z++)
	{
		if(ero2[z]>0)
			circle++;
	}

	fprintf(stderr, "circle=%d\n", circle);

	unsigned char dil[Size*Size] = {};
	erosion(Size, Size, Imagedata, dil, 0);
	write_pgm_image("dilation.pgm", Size, Size, dil);

	unsigned char edge[Size*Size] = {};
	for(i=0; i<Size*Size; i++)
		edge[i] = Imagedata[i]-ero[i];
	write_pgm_image("edge.pgm", Size, Size, edge);

	unsigned char conn[Size*Size]={};
	memcpy(conn, Imagedata, sizeof(Imagedata));
	conn_label(Size, Size, conn);
	write_pgm_image("connected.pgm", Size, Size, conn);


	exit(0);
	return 0;
}




	
