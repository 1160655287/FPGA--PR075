#include <stdio.h>
#include <stdlib.h>

#include "options.h"
#include "loadimage.h"
#include "greyscaletosvg.h"

unsigned char *rgbatogrey(unsigned char *rgba, int width, int height);

void usage()
{
	fprintf(stderr, "Raster to SVG converter (greyscale)\n");
	fprintf(stderr, "Usage imagetosvg [options] <filename>\n");
	fprintf(stderr, "<filename> - jpeg, bmp or gif file (converted to greyscale)\n");
	fprintf(stderr, "options:\n");
	fprintf(stderr, "   -Nlevels <N> - number of levels (default 64)\n");
	fprintf(stderr, "   -fit <pixels> - contour fit accuracy (default 3.0)\n");
}
int main(int argc, char **argv)
{
	int width = 100;
	int height = 62;
	unsigned char *rgba;
	unsigned char *grey;
	int err;
	FILE *fp;
	OPTIONS *opts;
	char *filename;
	int Nlevels = 64;
	double fit = 3.0;

	//goto test;
	opts = options(argc, argv, "");
	
	opt_get(opts, "-Nlevels", "%d", &Nlevels);
	opt_get(opts, "-fit", "%f", &fit);
	if (opt_Nargs(opts) != 1)
	{
		usage();
		exit(EXIT_FAILURE);
	}
	filename = opt_arg(opts, 0);
	if (opt_error(opts, stdout))
		exit(EXIT_FAILURE);
	killoptions(opts);

	if (Nlevels < 2 || Nlevels > 256)
	{
		fprintf(stderr, "Nlevels between 2 and 256\n");
		exit(EXIT_FAILURE);
	}
	if (fit < 0.1)
	{
		fprintf(stderr, "fit accuracy > 0.1 pixel\n");
		exit(EXIT_FAILURE);
	}
//test:
	//rgba = loadrgba("C:\\Users\\Malcolm\\Pictures\\lena.bmp", &width, &height, &err);
	rgba = loadrgba(filename, &width, &height, &err);
	if (!rgba)
	{
		fprintf(stderr, "Can't open input\n");
		exit(EXIT_FAILURE);
	}
	grey = rgbatogrey(rgba, width, height);
	if (!grey)
	{
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}
	fp = stdout;
	//fp = fopen("temp2.svg", "w");
	//if (!fp)
		//fprintf(stderr, "Can't open output file\n");


	free(rgba);
	rgba = 0;
	err = greyscaletosvg(grey, width, height, Nlevels, fit, fp);
	if (err < 0)
	{
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}
	free(grey);
	fclose(fp);
	//getchar();
	return 0;
}

unsigned char *rgbatogrey(unsigned char *rgba, int width, int height)
{
	unsigned char *answer;
	int ii;
	double Y;

	answer = malloc(width * height);
	if (!answer)
		return 0;
	for (ii = 0; ii < width*height; ii++)
	{
		Y = 0.2126 *rgba[ii * 4] + 0.7152 *rgba[ii * 4 + 1] + 0.0722 *rgba[ii * 4 + 2];
		answer[ii] = (unsigned char)Y;
	}

	return answer;
}