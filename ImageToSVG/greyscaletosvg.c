//
//  greyscaletosvg.c
//  ansiscratch
//
//  Created by Malcolm McLean on 13/07/2015.
//  Copyright (c) 2015 Malcolm McLean. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "binaryutils.h"
#include "bezierfit.h"
#include "greyscaletosvg.h"

static void markcontours(unsigned char *out, unsigned char *grey, int width, int height, int low, int high);
static void svgclosedpath(FILE *fp, double *bez, int Nsegs, int r, int g, int b, double opacity);
static void filtervalsprecious(double *x, int N, unsigned char *precious);

static const char *svgtag  = "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">";

int floodfill4(unsigned char *grey, int width, int height, int x, int y, unsigned char target, unsigned char dest);

/*
   Convert a greyscale image to an svg file
   Params: grey - the greyscale iamge
           width - image width
		   height - image height
		   Nlevels - number of lvels to use (default 64) 
		   fitting - bezier curve fit accuracy (default 3.0 pixels)
		   svg - output file for svg result
   Notes:
        reasonably good results on natural photographs.
*/
int greyscaletosvg(unsigned char *grey, int width, int height, int Nlevels, double fitting, FILE *svg)
{
	unsigned char *buff;
	int i, ii, iii;
	int cutoff;
	int Ncontours;
	int *cN;
	int Nsegs;
	double **cx, **cy;
	double *path;
	double opacity;
	unsigned char *precious;

	buff = malloc(width * height);
	if (!buff)
		goto out_of_memory;


	fprintf(svg, "<svg width = \"%d\" height = \"%d\">\n", width, height);
	fprintf(svg, "<rect x = \"0\" y = \"0\" width = \"%d\" height = \"%d\" stroke = \"none\" fill = \"#000000\" />",
		width, height);
	for (i = 1; i < Nlevels; i++)
	{
		cutoff = (i * 256) / (Nlevels);
		opacity = 1.0 / (Nlevels - i);
		markcontours(buff, grey, width, height, cutoff, ((i + 1) * 256) / Nlevels);
		Ncontours = getcontours(buff, width, height, &cx, &cy, &cN);

		for (ii = 0; ii < Ncontours; ii++)
		{
			if (cN[ii] > 8)
			{
				precious = malloc(cN[ii]);
				if (!precious)
					goto out_of_memory;
				for (iii = 0; iii < cN[ii]; iii++)
				{
					if (cx[ii][iii] == 0 || cx[ii][iii] == width ||
						cy[ii][iii] == 0 || cy[ii][iii] == height)
						precious[iii] = 1;
					else
						precious[iii] = 0;
				}
				filtervalsprecious(cx[ii], cN[ii], precious);
				filtervalsprecious(cy[ii], cN[ii], precious);
				free(precious);
			}
			path = bez_fitpoints(cx[ii], cy[ii], cN[ii], fitting, &Nsegs);

			svgclosedpath(svg, path, Nsegs, cutoff, cutoff, cutoff, 1.0);

			free(path);
		}

		for (ii = 0; ii < Ncontours; ii++)
		{
			free(cx[ii]);
			free(cy[ii]);
		}
		free(cx);
		free(cy);
		free(cN);


	}
	fprintf(svg, "</svg>\n");

	return 0;
out_of_memory:
	return -1;
}

static void markcontours(unsigned char *out, unsigned char *grey, int width, int height, int low, int high)
{
	int i;

	for (i = 0; i < width*height; i++)
		out[i] = (grey[i] >= low) ? 2 : 0;
	for (i = 0; i < width*height; i++)
	{
		if (out[i] == 2 && grey[i] >= low && grey[i] < high)
			floodfill4(out, width, height, i%width, i / width, 2, 1);
	}
	for (i = 0; i < width*height; i++)
	{
		if (out[i] == 2)
			out[i] = 0;
	}
}


static void svgclosedpath(FILE *fp, double *bez, int Nsegs, int r, int g, int b, double opacity)
{
    int i;
    
    if(Nsegs <= 0)
        return;
    fprintf(fp, "<path d=\"M %d %d ", (int) bez[0], (int) bez[1]);
	for (i = 0; i < Nsegs; i++)
	{
		fprintf(fp, "        C %d, %d %d, %d %d, %d \n",
			(int)bez[i * 8 + 2], (int)bez[i * 8 + 3],
			(int)bez[i * 8 + 4], (int)bez[i * 8 + 5],
			(int)bez[i * 8 + 6], (int)bez[i * 8 + 7]);
	}
    fprintf(fp, "        Z\"\n");
	fprintf(fp, "fill = \"#%02x%02x%02x\" ", r, g, b);
	fprintf(fp, "stroke = \"none\" ");
	fprintf(fp, "opacity = \"%3.3f\"/>\n", opacity);
    
}


static void filtervalsprecious(double *x, int N, unsigned char *precious)
{
    int i;
    double *res;
    int j, k;
    
    res = (double *) malloc(N * sizeof(double));
    for(i=0;i<N;i++)
    {
        j = i == 0 ? N- 1: i -1;
        k = i == N-1 ? 0 : i+1;
        res[j] = (x[j]+x[k])/2;
    }
	for (i = 0; i < N; i++)
	{
		if (!precious[i])
			x[i] = res[i];
	}
    free(res);
}


/*
Algorithm:
Flood-fill (node, target-color, replacement-color):
1. Set Q to the empty queue.
2. If the color of node is not equal to target-color, return.
3. Add node to Q.
4. For each element n of Q:
5.     If the color of n is equal to target-color:
6.         Set w and e equal to n.
7.         Move w to the west until the color of the node to the west of w no longer matches target-color.
8.         Move e to the east until the color of the node to the east of e no longer matches target-color.
9.         Set the color of nodes between w and e to replacement-color.
10.         For each node n between w and e:
11.             If the color of the node to the north of n is target-color, add that node to Q.
12.             If the color of the node to the south of n is target-color, add that node to Q.
13. Continue looping until Q is exhausted.
14. Return.
*/

/*
floodfill4 - floodfill, 4 connectivity
Params: grey - the image (formally it's greyscale but it coul dbe binary or indexed)
width - image width
heoght - image height
x, y - seed point
target - the colour to flood
dest - the colur to replace it by.
Returns: number of pixels flooded
*/
int floodfill4(unsigned char *grey, int width, int height, int x, int y, unsigned char target, unsigned char dest)
{
	int *qx = 0;
	int *qy = 0;
	int qN = 0;
	int qpos = 0;
	int qcapacity = 0;
	int wx, wy;
	int ex, ey;
	int tx, ty;
	int ix;
	int *temp;
	int answer = 0;

	if (grey[y * width + x] != target)
		return 0;
	qx = malloc(width * sizeof(int));
	qy = malloc(width * sizeof(int));
	if (qx == 0 || qy == 0)
		goto error_exit;
	qcapacity = width;
	qx[qpos] = x;
	qy[qpos] = y;
	qN = 1;

	while (qN != 0)
	{
		tx = qx[qpos];
		ty = qy[qpos];
		qpos++;
		qN--;

		if (qpos == 256)
		{
			memmove(qx, qx + 256, qN*sizeof(int));
			memmove(qy, qy + 256, qN*sizeof(int));
			qpos = 0;
		}
		if (grey[ty*width + tx] != target)
			continue;
		wx = tx;
		wy = ty;
		while (wx >= 0 && grey[wy*width + wx] == target)
			wx--;
		wx++;
		ex = tx;
		ey = ty;
		while (ex < width && grey[ey*width + ex] == target)
			ex++;
		ex--;


		for (ix = wx; ix <= ex; ix++)
		{
			grey[ty*width + ix] = dest;
			answer++;
		}

		if (ty > 0)
			for (ix = wx; ix <= ex; ix++)
			{
			if (grey[(ty - 1)*width + ix] == target)
			{
				if (qpos + qN == qcapacity)
				{
					temp = realloc(qx, (qcapacity + width) * sizeof(int));
					if (temp == 0)
						goto error_exit;
					qx = temp;
					temp = realloc(qy, (qcapacity + width) * sizeof(int));
					if (temp == 0)
						goto error_exit;
					qy = temp;
					qcapacity += width;
				}
				qx[qpos + qN] = ix;
				qy[qpos + qN] = ty - 1;
				qN++;
			}
			}
		if (ty < height - 1)
			for (ix = wx; ix <= ex; ix++)
			{
			if (grey[(ty + 1)*width + ix] == target)
			{
				if (qpos + qN == qcapacity)
				{
					temp = realloc(qx, (qcapacity + width) * sizeof(int));
					if (temp == 0)
						goto error_exit;
					qx = temp;
					temp = realloc(qy, (qcapacity + width) * sizeof(int));
					if (temp == 0)
						goto error_exit;
					qy = temp;
					qcapacity += width;
				}
				qx[qpos + qN] = ix;
				qy[qpos + qN] = ty + 1;
				qN++;
			}
			}
	}

	free(qx);
	free(qy);

	return answer;
error_exit:
	free(qx);
	free(qy);
	return -1;
}

/*
floodfill, 8 connectivity
*/
int floodfill8(unsigned char *grey, int width, int height, int x, int y, unsigned char target, unsigned char dest)
{
	int *qx = 0;
	int *qy = 0;
	int qN = 0;
	int qpos = 0;
	int qcapacity = 0;
	int wx, wy;
	int ex, ey;
	int tx, ty;
	int ix;
	int *temp;
	int answer = 0;

	if (grey[y * width + x] != target)
		return 0;
	qx = malloc(width * sizeof(int));
	qy = malloc(width * sizeof(int));
	if (qx == 0 || qy == 0)
		goto error_exit;
	qcapacity = width;
	qx[qpos] = x;
	qy[qpos] = y;
	qN = 1;

	while (qN != 0)
	{
		tx = qx[qpos];
		ty = qy[qpos];
		qpos++;
		qN--;

		if (qpos == 256)
		{
			memmove(qx, qx + 256, qN*sizeof(int));
			memmove(qy, qy + 256, qN*sizeof(int));
			qpos = 0;
		}
		if (grey[ty*width + tx] != target)
			continue;
		wx = tx;
		wy = ty;
		while (wx >= 0 && grey[wy*width + wx] == target)
			wx--;
		wx++;
		ex = tx;
		ey = ty;
		while (ex < width && grey[ey*width + ex] == target)
			ex++;
		ex--;


		for (ix = wx; ix <= ex; ix++)
		{
			grey[ty*width + ix] = dest;
			answer++;
		}

		if (ty > 0)
			for (ix = wx - 1; ix <= ex + 1; ix++)
			{
			if (ix < 0 || ix >= width)
				continue;
			if (grey[(ty - 1)*width + ix] == target)
			{
				if (qpos + qN == qcapacity)
				{
					temp = realloc(qx, (qcapacity + width) * sizeof(int));
					if (temp == 0)
						goto error_exit;
					qx = temp;
					temp = realloc(qy, (qcapacity + width) * sizeof(int));
					if (temp == 0)
						goto error_exit;
					qy = temp;
					qcapacity += width;
				}
				qx[qpos + qN] = ix;
				qy[qpos + qN] = ty - 1;
				qN++;

			}
			}
		if (ty < height - 1)
			for (ix = wx - 1; ix <= ex + 1; ix++)
			{
			if (ix < 0 || ix >= width)
				continue;

			if (grey[(ty + 1)*width + ix] == target)
			{
				if (qpos + qN == qcapacity)
				{
					temp = realloc(qx, (qcapacity + width) * sizeof(int));
					if (temp == 0)
						goto error_exit;
					qx = temp;
					temp = realloc(qy, (qcapacity + width) * sizeof(int));
					if (temp == 0)
						goto error_exit;
					qy = temp;
					qcapacity += width;
				}
				qx[qpos + qN] = ix;
				qy[qpos + qN] = ty + 1;
				qN++;

			}
			}
	}

	free(qx);
	free(qy);

	return answer;
error_exit:
	free(qx);
	free(qy);
	return -1;
}
