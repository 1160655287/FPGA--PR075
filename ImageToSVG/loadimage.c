#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bmp.h"
#include "jpeg.h"
#include "gif.h"
#include "lodepng.h"

#define FMT_UNKNOWN 0
#define FMT_JPEG 1
#define FMT_PNG 2
#define FMT_BMP 3
#define FMT_GIF 4

static int getformat(char *fname);
static char *getextension(char *fname);
static void makelower(char *str);
static char *mystrdup(const char *str);

unsigned char *loadrgba(char *fname, int *width, int *height, int *err);

unsigned char *loadasjpeg(char *fname, int *width, int *height);
unsigned char *loadasbmp(char *fname, int *width, int *height);
unsigned char *loadaspng(char *fname, int *width, int *height);
unsigned char *loadasgif(char *fname, int *width, int *height);

unsigned char *loadrgba(char *fname, int *width, int *height, int *err)
{
  int fmt;
  unsigned char *answer = 0;

  if(err)
    *err = 0;

  fmt = getformat(fname);

  switch(fmt)
  {
  case FMT_UNKNOWN:
    if(err)
      *err = -3;
    return 0;
  case FMT_JPEG:
    answer = loadasjpeg(fname, width, height);
    break;
  case FMT_PNG:
    answer =  loadaspng(fname, width, height);
    break;
  case FMT_BMP:
    answer =  loadasbmp(fname, width, height);
    break;
  case FMT_GIF:
    answer =  loadasgif(fname, width, height);
    break;
  }
  if(!answer)
    if(err)
      *err = -1;

  return answer;
}

unsigned char *loadasjpeg(char *fname, int *width, int *height)
{
  unsigned char *rgb;
  int w, h;
  unsigned char *answer;
  int i;

  rgb = loadjpeg(fname, &w, &h);
  if(!rgb)
    return 0;
  answer = malloc(w * h * 4);
  if(!answer)
  {
    free(rgb);
    return 0;
  }
  for(i=0;i<w*h;i++)
  {
    answer[i*4] = rgb[i*3];
    answer[i*4+1] = rgb[i*3+1];
    answer[i*4+2] = rgb[i*3+2];
    answer[i*4+3] = 0xFF;
  }
  free(rgb);
  *width = w;
  *height = h;
  return answer;
}



unsigned char *loadasbmp(char *fname, int *width, int *height)
{
  unsigned char *rgb;
  int w, h;
  unsigned char *answer;
  int i;

  rgb = loadbmp(fname, &w, &h);
  if(!rgb)
    return 0;
  answer = malloc(w * h * 4);
  if(!answer)
  {
    free(rgb);
    return 0;
  }
  for(i=0;i<w*h;i++)
  {
    answer[i*4] = rgb[i*3];
    answer[i*4+1] = rgb[i*3+1];
    answer[i*4+2] = rgb[i*3+2];
    answer[i*4+3] = 0xFF;
  }
  free(rgb);
  *width = w;
  *height = h;
  return answer;
}

unsigned char *loadaspng(char *fname, int *width, int *height)
{
  int err;
  unsigned char *argb;
  int w, h;
  unsigned char red,green, blue, alpha;
  int i;

  argb = 0;
  err = lodepng_decode32_file(&argb, &w, &h, fname); 
  if(err)
  {
    free(argb);
    return 0;
  }
  /*
  for(i=0;i<w*h;i++)
  {
    red = argb[i*4+1];
    green = argb[i*4+2];
    blue = argb[i*4+3];
    alpha = argb[i*4];
    argb[i*4] = red;
    argb[i*4+1] = green;
    argb[i*4+2] = blue;
    argb[i*4+3] = alpha;
  }
  */
  *width = w;
  *height = h;

  return argb;
}

unsigned char *loadasgif(char *fname, int *width, int *height)
{
  unsigned char *index;
  unsigned char pal[256*3];
  int w, h;
  int i;
  int transparent;
  unsigned char *answer;

  index = loadgif(fname, &w, &h, pal, &transparent);
  if(!index)
    return 0;
  answer = malloc(w * h * 4);
  if(!answer)
  {
    free(index);
    return 0;
  }

  for(i=0;i<w*h;i++)
  {
    answer[i*4] = pal[ index[i] * 3];
    answer[i*4+1] = pal[ index[i] * 3 + 1];
    answer[i*4+2] = pal[ index[i] * 3 + 2];
    if(index[i] == transparent)
      answer[i*4+3] = 0;
    else
      answer[i*4+3] = 0xFF;
  }

  free(index);
  *width = w;
  *height = h;
  return answer;
}

static int getformat(char *fname)
{
  char *ext;
  int answer = FMT_UNKNOWN;

  ext = getextension(fname);
  if(!ext)
    return FMT_UNKNOWN;
  makelower(ext);
  if(!strcmp(ext, "jpeg") || !strcmp(ext, "jpg"))
    answer = FMT_JPEG;
  if(!strcmp(ext, "png"))
    answer = FMT_PNG;
  if(!strcmp(ext, "bmp"))
    answer = FMT_BMP;
  if(!strcmp(ext, "gif"))
    answer = FMT_GIF;
  
  free(ext);
  return answer;
}

static char *getextension(char *fname)
{
  char *answer;
  
  answer = strrchr(fname, '.');
  if(!answer)
    return 0;
  if(strchr(answer, '/'))
    return 0;
  if(strchr(answer, '\\'))
    return 0;
  
  if(answer == fname)
    return 0;

  return mystrdup(answer+1);  
}

static void makelower(char *str)
{
  int i;

  for(i=0;str[i];i++)
    str[i] = tolower(str[i]);
}
  
static char *mystrdup(const char *str)
{
  char *answer;
 
  answer = malloc(strlen(str) + 1);
  if(answer)
    strcpy(answer, str);

  return answer;
}
