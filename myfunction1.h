
#ifndef __myfunc1_h
#define __myfunc1_h

#include "myutil.h"
#include "readBMP.h"

void initialize_pixel_sum(int *sum);
void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);


void rowBlurKernel_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);
void rowSharpKernel_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);
void blurKernel_withou_filter_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);
void blurKernel_with_filter_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);
void sharpKernel_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);




#endif

