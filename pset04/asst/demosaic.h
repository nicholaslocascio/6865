// demosaic.h
// Assignment 4

#ifndef __demosaic__h
#define __demosaic__h

using namespace std;

#include "Image.h"
#include "basicImageManipulation.h"
#include <iostream>
#include <math.h>


Image basicGreen(const Image &raw, int offset=1);
Image basicRorB(const Image &raw, int offsetX, int offsetY);
Image basicDemosaic(const Image &raw, int offsetGreen=1, int offsetRedX=1, int offsetRedY=1, int offsetBlueX=0, int offsetBlueY=0);
Image edgeBasedGreen(const Image &raw, int offset=1);
Image edgeBasedGreenDemosaic(const Image &raw, int offsetGreen=1, int offsetRedX=1, int offsetRedY=1, int offsetBlueX=0, int offsetBlueY=0);
Image greenBasedRorB(const Image &raw, Image &green, int offsetX, int offsetY);
Image improvedDemosaic(const Image &raw, int offsetGreen=1, int offsetRedX=1, int offsetRedY=1, int offsetBlueX=0, int offsetBlueY=0);

 
#endif
