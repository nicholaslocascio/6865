// basicImageManipulation.cpp
// Assignment 2

#include "basicImageManipulation.h"
using namespace std;

// PS02 - 3.1.1: create a new image that is k times bigger than the input by using nearest neighbor interpolation
Image scaleNN(const Image &im, const float &factor){
    return Image(0); // change this
}

// PS02 - 3.2.1: using bilinear interpolation to assign the value of a location from its neighboring pixel values
float interpolateLin(const Image &im, float x, float y, int z, bool clamp){
    return 0; // change this
}

// PS02 - 3.2.2: create a new image that is k times bigger than the input by using bilinear interpolation
Image scaleLin(const Image &im, const float &factor){
    return Image(0); // change this
}


// PS02 - 3.3.1 (6.865 required, 6.815 extra credit): rotate an image around its center by theta
Image rotate(const Image &im, const float &theta) {
    // center around which to rotate. Use this center coordinate.
    // float centerX = (im.width()-1.0)/2.0;
    // float centerY = (im.height()-1.0)/2.0;


    return Image(0); // change this
}


/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 //          Functions from previous problem sets           //
 *************************************************************/

Image brightness(const Image &im, const float &factor) {
  return im * factor;
}

Image contrast(const Image &im, const float &factor, const float &midpoint) {
  return (im - midpoint)*factor+midpoint;
}
