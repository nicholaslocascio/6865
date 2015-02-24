// basicImageManipulation.cpp
// Assignment 2

#include "basicImageManipulation.h"
using namespace std;

// PS02 - 3.1.1: create a new image that is k times bigger than the input by using nearest neighbor interpolation
Image scaleNN(const Image &im, const float &factor){
  cout << "scle NN" << "\n";
  int width = floor(factor*im.width());
  int height = floor(factor*im.height());

  Image output = Image(width, height, im.channels());

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        int xx = floor(x/factor);
        int yy = floor(y/factor);
        float value = im(xx, yy, z);
        output(x, y, z) = value;
      }
    }
  }

  return output;
}

// PS02 - 3.2.1: using bilinear interpolation to assign the value of a location from its neighboring pixel values
float interpolateLin(const Image &im, float x, float y, int z, bool clamp){
  int x_low = floor(x);
  int x_high = x_low + 1;

  float weight_x_low = 1.0 - (x - x_low);
  float weight_x_high = 1.0 - weight_x_low;

  int y_low = floor(y);
  int y_high = y_low + 1;
  float weight_y_low = 1.0 - (y - y_low);
  float weight_y_high = 1.0 - weight_y_low;

  float p00 = im.smartAccessor(x_low, y_low, z, clamp);
  float p10 = im.smartAccessor(x_high, y_low, z, clamp);
  float p01 = im.smartAccessor(x_low, y_high, z, clamp);
  float p11 = im.smartAccessor(x_high, y_high, z, clamp);

  float interp_x_with_y_low = weight_x_low*p00 + weight_x_high*p10;
  float interp_x_with_y_high = weight_x_low*p01 + weight_x_high*p11;

  float interp_final = weight_y_low*interp_x_with_y_low + weight_y_high*interp_x_with_y_high;
  return interp_final;
}

// PS02 - 3.2.2: create a new image that is k times bigger than the input by using bilinear interpolation
Image scaleLin(const Image &im, const float &factor){
  int width = factor*im.width();
  int height = factor*im.height();

  Image output = Image(width, height, im.channels());

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float xx = x/factor;
        float yy = y/factor;
        float value = interpolateLin(im, xx, yy, z, true);
        output(x, y, z) = value;
      }
    }
  }

  return output;
}


// PS02 - 3.3.1 (6.865 required, 6.815 extra credit): rotate an image around its center by theta
Image rotate(const Image &im, const float &theta) {
    // center around which to rotate. Use this center coordinate.
  float centerX = (im.width()-1.0)/2.0;
  float centerY = (im.height()-1.0)/2.0;

  Image output = Image(im.width(), im.height(), im.channels());

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float vector_x = (x - centerX);
        float vector_y = (y - centerY);
        float vector_theta = atan(vector_y/vector_x);

        float xx = vector_x*cos(-theta) - vector_y*sin(-theta) + centerX;
        float yy = vector_x*sin(-theta) + vector_y*cos(-theta) + centerY;

        float value = interpolateLin(im, xx, yy, z, false);
        output(x, y, z) = value;
      }
    }
  }
  return output;
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
