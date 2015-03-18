// demosaic.cpp
// Assignment 4

#include "demosaic.h"
#include <math.h>

using namespace std;


// PS04 - 3.1.1 - takes as input a raw image and returns a single-channel
// 2D image corresponding to the green channel using simple interpolation
Image basicGreen(const Image &raw, int offset){
  Image output = Image(raw.width(), raw.height(), 1);

  for (int y = 0; y < output.height(); y++) {
    for (int x = 0; x < output.width(); x++) {
      bool onEdge = (x == 0 || y == 0 || x == output.width()-1 || y == output.height()-1);
      if (onEdge) {
        output(x, y) = raw(x,y);
        continue;
      }
      bool hasValue = (x-offset + y)%2 == 0;
      if (hasValue) {
        output(x, y) = raw(x, y);
      } else {
        output(x, y) = (raw(x-1, y) + raw(x, y-1) + raw(x+1, y) + raw(x, y+1))/4.0;
      }
    }
  }
  return output;
}

// PS04 - 3.2.1 - takes as input a raw image and returns a single-channel
// 2D image corresponding to the red or blue channel using simple interpolation
Image basicRorB(const Image &raw, int offsetX, int offsetY){
  Image output = Image(raw.width(), raw.height(), 1);

  for (int y = 0; y < output.height(); y++) {
    for (int x = 0; x < output.width(); x++) {
      bool onEdge = (x == 0 || y == 0 || x == output.width()-1 || y == output.height()-1);
      if (onEdge) {
        output(x, y) = raw(x,y);
        continue;
      }
      bool hasValue = ((x-offsetX)%2 + (y-offsetY)%2) == 0;
      if (hasValue) {
        output(x, y) = raw(x, y);
      } else {
        if (((y+offsetY)%2 == 0) && ((x+offsetX)%2 == 1)) {
          output(x, y) = (raw(x-1, y) + raw(x+1, y))/2.0;
        } else if (((x+offsetX)%2 == 0) && ((y+offsetY)%2 == 1)) {
            output(x, y) = (raw(x, y-1) + raw(x, y+1))/2.0;
        } else {
          output(x, y) = (raw(x-1, y-1) + raw(x+1, y-1) + raw(x-1, y+1) + raw(x+1, y+1))/4.0;
        }
      }
    }
  }
  return output;
}

// PS04 - 3.2.2 - takes as input a raw image and returns an rgb image
// using simple interpolation to demosaic each of the channels
Image basicDemosaic(const Image &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY){
  Image green = basicGreen(raw, offsetGreen);
  Image red = basicRorB(raw, offsetRedX, offsetRedY);
  Image blue = basicRorB(raw, offsetBlueX, offsetBlueY);

  Image output = Image(blue.width(), blue.height(), 3);

  for (int y = 0; y < output.height(); y++) {
    for (int x = 0; x < output.width(); x++) {
      output(x, y, 0) = red(x, y);
      output(x, y, 1) = green(x, y);
      output(x, y, 2) = blue(x, y);
    }
  }

  return output;
}

// PS04 - 4.1 - takes a raw image and outputs a single-channel
// image corresponding to the green channel taking into account edges
Image edgeBasedGreen(const Image &raw, int offset){
  Image output = Image(raw.width(), raw.height(), 1);

  for (int y = 0; y < output.height(); y++) {
    for (int x = 0; x < output.width(); x++) {
      bool onEdge = (x == 0 || y == 0 || x == output.width()-1 || y == output.height()-1);
      if (onEdge) {
        output(x, y) = raw(x,y);
        continue;
      }
      bool hasValue = (x-offset + y)%2 == 0;
      if (hasValue) {
        output(x, y) = raw(x, y);
      } else {
        float varianceY = fabs(raw(x, y-1) - raw(x, y+1));
        float varianceX = fabs(raw(x-1, y) - raw(x+1, y));

        bool interpolateUp = varianceY < varianceX;
        if (interpolateUp) {
          output(x, y) = (raw(x, y-1) + raw(x, y+1))/2.0;
        } else {
          output(x, y) = (raw(x-1, y) + raw(x+1, y))/2.0;
        }
      }
    }
  }
  cout << "happy" << endl;
  return output;
}

// PS04 - 4.2 - takes as input a raw image and returns an rgb image
// using edge-based green demosaicing for the green channel and
// simple interpolation to demosaic the red and blue channels
Image edgeBasedGreenDemosaic(const Image &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY){
  Image green = edgeBasedGreen(raw, offsetGreen);
  Image red = basicRorB(raw, offsetRedX, offsetRedY);
  Image blue = basicRorB(raw, offsetBlueX, offsetBlueY);

  Image output = Image(blue.width(), blue.height(), 3);

  for (int y = 0; y < output.height(); y++) {
    for (int x = 0; x < output.width(); x++) {
      output(x, y, 0) = red(x, y);
      output(x, y, 1) = green(x, y);
      output(x, y, 2) = blue(x, y);
    }
  }

  return output;
}


// PS04 - 5.1 - takes as input a raw image and returns a single-channel
// 2D image corresponding to the red or blue channel using green based interpolation
Image greenBasedRorB(const Image &raw, Image &green, int offsetX, int offsetY){
  Image output = Image(raw.width(), raw.height(), 1);

  for (int y = 0; y < output.height(); y++) {
    for (int x = 0; x < output.width(); x++) {
      bool onEdge = (x == 0 || y == 0 || x == output.width()-1 || y == output.height()-1);
      if (onEdge) {
        output(x, y) = raw(x,y);
        continue;
      }
      bool hasValue = ((x-offsetX)%2 + (y-offsetY)%2) == 0;
      if (hasValue) {
        output(x, y) = raw(x, y);
      } else {
        if (((y+offsetY)%2 == 0) && ((x+offsetX)%2 == 1)) {
          output(x, y) = green(x, y) + (raw(x-1, y)-green(x-1, y) + raw(x+1, y)-green(x+1, y))/2.0;
        } else if (((x+offsetX)%2 == 0) && ((y+offsetY)%2 == 1)) {
            output(x, y) = green(x, y) + (raw(x, y-1)-green(x, y-1) + raw(x, y+1)-green(x, y+1))/2.0;
        } else {
          output(x, y) = green(x, y) + (raw(x-1, y-1)-green(x-1, y-1) + raw(x+1, y-1)-green(x+1, y-1) + raw(x-1, y+1)-green(x-1, y+1) + raw(x+1, y+1)-green(x+1, y+1))/4.0;
        }
      }
    }
  }
  return output;
}

// PS04 - 5.2 - takes as input a raw image and returns an rgb image
// using edge-based green demosaicing for the green channel and
// simple green based demosaicing of the red and blue channels
Image improvedDemosaic(const Image &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY){

  Image green = edgeBasedGreen(raw, offsetGreen);
  Image red = greenBasedRorB(raw, green, offsetRedX, offsetRedY);
  Image blue = greenBasedRorB(raw, green, offsetBlueX, offsetBlueY);
  Image output = Image(blue.width(), blue.height(), 3);

  for (int y = 0; y < output.height(); y++) {
    for (int x = 0; x < output.width(); x++) {
      output(x, y, 0) = red(x, y);
      output(x, y, 1) = green(x, y);
      output(x, y, 2) = blue(x, y);
    }
  }

  return output;
}
