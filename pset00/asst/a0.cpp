#include "a0.h"
#include "math.h"

using namespace std;

// Print the sum of two numbers
void helloworld(float a, float b) {
  float c = a + b;
  cout << "The value of a is " << a << "." << endl;
  cout << "The value of b is " << b << "." << endl;
  cout << "The sum of a and b is " << c << "." << endl;
}

// Create an image and return it
// You can test your function by calling my_im.write(filename) to write the output somewhere
Image readAnImage(const std::string &filename) {
  Image my_im = Image(filename);
  return my_im;
}


// Change the brightness of the image
// const Image & means a reference to im wil get passed to you, 
// but the compiler won't let you modify it
Image brightness(const Image &im, const float &factor) {
  // Image output(im.width(), im.height(), im.channels());
  // Modify image brightness
  // return output;
  if (factor < 0) {
    throw OutOfBoundsException();
  }

  Image brightened_image = Image(im.width(), im.height(), im.channels(), "brightened");

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      for (int z = 0; z < im.channels(); z++) {
        float old_value = im(x, y, z);
        float new_value = old_value*factor;
        new_value = fmin(fmax(new_value, 0), 1);
        brightened_image(x, y, z) = new_value;
      }
    }
  }

  return brightened_image;
}

Image contrast(const Image &im, const float &factor, const float &midpoint) {
  // Image output(im.width(), im.height(), im.channels());
  // Modify image contrast
  // return output;
  if (factor < 0) {
    throw OutOfBoundsException();
  }
  if (midpoint < 0 || midpoint > 1) {
    throw OutOfBoundsException();
  }

  Image contrasted_image = Image(im.width(), im.height(), im.channels(), "contrasted");

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      for (int z = 0; z < im.channels(); z++) {
        float old_value = im(x, y, z);
        float new_value = factor*(old_value-midpoint)+midpoint;
        new_value = fmin(fmax(new_value, 0), 1);
        contrasted_image(x, y, z) = new_value;
      }
    }
  }

  return contrasted_image;
}
