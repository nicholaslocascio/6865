#include "a0.h"

using namespace std;

// Print the sum of two numbers
void helloworld(float a, float b) {
  // Create a variable c that is the sum of a and b
  // Use cout to display the values
  
}

// Create an image and return it
// You can test your function by calling my_im.write(filename) to write the output somewhere
Image readAnImage(const std::string &filename) {
  // Use the constructor in Image.h
  return Image(1,1,1); // Change this
}



// Change the brightness of the image
// const Image & means a reference to im wil get passed to you, 
// but the compiler won't let you modify it
Image brightness(const Image &im, const float &factor) {
  // Image output(im.width(), im.height(), im.channels());
  // Modify image brightness
  // return output;
  return Image(1,1,1); //Change this
}

Image contrast(const Image &im, const float &factor, const float &midpoint) {
  // Image output(im.width(), im.height(), im.channels());
  // Modify image contrast
  // return output;
  return Image(1,1,1); //Change this
}
