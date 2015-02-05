#include "a0.h"

using namespace std;

// This is a way for you to test your functions. 
// We will only grade the contents of a0.cpp and Image.cpp
int main() {
  helloworld(-10.0, 4.2);
  Image boston = readAnImage("Input/Boston_low_contrast.png");
  // my_im.write("test-read.png");
  Image im_2 =  Image(100, 100, 3, "test size");
  cout << "size of data is " << im_2.number_of_elements() << endl;
  cout << "image data at 4 is " << boston(4) << endl;
  cout << "image data at 4,3 is " << boston(4,3) << endl;
  cout << "image data at 4,3,1 is " << boston(4,3,1) << endl;
  Image brightened_image = brightness(boston, 100.0);
  brightened_image.write("test-brightened.png");

  Image contrasted_image = contrast(boston, 7.0, 0.5);
  contrasted_image.write("contrasted_image.png");
}
