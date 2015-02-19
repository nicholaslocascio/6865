// a1.h
// Assignment 1


#include "a1.h"
using namespace std;

// The implementation of these functions will be given to you after the pset 0 deadline
Image brightness(const Image &im, const float &factor) {
  return im * factor;
}

Image contrast(const Image &im, const float &factor, const float &midpoint) {
  return (im - midpoint)*factor+midpoint;
}



Image changeGamma(const Image & im, const double & old_gamma, const double & new_gamma) {
  // Image output(im.width(), im.height(), im.channels());
  // Figure out what power to take the values of im, to get the values of output
  // return output;

  if (new_gamma < 0 ) {
    throw OutOfBoundsException();
  }

  Image gammad_image = Image(im.width(), im.height(), im.channels());

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      for (int z = 0; z < im.channels(); z++) {
        float old_value = im(x, y, z);
        float new_value = (float)(pow(old_value, new_gamma/old_gamma));
        gammad_image(x, y, z) = new_value;
      }
    }
  }

  return gammad_image;
}

// Change the exposure of the image. This is different than brightness because
// it means you multiply an image's intensity in the linear domain.
Image exposure(const Image &im, const float &factor) {
  if (factor < 0) {
    throw OutOfBoundsException();
  }
  return changeGamma(brightness(changeGamma(im, 1.0/2.2, 1.0), factor), 1.0, 1.0/2.2);
}


Image color2gray(const Image &im, const std::vector<float> &weights) {
  Image output(im.width(), im.height(), 1);

  Image gammad_image = changeGamma(im, 1.0/2.2, 1.0);

  for (int x = 0; x < gammad_image.width(); x++) {
    for (int y = 0; y < gammad_image.height(); y++) {
      float total = 0;
      for (int z = 0; z < 3; z++) {
        float weight = weights[z];
        float old_value = gammad_image(x, y, z);
        total += weight*old_value;
      }
      output(x, y) = total;
    }
  }

  output = changeGamma(output, 1.0, 1.0/2.2);
  return output;
}


// For this function, we want two outputs, a single channel luminance image 
// and a three channel chrominance image. Return them in a vector with luminance first
std::vector<Image> lumiChromi(const Image &im) {
  Image im_chrominance = Image(im.width(), im.height(), im.channels());

  Image im_luminance = color2gray(im);

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      float luminance_value = im_luminance(x, y);

      for (int z = 0; z < 3; z++) {
        float old_value = im(x, y, z);
        float new_value = old_value/luminance_value;
        im_chrominance(x, y, z) = new_value;
      }
    }
  }

  std::vector<Image> output;
  output.push_back(im_luminance);
  output.push_back(im_chrominance);
  return output;
}

// Modify brightness then contrast
Image brightnessContrastLumi(const Image &im, const float &brightF, const float & contrastF, const float & midpoint) {
  std::vector<Image> lumi_and_chromi = lumiChromi(im);
  Image im_luminance = lumi_and_chromi[0];
  Image im_chrominance = lumi_and_chromi[1];
  Image output = Image(im.width(), im.height(), im.channels());

  Image modified_luminance = contrast(brightness(im_luminance, brightF), contrastF, midpoint);

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      float luminance_value = modified_luminance(x, y);
      for (int z = 0; z < 3; z++) {
        float chrom_value = im_chrominance(x, y, z);
        float new_value = chrom_value*luminance_value;
        new_value = fmin(fmax(new_value, 0), 1);
        output(x, y, z) = new_value;
      }
    }
  }

  return output;
}

Image rgb2yuv(const Image &im) {
  // Create output image of appropriate size
  // Change colorspae
  Image output = Image(im.width(), im.height(), im.channels());

  float a_1 = 0.299;
  float a_2 = 0.587;
  float a_3 = 0.114;

  float b_1 = -0.147;
  float b_2 = -0.289;
  float b_3 = 0.436;

  float c_1 = 0.615;
  float c_2 = -0.515;
  float c_3 = -0.100;

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      float r = im(x, y, 0);
      float g = im(x, y, 1);
      float b = im(x, y, 2);

      float y_val = a_1*r + a_2*g + a_3*b;
      float u_val = b_1*r + b_2*g + b_3*b;
      float v_val = c_1*r + c_2*g + c_3*b;

      output(x, y, 0) = y_val;
      output(x, y, 1) = u_val;
      output(x, y, 2) = v_val;
    }
  }

  return output; 
}

Image yuv2rgb(const Image &im) {
  Image output = Image(im.width(), im.height(), im.channels());

  float a_1 = 1.0;
  float a_2 = 0;
  float a_3 = 1.14;

  float b_1 = 1.0;
  float b_2 = -0.395;
  float b_3 = -0.581;

  float c_1 = 1.0;
  float c_2 = 2.032;
  float c_3 = 0.0;

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      float y_val = im(x, y, 0);
      float u_val = im(x, y, 1);
      float v_val = im(x, y, 2);

      float r = a_1*y_val + a_2*u_val + a_3*v_val;
      float g = b_1*y_val + b_2*u_val + b_3*v_val;
      float b = c_1*y_val + c_2*u_val + c_3*v_val;

      output(x, y, 0) = r;
      output(x, y, 1) = g;
      output(x, y, 2) = b;
    }
  }

  return output;
}

Image saturate(const Image &im, const float &factor) {
  Image output = Image(im.width(), im.height(), im.channels());
  Image yuv = rgb2yuv(im);

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      float y_val = yuv(x, y, 0);
      float u_val = yuv(x, y, 1);
      float v_val = yuv(x, y, 2);

      output(x, y, 0) = y_val;
      output(x, y, 1) = u_val*factor;
      output(x, y, 2) = v_val*factor;
    }
  }

  output = yuv2rgb(output);
  return output;
}

// Return two images in a C++ vector
std::vector<Image> spanish(const Image &im) {

  std::vector<Image> lumi_chromi_images = lumiChromi(im);
  Image im_luminance = lumi_chromi_images[0];
  Image im_chrominance = lumi_chromi_images[1];

  Image im_2 = im_luminance;
  Image im_1 = Image(im.width(), im.height(), im.channels());
  Image yuv = rgb2yuv(im_chrominance);

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      float u_val = yuv(x, y, 1);
      float v_val = yuv(x, y, 2);

      im_1(x, y, 0) = 0.5;
      im_1(x, y, 1) = -1.0*u_val;
      im_1(x, y, 2) = -1.0*v_val;
    }
  }

  im_1 = yuv2rgb(im_1);

  int middle_x = floor(im_1.width()/2);
  int middle_y = floor(im_1.height()/2);
  im_1(middle_x, middle_y, 0) = 0.0;
  im_1(middle_x, middle_y, 1) = 0.0;
  im_1(middle_x, middle_y, 2) = 0.0;

  std::vector<Image> output;
  output.push_back(im_1);
  output.push_back(im_2);
  return output;
}


// White balances an image using the gray world assumption
Image grayworld(const Image & im) {
  Image output = Image(im.width(), im.height(), im.channels());

  float r_sum = 0.0;
  float g_sum = 0.0;
  float b_sum = 0.0;

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      r_sum += im(x, y, 0);
      g_sum += im(x, y, 1);
      b_sum += im(x, y, 2);
    }
  }

  int num_pix = im.width()*im.height();
  float r_ave = r_sum/num_pix;
  float g_ave = g_sum/num_pix;
  float b_ave = b_sum/num_pix;

  float r_fact = g_ave/r_ave;
  float g_fact = 1.0;
  float b_fact = g_ave/b_ave;

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      output(x, y, 0) = r_fact*im(x, y, 0);
      output(x, y, 1) = g_fact*im(x, y, 1);
      output(x, y, 2) = b_fact*im(x, y, 2);
    }
  }

  return output;
}
