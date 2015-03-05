// filtering.cpp
// Assignment 3


#include "filtering.h"
#include <math.h>

using namespace std;

/**************************************************************
 //            IMAGE CONVOLUTION AND FILTERING               //
 *************************************************************/


// PS03 - 2.1.1 - convolve an image with a box filter of size k by k
Image boxBlur(const Image &im, const int &k, bool clamp) {
Image output = Image(im.width(), im.height(), im.channels());

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float value = 0.0;
        int half_box_width = k/2;
        float numPixels = k*k;
          for (int xx = x-half_box_width; xx < x+half_box_width + 1; xx++) {
            for (int yy = y-half_box_width; yy < y+half_box_width + 1; yy++) {
                value += im.smartAccessor(xx, yy, z, clamp)/numPixels;
            }
        }
        output(x, y, z) = value;
      }
    }
  }

    return output;
}


// PS03 - 2.2.2 - reimeplement the box filter using the filter class.
// check that your results math those in the previous function "boxBlur"
Image boxBlur_filterClass(const Image &im, const int &k, bool clamp) {

    int filterWidth = k;
    int filterHeight = k;
    int numBoxes = filterWidth*filterHeight;
    vector<float> fData (numBoxes, 0.0);
    for(int i = 0; i < numBoxes; i++) {
       fData[i] = 1.0/numBoxes;
   }

    Filter f = Filter(fData, filterWidth, filterHeight);
    Image output = f.Convolve(im, clamp);
  
    return output;
}


// PS03 - 2.3.1 - uses a Sobel kernel to compute the horizontal and vertical
// components of the gradient of an image and returns the gradient magnitude.
Image gradientMagnitude(const Image &im, bool clamp){

    float fDataXArray[] = { -1.0, 0.0, 1.0, -2.0, 0.0, 2.0, -1.0, 0.0, 1.0 };
    vector<float> fDataX (fDataXArray, fDataXArray + sizeof(fDataXArray) / sizeof(float) );
    Filter sobelX(fDataX, 3, 3);

    float fDataYArray[] = { -1.0, -2.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
    vector<float> fDataY (fDataYArray, fDataYArray + sizeof(fDataYArray) / sizeof(float) );
    Filter sobelY(fDataY, 3, 3);

    Image sobelXOutput = sobelX.Convolve(im, clamp);
    Image sobelYOutput = sobelX.Convolve(im, clamp);

    Image output = Image(im.width(), im.height(), im.channels());

    for (int x = 0; x < output.width(); x++) {
        for (int y = 0; y < output.height(); y++) {
            for (int z = 0; z < output.channels(); z++) {
                output(x, y, z) = pow(pow(sobelXOutput(x, y), 2.0) + pow(sobelYOutput(x, y), 2.0), 0.5);
            }
        }
    }

    return output;
    
}

float compute_gaussian(float r, float sigma) {
    float e = 2.718281828459045235360;
    float pi = 3.1415926535;
    float a = 1.0/(sigma*pow(2.0*pi, 0.5));
    float c = sigma;

    float r2 = pow(r, 2.0);
    float c2_2 = 2.0*pow(c, 2.0);
    float value = a*pow(e, -1.0*r2/c2_2);
    return value;
}

// PS03 - 2.4.1 - create a vector containing the normalized values in a 1D Gaussian filter
vector<float> gauss1DFilterValues(float sigma, float truncate){

    int numBoxes = 1+2*ceil(sigma * truncate);

    vector<float> fData (numBoxes, 0.0);

    float e = 2.718281828459045235360;
    float pi = 3.1415926535;
    float a = 1.0/(sigma*pow(2.0*pi, 0.5));
    float c = sigma;
    float totalWeight = 0.0;

    for(int i = 0; i < numBoxes; i++) {
        float dx = floor(numBoxes/2) - i;
        float r = dx;
        float r2 = pow(r, 2.0);
        float c2_2 = 2.0*pow(c, 2.0);
        float value = a*pow(e, -1.0*r2/c2_2);
        fData[i] = value;
        totalWeight += value;
   }
   for(int i = 0; i < numBoxes; i++) {
        fData[i] = fData[i]/totalWeight;
    }
    
    return fData;
}

// PS03 - 2.4.2 - blur across the rows of an image
Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate, bool clamp){

    vector<float> fData = gauss1DFilterValues(sigma, truncate);
    int width = fData.size();
    int height = 1;

    Filter horizontal_gaussian(fData, width, height);

    Image output = horizontal_gaussian.Convolve(im, clamp);
  
    return output;
}

// PS03 - 2.4.3 - create a vector containing the normalized values in a 2D Gaussian filter
vector<float> gauss2DFilterValues(float sigma, float truncate){
    
    int filterWidth = 1+2*ceil(sigma * truncate);
    int filterHeight = 1+2*ceil(sigma * truncate);
    int numBoxes = filterWidth*filterHeight;

    vector<float> fData (numBoxes, 0.0);

    float e = 2.718281828459045235360;
    float pi = 3.1417;
    float a = 1.0/(sigma*pow(2.0*pi, 0.5));
    float c = sigma;
    float totalWeight = 0.0;

    for(int x = 0; x < filterWidth; x++) {
        for(int y = 0; y < filterWidth; y++) {
        float dx = floor(filterWidth/2.0) - x;
        float dy = floor(filterWidth/2.0) - y;
        float r = pow(pow(dx, 2.0) + pow(dy, 2.0), 0.5);
        float r2 = pow(r, 2.0);
        float c2_2 = 2.0*pow(c, 2.0);
        float value = a*pow(e, -1.0*r2/c2_2);
        int i = x*filterWidth+y;
        fData[i] = value;
        totalWeight += value;
    }
   }

    for(int x = 0; x < filterWidth; x++) {
        for(int y = 0; y < filterWidth; y++) {
            int i = x*filterWidth+y;
            fData[i] = fData[i]/totalWeight;
        }
    }
    
    return fData;}


// PS03 - 2.4.4 - Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
Image gaussianBlur_2D(const Image &im, float sigma, float truncate, bool clamp){
    vector<float> fData = gauss2DFilterValues(sigma, truncate);

    int width = sqrt(fData.size());
    int height = width;

    Filter gauss2DFilter(fData, width, height);

    Image output = gauss2DFilter.Convolve(im, clamp);
  
    return output;
}

// PS03 - 2.4.5 - Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
Image gaussianBlur_seperable(const Image &im, float sigma, float truncate, bool clamp){
    vector<float> fData = gauss1DFilterValues(sigma, truncate);
    int width = fData.size();
    int height = 1;

    Filter horizontal_gaussian(fData, width, height);
    Filter vertical_gaussian(fData, height, width);

    Image horizontal_output = horizontal_gaussian.Convolve(im, clamp);
    Image output = vertical_gaussian.Convolve(horizontal_output, clamp);
  
    return output;
}


// PS03 - 2.5.1 - sharpen an image
Image unsharpMask(const Image &im, float sigma, float truncate, float strength, bool clamp){
    Image blurred2D = gaussianBlur_seperable(im, sigma, truncate, clamp);
    Image high_pass = im - blurred2D;

    Image output = im + high_pass*strength;
    return output;
    
}


// PS03 - 3.0.1 -  Denoise an image using bilateral filtering
Image bilateral(const Image &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp){
Image output = Image(im.width(), im.height(), im.channels());

int filterWidth = 1+2*ceil(sigmaDomain * truncateDomain);
int filterHeight = 1+2*ceil(sigmaDomain * truncateDomain);

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
        float k = 0.0;
        float r_out = 0.0;
        float g_out = 0.0;
        float b_out = 0.0;

        float r_center = im.smartAccessor(x, y, 0);
        float g_center = im.smartAccessor(x, y, 1);
        float b_center = im.smartAccessor(x, y, 2);

        int half_box_width = filterWidth/2;
        int half_box_height = filterHeight/2;
        for (int xx = 0; xx < filterWidth; xx++) {
            for (int yy = 0; yy < filterHeight; yy++) {
                float xxx = x - xx + half_box_width;
                float yyy = y - yy + half_box_height;

                float r_in = im.smartAccessor(xxx, yyy, 0, clamp);
                float g_in = im.smartAccessor(xxx, yyy, 1, clamp);
                float b_in = im.smartAccessor(xxx, yyy, 2, clamp);

                float dx2 = pow(xxx - x, 2.0);
                float dy2 = pow(yyy - y, 2.0);

                float distance = pow(dx2 + dy2, 0.5);

                float dr2 = pow(r_center - r_in, 2.0);
                float dg2 = pow(g_center - g_in, 2.0);
                float db2 = pow(b_center - b_in, 2.0);

                float color_distance = pow(dr2 + dg2 + db2, 0.5);

                float gaussian_distance_weight = compute_gaussian(distance, sigmaDomain);
                float gaussian_color_weight = compute_gaussian(color_distance, sigmaRange);
                float gaussians_multiplied = gaussian_distance_weight*gaussian_color_weight;

                k += gaussians_multiplied;

                float r_val = gaussians_multiplied*r_in;
                float g_val = gaussians_multiplied*g_in;
                float b_val = gaussians_multiplied*b_in;

                r_out += r_val;
                g_out += g_val;
                b_out += b_val;

            }
        }
        r_out = r_out/k;
        g_out = g_out/k;
        b_out = b_out/k;

        output(x, y, 0) = r_out;
        output(x, y, 1) = g_out;
        output(x, y, 2) = b_out;

    }
  }

    return output;
}


// PS03 - 3.1.1 - 6.865 only: Bilaterial Filter an image seperatly for
// the Y and UV components of an image
Image bilaYUV(const Image &im, float sigmaRange, float sigmaY, float sigmaUV, float truncateDomain, bool clamp){

Image im2 = rgb2yuv(im);
Image output = Image(im.width(), im.height(), im.channels());

// int filterWidth = 1+2*ceil(sigmaRange * truncateDomain);
// int filterHeight = 1+2*ceil(sigmaRange * truncateDomain);

int Y_kernelRadius = ceil(sigmaY * truncateDomain);
int Y_kernelSize   = 1 + 2 * Y_kernelRadius;

int UV_kernelRadius = ceil(sigmaUV * truncateDomain);
int UV_kernelSize   = 1 + 2 * UV_kernelRadius;

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
        float ky = 0.0;
        float ku = 0.0;
        float kv = 0.0;

        float y_out = 0.0;
        float u_out = 0.0;
        float v_out = 0.0;

        float y_center = im2.smartAccessor(x, y, 0);
        float u_center = im2.smartAccessor(x, y, 1);
        float v_center = im2.smartAccessor(x, y, 2);

        int half_box_width = Y_kernelSize/2;
        int half_box_height = Y_kernelSize/2;
        for (int xx = 0; xx < Y_kernelSize; xx++) {
            for (int yy = 0; yy < Y_kernelSize; yy++) {
                float xxx = x - xx + half_box_width;
                float yyy = y - yy + half_box_height;

                float y_in = im2.smartAccessor(xxx, yyy, 0, clamp);
                float u_in = im2.smartAccessor(xxx, yyy, 1, clamp);
                float v_in = im2.smartAccessor(xxx, yyy, 2, clamp);

                float dx2 = pow(xxx - x, 2.0);
                float dy2 = pow(yyy - y, 2.0);

                float distance = pow(dx2 + dy2, 0.5);

                float dyy2 = pow(y_center - y_in, 2.0);

                float y_distance = pow(dyy2, 0.5);

                float gaussian_distance_weight = compute_gaussian(distance, sigmaRange);

                float gaussian_y_weight = compute_gaussian(y_distance, sigmaY);

                ky += gaussian_y_weight*gaussian_distance_weight;

                float y_val = gaussian_y_weight*gaussian_distance_weight*y_in;

                y_out += y_val;

            }
        }
        half_box_width = UV_kernelSize/2;
        half_box_height = UV_kernelSize/2;
        for (int xx = 0; xx < UV_kernelSize; xx++) {
            for (int yy = 0; yy < UV_kernelSize; yy++) {
                float xxx = x - xx + half_box_width;
                float yyy = y - yy + half_box_height;

                float u_in = im2.smartAccessor(xxx, yyy, 1, clamp);
                float v_in = im2.smartAccessor(xxx, yyy, 2, clamp);

                float dx2 = pow(xxx - x, 2.0);
                float dy2 = pow(yyy - y, 2.0);

                float distance = pow(dx2 + dy2, 0.5);

                float du2 = pow(u_center - u_in, 2.0);
                float dv2 = pow(v_center - v_in, 2.0);

                float u_distance = pow(du2, 0.5);
                float v_distance = pow(dv2, 0.5);
                float gaussian_u_weight = compute_gaussian(u_distance, sigmaUV);
                float gaussian_v_weight = compute_gaussian(v_distance, sigmaUV);

                float gaussian_distance_weight = compute_gaussian(distance, sigmaRange);

                ku += gaussian_u_weight*gaussian_distance_weight;
                kv += gaussian_v_weight*gaussian_distance_weight;

                float u_val = gaussian_u_weight*gaussian_distance_weight*u_in;
                float v_val = gaussian_v_weight*gaussian_distance_weight*v_in;

                u_out += u_val;
                v_out += v_val;
            }
        }

        output(x, y, 0) = y_out;
        output(x, y, 1) = u_out;
        output(x, y, 2) = v_out;

    }
  }

  output = yuv2rgb(output);

  return output;
}


/**************************************************************
 //                 FILTER CLASS FUNCTIONS                  //
 *************************************************************/


// PS03 - 2.2.1 - write a convolution function for the filter class
Image Filter::Convolve(const Image &im, bool clamp){
    
Image output = Image(im.width(), im.height(), im.channels());

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float value = 0.0;
        int half_box_width = width/2;
        int half_box_height = height/2;
          for (int xx = 0; xx < width; xx++) {
            for (int yy = 0; yy < height; yy++) {
                float kernelValue = (*this)(xx, yy);
                float xxx = x - xx + half_box_width;
                float yyy = y - yy + half_box_height;
                float pixelValue = im.smartAccessor(xxx, yyy, z, clamp);
                value += pixelValue*kernelValue;
            }
        }
        output(x, y, z) = value;
      }
    }
  }

    return output;
}


/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
Image impulseImg(const int &k){
    
    // initlize a kxkx1 image of all 0's
    Image impulse(k, k, 1);
    
    // set the center pixel to have intensity 1
    int center = floor(k/2);
    impulse(center,center,0) = 1;
    
    return impulse;
}


// Filter class constructor
Filter::Filter(const vector<float> &fData, const int &fWidth, const int &fHeight) {
    
    kernel = fData;
    width = fWidth;
    height = fHeight;
    
}

Filter::Filter(const int &fWidth, const int &fHeight) {
  width = fWidth;
  height = fHeight;
  kernel = std::vector<float>(width*height,0);
}

const float & Filter::operator()(int x, int y) const {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();
    
    return kernel[x +y*width];
    
}


float & Filter::operator()(int x, int y) {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();
    
    return kernel[x +y*width];
}


Filter::~Filter() {}
