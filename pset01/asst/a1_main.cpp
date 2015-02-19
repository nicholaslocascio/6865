#include "a1.h"

using namespace std;

// This is a way for you to test your functions. 
// We will only grade the contents of a1.cpp
int main() {
  Image original = Image("Input/flower.png");
  Image castle = Image("Input/castle_small.png");

  Image im_2 =  Image(100, 100, 3, "test size");

  Image gammad_image = changeGamma(original, 2.2, 4.0);
  gammad_image.write("Output/gammad_image.png");

  Image brightened_image = brightness(original, 5.0);
  brightened_image.write("Output/brightened_image.png");

  Image exposure_image = exposure(original, 5.0);
  exposure_image.write("Output/exposure_image.png");

  float w[3] = {1.0/3, 1.0/3, 1.0/3};
  std::vector<float> weights = std::vector<float>(w, w+3);

  Image gray_image = color2gray(original, weights);
  gray_image.write("Output/gray_image_weight1.png");


  Image gray_image_default = color2gray(original);
  gray_image_default.write("Output/gray_image_default.png");


  std::vector<Image> lumi_and_chromi = lumiChromi(original);
  Image lumi = lumi_and_chromi[0];
  Image chromi = lumi_and_chromi[1];
  lumi.write("Output/lumi.png");
  chromi.write("Output/chromi.png");

  Image brightness_and_contrast = brightnessContrastLumi(original, 1.0, 1.0, 0.5);
  brightness_and_contrast.write("Output/brightness_and_contrast.png");

  Image rgb_to_yuv_back = yuv2rgb(rgb2yuv(original));
  rgb_to_yuv_back.write("Output/rgb_to_yuv_back.png");

  Image saturate_image = saturate(original, 5.0);
  saturate_image.write("Output/saturate_image.png");

  std::vector<Image> illusion_ims = spanish(castle);
  Image castle_lum = illusion_ims[0];
  Image castle_chrom = illusion_ims[1];
  castle_lum.write("Output/castle_lum.png");
  castle_chrom.write("Output/castle_chrom.png");


  Image white_balance = grayworld(original);
  white_balance.write("Output/white_balance.png");
}
