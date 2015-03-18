// hdr.cpp
// Assignment 5


#include "hdr.h"
#include "filtering.h"
#include <math.h>
#include <algorithm>


using namespace std;

/**************************************************************
 //                       HDR MERGING                        //
 *************************************************************/

// PS05 - 2.1 - generate a weight image that indicates which pixels are good to use in hdr
Image computeWeight(const Image &im, float epsilonMini, float epsilonMaxi){
Image output = Image(im.width(), im.height(), im.channels());

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float value = im(x,y,z);
        if (value >= epsilonMini && value <= epsilonMaxi) {
          output(x,y,z) = 1;
        } else {
          output(x,y,z) = 0;
        }
      }
    }
  }
     return output; // change this
}

// PS05 - 2.2 - compute the multiplication factor between a pair of images
float computeFactor(const Image &im1, const Image &w1, const Image &im2, const Image &w2){
  Image ratioImage = Image(im1.width(), im1.height(), im1.channels());

  vector<float> ratioVect;

  for (int x = 0; x < ratioImage.width(); x++) {
    for (int y = 0; y < ratioImage.height(); y++) {
      for (int z = 0; z < ratioImage.channels(); z++) {
        float v1 = im1(x,y,z) + 10e-10;
        float v2 = im2(x,y,z);
        float ratio = v2/v1;
        if (w1(x,y,z) == 1 && w2(x,y,z) == 1) {
          ratioVect.push_back(ratio);
        }
      }
    }
  }

  float median = 1.0;
  if (ratioVect.size() > 0) {
    sort(ratioVect.begin(), ratioVect.end());
    int middleIndex = ratioVect.size()/2;
    median = ratioVect[middleIndex];
  } else {
    cout << "f" << endl;
  }

  if (median==0) {
    cout << "very bad" << endl;
  }
  cout << "med" << median << endl;
  cout << ratioVect.size() << endl;
  return median;
}

// PS05 - 2.3 - merge images to make a single hdr image
Image makeHDR(vector<Image> &imSeq, float epsilonMini, float epsilonMaxi){
  Image im1 = imSeq[0];
  Image output = Image(im1.width(), im1.height(), im1.channels());
  vector<Image> weightSeq;
  vector<float> factorSeq;

  for (int i=0; i < ((int) imSeq.size()); i++) {
    Image im1 = imSeq[i];
    Image weightImage = Image(im1.width(), im1.height(), im1.channels());

    if (i == 0) {
      weightImage = computeWeight(im1, epsilonMini, 1.0);
    } else if (i == imSeq.size()-1) {
      weightImage = computeWeight(im1, 0.0, epsilonMaxi);
    } else {
      weightImage = computeWeight(im1, epsilonMini, epsilonMaxi);
    }

    weightSeq.push_back(weightImage);
  }

  for (int i=0; i < ((int) imSeq.size())-1; i++) {
    Image im1 = imSeq[i];
    Image im2 = imSeq[i+1];
    Image weightImage1 = weightSeq[i];
    Image weightImage2 = weightSeq[i+1];
    float factor = computeFactor(im1, weightImage1, im2, weightImage2);
    factorSeq.push_back(factor);
  }

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float sumWeight2 = 0.0;
        float kj = 1.0;
        float sumValue = 0.0;
          for (int i=0; i < ((int) imSeq.size())-1; i++) {
            float weight1 = weightSeq[i](x, y, z);
            float weight2 = weightSeq[i+1](x, y, z);
            float factor = factorSeq[i+1];
            float ki = factor*kj;
            kj = ki;
            sumWeight2 += weight2;
            sumValue += (weight2/ki)*imSeq[i+1](x, y, z);
          }
          if(sumWeight2 == 0 || sumWeight2 == (imSeq.size()-1)) {
            output(x, y, z) = imSeq[0](x,y,z);
          } else {
            output(x, y, z) = sumValue/sumWeight2;
          }
      }
    }
  }

   return output;
}

/**************************************************************
 //                      TONE MAPPING                        //
 *************************************************************/


// PS05 - 3 - tone map an hdr image
Image toneMap(const Image &im, float targetBase, float detailAmp, bool useBila, float sigmaRange) {
    vector<Image> lumiChromImages = lumiChromi(im);
    Image lumi = lumiChromImages[0];
    Image chromi = lumiChromImages[1];
    Image logLumi = log10Image(lumi);
    float sigma = max(im.width(), im.height())/50.0;
    float sigmaDomain = sigma;
    int truncate = 3;

    Image blurred = Image(im.width(), im.height(), im.channels());
    if (useBila) {
      blurred = bilateral(logLumi, sigmaRange, sigmaDomain, truncate, false);
    } else {
      blurred = gaussianBlur_seperable(logLumi, sigma, truncate, false);
    }

    float blurMaxVal = blurred.max();
    float blurMinVal = image_minnonzero(blurred);
    float k = log10(targetBase)/(blurMaxVal-blurMinVal);

    Image detailLog = logLumi - blurred;
    Image lumiOut = k*(blurred-blurMaxVal) + detailAmp*detailLog;
    Image lumiOutExp = exp10Image(lumiOut);
    Image output = lumiChromi2rgb(lumiOutExp, chromi);

    // add gamma correction back into the image right before returning
    output = changeGamma(output, 1.0f, 1/2.2);
    return detailLog;
}


// Tone Mapping Helpful Functions

// image --> log10Image
Image log10Image(const Image &im) {
  Image output = Image(im.width(), im.height(), im.channels());
  // Taking a linear image im, transform to log10 scale.
  // To avoid infinity issues, make any 0-valued pixel be equal the the minimum
  // non-zero value. See image_minnonzero(im).
  float minvalue = image_minnonzero(im);
  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float value = im(x,y,z);
        if (value == 0.0) {
          output(x, y, z) = log10(minvalue);
        } else {
          output(x, y, z) = log10(value);
        }
      }
    }
  }

    
  return output;
}

// Image --> 10^Image
Image exp10Image(const Image &im) {
  Image output = Image(im.width(), im.height(), im.channels());
  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float value = im(x,y,z);
        output(x, y, z) = (float)(pow(10.0, value));
      }
    }
  }
  return output;
}

// min non-zero pixel value of image
float image_minnonzero(const Image &im) {
  float min = 1e15;

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      for (int z = 0; z < im.channels(); z++) {
        float value = im(x,y,z);
        if (value < min && value != 0.0) {
          min = value;
        }
      }
    }
  }
  return min;
}
