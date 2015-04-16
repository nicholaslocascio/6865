#include "blending.h"
#include "matrix.h"
#include <math.h>
#include <ctime>

using namespace std;

/*****************************************************************************
 * blending related functions re-written from previous asasignments
 *****************************************************************************/

// PSet 08 instead of writing source in out, *add* the source to out based on the weight
// so out(x,y) = out(x, y) + weight * image
void applyhomographyBlend(const Image &source, const Image &weight, Image &out, Matrix &H, bool bilinear) {
  Matrix H_inverse = H.inverse();

  int x_start = 0;
  int x_end = out.width();
  int y_start = 0;
  int y_end = out.height();

  for (int x = x_start; x < x_end; x++) {
    for (int y = y_start; y < y_end; y++) {
      float out_data[]  = {0, 0, 0};
      out_data[0] = x;
      out_data[1] = y;
      out_data[2] = 1;

      Matrix out_matrix(1, 3, out_data);
      Matrix in_matrix = H_inverse.multiply(out_matrix);
      in_matrix = homogenize(in_matrix);
      float x_pos = in_matrix(0, 0);
      float y_pos = in_matrix(0, 1);
      bool valid_pixel = (x_pos >= 0 && x_pos < source.width() && y_pos >= 0 && y_pos < source.height());
      if (valid_pixel) {
          float weight_value = interpolateLin(weight, x_pos, y_pos, 0, bilinear);
          for (int c = 0; c < source.channels(); c++) {
            float old_value = out(x, y, c);
            float addition = interpolateLin(source, x_pos, y_pos, c, bilinear);
            float weighted_addition = weight_value*addition;
            float new_value = old_value + weighted_addition;
            out(x, y, c) = new_value;
          }
      }
    }
  }

}

// PSet 08 stitch using image weights.
// note there is no weight normalization.
Image stitchLinearBlending(const Image &im1, const Image &im2, const Image &we1, const Image &we2, Matrix H) {
  vector<float> bbox_1 = computeTransformedBBox(im1.width(), im1.height(), H);
  vector<float> bbox_2 = vector<float>(4, 0.0);

  bbox_2[0] = 0.0;
  bbox_2[1] = im2.width() - 1.0;
  bbox_2[2] = 0.0;
  bbox_2[3] = im2.height() - 1.0;

  vector<float> union_box = bboxUnion(bbox_1, bbox_2);
  Matrix bbox_translated = translate(union_box);

  int out_width = union_box[1] + bbox_translated(2, 0);
  int out_height = union_box[3] + bbox_translated(2, 1);
  Image out(out_width, out_height, im1.channels());

  Matrix bbt_H = bbox_translated.multiply(H);
  applyhomographyBlend(im2, we2, out, bbox_translated, true);
  applyhomographyBlend(im1, we1, out, bbt_H, true);

  return out;
}

/*****************************************************************************
 * blending functions Pset 08
 *****************************************************************************/

// PSet 08 weight image
Image blendingweight(int imwidth, int imheight) {
  Image weights(imwidth, imheight, 1);

  for (int x = 0; x < imwidth; x++) {
    for (int y = 0; y < imheight; y++) {
      float x_center = imwidth/2.0;
      float y_center = imheight/2.0;

      float dx = x_center - fabs(x - x_center);
      float dy = y_center - fabs(y - y_center);

      float x_val = dx / x_center;
      float y_val = dy / y_center;

      float w = x_val*y_val;
      weights(x, y, 0) = w;
    }
  }
  return weights;
}

// Optional: low freq and high freq (2-scale) decomposition
vector<Image> scaledecomp(const Image &im, float sigma) {
  vector <Image> ims;

  Image low_freq_dec = gaussianBlur_seperable(im, sigma);
  ims.push_back(low_freq_dec);

  Image high_freq_dec = (im - low_freq_dec);
  ims.push_back(high_freq_dec);

  return ims;
}

// PSet 08 stitch using different blending models
// blend can be 0 (none), 1 (linear) or 2 (2-layer)
// use more helper functions if necessary
Image stitchBlending(Image &im1, Image &im2, Matrix H, int blend) {
  Image output(0);
  if (blend == 0) {
      vector<float> bbox_1 = computeTransformedBBox(im1.width(), im1.height(), H);
      vector<float> bbox_2 = vector<float>(4, 0.0);

      bbox_2[0] = 0.0;
      bbox_2[1] = im2.width() - 1.0;
      bbox_2[2] = 0.0;
      bbox_2[3] = im2.height() - 1.0;

      vector<float> union_box = bboxUnion(bbox_1, bbox_2);
      Matrix bbox_translated = translate(union_box);

      int out_width = union_box[1] + bbox_translated(2, 0);
      int out_height = union_box[3] + bbox_translated(2, 1);
      Image out(out_width, out_height, im1.channels());

      Matrix bbt_H = bbox_translated.multiply(H);
      applyhomography(im2, out, bbox_translated, true);
      applyhomography(im1, out, bbt_H, true);
      output = out;
  } else if (blend == 1) {
      Image w1 = blendingweight(im1.width(), im1.height());
      Image w2 = blendingweight(im2.width(), im2.height());

      Image empty1 = Image(w1.width(), w2.height(), 1);
      Image white1 = empty1 + 1;
      Image empty2 = Image(w1.width(), w2.height(), 1);
      Image white2 = empty2 + 1;

      Image sum_of_weights = stitchLinearBlending(w1, w2, white1, white2, H);

      Image blended_images_with_weights = stitchLinearBlending(im1, im2, w1, w2, H);
      Image out = comp_division(blended_images_with_weights, sum_of_weights);
      output = out;
  } else if (blend == 2) {
      Image w1 = blendingweight(im1.width(), im1.height());
      Image w2 = blendingweight(im2.width(), im2.height());

      Image empty1 = Image(w1.width(), w2.height(), 1);
      Image white1 = empty1 + 1;
      Image empty2 = Image(w1.width(), w2.height(), 1);
      Image white2 = empty2 + 1;

      Image sum_of_weights = stitchLinearBlending(w1, w2, white1, white2, H);

      vector<Image> im_1_decomp = scaledecomp(im1);
      vector<Image> im_2_decomp = scaledecomp(im2);

      Image low_freq_im1 = im_1_decomp[0];
      Image high_freq_im1 = im_1_decomp[1];

      Image low_freq_im2 = im_2_decomp[0];
      Image high_freq_im2 = im_2_decomp[1];

      Image high_1 = stitchLinearBlending(high_freq_im1, high_freq_im2, white1, empty2, H);
      Image high_2 = stitchLinearBlending(high_freq_im1, high_freq_im2, empty1, white2, H);

      Image w_1O = stitchLinearBlending(w1, w2, white1, empty2, H);
      Image w_2O = stitchLinearBlending(w1, w2, empty1, white2, H);

      Image blended_lows = stitchLinearBlending(low_freq_im1, low_freq_im2, w1, w2, H);
      Image blended_low_normalized = comp_division(blended_lows, sum_of_weights);

      for (int x = 0; x < high_2.width(); x++) {
        for (int y = 0; y < high_2.height(); y++) {
          for (int c = 0; c < high_2.channels(); c++) {
            float add_value = 0.0;
            if (w_1O(x, y, 0) > w_2O(x, y,0)) {
              add_value = high_1(x, y, c);
            } else {
              add_value = high_2(x, y, c);
            }
            blended_low_normalized(x, y, c) = blended_low_normalized(x, y, c) + add_value;
          }
        }
      }
      output = blended_low_normalized;
    }
  return output;
}

Image comp_division(Image &a, Image &b) {
  Image out(a.width(), a.height(), a.channels());
  for (int x = 0; x < out.width(); x++) {
        for (int y = 0; y < out.height(); y++) {
            for (int c = 0; c < out.channels(); c++) {
              float a_val = a(x, y, c);
              float b_val = b(x, y, 0);
              if (b_val == 0) {
                out(x, y, c) = 0;
              } else {
                out(x, y, c) = a_val/b_val;
              }
            }
        }
  }

    return out;
}

// PSet 08 auto stitch
Image autostitch(Image &im1, Image &im2, int blend, float blurDescriptor, float radiusDescriptor) {

  vector<Point> corners_1 = HarrisCorners(im1);
  vector<Feature> feature_1 = computeFeatures(im1, corners_1);

  vector<Point> corners_2 = HarrisCorners(im2);
  vector<Feature> feature_2 = computeFeatures(im2, corners_2);

  vector<Correspondance> correspondances = findCorrespondences(feature_1, feature_2);
  Matrix H = RANSAC(correspondances, 120);

  return stitchBlending(im1, im2, H, blend);
}

/************************************************************************
 * Coolness: mini planets.
 ************************************************************************/

Image pano2planet(const Image &pano, int newImSize, bool clamp) {
  float new_mid = newImSize/2.0;
  float R_boundary = (newImSize/2.0)*sqrt(2.0);
  Image out(newImSize, newImSize, pano.channels());

  for (int x = 0; x < newImSize; x++) {
    for (int y = 0; y < newImSize; y++) {
      float dx = (x - new_mid);
      float dy = (y - new_mid);

      float r = sqrt(pow(dy, 2.0) + pow(dx,2.0)); 
      float theta = atan2(dx, dy);

      int old_width = pano.width();
      int old_height = pano.height();
      float old_mid = old_width / 2.0;

      float xx = 0.0;
      if (theta > 0.0) {
        xx = old_mid + (1.0 - theta/pi)*old_mid;
      } else {
        xx =  -1.0*(theta/pi)*old_mid;
      }
      float r_perc = r/R_boundary;
      float yy = (1.0 - r_perc)*old_height;
  
      for (int c = 0; c < pano.channels(); c++) {
        out(x, y, c) = interpolateLin(pano, xx, yy, c, true);
      }
    }
  }
  return out;
}


/************************************************************************
 * 6.865: Stitch N images into a panorama
 ************************************************************************/

// Pset08-865. Compute sequence of N-1 homographies going from Im_i to Im_{i+1}
// Implement me!
vector<Matrix> sequenceHs(vector<Image> ims, float blurDescriptor, float radiusDescriptor) {
  vector<Matrix> Hs;

  for (int i = 0; i < ims.size()-1; i++) {
    Image im1 = ims[i];
    Image im2 = ims[i+1];
    vector<Point> corners_1 = HarrisCorners(im1);
    vector<Feature> feature_1 = computeFeatures(im1, corners_1);

    vector<Point> corners_2 = HarrisCorners(im2);
    vector<Feature> feature_2 = computeFeatures(im2, corners_2);

    vector<Correspondance> correspondances = findCorrespondences(feature_1, feature_2);
    Matrix H = RANSAC(correspondances, 120);
    Hs.push_back(H);
  }
  return Hs;

}

// Pset08-865. stack homographies:
//   transform a list of (N-1) homographies that go from I_i to I_i+1
//   to a list of N homographies going from I_i to I_refIndex.
vector <Matrix> stackHomographies(vector <Matrix> Hs, int refIndex) {
  vector <Matrix> gHs;

  int N = Hs.size() + 1;
  Matrix data(3, 3);

  gHs = vector<Matrix>(N, data);
  Matrix eye_matrix = eye(3);
  gHs[refIndex] = eye_matrix;

  int start_i = refIndex-1;
  int end_i = -1;

  for (int i = start_i; i > end_i; i--) {
    Matrix M = gHs[i+1].multiply(Hs[i]);
    gHs[i] = M;
  }

  int start_i_2 = refIndex + 1;
  int end_i_2 = gHs.size();

  for (int i2 = start_i_2; i2 < end_i_2; i2++) {
    Matrix gg = gHs[i2-1].multiply(Hs[i2-1]);
    Matrix gginv = gg.inverse();
    gHs[i2] = gginv;
  }

  return gHs;
}

// Pset08-865: compute bbox around N images given one main reference.
vector<float> bboxN(const vector<Matrix> &Hs, const vector<Image> &ims) {
  vector<float> ubbox;

  for (int i = 0; i < ims.size()-1; i++) {
    Matrix H = Hs[i];
    Image im2 = ims[i];
    Image im1 = ims[i+1];

    vector<float> bbox_1 = computeTransformedBBox(im1.width(), im1.height(), H);
    vector<float> bbox_2 = vector<float>(4, 0.0);

    bbox_2[0] = 0.0;
    bbox_2[1] = im2.width() - 1.0;
    bbox_2[2] = 0.0;
    bbox_2[3] = im2.height() - 1.0;

    vector<float> union_box = bboxUnion(bbox_1, bbox_2);
    if (i == 0) {
      ubbox = union_box;
    } else {
      ubbox = bboxUnion(union_box, ubbox);
    }
  }

  return ubbox;
}

// Pset08-865.
// Implement me!
Image autostitchN(vector<Image> ims, int refIndex, float blurDescriptor, float radiusDescriptor) {
  
  vector<Matrix> Hs = sequenceHs(ims, blurDescriptor, radiusDescriptor);
  vector <Matrix> gHs = stackHomographies(Hs, refIndex);

  vector<float> bbox_all = bboxN(gHs, ims);
  Matrix bbox_all_translated = translate(bbox_all);

  int out_width = bbox_all[1] + bbox_all_translated(2, 0);
  int out_height = bbox_all[3] + bbox_all_translated(2, 1);

  Image output(out_width, out_height, 3);

  for (int i = 0; i < ims.size()-1; i++) {
    Image im1(0);
    Image im2(0);

    if (i < refIndex) {
      im1 = ims[i];
      im2 = ims[i+1];
    } else {
      im2 = ims[i];
      im1 = ims[i+1];
    }

    Matrix H = gHs[i];

    vector<float> bbox_1 = computeTransformedBBox(im1.width(), im1.height(), H);
    vector<float> bbox_2 = vector<float>(4, 0.0);

    bbox_2[0] = 0.0;
    bbox_2[1] = im2.width() - 1.0;
    bbox_2[2] = 0.0;
    bbox_2[3] = im2.height() - 1.0;

    vector<float> union_box = bboxUnion(bbox_1, bbox_2);
    Matrix bbox_translated = translate(union_box);

    Matrix bbt_H = bbox_translated.multiply(H);
    applyhomography(im2, output, bbox_translated, true);
    applyhomography(im1, output, bbt_H, true);
  }

  return output;
}


/******************************************************************************
 * Helpful functions
 *****************************************************************************/

 // copy a single-channeled image to several channels
Image copychannels(const Image &im, int nChannels) {
 assert(im.channels() == 1, "image must have one channel");
 Image oim(im.width(), im.height(), nChannels);

 for (int i = 0; i < im.width(); i++) {
   for (int j = 0; j < im.height(); j++) {
     for (int c = 0; c < nChannels; c++) {
       oim(i, j, c) = im(i, j);
     }
   }
 }
 return oim;
}

// create an n x n identity matrix
Matrix eye(int n) {
  Matrix m(n, n);
  for (int i = 0; i < n; i++) m(i, i) = 1;
  return m;
}
