#include "homography.h"
#include "matrix.h"

using namespace std;


// PS06: apply homography
void applyhomography(const Image &source, Image &out, Matrix &H, bool bilinear) {
  int width = out.width();
  int height = out.height();

  applyhomographyGeneral(source, out, H, bilinear, 0, width, 0, height);
}

// PS06: 6.865 or extra credit: apply homography fast
void applyhomographyFast(const Image &source, Image &out, Matrix &H, bool bilinear) {
  // do something :)
  vector<float> bbox = computeTransformedBBox(source.width(), source.height(), H);

  int x_min = bbox[0];
  int x_max = bbox[1];
  int y_min = bbox[2];
  int y_max = bbox[3];

  applyhomographyGeneral(source, out, H, bilinear, x_min, x_max, y_min, y_max);
}

void applyhomographyGeneral(const Image &source, Image &out, Matrix &H, bool bilinear, int x_start, int x_end, int y_start, int y_end) {
  Matrix H_inverse = H.inverse();

  float x_pos = 0;
  float y_pos = 0;

  float m[]  = {0, 0, 0};

  for (int x = x_start; x < x_end; x++) {
    for (int y = y_start; y < y_end; y++) {
      m[0] = x;
      m[1] = y;
      m[2] = 1;

      Matrix out_matrix(1, 3, m);
      Matrix in_matrix = H_inverse.multiply(out_matrix);
      x_pos = in_matrix(0, 0) / in_matrix(0, 2);
      y_pos = in_matrix(0, 1) / in_matrix(0, 2);
      bool valid_pixel = (x_pos >= 0 && x_pos < source.width() && y_pos >= 0 && y_pos < source.height());
      if (valid_pixel) {
          for (int c = 0; c < source.channels(); c++) {
            out(x, y, c) = source.smartAccessor(x_pos, y_pos, c, bilinear);
          }
      }
    }
  }
  
}

// PS06: compute homography given a list of point pairs
Matrix computeHomography(const float listOfPairs[4][2][3]) {

  Matrix A(8,8);

  for (int i = 0; i < 4; i++) {
    addConstraint(A, 2*i, listOfPairs[i]);
  }

  Matrix B (1, 8);

  for (int i = 0; i < 4; i++) {
    B(0, 2*i) = listOfPairs[i][1][0];
    B(0, 2*i+1) = listOfPairs[i][1][1];
  }

  Matrix X = A.inverse().multiply(B);

  Matrix Y (3,3);
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      int index = i+3*j;
      if (i+3*j < 8) {
        Y(i, j) = X(0, i+3*j);
      }
    }
  }
  Y(2,2) = 1;

  return Y;
}

// PS06: optional function that might help in computeHomography
void addConstraint(Matrix &systm,  int i, const float constr[2][3]) {
  float x = constr[0][0];
  float y = constr[0][1];

  float xp = constr[1][0];
  float yp = constr[1][1];

  systm(0, i) = x;
  systm(1, i) = y;
  systm(2, i) = 1;
  systm(3, i) = 0;
  systm(4, i) = 0;
  systm(5, i) = 0;
  systm(6, i) = -x*xp;
  systm(7, i) = -y*xp;

  systm(0, i+1) = 0;
  systm(1, i+1) = 0;
  systm(2, i+1) = 0;
  systm(3, i+1) = x;
  systm(4, i+1) = y;
  systm(5, i+1) = 1;
  systm(6, i+1) = -x*yp;
  systm(7, i+1) = -y*yp;
}

// PS06: compute a transformed bounding box
// returns [xmin xmax ymin ymax]
vector<float> computeTransformedBBox(int imwidth, int imheight, Matrix H) {

  float top_left[3] = {0.0, 0.0, 1.0};
  float bottom_left[3] = {imwidth-1.0, 0.0, 1.0};

  float top_right[3] = {0.0, imheight-1.0, 1.0};
  float bottom_right[3] = {imwidth-1.0, imheight-1.0, 1.0};

  Matrix m_top_left(1, 3, top_left);
  Matrix m_bottom_left(1, 3, bottom_left);

  Matrix m_top_right(1, 3, top_right);
  Matrix m_bottom_right(1, 3, bottom_right);

  Matrix m_top_left_T = H.multiply(m_top_left);
  Matrix m_bottom_left_T = H.multiply(m_bottom_left);
  Matrix m_top_right_T = H.multiply(m_top_right);
  Matrix m_bottom_right_T = H.multiply(m_bottom_right);

  Matrix m_h_top_left_T = homogenize(m_top_left_T);
  Matrix m_h_top_right_T = homogenize(m_top_right_T);
  Matrix m_h_bottom_left_T = homogenize(m_bottom_left_T);
  Matrix m_h_bottom_right_T = homogenize(m_bottom_right_T);

  vector<float> y_values = vector<float>(4, 0.0);
  y_values[0] = m_h_top_left_T(0, 1);
  y_values[1] = m_h_top_right_T(0, 1);
  y_values[2] = m_h_bottom_left_T(0, 1);
  y_values[3] = m_h_bottom_right_T(0, 1);

  vector<float> x_values = vector<float>(4, 0.0);
  x_values[0] = m_h_top_left_T(0, 0);
  x_values[1] = m_h_top_right_T(0, 0);
  x_values[2] = m_h_bottom_left_T(0, 0);
  x_values[3] = m_h_bottom_right_T(0, 0);

  float y_min = min_vec_elem(y_values);
  float y_max = max_vec_elem(y_values);
  float x_min = min_vec_elem(x_values);
  float x_max = max_vec_elem(x_values);

  vector<float> minmax = vector<float>(4, 0.0);
  minmax[0] = x_min;
  minmax[1] = x_max;
  minmax[2] = y_min;
  minmax[3] = y_max;

  return minmax;
}

// PS06: homogenize vector v.
// this function is not required, but would be useful for use in
// computeTransformedBBox()
Matrix homogenize(Matrix &v) {
  Matrix hom(1,3);

  hom(0,0) = v(0,0)/v(0,2);
  hom(0,1) = v(0,1)/v(0,2);
  hom(0,2) = 1;

  return hom;
}

// PS06: compute a 3x3 translation Matrix
Matrix translate(vector<float> B) {
  float yy = -1*B[2];
  float xx = -1*B[0];

  float matrix_data[9] = {
    1, 0, xx,
    0, 1, yy,
    0, 0, 1
  };

  Matrix trans(3, 3, matrix_data);
  return trans;
}

// PS06: compute the union of two bounding boxes
vector <float> bboxUnion(vector<float> B1, vector<float> B2) {

  float xmin = min(B1[0], B2[0]);
  float xmax = max(B1[1], B2[1]);
  float ymin = min(B1[2], B2[2]);
  float ymax = max(B1[3], B2[3]);

  vector<float> B = vector<float>(4, 0.0);
  B[0] = xmin;
  B[1] = xmax;
  B[2] = ymin;
  B[3] = ymax;

  return B;
}

// PS06: stitch two images given a list or 4 point pairs
Image stitch(const Image &im1, const Image &im2, const float listOfPairs[4][2][3]) {
  Matrix H = computeHomography(listOfPairs);

  vector<float> bbox1 = computeTransformedBBox(im1.width(), im1.height(), H);

  vector<float> bbox2 = vector<float>(4, 0.0);
  bbox2[0] = 0.0;
  bbox2[1] = im2.width() - 1.0;
  bbox2[2] = 0.0;
  bbox2[3] = im2.height() - 1.0;

  vector<float> union_box = bboxUnion(bbox1, bbox2);

  Matrix translated_bbox = translate(union_box);

  float tx = translated_bbox(2, 0);
  float ty = translated_bbox(2, 1);

  int out_width = union_box[1] + tx;
  int out_height = union_box[3] + ty;

  Image output(out_width, out_height, im1.channels());

  Matrix tb_H = translated_bbox.multiply(H);

  applyhomography(im2, output, translated_bbox, true);
  applyhomography(im1, output, tb_H, true);

  return output;
}

// PS06: useful for bounding-box debugging.
Image drawBoundingBox(const Image &im, vector<float> minmax) {
  int x_min = minmax[0];
  int x_max = minmax[1];
  int y_min = minmax[2];
  int y_max = minmax[3];

  Image output(im.width(), im.height(), im.channels());

  for (int x = 0; x < im.width(); x++) {
    for (int y = 0; y < im.height(); y++) {
      for (int c = 0; c < im.channels(); c++) {
        bool in_x_bounds = (x <= x_max && x >= x_min);
        bool in_y_bounds = (y <= y_max && y >= y_min);

        bool onTopLine = (y == y_max) && in_x_bounds;
        bool onBottomLine = (y == y_min) && in_x_bounds;
        bool onLeftLine = (x == x_min) && in_y_bounds;
        bool onRightLine = (x == x_max) && in_y_bounds;

        if (onTopLine || onBottomLine || onLeftLine || onRightLine) {
          output(x, y, c) = 0;
        } else {
          output(x, y, c) = im(x, y, c);
        }
      }
    }
  }
  return output;
}


/***********************************************************************
 * Do not edit code  below this line
 **********************************************************************/
// get the minimum vector element
float min_vec_elem(vector<float> v) {
  float mmin = FLT_MAX;
  for (int i = 0; i < (int)v.size(); i++) mmin = min(mmin, v[i]);
  return mmin;
}

// get the maximum vector element
float max_vec_elem(vector<float> v) {
  float mmax = -FLT_MAX;
  for (int i = 0; i < (int)v.size(); i++) mmax = max(mmax, v[i]);
  return mmax;
}
