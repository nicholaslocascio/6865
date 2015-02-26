// morphing.cpp
// Assignment 2


#include "morphing.h"
using namespace std;

/**************************************************************
 *            IMAGE WARPING/MORPHING FUNCTIONS                *
 *************************************************************/

// PS02 - 4.3.1: warp an entire image according to a pair of segments.
Image warpBy1(const Image &im, Segment &segBefore, Segment &segAfter){

  Image output = Image(im.width(), im.height(), im.channels());

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float xxx = (float) x;
        float yyy = (float) y;

        float u = segAfter.getU(xxx, yyy);
        float v = segAfter.getV(xxx, yyy);

        vector<float> X = segBefore.UVtoX(u, v);
        float xx = X[0];
        float yy = X[1];

        float value = interpolateLin(im, xx, yy, z, true);
        output(x, y, z) = value;
      }
    }
  }
  return output;
}

// PS02 - 4.4.2: warp an image according to a vector of before and after segments using segment weighting
Image warp(const Image &im, vector<Segment> &segsBefore, vector<Segment> &segsAfter, float a, float b, float p){
  Image output = Image(im.width(), im.height(), im.channels());

  for (int x = 0; x < output.width(); x++) {
    for (int y = 0; y < output.height(); y++) {
      for (int z = 0; z < output.channels(); z++) {
        float xValue = 0.0;
        float yValue = 0.0;
        float totalWeight = 0.0;
        for(int s = 0; s < segsBefore.size(); s++) {

            Segment segAfter = segsAfter[s];
            Segment segBefore = segsBefore[s];

            float xxx = (float) x;
            float yyy = (float) y;

            float w = segBefore.weight(xxx, yyy, a, b, p);
            totalWeight += w;

            float u = segAfter.getU(xxx, yyy);
            float v = segAfter.getV(xxx, yyy);

            vector<float> X = segBefore.UVtoX(u, v);
            float xx = X[0];
            float yy = X[1];

            xValue += w*xx;
            yValue += w*yy;
        }
        xValue /= totalWeight;
        yValue /= totalWeight;

        output(x, y, z) = interpolateLin(im, xValue, yValue, z, true);;
      }
    }
  }
  return output;
}

//PS02 - 4.5.1: return a vector of N images in addition to the two inputs that morphs going between im1 and im2 for the corresponding segments
vector<Image> morph(const Image &im1, const Image &im2, vector<Segment> &segsBefore, vector<Segment> &segsAfter, int N, float a, float b, float p){
    vector<Image> images =  vector<Image>();
        cout << "N = " << N << endl;

    for (int i = 0; i < N+1; i++) {
        cout << "i = " << i << endl;
        Image imageA = warp(im1, segsBefore, segsAfter, a, b, p);
        Image imageB = warp(im2, segsBefore, segsAfter, a, b, p);
        float alpha = ((float)(i))/(N);
        Image imageC = imageA*(1.0-alpha) + imageB*(alpha);

        cout << alpha << endl;

        images.push_back(imageC);
    }
    return images;
}

/**************************************************************
 *                 SEGMENT CLASS FUNCTIONS                    *
 *************************************************************/

// PS02 - 4.2.1: Segment constructor takes in 2 points (x1,y1) and (x2,y2) correspoding to the ends of a segment and computes:
// P - 2-element vector to point (x1, y1)
// Q - 2-element vector to pont (x2, y2)
// PQ - 2-element vector from P to Q
// PQ2 - float containing squared distance between P and Q
// PQlength - float containing distance between P and Q
// PQDivByPQlength2 - 2-element vector PQ normalized by PQ2
// perpPQ - 2-element vector perpendicular to PQ
// perpPQDivByPQlength - 2-element vector perpPQ normalized by PQlength
Segment::Segment(const float &x1, const float &y1, const float &x2, const float &y2) {

    P       = vector<float>(2,0);
    Q       = vector<float>(2,0);
    perpPQ  = vector<float>(2,0);

    P[0] = x1;
    P[1] = y1;
    Q[0] = x2;
    Q[1] = y2;

    PQ = Segment::subtract(P, Q);
    perpPQ[0] = -1*PQ[1];
    perpPQ[1] = PQ[0];

    PQ2 = dot(PQ, PQ);
    PQlength = pow(PQ2, 0.5);

    float tmp = 1.0 / PQ2;
    PQDivByPQlength2 = Segment::scalarMult(PQ, tmp);
    tmp = 1.0 / PQlength;
    perpPQDivByPQlength = Segment::scalarMult(perpPQ, tmp);
}


// PS02 - 4.2.2: Implement the computation of the u coordinate of a point (x,y) with respect to a segment
float Segment::getU(float &x, float &y){
    vector<float> X = vector<float>(2);
    X[0] = x;
    X[1] = y;
    vector<float> XP = Segment::subtract(P, X);
    float u = Segment::dot(XP, PQ)/PQ2;
    return u;
}


// PS02 - 4.2.2: Implement the computation of the v coordinate of a point (x,y) with respect to a segment
float Segment::getV(float &x, float &y){
    vector<float> X = vector<float>(2);
    X[0] = x;
    X[1] = y;
    vector<float> XP = Segment::subtract(P, X);
    float v = Segment::dot(XP, perpPQ)/PQlength;
    return v;
}

// PS02 - 4.2.2: compute the new (x, y) position of a point given by the (u,v) location relative to another segment.
// return the point (x,y) in a 2-element vector
vector<float> Segment::UVtoX(float &u, float &v){
    vector<float> newXY = vector<float>(2);
    newXY[0] =  P[0] + PQ[0]*u + perpPQ[0]*v/PQlength;
    newXY[1] =  P[1] + PQ[1]*u + perpPQ[1]*v/PQlength;
    return newXY;
}

// PS02 - 4.2.3: Implement distance from a point (x,y) to the segment. Remember the 3 cases from class
float Segment::dist(float &x, float &y){
    float u = getU(x, y); 
    float v = getV(x, y);
    vector<float> X = vector<float>(2);
    X[0] = x;
    X[1] = y;

    if (u<1 && u>0) {
        return fabs(v);
    } else if (u<0) {
        vector<float> px = Segment::subtract(X, P);
        return pow(px[0]*px[0] + px[1]*px[1], 0.5);
    } else {
        vector<float> qx = Segment::subtract(X, Q);
        return pow(qx[0]*qx[0] + qx[1]*qx[1], 0.5);
    }
}

// PS02 - 4.4.1: compute the weight of a segment to a point (x,y) given the weight parameters a,b, and p
float Segment::weight(float &x, float &y, float &a, float &b, float &p){
    float length = PQlength;
    float d = dist(x, y);
    float w = pow(pow(length, p)/(a + d), b);
    return w;
}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// subtracts 2 vectors of the same length.
vector<float> Segment::subtract(vector<float> &vec1, vector<float> &vec2){
// create vector from vec1 to vec2

    vector<float> vec_1_to_2 (vec1.size(), 0);

    if(vec1.size() == vec2.size()){

        for (int i=0; i < (int) vec1.size(); i++){
            vec_1_to_2[i] = vec2[i] - vec1[i];
        }

    }else{
        throw MismatchedSizeException();
    }

    return vec_1_to_2;
}

// takes the dot product between 2 vectors of the same length
float Segment::dot(vector<float> &vec1, vector<float> &vec2){

    float dotProd = 0;

    if(vec1.size() == vec2.size()){

        for (int i=0; i< (int) vec1.size(); i++){
            dotProd += vec2[i]*vec1[i];
        }

    }else{
        throw MismatchedSizeException();
    }

    return dotProd;
}

// mutliplies an entire vector by a scalor value
vector<float> Segment::scalarMult(vector<float> &vec, float &factor){

    vector<float> nVec (vec.size(), 0);
    for(int i=0; i< (int) vec.size(); i++){
        nVec[i] = vec[i]*factor;
    }
    return nVec;
}

// destructor
Segment::~Segment() { } // Nothing to clean up
