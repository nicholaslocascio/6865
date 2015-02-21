// morphing.cpp
// Assignment 2


#include "morphing.h"
using namespace std;

/**************************************************************
 *            IMAGE WARPING/MORPHING FUNCTIONS                *
 *************************************************************/

// PS02 - 4.3.1: warp an entire image according to a pair of segments.
Image warpBy1(const Image &im, Segment &segBefore, Segment &segAfter){
    return Image(0); // change this
}

// PS02 - 4.4.2: warp an image according to a vector of before and after segments using segment weighting
Image warp(const Image &im, vector<Segment> &segsBefore, vector<Segment> &segsAfter, float a, float b, float p){
    return Image(0); // change this
}

//PS02 - 4.5.1: return a vector of N images in addition to the two inputs that morphs going between im1 and im2 for the corresponding segments
vector<Image> morph(const Image &im1, const Image &im2, vector<Segment> &segsBefore, vector<Segment> &segsAfter, int N, float a, float b, float p){
    return vector<Image>(); // change this
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

    // finish implementing this
}


// PS02 - 4.2.2: Implement the computation of the u coordinate of a point (x,y) with respect to a segment
float Segment::getU(float &x, float &y){
    return 0; // change this
}


// PS02 - 4.2.2: Implement the computation of the v coordinate of a point (x,y) with respect to a segment
float Segment::getV(float &x, float &y){
    return 0; // change this
}


// PS02 - 4.2.2: compute the new (x, y) position of a point given by the (u,v) location relative to another segment.
// return the point (x,y) in a 2-element vector
vector<float> Segment::UVtoX(float &u, float &v){
    return vector<float>(); // change this
}

// PS02 - 4.2.3: Implement distance from a point (x,y) to the segment. Remember the 3 cases from class
float Segment::dist(float &x, float &y){
    return 0; // change this
}

// PS02 - 4.4.1: compute the weight of a segment to a point (x,y) given the weight parameters a,b, and p
float Segment::weight(float &x, float &y, float &a, float &b, float &p){
    return 0; // change this
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
