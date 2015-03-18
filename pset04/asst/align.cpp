// align.cpp
// Assignment 4


#include "align.h"
#include <math.h>

using namespace std;

// PS04 - 2.1.1 - Basic denoising by computing the average of a sequence of images
Image denoiseSeq(const vector<Image> &imSeq){
    float length = imSeq.size();
    Image firstImage = imSeq[0];
    Image output = Image(firstImage.width(), firstImage.height(), firstImage.channels());

    for (int i = 0; i < length; i++) {
        Image im = imSeq[i];
        output = output + im/length;
    }

    return output;
}


// PS04 - 2.2.1 - returns an image visualizing the per-pixel and
// per-channel log of the signal-to-noise ratio scaled by scale.
Image logSNR(const vector<Image> &imSeq, float scale){
    int numImages = imSeq.size();
    Image firstImage = imSeq[0];
    Image meanImage = denoiseSeq(imSeq);
    Image ex2Image = Image(firstImage.width(), firstImage.height(), firstImage.channels());
    Image e2xImage = meanImage*meanImage;
    Image logsnrImage = Image(firstImage.width(), firstImage.height(), firstImage.channels());
    Image varianceImage = Image(firstImage.width(), firstImage.height(), firstImage.channels());

    vector<Image> diffImages;

    for (int i = 0; i < numImages; i++) {
        Image image = imSeq[0];
        Image diff = image-meanImage;
        diffImages.push_back(diff);
        ex2Image = ex2Image + (image*image);
    }
    ex2Image = ex2Image/numImages;


    for (int i = 0; i < numImages; i++) {
        Image diff = diffImages[i];
        varianceImage = varianceImage + diff*diff;
    }

    varianceImage = varianceImage/numImages;

    for (int y = 0; y < varianceImage.height(); y++) {
        for (int x = 0; x < varianceImage.width(); x++) {
            for (int z = 0; z < varianceImage.channels(); z++) {
                float varValue = varianceImage(x, y, z);
                float epsilon = 1.0e-06;
                varianceImage(x, y, z) = max(varValue, epsilon);
            }
        }
    }

    Image snrImage = ex2Image/varianceImage;

    for (int y = 0; y < snrImage.height(); y++) {
        for (int x = 0; x < snrImage.width(); x++) {
            for (int z = 0; z < snrImage.channels(); z++) {
                logsnrImage(x, y, z) = 10.0*log10(snrImage(x, y, z))*scale;
            }
        }
    }

    return logsnrImage;
}


// PS04 - 2.3.1 -  returns the (x,y) offset that best aligns im2 to match im1.
vector<int> align(const Image &im1, const Image &im2, int maxOffset){
    vector<int> bestOffset = vector<int>(2,0);
    float smallestError = 1e25;
    int sideLength = maxOffset;
    for (int dx = -sideLength; dx < sideLength; dx++) {
        for (int dy = -sideLength; dy < sideLength; dy++) {
            float totalError = 0.0;
            for (int y = sideLength; y < im1.height() - sideLength; y++) {
                for (int x = sideLength; x < im1.width() - sideLength; x++) {
                    for (int z = 0; z < im1.channels(); z++) {
                        int xx = x + dx;
                        int yy = y + dy;

                        float value1 = im1(x, y, z);
                        float value2 = im2(xx, yy, z);
                        float error2 = pow(value2 - value1, 2.0);
                        totalError += error2;
                    }
                }
            }
            if (totalError < smallestError) {
                smallestError = totalError;
                bestOffset[0] = -1*dx;
                bestOffset[1] = -1*dy;
            }
        }
    }

    return bestOffset; // change this
}

// PS04 - 2.3.2 -  registers all images to the first one in a sequence and outputs
// a denoised image even when the input sequence is not perfectly aligned.
Image alignAndDenoise(const vector<Image> &imSeq, int maxOffset){
    float length = imSeq.size();
    Image firstImage = imSeq[0];
    vector<Image> aligned;
    aligned.push_back(firstImage);
    for (int i = 0; i < length; i++) {
        Image image = imSeq[i];
        vector<int> bestOffset = align(firstImage, image, maxOffset);
        float dx = bestOffset[0];
        float dy = bestOffset[1];
        Image rolled = roll(image, dx, dy);
        aligned.push_back(rolled);
    }

    return denoiseSeq(aligned);
}

// PS04 - 6.1.1 -  6.865 only:
// split a Sergey images to turn it into one 3-channel image.
Image split(const Image &sergeyImg){
    int imageHeight = floor(sergeyImg.height()/3);
    Image output = Image(sergeyImg.width(), imageHeight, 3);

    for (int y = 0; y < imageHeight; y++) {
        for (int x = 0; x < sergeyImg.width(); x++) {
            output(x,y,2) = sergeyImg(x, y);
            output(x,y,1) = sergeyImg(x, y+imageHeight);
            output(x,y,0) = sergeyImg(x, y+imageHeight*2);
        }
    }

    return output;
}

// PS04 - 6.1.2 -  6.865 only:
// aligns the green and blue channels of your rgb channel of a sergey
// image to the red channel. This should return an aligned RGB image
Image sergeyRGB(const Image &sergeyImg, int maxOffset){
    Image colorImage = split(sergeyImg);

    Image output = Image(colorImage.width(), colorImage.height(), 3);

    Image redImage = Image(colorImage.width(), colorImage.height(), 3);
    Image greenImage = Image(colorImage.width(), colorImage.height(), 3);
    Image blueImage = Image(colorImage.width(), colorImage.height(), 3);

    for (int y = 0; y < colorImage.height(); y++) {
        for (int x = 0; x < colorImage.width(); x++) {
            redImage(x,y,0) = colorImage(x, y, 0);
            greenImage(x,y,0) = colorImage(x, y, 1);
            blueImage(x,y,0) = colorImage(x, y, 2);
        }
    }

    vector<int> bestGreenOffset = align(redImage, greenImage, maxOffset);
    vector<int> bestBlueOffset = align(redImage, blueImage, maxOffset);
    int g_dx = bestGreenOffset[0];
    int g_dy = bestGreenOffset[1];
    int b_dx = bestBlueOffset[0];
    int b_dy = bestBlueOffset[1];

    Image greenRolled = roll(greenImage, g_dx, g_dy);
    Image blueRolled = roll(blueImage, b_dx, b_dy);

    for (int y = 0; y < colorImage.height(); y++) {
        for (int x = 0; x < colorImage.width(); x++) {
            output(x, y, 0) = redImage(x, y, 0);
            output(x, y, 1) = greenRolled(x, y, 0);
            output(x, y, 2) = blueRolled(x, y, 0);
        }
    }

    return output;
}


/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// This circularly shifts an image by xRoll in the x direction and
// yRoll in the y direction. xRoll and yRoll can be negative or postive
Image roll(const Image &im, int xRoll, int yRoll){
    
    int xNew, yNew;
    Image imRoll(im.width(), im.height(), im.channels());
    
    // for each pixel in the original image find where its corresponding
    // location is in the rolled image
    for (int x=0; x<im.width(); x++){
        for (int y=0; y<im.height(); y++){
            
            // use modulo to figure out where the new location is in the
            // rolled image. Then take care of when this returns a negative number
            xNew = (x + xRoll) % im.width();
            yNew = (y + yRoll) % im.height();
            xNew = (xNew<0)*(imRoll.width() + xNew) + (xNew>=0)*xNew;
            yNew = (yNew<0)*(imRoll.height() + yNew) + (yNew>=0)*yNew;
            
            // assign the rgb values for each pixel in the original image to
            // the location in the new image
            for (int z=0; z<im.channels(); z++){
                imRoll(xNew, yNew, z) = im(x,y,z);
            }
        }
    }
    
    // return the rolled image
    return imRoll;
}
