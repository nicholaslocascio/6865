#include "morphing.h"


using namespace std;


// test the smart accessor function
void testSmartAccessor(){
    
    // load an image and create 2 images that will test the smart accessor
    const Image input("./Input/bear.png");
    Image clampTrue(input.width(), input.height(), input.channels());
    Image clampFalse(input.width(), input.height(), input.channels());
    
    for (int z=0; z<input.channels(); z++){
        for (int x=0; x<input.width(); x++){
            for (int y=0; y<input.height(); y++){
                
                // replace non-valid pixel values with the value of the nearest pixel
                clampTrue(x,y,z) = input.smartAccessor(x-10, y-10, z, true);
                // replace non-valid pixel values with black (value=0)
                clampFalse(x,y,z) = input.smartAccessor(x-10, y-10, z);
            }
        }
    }
    
    clampTrue.write("./Output/smartAccessor_clampTrue.png");
    clampFalse.write("./Output/smartAccessor_clampFalse.png");
    
}

// a function to test scaling
void testScaling(){
    
    // load in the image and print out original size
    float fact = 2.253;
    const Image bostonim("./Input/BostonRainbow-crop-400.png");
    printf("Boston image is %dx%dx%d\n", bostonim.width(), bostonim.height(), bostonim.channels());
    
    // scale using NN interpolation and print the size of the new image
    Image scaledNN = scaleNN(bostonim, fact);
    scaledNN.write("./Output/BostonRainbow-scaled-NN.png");
    printf("Scaled-NN image is %dx%dx%d\n", scaledNN.width(), scaledNN.height(), scaledNN.channels());
    
    // scale using bilinear interpolation and print the size of the new image
    Image scaledLin = scaleLin(bostonim, fact);
    scaledLin.write("./Output/BostonRainbow-scaled-Lin.png");
    printf("Scaled-Lin image is %dx%dx%d\n", scaledLin.width(), scaledLin.height(), scaledLin.channels());
    
}

// a function to test rotation for those who have done it
void testRotation(){
    
    float theta = 3.141/4;
    
    const Image bostonim("./Input/BostonRainbow-crop-400.png");
    
    Image rot = rotate(bostonim, theta);
    rot.write("./Output/BostonRainbow-rotated.png");
    printf("Scaled-Lin image is %dx%dx%d\n", rot.width(), rot.height(), rot.channels());
    
}


// test warp by 1
void testWarpBy1(){
    
    Image bearim("./Input/bear.png");
    
    // define before and after segments
    Segment segBefore(0,0, 10,0);
    Segment segAfter(10, 10, 30, 15);
    
    Image warp1 = warpBy1(bearim, segBefore, segAfter);
    
    warp1.write("./Output/warpBy1.png");
    
}


// a function to test your morphing function with the fredo and werewolf images
void testMorph(){
    
    // load the images
    Image fredo("./Input/fredo.png");
    Image werewolf("./Input/werewolf.png");
    
    // print out the size of the two images
    printf("Fredo image is %dx%dx%d\n", fredo.width(), fredo.height(), fredo.channels());
    printf("Werewolf image is %dx%dx%d\n", werewolf.width(), werewolf.height(), werewolf.channels());
    
    
    // define the segments corresponding to fredo (segsBefore) and the werewolf (segsAfter)
    vector<Segment> segsBefore, segsAfter;
    
    segsBefore.push_back(Segment(87, 128, 109, 127));
    segsBefore.push_back(Segment(143, 127, 162, 131));
    segsBefore.push_back(Segment(96, 197, 129, 190));
    segsBefore.push_back(Segment(118, 221, 132, 200));
    segsBefore.push_back(Segment(140, 238, 165, 170));
    segsBefore.push_back(Segment(71, 242, 44, 196));
    segsBefore.push_back(Segment(9, 46, 34, 14));
    segsBefore.push_back(Segment(175, 66, 154, 27));
    
    segsAfter.push_back(Segment(83, 114, 107, 109));
    segsAfter.push_back(Segment(139, 103, 157, 101));
    segsAfter.push_back(Segment(100, 170, 132, 151));
    segsAfter.push_back(Segment(125, 198, 145, 171));
    segsAfter.push_back(Segment(142, 196, 158, 139));
    segsAfter.push_back(Segment(96, 211, 75, 180));
    segsAfter.push_back(Segment(11, 41, 33, 7));
    segsAfter.push_back(Segment(175, 56, 155, 13));
    
    
    // create an image morphing between fredo and werewolf at time t=0.5
    vector<Image> imMorph = morph(fredo, werewolf, segsBefore, segsAfter, 8);
    
    // write out images
    char buffer [50];
    for ( int n=0; n<10; n++){
        
        Image im = imMorph[n];
        sprintf (buffer, "./Output/fredo_werewolf_morph_%d.png", n);
        
        im.write(buffer);
    }
}

void morphClass() {

    // load the images
    Image a("./Input/me.png");
    Image b("./Input/you.png");
    
    // print out the size of the two images
    printf("Image 1 is %dx%dx%d\n", a.width(), a.height(), a.channels());
    printf("Image 2 is %dx%dx%d\n", b.width(), b.height(), b.channels());
     
vector<Segment> segsBefore;
segsBefore.push_back(Segment(72, 126, 106, 126)); 
segsBefore.push_back(Segment(75, 138, 104, 135)); 
segsBefore.push_back(Segment(55, 86, 79, 85)); 
segsBefore.push_back(Segment(96, 83, 122, 85)); 
segsBefore.push_back(Segment(55, 144, 90, 163)); 
segsBefore.push_back(Segment(102, 162, 126, 126)); 
segsBefore.push_back(Segment(40, 96, 44, 120)); 
segsBefore.push_back(Segment(137, 95, 132, 120)); 
segsBefore.push_back(Segment(64, 32, 100, 30)); 
segsBefore.push_back(Segment(107, 36, 128, 57)); 
segsBefore.push_back(Segment(57, 39, 40, 66)); 
segsBefore.push_back(Segment(143, 166, 182, 185)); 
segsBefore.push_back(Segment(8, 188, 43, 174)); 

vector<Segment> segsAfter;
segsAfter.push_back(Segment(84, 126, 114, 124)); 
segsAfter.push_back(Segment(86, 134, 113, 133)); 
segsAfter.push_back(Segment(64, 73, 88, 74)); 
segsAfter.push_back(Segment(103, 70, 130, 74)); 
segsAfter.push_back(Segment(61, 135, 97, 157)); 
segsAfter.push_back(Segment(114, 154, 142, 123)); 
segsAfter.push_back(Segment(45, 93, 54, 123)); 
segsAfter.push_back(Segment(148, 85, 149, 115)); 
segsAfter.push_back(Segment(77, 19, 114, 19)); 
segsAfter.push_back(Segment(123, 20, 145, 46)); 
segsAfter.push_back(Segment(64, 25, 40, 60)); 
segsAfter.push_back(Segment(148, 163, 186, 178)); 
segsAfter.push_back(Segment(17, 185, 52, 169));
    
    
    vector<Image> imMorph = morph(a, b, segsBefore, segsAfter, 15);
    
    // write out images
    char buffer [50];
    for ( int n=0; n<15; n++){
        
        Image im = imMorph[n];
        sprintf (buffer, "./Output/class_morph_%d.png", n);
        
        im.write(buffer);
    }

}


// This is a way for you to test your functions.
// We will not grade the contents of the main function
int main() {
    
    // uncomment these test functions as you complete the assignment to test your code
    
    // testSmartAccessor();
    // testScaling();
    // testRotation();
    // testWarpBy1();
    //testMorph();
    morphClass();
    
    
}
