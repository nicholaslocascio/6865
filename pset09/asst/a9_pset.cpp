/**
 * Instructions.
 *
 * This file contains a number of empty function stubs with the problem
 * statements. You need to complete fill in all details. These problems
 * are designed as extensions to the tutorials in this Pset. Please run
 * each tutorial and go through the code and comments. You are welcome
 * to discuss details on Piazza if you do not understand the tutorials.
 *
 * The main function is in a9_main.cpp. You can modify the contents for
 * your debugging purposes but restore the original contents.
 *
 * Live long and propser.
 */

#include <a9.h>
#include <vector>

// This applied a compute_root() schedule to all the Func's that are consumed by
// the calling Func
void apply_auto_schedule(Func F) {
    map<string,Internal::Function> flist = Internal::find_transitive_calls(F.function());
    flist.insert(std::make_pair(F.name(), F.function()));
    map<string,Internal::Function>::iterator fit;
    for (fit=flist.begin(); fit!=flist.end(); fit++) {
        Func f(fit->second);
        f.compute_root();
        cout << "Warning: applying default schedule to " << f.name() << endl;
    }
    cout << endl;
}

// ---------------------------------- PART 1 -----------------------------------
// For all Halide Func's you write in Part 1, you do not need to worry about
// schedule. However, things will be faster if you just use:
// apply_auto_schedule(final_Result_Func);


// Use Halide to compute a 512x512 smooth gradient where pixel (x,y) is
// equal to x+y divided by 1024.
// Hint: tutorial 1 has a very similar function except that it does not
// normalize by 1024
// Schedule: use compute_root() on all the Func's you create.
// \return a single channel output image of size 512x512
Image<uint8_t> SmoothGradNormalized(void) {
    Func gradient;

    Var x("x");
    Var y("y");

    Expr e = 255.0f*cast<float>(x + y)/1024.0f;
    e = clamp(e, 0.0f, 255.0f);
    e = cast<uint8_t>(e);
    gradient(x, y) = e;

    Image<uint8_t> output = gradient.realize(512, 512);

    return output;
}

// Use a Halide Func to compute a wavy RGB image like that obtained
// by the formula: output(x,y,c) = (1-c)*cos(x)*cos(y)
// Schedule: use compute_root() on all the Func's you create.
// Hint: you need one more dimension than above
// \return 512x512 RGB image
Image<uint8_t> WavyRGB(void) {
    Func hal_func;

    Var x("x");
    Var y("y");
    Var c("c");

    Expr e = 255.0f*(1.0f-c)*cos(x)*cos(y);
    e = clamp(e, 0.0f, 255.0f);
    e = cast<uint8_t>(e);
    hal_func(x, y, c) = e;

    Image<uint8_t> output = hal_func.realize(512, 512, 3);

    return output;
}

// Use Halide to compute a 1-channel luminance image representing 0.3R+0.6G+0.1B
// Input: RGB image
// Return: single channel luminance image the same size as the input image
Image<uint8_t> Luminance(Image<uint8_t> input) {
    Func hal_func;
    Var x("x");
    Var y("y");

    Expr e = input(x, y, 0)*0.3f + input(x, y, 1)*0.6f + input(x, y, 2)*0.1f;
    e = clamp(e, 0.0f, 255.0f);
    e = cast<uint8_t>(e);

    hal_func(x, y) = e;

    Image<uint8_t> output = hal_func.realize(input.width(), input.height());

    return output;
}

// Use Halide to apply a Sobel filter and return the gradient magnitude.
// Use the following kernels to compute the horizontal and vertical gradients
//
//            | -1 0 1 |             | -1 -2 -1 |
// x kernel = | -2 0 2 |, y kernel = |  0  0  0 |
//            | -1 0 1 |             |  1  2  1 |
//
// Input: single channel luminance image
// Return: single channel image where each pixel is the magnitude of Sobel x and y filters
Image<uint8_t> Sobel(Image<uint8_t> input) {
    Func hal_func;

    Func grad_x("grad_x");
    Func grad_y("grad_y");

    Var x("x");
    Var y("y");

    Func clamped("clamped");
    clamped(x, y) = cast<float>(
        input(clamp(x, 0, input.width()-1), clamp(y, 0, input.height()-1)));

    grad_x(x, y) = (-1.0f*clamped(x-1, y-1) + 0.0f*clamped(x, y-1) + 1.0f*clamped(x+1, y-1) +
                        -2.0f*clamped(x-1, y) + 0.0f*clamped(x, y) + 2.0f*clamped(x+1, y) +
                        -1.0f*clamped(x-1, y+1) + 0.0f*clamped(x, y+1) + 1.0f*clamped(x+1, y+1))/4.0f;

    grad_y(x, y) = (-1.0f*clamped(x-1, y-1) + -2.0f*clamped(x, y-1) + -1.0f*clamped(x+1, y-1) +
                        0.0f*clamped(x-1, y) + 0.0f*clamped(x, y) + 0.0f*clamped(x+1, y) +
                        1.0f*clamped(x-1, y+1) + 2.0f*clamped(x, y+1) + 1.0f*clamped(x+1, y+1))/4.0f;


    Expr e = 255.0f*sqrt(pow(grad_x(x, y), 2.0f) + pow(grad_y(x, y), 2.0f));
    e = clamp(e, 0.0f, 255.0f);
    e = cast<uint8_t>(e);

    hal_func(x, y) = e;

    Image<uint8_t> output = hal_func.realize(input.width(), input.height());

    return output;
}

// Write the cpp nested loops corresponding to the 3x3 box schedule 5 in tutorial 6
// print the order of evaluation. Each time you perform a computation of blur_x or
// blur_y, put a triplet with the name of the function "blur_x" or "blur_y") and
// the output coordinates x and y. For example, when you compute (50,75) print:
// blur_x 50 75
// in a new line.
//
// Use cerr to print instead or cout. No empty lines in your output.
//
// Repeat the above procedute for Schedule 6 of tutorial 6
//
// Schedule 5:
//   blur_y.compute_root()
//   blur_x.compute_at(blur_y, x)
//
void cppCovoidorBoxSchedule5(Image<uint8_t> input) {
    cerr << "Started cppCovoidorBoxSchedule5" << endl;

    for (int y = 0; y < input.height()-2; y++) {
        for (int x = 0; x < input.width()-2; x++) {
            cerr << "blur_x " << x << " " <<  y << endl;
        }
        for (int x = 0; x < input.width()-2; x++) {
            cerr << "blur_y " << x << " " <<  y << endl;
        }
    }

    cerr << "Completed cppCovoidorBoxSchedule5" << endl;
}

// Write the cpp nested loops corresponding to the 3x3 box schedule 5 in tutorial 6
// print the order of evaluation. Each time you compute a pixel {x,y} of blur_x or
// blur_y, print the name of the function, x and y. For example, when you compute
// (50,75) print the following on a new line.
// blur_x 50 75
//
// Use cerr to print instead or cout.
//
// Repeat the above procedute for Schedule 6 of tutorial 6
//
// Schedule 6
//   blur_y.tile(x, y, xo, yo, xi, yi, 2, 2)
//   blur_x.compute_at(blur_y, yo)
//
void cppCovoidorBoxSchedule6(Image<uint8_t> input) {
    cerr << "Started   cppCovoidorBoxSchedule6" << endl;

    int max_x = input.width()-3;
    int max_y = input.height()-3;

    for (int x = 0; x < (input.width()-1)/2; x++) {
        for (int y = 0; y < input.height()-1; y++) {
            for (int b = 0; b < 2; b++) {
                string function_name = "blur_x";
                if (b > 0) {
                    function_name = "blur_y";
                }

                int x0 = fmin(x*2, max_x);
                int x1 = fmin(x*2+1, max_x);
                int x2 = fmin(x*2+2, max_x);
                int x3 = fmin(x*2+3, max_x);

                int y0 = fmin(y*2, max_y);
                int y1 = fmin(y*2+1, max_y);
                int y2 = fmin(y*2+2, max_y);
                int y3 = fmin(y*2+3, max_y);

                cerr << function_name << " " << x0 << " " << y0  << endl;
                cerr << function_name << " " << x0 << " " <<  y1 << endl;
                cerr << function_name << " " << x0 << " " <<  y2 << endl;
                cerr << function_name << " " << x0 << " " <<  y3 << endl;

                cerr << function_name << " " << x0 << " " << y0  << endl;
                cerr << function_name << " " << x1 << " " <<  y1 << endl;
                cerr << function_name << " " << x2 << " " <<  y2 << endl;
                cerr << function_name << " " << x3 << " " <<  y3 << endl;
            }
        }
    }

    cerr << "Completed cppCovoidorBoxSchedule6" << endl;
}

// Repeat the above procedute for Schedule 7 of tutorial 6
//
// Schedule 7
//   blur_y.split(x, xo, xi, 2)
//   blur_x.compute_at(blur_y, y)
void cppCovoidorBoxSchedule7(Image<uint8_t> input) {
    cerr << "Started   cppCovoidorBoxSchedule7" << endl;

    int max_x = input.width()-3;

    for (int x = 0; x < input.width()-2; x++) {
        for (int y = 0; y < input.height()-2; y++) {
            cerr << "blur_x " << x << " " <<  y << endl;
        }
    }

    for (int y = 0; y < input.height()-2; y++) {
        for (int x = 0; x < (input.width()-1)/2; x++) {
            for (int x2 = 0; x2 < (input.width()-1)/2; x2++) {
                int xx = fmin(x*2 + x2, max_x);
                cerr << "blur_y " << xx << " " <<  y << endl;
            }
        }
    }

    cerr << "Completed cppCovoidorBoxSchedule7" << endl;
}


// ---------------------------------- PART 2 -----------------------------------
// For all the problems in Part 2, yo are expected to write schedules. You will
// score points for faster schedules. You should also report the machine used
// for testing in the form accompanying this Pset. Each problem is supposed to
// return a Halide Func that computes the result in uint8_t. All intermediate
// operations should be performed in float for Gaussian and unsharp mask.
//
// First get the algorithm right and make sure it is generating the
// correct output. Use apply_auto_schedule() during this phase. Once you
// are confident of your result, remove apply_auto_schedule() and write
// your own schedule.
//
// Hint: your schedule should beat apply_auto_schedule()


// Compute the local maximum in a window of given size centered at the pixel.
// Return white color if it's a local maximum and black otherwise.
//
// Handle boundary conditions by clamping at boundaries (see tutorials)
//
// Input: single channel luminance image
// Return: Halide Func that computes the localMax
//
Func LocalMax(Image<uint8_t> lumi, int window) {
    Func LM("LM");
    Func MAX_VAL("MAX_VAL");

    Var x("x");
    Var y("y");

    Func clamped("clamped");
    clamped(x, y) = cast<float>(
        lumi(clamp(x, 0, lumi.width()-window/2), clamp(y, 0, lumi.height()-window/2)));

    RDom r_x(-window, window*2, -window, window*2, "r_x");

    MAX_VAL(x, y) += cast<uint8_t>(max(clamped(x,y), clamped(x + r_x.x, y + r_x.y)));

    LM(x, y) = cast<uint8_t>(select(clamped(x,y) == MAX_VAL(x, y), 255, 0));

    apply_auto_schedule(LM);

    return LM;
}

float gauss1DFilterTotalWeight(float sigma, float truncate){

    int radius = ceil(sigma * truncate);
    float c = sigma;
    float totalWeight = 1.0;

    for(int i = 1; i < radius; i++) {
        cout << i << endl;
        float c2_2 = 2.0*pow(c, 2.0);
        float value = exp(-1.0*pow(i, 2.0)/c2_2);
        totalWeight += 2*value;
   }
    
    return totalWeight;
}


// Perform 2D gaussian blur. use what you learnt in previous Pset about
// fast Gaussian blur algorithm. The Gaussians have variance sigma and the
// kernels are truncated at trunc*sigma on both sides
//
// Use the sum() sugar (see tutorials)
// First generate the kernel as a Halide Func
// Handle boundary conditions
//
// Input: single channel luminance image
// Return: Func that computes Gaussian blur
//
Func Gaussian(Image<uint8_t> input, float sigma, float truncate) {
    Func GB("Gaussian");
    Func blur_x("blur_x");
    Func blur_y("blur_y");
    Func gauss_kernal("gauss_kernal");

    Var x("x");
    Var y("y");

    int radius = ceil(sigma*truncate);
    int length = 2*radius + 1;

    Func clamped("clamped");
    clamped(x, y) = cast<float>(
        input(clamp(x, 0, input.width()-1), clamp(y, 0, input.height()-1)));

    float c2_2 = 2.0*pow(sigma, 2.0);
    gauss_kernal(x) = exp(-1.0f*pow(x, 2)/c2_2);

    RDom r_x(0, length, "r_x");
    float kernal_sum = gauss1DFilterTotalWeight(sigma, truncate);

    blur_x(x, y) = sum(gauss_kernal(r_x.x-radius)/kernal_sum*clamped(x+r_x.x-radius, y));
    blur_y(x, y) = sum(gauss_kernal(r_x.x-radius)/kernal_sum*blur_x(x, y+r_x.x-radius));

    GB(x, y) = cast<uint8_t>(clamp(blur_y(x, y), 0.0f, 255.0f));

    Var x_outer("x_outer");
    Var y_outer("y_outer");
    Var x_inner("x_inner");
    Var y_inner("y_inner");
    Var tile_index("tile_index");

    GB.tile(x, y, x_outer, y_outer, x_inner, y_inner, 2, 2)
            .fuse(x_outer, y_outer, tile_index)
            .parallel(tile_index);

    apply_auto_schedule(GB);

    profile(GB, input.width(), input.height(), input.channels());

    return GB;
}

// Implement Unsharp Mask. Use a Gaussian of standard deviation sigma to extract
// a lowpassed version of the image. Subtract that lowpassed version from the
// original image to produce a highpassed version of the image and then add the
// highpassed version back to it strength times.
//
// Handle boundary conditions
// Use the Gaussian above for blur
// You may want to write a new Gaussian with a new schedule because it may be
// possible to use a different schedule for Gaussian blur such that it works
// better with subsequent stages
//
// Input: single channel luminance image
// Return: Halide Func that computes the localMax
//
Func UnsharpMask(Image<uint8_t> input, float sigma, float truncate, float strength) {
    Func USM("USM");
    Func highpassed("highpassed");

    Image<uint8_t> lumi = input;

    Var x("x");
    Var y("y");

    Func gaussian = Gaussian(input, sigma, truncate);

    Func clamped("clamped");
    clamped(x, y) = cast<float>(
        input(clamp(x, 0, input.width()-1), clamp(y, 0, input.height()-1)));

    highpassed(x, y) = clamped(x, y) - gaussian(x, y);

    USM(x, y) = cast<uint8_t>(clamp(clamped(x, y) + strength*highpassed(x, y), 0.0f, 255.0f));

    apply_auto_schedule(USM);

    return USM;
}
