#ifndef _A9_H_
#define _A9_H_

#include <iostream>
#include <cmath>
#include <Halide.h>
#include <image_io.h>

using std::string;
using std::map;
using std::cout;
using std::cerr;
using std::endl;

void apply_auto_schedule(Func F);
Image<uint8_t> SmoothGradNormalized(void);
Image<uint8_t> WavyRGB(void);
Image<uint8_t> Luminance(Image<uint8_t> input);
Image<uint8_t> Sobel(Image<uint8_t> input);
void cppCovoidorBoxSchedule5(Image<uint8_t> input);
void cppCovoidorBoxSchedule6(Image<uint8_t> input);
void cppCovoidorBoxSchedule7(Image<uint8_t> input);
Func LocalMax(Image<uint8_t> lumi, int window=5);
Func Gaussian(Image<uint8_t> input, float sigma, float truncate=3.0f);
Func UnsharpMask(Image<uint8_t> input, float sigma, float truncate=3.0f, float strength=1.0f);

#endif // _A9_H_
