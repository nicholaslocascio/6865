#include <a9.h>
#include <timing.h>

int main(int argc, char** argv)
{
    Image<uint8_t> input = load<uint8_t>("images/rgb.png");
    Image<uint8_t> lumi  = load<uint8_t>("images/grey.png");

    // Part 1
    {
        Image<uint8_t> out1 = SmoothGradNormalized();
        Image<uint8_t> out2 = WavyRGB();
        Image<uint8_t> out3 = Luminance(input); // CHANGE: Fixed bug here
        Image<uint8_t> out4 = Sobel(lumi);

        save(out1, "Output/Smooth.png");
        save(out2, "Output/WavyRGB.png");
        save(out3, "Output/Luminance.png");
        save(out4, "Output/Sobel.png");

        Image<uint8_t> tiny(10, 10);            // CHANGE: next three lines
        cppCovoidorBoxSchedule5(tiny);          // to run these functions
        cppCovoidorBoxSchedule6(tiny);          // on a small image
        cppCovoidorBoxSchedule7(tiny);
    }

    // Part 2
    {
        Image<uint8_t> out1, out2, out3;

        int w = lumi.width();
        int h = lumi.height();

        Func r1 = LocalMax   (lumi);
        Func r2 = Gaussian   (lumi, 3.0f);
        Func r3 = UnsharpMask(lumi, 3.0f);

        profile(r1, w, h);
        profile(r2, w, h);
        profile(r2, w, h);

        out1 = r1.realize(w, h);
        out2 = r2.realize(w, h);
        out3 = r3.realize(w, h);

        save(out1, "Output/LocalMax.png");
        save(out2, "Output/Gaussian.png");
        save(out3, "Output/UnsharpMask.png");
    }

    return EXIT_SUCCESS;
}
