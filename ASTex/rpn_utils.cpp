#include "rpn_utils.h"

void saveFourierModulusPhaseGray(const std::string &out_path, const std::string& in_texture)
{
    try
    {
        create_directory(out_path);

        ImageGrayd image;
        IO::loadu8_in_01(image, in_texture);

        ImageSpectrald fftModulus, fftPhase;

        Fourier::fftForwardModulusAndPhase(image, fftModulus, fftPhase);

        IO::save_spectrum(fftModulus, out_path+"/modulus.png");
        IO::save_phase(fftPhase, out_path+"/phase.png");
    }
    catch(itk::ExceptionObject e)
    {
        std::cerr << "Error! " << e.what() << std::endl;
    }
}

void rpn_scalar(const ImageSpectrald& modulus, ImageSpectrald& phase, ImageGrayd& output)
{
    phase.initItk(modulus.width(), modulus.height());
    Fourier::randomPhase(phase, [](int,int){return true;}, false);
    //Fourier::randomPhase_normal(phase, [](int,int){return true;}, false, 0, M_PI/2);
    Fourier::fftInverseModulusAndPhase(modulus, phase, output);
    return;
}

void compute_transition(ImageGrayd& transition, float alpha)
{
    float *bumpx, *bumpy;       /* transition part along
                                 * x-axis and y-axis */
    float *ptr_bumpx, *ptr_bumpy; /* pointers towards bumpx and bumpy */
    int margx, margy;           /* width of margins */
    int x, y;                   /* indices for for loops */
    float tmp, tmpy, tmpx;      /* temp values */

    /* Sizes for memory allocation */
    margx = std::ceil(alpha * transition.width());
    margy = std::ceil(alpha * transition.height());

    /* Memory allocation */
    bumpx = new float[margx];
    bumpy = new float[margy];

    /* Computation of the radial transition bumpx: it is the
     * cumulative histogram of the molifier */
    ptr_bumpx = bumpx;
    tmp = 0.;
    for (x = 0; x < margx; x++) {
        tmpx = 2 * x / (margx - 1) - 1; /* affine change of variable
                                         * [0,margx-1] -> [-1,1] */
        tmp += std::exp(-1. / (1. - tmpx * tmpx));
        *ptr_bumpx = tmp;
        ptr_bumpx++;
    }
    /* Normalization of the cumulative histogram: the total sum
     * must equal 1 */
    ptr_bumpx = bumpx;
    for (x = 0; x < margx; x++) {
        *ptr_bumpx /= tmp;
        ptr_bumpx++;
    }

    /* Same procedure for bumpy */
    ptr_bumpy = bumpy;
    tmp = 0.;
    for (y = 0; y < margy; y++) {
        tmpy = 2 * y / (margy - 1) - 1;
        tmp += std::exp(-1. / (1. - tmpy * tmpy));
        *ptr_bumpy = tmp;
        ptr_bumpy++;
    }
    ptr_bumpy = bumpy;
    for (y = 0; y < margy; y++) {
        *ptr_bumpy /= tmp;
        ptr_bumpy++;
    }

    /* Computation of the float image transition */
    ptr_bumpy = bumpy;
    ptr_bumpx = bumpx;

    for (y = 0; y < (int) transition.height(); y++) {
        /* Test to determine the atenuation coefficient tmpy:
         * we go back and forth on the array bumpy */
        if (y < margy) {
            tmpy = *ptr_bumpy;
            ptr_bumpy++;
        }
        else if (y >= (int) transition.height() - margy) {
            ptr_bumpy--;
            tmpy = *ptr_bumpy;
        }
        else {
            tmpy = 1.;
        }
        for (x = 0; x < (int) transition.width(); x++) {
            /* Test to determine the atenuation coefficient
             * tmpx: we go back and forth on the array bumpx */
            if (x < margx) {
                tmpx = *ptr_bumpx;
                ptr_bumpx++;
            }
            else if (x >= (int) transition.width() - margx) {
                ptr_bumpx--;
                tmpx = *ptr_bumpx;
            }
            else {
                tmpx = 1.;
            }
            /* Computation of the coefficient of transition */
            transition.pixelAbsolute(x, y) = tmpy * tmpx;
        }
    }

    /* Free memory */
    delete[] bumpx;
    delete[] bumpy;

    return;
}

void gray_RPN(const ImageGrayd& in, ImageGrayd& out, unsigned int extendX, unsigned int extendY, bool crop, bool periodic_component, bool call_srand)
{
    if(!in.is_initialized())
        return;
    try
    {
        ImageGrayd pc;
        if(periodic_component)
            image2periodicComponent(in, pc);
        else
            pc.copy_pixels(in);

        //create an extended spot
        if(extendX>0 || extendY>0)
        {
            if(!out.is_initialized() || (unsigned)out.width()!=(pc.width()+extendX) || (unsigned)out.height()!=(pc.height()+extendY))
                out.initItk(pc.width()+extendX, pc.height()+extendY);

            //compute the mean value of the periodic component
            ImageGrayd::PixelType mean=0;
            pc.for_all_pixels([&mean] (ImageGrayd::PixelType& pix)
            {
                mean += pix;
            });
            mean /= (pc.width()*pc.height());

            //fill the extended spot with the mean
            out.for_all_pixels([&mean] (ImageGrayd::PixelType& pix)
            {
                pix=mean;
            });

            //find middle of extended spot
            int m_x=out.width()/2, m_y=out.height()/2;

            //find origin (top left) of extended spot
            int o_x=m_x-in.width()/2, o_y=m_y-in.height()/2;

            //pre-compute the constant for formula 6 of GGM11
            double sqrtSizes = std::sqrt(double(out.width()*out.height()) / (pc.width()*pc.height()));

            //compute the transition function image
            ImageGrayd transition;
            transition.initItk(pc.width(), pc.height());
            compute_transition(transition, 0.02);

//            IO::save01_in_u8(transition, "transition.png");

            for(int x=0; x<in.width(); ++x)
            {
                for(int y=0; y<in.height(); ++y)
                {
                    out.pixelAbsolute(o_x+x, o_y+y) += (pc.pixelAbsolute(x, y) - mean) * sqrtSizes * transition.pixelAbsolute(x, y);
                }
            }
        }
        else
        {
            out=pc;
        }

        //compute a random phase
        ImageSpectrald randomPhase;
        randomPhase.initItk(out.width(), out.height(), true);
        if(call_srand)
            srand(time(NULL));

        Fourier::randomPhase(randomPhase, [](int,int){return true;}, false);
        //Fourier::randomPhase_normal(randomPhase, [](int,int){return true;}, false, M_PI, M_PI);

        auto addRandomPhase = [&randomPhase] (ImageSpectrald::PixelType& pix, int x, int y)
        {
            pix += randomPhase.pixelAbsolute(x, y);
            pix = pix>=2*M_PI ? pix - 2*M_PI : (pix <= -2*M_PI ? pix + 2*M_PI : pix);
        };

        ImageSpectrald modulus, phase;
        Fourier::fftForwardModulusAndPhase(out, modulus, phase);

        phase.for_all_pixels(addRandomPhase);

        Fourier::fftInverseModulusAndPhase(modulus, phase, out);
    }
    catch(itk::ExceptionObject e)
    {
        std::cerr << "Error! " << e.what() << std::endl;
    }
    return;
}

void colored_RPN(const ImageRGBd& in, ImageRGBd& out, color_space_t colorSpace, color_dephasing_mode_t mode,
                 unsigned int extendX, unsigned int extendY, bool crop, bool periodic_component, bool call_srand, double scale_randomPhase, const ImageSpectrald *phase)
{
    if(!in.is_initialized() || mode==LUMINANCE || mode==ALL)
        return;

    try
    {
        ImageRGBd pc;
        pc.initItk(in.width(), in.height());
        if(periodic_component)
            image2periodicComponent(in, pc);
        else
            pc.copy_pixels(in);

        //create an extended spot
        if(extendX>0 || extendY>0)
        {
            if(!out.is_initialized() || (unsigned)out.width()!=(pc.width()+extendX) || (unsigned)out.height()!=(pc.height()+extendY))
                out.initItk(pc.width()+extendX, pc.height()+extendY);

            //compute the mean value of the periodic component
            ImageRGBd::PixelType mean;
            for(int i=0; i<3; ++i) mean[i]=0;
            pc.for_all_pixels([&mean] (ImageRGBd::PixelType& pix)
            {
                mean += pix;
            });
            for(int i=0; i<3; ++i) mean[i] /= (pc.width()*pc.height());

            //fill the extended spot with the mean
            out.for_all_pixels([&mean] (ImageRGBd::PixelType& pix)
            {
                pix=mean;
            });

            //find middle of extended spot
            int m_x=out.width()/2, m_y=out.height()/2;

            //find origin (top left) of extended spot
            int o_x=m_x-in.width()/2, o_y=m_y-in.height()/2;

            //pre-compute the constant for formula 6 of GGM11
            double sqrtSizes = std::sqrt(double(out.width()*out.height()) / (pc.width()*pc.height()));

            //compute the transition function image
            ImageGrayd transition;
            transition.initItk(pc.width(), pc.height());
            compute_transition(transition, 0.5);

//            IO::save01_in_u8(transition, "transition.png");

            for(int x=0; x<pc.width(); ++x)
            {
                for(int y=0; y<pc.height(); ++y)
                {
                    for(int i=0; i<3; ++i) out.pixelAbsolute(o_x+x, o_y+y)[i] += (pc.pixelAbsolute(x, y)[i] - mean[i]) * sqrtSizes * transition.pixelAbsolute(x, y);
                }
            }
        }
        else
        {
            out=pc;
        }

        //compute a random phase
        ImageSpectrald randomPhase;
        randomPhase.initItk(out.width(), out.height(), true);
        if(call_srand)
            srand(time(NULL));

        if(phase!=NULL)
            randomPhase.copy_pixels(*phase);
        else
            Fourier::randomPhase(randomPhase, [](int,int){return true;}, false);

        if(scale_randomPhase!=1.0)
            randomPhase.for_all_pixels([&] (ImageSpectrald::PixelType &pix)
            {
               pix*=scale_randomPhase;
            });

        //Fourier::randomPhase_normal(randomPhase, [](int,int){return true;}, false, 0, M_PI_2/2);
        //extract
        ImageGrayd redChannel, greenChannel, blueChannel;

        if(colorSpace==LUV_SPACE)
            extractLuv(out, redChannel, greenChannel, blueChannel);
        else if(colorSpace==LAB_SPACE)
            extractLab(out, redChannel, greenChannel, blueChannel);
        else
            extract3Channels(out, redChannel, greenChannel, blueChannel);

        if(mode==PCA_SYNTHESIS)
            fold3Channels(out, redChannel, greenChannel, blueChannel);
        PCA pca(out);

        if(mode==PCA_SYNTHESIS)
        {
            //as long as the distribution is a real gaussian, this works.

            MaskBool mb(out.width(), out.height());
            mb |= [] (int, int) {return true;};

            pca.computePCA(mb);

            //channels are projected into principal components
            pca.project(redChannel, greenChannel, blueChannel);
        }

        ImageSpectrald redFftInput, greenFftInput, blueFftInput;
        ImageSpectrald redPhase, greenPhase, bluePhase;

        Fourier::fftForwardModulusAndPhase(redChannel, redFftInput, redPhase);
        Fourier::fftForwardModulusAndPhase(greenChannel, greenFftInput, greenPhase);
        Fourier::fftForwardModulusAndPhase(blueChannel, blueFftInput, bluePhase);

        if(mode==ONE_PHASE)
        {
            //this will compute the wrong RPN, keeping the normal distributions for each channel, but modifying the covariance.
            //this actually works if in is monochromatic, that is, if its PCA has only one major component.

            Fourier::fftInverseModulusAndPhase(redFftInput, randomPhase, redChannel);
            Fourier::fftInverseModulusAndPhase(greenFftInput, randomPhase, greenChannel);
            Fourier::fftInverseModulusAndPhase(blueFftInput, randomPhase, blueChannel);
        }
        else if(mode==SAME_PHASE)
        {
            //this is a debug mode. It tests the fourier transform and should return the same image as in.

            Fourier::fftInverseModulusAndPhase(redFftInput, redPhase, redChannel);
            Fourier::fftInverseModulusAndPhase(greenFftInput, greenPhase, greenChannel);
            Fourier::fftInverseModulusAndPhase(blueFftInput, bluePhase, blueChannel);
        }
        else if(mode==NULL_PHASE)
        {
            //This will compute a false texton, with wrong colors, unless the phases are highly correlated.

            const Offset off = randomPhase.getCenter();
            randomPhase.parallel_for_all_pixels([&off] (ImageSpectrald::PixelType &pix, int x, int y)
            {
                if(x!=off[0] && y!=off[1])
                    pix = 0;
            });

            Fourier::fftInverseModulusAndPhase(redFftInput, randomPhase, redChannel);
            Fourier::fftInverseModulusAndPhase(greenFftInput, randomPhase, greenChannel);
            Fourier::fftInverseModulusAndPhase(blueFftInput, randomPhase, blueChannel);
        }
        else if(mode==TEXTON)
        {
            //This will compute the good texton.

            const Offset off = randomPhase.getCenter();
            auto textonPhase = [&off, &redPhase] (ImageSpectrald::PixelType &pix, int x, int y)
            {
                if(x!=off[0] && y!=off[1])
                    pix -= redPhase.pixelAbsolute(x, y);
            };

            greenPhase.parallel_for_all_pixels(textonPhase);
            bluePhase.parallel_for_all_pixels(textonPhase);
            //negate redPhase
            redPhase.parallel_for_all_pixels(textonPhase);

            Fourier::fftInverseModulusAndPhase(redFftInput, redPhase, redChannel);
            Fourier::fftInverseModulusAndPhase(greenFftInput, greenPhase, greenChannel);
            Fourier::fftInverseModulusAndPhase(blueFftInput, bluePhase, blueChannel);

            ImageGrayd *ptr_channel;
            auto shiftSpatialX = [&ptr_channel, &off] (ImageGrayd::PixelType &pix, int x, int y)
            {
                if(x<off[0])
                {
                    ImageGrayd::PixelType tmpPix=pix;
                    pix=ptr_channel->pixelAbsolute(x+off[0], y);
                    ptr_channel->pixelAbsolute(x+off[0], y)=tmpPix;
                }
            };

            auto shiftSpatialY = [&ptr_channel, &off] (ImageGrayd::PixelType &pix, int x, int y)
            {
                if(y<off[1])
                {
                    ImageGrayd::PixelType tmpPix=pix;
                    pix=ptr_channel->pixelAbsolute(x, y+off[1]);
                    ptr_channel->pixelAbsolute(x, y+off[1])=tmpPix;
                }
            };

            ptr_channel=&redChannel;
            ptr_channel->for_all_pixels(shiftSpatialX);
            ptr_channel->for_all_pixels(shiftSpatialY);
            ptr_channel=&greenChannel;
            ptr_channel->for_all_pixels(shiftSpatialX);
            ptr_channel->for_all_pixels(shiftSpatialY);
            ptr_channel=&blueChannel;
            ptr_channel->for_all_pixels(shiftSpatialX);
            ptr_channel->for_all_pixels(shiftSpatialY);
        }
        else if(mode==SEVERAL_PHASES || mode==PCA_SYNTHESIS)
        {
            //Wrong colors, but normal distributions for each channel.
            //It can produce acceptable colors if the covariance of each channel is close to zero and the texture in big enough,
            //however you shouldn't rely on that for most images.
            //Perfectly fine within a PCA though.

            rpn_scalar(redFftInput, randomPhase, redChannel);
            rpn_scalar(greenFftInput, randomPhase, greenChannel);
            rpn_scalar(blueFftInput, randomPhase, blueChannel);
        }
        else //mode=normal or pca synthesis
        {
            //good colors, and a multivariate normal distribution of intensities that fits in's means and covariance matrix.
            //considered the best method yet, but it only yields good colors if the input texture is Gaussian.
            //the input is Gaussian if its phase is a uniform random phase. This can be approximated by having a multivariate normal distribution
            //as the histogram of intensities, representing the first order statistics of the input exemplar.

            auto addRandomPhase = [&] (ImageSpectrald::PixelType& pix, int x, int y)
            {
                pix += randomPhase.pixelAbsolute(x, y);
                if(pix>2*M_PI)
                    pix -= 2*M_PI;
                else if(pix<-2*M_PI)
                    pix += 2*M_PI;
            };

            redPhase.for_all_pixels(addRandomPhase);
            greenPhase.for_all_pixels(addRandomPhase);
            bluePhase.for_all_pixels(addRandomPhase);

            Fourier::fftInverseModulusAndPhase(redFftInput, redPhase, redChannel);
            Fourier::fftInverseModulusAndPhase(greenFftInput, greenPhase, greenChannel);
            Fourier::fftInverseModulusAndPhase(blueFftInput, bluePhase, blueChannel);
        }

        if(mode==PCA_SYNTHESIS)
            pca.back_project(redChannel, greenChannel, blueChannel, out);
        else
            fold3Channels(out, redChannel, greenChannel, blueChannel);

        if(colorSpace==LUV_SPACE)
            foldRGBfromLuv(out, out);
        else if(colorSpace==LAB_SPACE)
            foldRGBfromLab(out, out);

        if(crop)
        {   //that simply means cutting pixels off to their closest boundary.
            out.for_all_pixels([] (ImageRGBd::PixelType& pix)
            {
                pix[0]=std::max(0.0, std::min(pix[0], 1.0));
                pix[1]=std::max(0.0, std::min(pix[1], 1.0));
                pix[2]=std::max(0.0, std::min(pix[2], 1.0));
            });
        }
    }
    catch(itk::ExceptionObject e)
    {
        std::cerr << "Error! " << e.what() << std::endl;
    }

    return;
}

double compute_crossCorrelation_diff(const ImageRGBd& in1, const ImageRGBd& in2, int channel1, int channel2, ImageGrayd& diff)
{
    assert(in1.width() == in2.width() && in1.height() == in2.height());

    if(!diff.is_initialized() || diff.width() != in1.width() || diff.height() != in1.height())
        diff.initItk(in1.width(), in1.height());

    ImageGrayd crossCorrelation1, crossCorrelation2;
    crossCorrelation1.initItk(in1.width(), in1.height());
    crossCorrelation2.initItk(in1.width(), in1.height());

    MaskBool mb(in1.width(), in1.height());
    mb |= [] (int, int) {return true;};

    Fourier::crossCorrelation_full_size(in1, crossCorrelation1, mb, channel1, channel2);
    Fourier::crossCorrelation_full_size(in2, crossCorrelation2, mb, channel1, channel2);

    double diffsum=0;

    diff.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
    {
        pix = std::abs(crossCorrelation1.pixelAbsolute(x, y) - crossCorrelation2.pixelAbsolute(x, y));
        diffsum += pix;
    });

    return diffsum;
}
