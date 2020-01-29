#include "gaussian_transfer.h"

namespace ASTex {


float erfinv(float x) //ripped off from stackoverflow. I do not own this. Seems to have a maximum error of about 1.0e-3.
{
   float tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = std::log(x);

   tt1 = 2/(M_PI*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*std::sqrt(-tt1 + std::sqrt(tt1*tt1 - tt2)));
}

GaussianTransfer::GaussianTransfer() : m_boundaries(), m_pixels()
{}

//reminder : bins must be sorted from smaller to bigger
void GaussianTransfer::transform(ImageGrayd& image)
{
    HistogramGrayd h(image);

    m_boundaries.resize(h.binsNumber());
    m_pixels.resize(h.binsNumber());

    //mean and variance could literally be anything and is even handier if they are 0 and 1,
    //but let's stick to the image's mean and variance to test it
    double mean=h.mean();
    double variance=h.variance();
    double totalFrequency = 0;
    int i=0;
    double x2;

    for(const auto& bin : h)
    {
        //reminder: bin.first represents the bin value, and bin.second the occurences of this bin.

        int occ = bin.second;
        //if(occ==0) //this would be a major inconvenience : useless boundaries and erfinv(+-1) is (supposed to be) +-infinite.
        if(i!=h.binsNumber()-1)
        {
            //compute the frequency of this bin.
            double frequency = double(occ);
            totalFrequency += frequency;
            //compute the equation which gives the value x2 such that for a variable X~N(mean, variance), P(x1 < X < x2)=frequency. x1 corresponds to the variable previousBoundary.
            //Thus, we know the numbers x1 and x2 such that a normal distribution is likely to produce a value that occurs as many times as bin.first in the histogram !
            x2 = std::sqrt(2) * std::sqrt(variance) * erfinv(2*totalFrequency/h.size() - 1) + mean;
            m_boundaries[i]=x2;
            m_pixels[i++]=bin.first;
        }
        else
        {
            m_boundaries[i]=x2+0.02;
            m_pixels[i]=bin.first;
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    //both boundaries and pixels are sorted, making the search easier.
    image.for_all_pixels([&] (ImageGrayd::PixelType &pix)
    {
        std::vector<ImageGrayd::PixelType>::const_iterator c_it=std::lower_bound(m_pixels.begin(), m_pixels.end(), pix);
        int index=c_it - m_pixels.begin(); //index is necessarily inside both m_pixels and m_boundaries (by semantic), but index-1 might not be
        double upper_bound = m_boundaries[index];
        double lower_bound = index > 0 ? m_boundaries[index-1] : m_boundaries[index]-0.02; //what can I say. An intern knows where their limits are

        std::uniform_real_distribution<> u(lower_bound, upper_bound);
        pix=u(gen);
    });

    m_boundaries[m_boundaries.size()-1]=std::numeric_limits<double>::infinity();

    return;
}

void GaussianTransfer::transform_inv(ImageGrayd::PixelType& pix) const
{
        std::vector<double>::const_iterator c_it=std::upper_bound(m_boundaries.begin(), m_boundaries.end(), (double)pix);
        int index=c_it - m_boundaries.begin();
        pix = m_pixels[index];
}

void GaussianTransfer::transform_inv(ImageGrayd& image) const
{
    image.for_all_pixels([&] (ImageGrayd::PixelType &pix)
    {
        transform_inv(pix);
    });
}

GaussianTransfer_Color_Recursive::GaussianTransfer_Color_Recursive() : m_BnPAssociation()
{}

void GaussianTransfer_Color_Recursive::transform(ImageRGBd &image)
{
    //initializations

    HistogramRGBd h(image);

    Eigen::Vector3d mean;
    for(int i =0; i < 3 ;i++)
        mean(i) = h.mean(i);
    Eigen::Matrix3d cov;
    for(int i = 0 ; i < 3 ; i++)
        for(int j = 0; j < 3; j++)
            cov(i,j) = h.covariance(i,j);

    ImageRGBd::PixelType old_x;
    for(int i=0; i<3; ++i)
        old_x[i] = std::numeric_limits<double>::min();

    //intermediate lambdas

    const auto push_blueLink = [&] (GreenLink &greenLink, const ImageRGBd::PixelType &pix)
    {
        greenLink.blueList.push_back(BlueLink());
        BlueLink &blueLink=greenLink.blueList.back();
        blueLink.pixel=pix[2];
        old_x[2]=pix[2];
    };

    const auto push_greenLink = [&] (RedLink &redLink, const ImageRGBd::PixelType &pix)
    {
        redLink.greenList.push_back(GreenLink());
        GreenLink &greenLink=redLink.greenList.back();
        greenLink.pixel=pix[1];
        old_x[1]=pix[1];
        push_blueLink(greenLink, pix);
    };

    const auto push_redLink = [&] (const ImageRGBd::PixelType &pix)
    {
        m_BnPAssociation.push_back(RedLink());
        RedLink &redLink=m_BnPAssociation.back();
        redLink.pixel=pix[0]; //pixel is known, boundary is unknown
        old_x[0]=pix[0];
        push_greenLink(redLink, pix);
    };

    //main iteration

    int occ[3] = {};

    for(const auto& bin : h)
    {
        ImageRGBd::PixelType pix=bin.first; //pixel

        //computing

        if(m_BnPAssociation.size()>0)
        {
            if(old_x[0]!=pix[0]) //check if this is a new red bin
            {
                //push a fresh red link, green link, and blue link
                occ[1]=0;
                occ[2]=0;

                push_redLink(pix);
            }
            else if(old_x[1]!=pix[1]) //not a new red, but a new green
            {
                //push a fresh green link and blue link
                occ[2]=0;

                push_greenLink(m_BnPAssociation.back(), pix);
            }
            else //has to be a new blue
            {
                //push a fresh blue link
                push_blueLink(m_BnPAssociation.back().greenList.back(), pix);
            }
        }
        else
        {
            push_redLink(pix);
        }

        for(int i=0; i<3; ++i)
            occ[i]+=bin.second; //number of times pixel appear

        //update occurences for any last link
        RedLink &redLink=m_BnPAssociation.back();
        redLink.boundary=(double)occ[0];

        GreenLink &greenLink=redLink.greenList.back();
        greenLink.boundary=(double)occ[1];

        BlueLink &blueLink=greenLink.blueList.back();
        blueLink.boundary=(double)occ[2];
    }

    m_mean = mean;
    for(int i=0; i<3; ++i)
        m_variance(i) =cov(i,i);

//    m_mean = glm::vec3(0.5, 0.5, 0.5);
//    for(int i=0; i<3; ++i)
//    {
//        m_variance[i]=0.02;
//    }


    closeList();

    std::random_device rd;
    std::mt19937 gen(rd());

    image.for_all_pixels([&] (ImageRGBd::PixelType &pix)
    {
       turn_pixel_into_a_more_gaussian_pixel(pix, gen);
    });

    return;
}

void GaussianTransfer_Color_Recursive::closeList()
{
    double red_occ=0;

    RedLink &last_red_link = m_BnPAssociation.back();

    for(RedLink &redLink : m_BnPAssociation)
    {
        red_occ = redLink.boundary;

        double green_occ=0;

        GreenLink &last_green_link = redLink.greenList.back();

        for(GreenLink &greenLink : redLink.greenList)
        {
            green_occ = greenLink.boundary;

            double blue_occ=0;

            BlueLink &last_blue_link = greenLink.blueList.back();

            for(BlueLink &blueLink : greenLink.blueList)
            {
                blue_occ = blueLink.boundary;

                if(&blueLink == &last_blue_link)
                    blueLink.boundary = std::numeric_limits<double>::infinity();
                else //last link boundary is the total number of occ for this list
                    blueLink.boundary = std::sqrt(2) * std::sqrt(m_variance[2]) * erfinv(2*blue_occ/last_blue_link.boundary - 1) + m_mean[2];
            }
            if(&greenLink == &last_green_link)
                greenLink.boundary = std::numeric_limits<double>::infinity();
            else
                greenLink.boundary = std::sqrt(2) * std::sqrt(m_variance[1]) * erfinv(2*green_occ/last_green_link.boundary - 1) + m_mean[1];
        }
        if(&redLink == &last_red_link)
            redLink.boundary = std::numeric_limits<double>::infinity();
        else
            redLink.boundary = std::sqrt(2) * std::sqrt(m_variance[0]) * erfinv(2*red_occ/last_red_link.boundary - 1) + m_mean[0];
    }
}

void GaussianTransfer_Color_Recursive::turn_pixel_into_a_more_gaussian_pixel(ImageRGBd::PixelType &pix, std::mt19937 &gen) const
{
    RedLink pRed;
    pRed.pixel = pix[0];

    GreenLink pGreen;
    pGreen.pixel = pix[1];

    BlueLink pBlue;
    pBlue.pixel = pix[2];

    std::list<RedLink>::const_iterator c_it_red = std::lower_bound(m_BnPAssociation.begin(), m_BnPAssociation.end(), pRed, Compare_Red(false));
    const RedLink &redLink = (*c_it_red);

    std::list<GreenLink>::const_iterator c_it_green = std::lower_bound(redLink.greenList.begin(), redLink.greenList.end(), pGreen, Compare_Green(false));
    const GreenLink &greenLink = (*c_it_green);

    std::list<BlueLink>::const_iterator c_it_blue = std::lower_bound(greenLink.blueList.begin(), greenLink.blueList.end(), pBlue, Compare_Blue(false));
    const BlueLink &blueLink = (*c_it_blue);

    //now get them boundaries and throw them pixels yarr

    if(redLink.boundary != std::numeric_limits<double>::infinity())
    {
        if(c_it_red!=m_BnPAssociation.begin())
        {
            const RedLink &redLink_prev = *(--c_it_red);

            std::uniform_real_distribution<> u(redLink_prev.boundary, redLink.boundary);

            pix[0]=u(gen);
        }
        else
        {
            std::uniform_real_distribution<> u(redLink.boundary-0.1*m_variance[0], redLink.boundary);

            pix[0]=u(gen);
        }
    }
    else
    {
        if(c_it_red!=m_BnPAssociation.begin())
        {
            const RedLink &redLink_prev = *(--c_it_red);

            std::uniform_real_distribution<> u(redLink_prev.boundary, redLink_prev.boundary+0.1*m_variance[0]);

            pix[0]=u(gen);
        }
        else
        {
//            std::uniform_real_distribution<> u(m_mean[0]-0.1*m_variance[0], m_mean[0]+0.1*m_variance[0]);
            std::normal_distribution<> n(m_mean[0], m_variance[0]);

            pix[0]=n(gen);
        }
    }

    //green, 1

    if(greenLink.boundary != std::numeric_limits<double>::infinity())
    {
        if(c_it_green!=redLink.greenList.begin())
        {
            const GreenLink &greenLink_prev = *(--c_it_green);

            std::uniform_real_distribution<> u(greenLink_prev.boundary, greenLink.boundary);

            pix[1]=u(gen);
        }
        else
        {
            std::uniform_real_distribution<> u(greenLink.boundary-0.1*m_variance[1], greenLink.boundary);

            pix[1]=u(gen);
        }
    }
    else
    {
        if(c_it_green!=redLink.greenList.begin())
        {
            const GreenLink &greenLink_prev = *(--c_it_green);

            std::uniform_real_distribution<> u(greenLink_prev.boundary, greenLink_prev.boundary+0.1*m_variance[1]);

            pix[1]=u(gen);
        }
        else
        {
//            std::uniform_real_distribution<> u(m_mean[1]-0.1*m_variance[1], m_mean[1]+0.1*m_variance[1]);

//            pix[1]=u(gen);

            std::normal_distribution<> n(m_mean[1], m_variance[1]);

            pix[1]=n(gen);
        }
    }

    //blue, 2

    if(blueLink.boundary != std::numeric_limits<double>::infinity())
    {
        if(c_it_blue!=greenLink.blueList.begin())
        {
            const BlueLink &blueLink_prev = *(--c_it_blue);

            std::uniform_real_distribution<> u(blueLink_prev.boundary, blueLink.boundary);

            pix[2]=u(gen);
        }
        else
        {
            std::uniform_real_distribution<> u(blueLink.boundary-0.1*m_variance[2], blueLink.boundary);

            pix[2]=u(gen);
        }
    }
    else
    {
        if(c_it_blue!=greenLink.blueList.begin())
        {
            const BlueLink &blueLink_prev = *(--c_it_blue);

            std::uniform_real_distribution<> u(blueLink_prev.boundary, blueLink_prev.boundary+0.1*m_variance[2]);

            pix[2]=u(gen);
        }
        else
        {
//            std::uniform_real_distribution<> u(m_mean[2]-0.1*m_variance[2], m_mean[2]+0.1*m_variance[2]);

//            pix[2]=u(gen);

            std::normal_distribution<> n(m_mean[2], m_variance[2]);

            pix[2]=n(gen);
        }
    }
}

void GaussianTransfer_Color_Recursive::transform_inv(ImageRGBd::PixelType& pix) const
{
    RedLink pRed;
    pRed.boundary = pix[0];

    GreenLink pGreen;
    pGreen.boundary = pix[1];

    BlueLink pBlue;
    pBlue.boundary = pix[2];

    std::list<RedLink>::const_iterator c_it_red=std::lower_bound(m_BnPAssociation.begin(), m_BnPAssociation.end(), pRed, Compare_Red(true));

    const RedLink &redLink = (*c_it_red);

    std::list<GreenLink>::const_iterator c_it_green=std::lower_bound(redLink.greenList.begin(), redLink.greenList.end(), pGreen, Compare_Green(true));

    const GreenLink &greenLink = (*c_it_green);

    std::list<BlueLink>::const_iterator c_it_blue =std::lower_bound(greenLink.blueList.begin(), greenLink.blueList.end(), pBlue, Compare_Blue(true));

    const BlueLink &blueLink = (*c_it_blue);

    pix[0]=redLink.pixel;
    pix[1]=greenLink.pixel;
    pix[2]=blueLink.pixel;

}

void GaussianTransfer_Color_Recursive::transform_inv(ImageRGBd& image) const
{
    image.for_all_pixels([this] (ImageRGBd::PixelType &pix)
    {
        transform_inv(pix);
    });
}

}
