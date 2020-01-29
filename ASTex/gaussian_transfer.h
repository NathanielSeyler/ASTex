#ifndef __GAUSSIAN_TRANSFER__H__
#define __GAUSSIAN_TRANSFER__H__

#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <ASTex/mask.h>
#include <type_traits>
#include <Eigen/Eigen>
#include <ASTex/histogram.h>
#include <random>
#include <bitset>

///About this class: to use, create an instance of either class (GaussianTransfer for grayd, and GaussianTransfer_Color_Recursive for rgbd),
///then call transform with an image and it'll become gaussian in order 1. Call transform_inv of the same instance and the histogram will be changed back.
///it distributes each pixel randomly into an interval, automatically determined into a gaussian.
///The distribution is uniform so far, thus the result only approximates a gaussian in order 1, but is sufficient for enough different pixels in one image.
///The glm dependency is to be removed.

namespace ASTex{

/**
 * \brief erfinv computes a single-precision floating point approximation of the inverse erf (std::erf),
 * ripped from on https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf ; this is needed for the GaussianTransfer main process
 * \param x parameter to be computed
 * \return approximation of erfinv(x)
 */
float erfinv(float x);

/**
 * \brief GaussianTransfer_Color_Recursive a fairly clean class written by a fairly good programmer.
 * It's also used to transform a ImageGrayd into another, except it becomes a gaussian texture in order 1 (or more if the input was a micro-texture).
 */
class GaussianTransfer
{
public:

    GaussianTransfer();

    /**
     * @brief transform determines a one-dimensional cut of a Gaussian of same mean and same variance than image and transforms the image so that it becomes a Gaussian texture in order 1.
     * The cut represents the unique combination of intervals such that the probability of a pixel being in one interval is equal to the probability in the histogram for this pixel,
     * and that these intervals are sorted according to the bins.
     * @pre image has a quantized number of bins, such as between 0 and 255.
     * @param image image to be transformed.
     */
    void transform(ImageGrayd &image);
    void transform_inv(ImageGrayd::PixelType& pix) const;
    void transform_inv(ImageGrayd& image) const;


private:

    std::vector<double> m_boundaries;
    std::vector<ImageGrayd::PixelType> m_pixels;

};

///pre-condition: used histograms have a strict order

/**
 * \brief GaussianTransfer_Color_Recursive a horrible class written by horrible people who did not want to bother.
 * It's also used to transform a ImageRGBd into another, except it's a gaussian texture in order 1 (or more if the input was a micro-texture).
 */
class GaussianTransfer_Color_Recursive
{
public:

    GaussianTransfer_Color_Recursive();

    /**
     * @brief transform determines a multidimensionnal cut of a Gaussian of same mean and same variance than image and transforms the image so that it becomes a Gaussian texture in order 1.
     * The cut represents the unique combination of intervals such that the probability of a pixel being in one interval is equal to the probability in the histogram for this pixel,
     * that these intervals are sorted according to the bins,
     * and that the blue channel's cut is constrained according to the green channel's cut, which is constrained by the red channel's cut.
     * @pre image has a quantized number of bins, such as between 0,0,0 and 255,255,255.
     * @param image image to be transformed.
     */
    void transform(ImageRGBd &image);

    /**
     * @brief inverse transformation
     * @pre transform was called earlier
     * @param image image to be transformed
     */
    void transform_inv(ImageRGBd& image) const;
    void transform_inv(ImageRGBd::PixelType& pix) const;

private:

    class BlueLink
    {
    public:
        BlueLink():pixel(0), boundary(0){}

        double pixel;
        double boundary;
    };

    class GreenLink
    {
    public:
        GreenLink():pixel(0), boundary(0), blueList(){}

        double pixel;
        double boundary;
        std::list<BlueLink> blueList;
    };

    class RedLink
    {
    public:
        RedLink():pixel(0), boundary(0), greenList(){}

        double pixel;
        double boundary;
        std::list<GreenLink> greenList;
    };

    class Compare_Red
    {
    public:
        Compare_Red(bool compare_boundary):m_cb(compare_boundary){}

        bool operator()(const RedLink &object, const RedLink &other) const {return m_cb ? object.boundary < other.boundary : object.pixel < other.pixel;}

    private:

        bool m_cb;
    };

    class Compare_Green
    {
    public:
        Compare_Green(bool compare_boundary):m_cb(compare_boundary){}

        bool operator()(const GreenLink &object, const GreenLink &other) const {return m_cb ? object.boundary < other.boundary : object.pixel < other.pixel;}

    private:

        bool m_cb;
    };

    class Compare_Blue
    {
    public:
        Compare_Blue(bool compare_boundary):m_cb(compare_boundary){}

        bool operator()(const BlueLink &object, const BlueLink &other) const {return m_cb ? object.boundary < other.boundary : object.pixel < other.pixel;}

    private:

        bool m_cb;
    };

    //closeList is used to compute the probabilities in the boundary and pixels association (BnP) structure, after all bins have been inserted.
    void closeList();

    //used on each pixels of the image to distribute it randomly inside their interval.
    void turn_pixel_into_a_more_gaussian_pixel(ImageRGBd::PixelType &pix, std::mt19937 &gen) const;

    //boundary and pixels association. "recursive" structure used to store the transformation.
    std::list<RedLink> m_BnPAssociation;

    Eigen::Vector3d m_mean;
    Eigen::Vector3d m_variance;

};

}

#endif //__GAUSSIAN_TRANSFER__H__
