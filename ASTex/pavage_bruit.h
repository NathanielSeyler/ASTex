#ifndef PAVAGE_BRUIT_H
#define PAVAGE_BRUIT_H

#include <ASTex/simplex.h>
#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/gaussian_transfer.h>
#include <ASTex/easy_io.h>
//#include <map>
#include <iostream>
#include <ctime>

namespace ASTex {

template<typename IMG>
class PavageBruit
{
private:

    static inline void boundMod(Vec2 & v,const int &w,const int &h)
    {
        v(0) = (static_cast<int> (v(0)) % w + w) % w;
        v(1) = (static_cast<int> (v(1)) % h + h) % h;
    }

    static inline void boundMirror(Vec2 & v,const int &w,const int &h)
    {
        int x = static_cast<int> (v(0));
        int y = static_cast<int> (v(1));
        if( x < 0)
            x = -x;
        else if(x >= w )
            x = 2 * (w - 1) - x;

        if(y < 0)
            y = -y;
        else if(y >= h)
            y = 2 * (h - 1) - y;

        v(0) = x;
        v(1) = y;
    }

    static inline void boundCenter(Vec2 & v,const int &w,const int &h,const int &l)
    {
        int x = static_cast<int> (v(0));
        int y = static_cast<int> (v(1));
        x = x %(w - 2*l )+l;
        y = y %(h - 2*l )+l;
        v(0) = x;
        v(1) = y;
    }

    static inline unsigned int hash(const Vec2 &v,const int &n,const int &random_offset)
    {
        unsigned int x = static_cast<unsigned int>(v(0));
        unsigned int y = static_cast<unsigned int>(v(1));
        unsigned int seed = (y%n) * n + (x%n) + random_offset;
        //int seed = morton(v[0],v[1]) + random_offset;
        unsigned int next  = 3039177861 * seed;
        return next;
    }

    static inline void clamp(ImageGrayd::DoublePixelEigen &c)
    {
        if(c < 0.)
            c = 0.;
        if(c > 1.)
            c = 1.;
    }

    static inline void clamp(ImageRGBd::DoublePixelEigen &c)
    {
        for(int i = 0 ; i<3 ;i++)
        {
            if(c(i) < 0.)
                c(i) = 0.;
            if(c(i) > 1.)
                c(i) = 1.;
        }
    }

    static inline GaussianTransfer transform(ImageGrayd & im)
    {
        GaussianTransfer gt;
        gt.transform(im);
        return gt;
    }

    static inline void transform_inv(const GaussianTransfer &gt,ImageGrayd & im)
    {
        gt.transform_inv(im);
    }

    static inline GaussianTransfer_Color_Recursive transform(ImageRGBd & im)
    {
        GaussianTransfer_Color_Recursive gt;
        gt.transform(im);
        return gt;
    }

    static inline void transform_inv(const GaussianTransfer_Color_Recursive &gt,ImageRGBd & im)
    {
        gt.transform_inv(im);
    }


//    static inline ImageGrayd::DoublePixelEigen erf(const ImageGrayd::DoublePixelEigen & p)
//    {
//        return 0.5 + 0.5 * std::erf((p - 0.5) / (6 * std::sqrt(2.)));
//    }

//    static inline ImageRGBd::DoublePixelEigen erf(const ImageRGBd::DoublePixelEigen & p)
//    {
//        ImageRGBu8::DoublePixelEigen res;
//        for(int i = 0; i <3 ;i++)
//            res(i) = 0.5 + 0.5 * std::erf((p(i) - 0.5) / (6 * std::sqrt(2.)));
//        return res;
//    }

    static inline typename IMG::PixelType linear_blending(
            const typename IMG::DoublePixelEigen &p0,double &w0,
            const typename IMG::DoublePixelEigen &p1,double &w1,
            const typename IMG::DoublePixelEigen &p2,double &w2)
    {
        typename IMG::DoublePixelEigen color = (w0 * p0 + w1 * p1 + w2 * p2 ) / (w0 + w1 + w2);

        clamp(color);
        return IMG::itkPixel(color);
    }

    static inline typename IMG::PixelType cov_blending(
            const typename IMG::DoublePixelEigen &p0,double &w0,
            const typename IMG::DoublePixelEigen &p1,double &w1,
            const typename IMG::DoublePixelEigen &p2,double &w2,
            const typename IMG::DoublePixelEigen &mean)
    {
        double sum = w0 + w1 + w2 ;
        w0 /= sum;
        w1 /= sum;
        w2 /= sum;
        double W = std::sqrt(w0 * w0 + w1 * w1 + w2 * w2);

        typename IMG::DoublePixelEigen color =
                ( w0 * (p0 - mean) + w1 * (p1 - mean) + w2 * (p2 - mean)) / W + mean;

        clamp(color);
        return IMG::itkPixel(color);
    }

//    static inline typename IMG::PixelType lookUpTable(typename IMG::DoublePixelEigen &cov)
//    {
//        //clamp(cov);
//        //typename IMG::DoublePixelEigen U = erf(cov);
//        return IMG::itkPixel(cov);
//    }
public:
    PavageBruit(){}

    static inline IMG tile( IMG &in,const int &width,const int &height,const int &l)
    {
        std::srand(std::time(0));
        int random_offset(std::rand());

        auto gt = transform(in);
        //IO::save01_in_u8(in,"in_gaussian.png");

        int in_width = in.width();
        int in_height = in.height();
        int nb_pixels = in_width * in_height;

        typename IMG::DoublePixelEigen mean = IMG::eigenDoublePixel(0.);
        in.for_all_pixels([&](const typename IMG::PixelType &pix)
        {
            mean += IMG::eigenDoublePixel(pix);
        });

        mean = mean / double(nb_pixels);
//        typename IMG::DoublePixelEigen mean= IMG::eigenDoublePixel(0.5);

        IMG out(width,height);
        ImageGrayu8 grid(width,height);
        ImageGrayu8 center(in_width,in_height);

        out.parallel_for_all_pixels([&](typename IMG::PixelType &pix,int x,int y)
        {
            Vec2 p_l(double(x)/l,double(y)/l);

            std::vector<Vec2> vertices = Simplex2D::vertices(p_l,grid,x,y);
            Vec2 v0_prime = vertices[0];
            Vec2 v1_prime = vertices[1];
            Vec2 v2_prime = vertices[2];
            //weights

            Vec2 d0 = p_l - Simplex2D::unskew(v0_prime);
            Vec2 d1 = p_l - Simplex2D::unskew(v1_prime);
            Vec2 d2 = p_l - Simplex2D::unskew(v2_prime);

            double w0 = Simplex2D::contribution(d0);
            double w1 = Simplex2D::contribution(d1);
            double w2 = Simplex2D::contribution(d2);

            //trouver les pixels dans l'input image a partir des tuiles

            int c0 = hash(v0_prime,width * height /(l*l),random_offset) % nb_pixels;
            Vec2 c00(c0%in_width,c0/in_width);
            boundCenter(c00,in_width,in_height,l);

            int c1 = hash(v1_prime,width * height/(l*l),random_offset) % nb_pixels;
            Vec2 c11(c1%in_width,c1/in_width);
            boundCenter(c11,in_width,in_height,l);

            int c2 = hash(v2_prime,width * height/(l*l),random_offset) % nb_pixels;
            Vec2 c22(c2%in_width,c2/in_width);
            boundCenter(c22,in_width,in_height,l);

            center.pixelAbsolute(c00(0),c00(1)) = 255;
            center.pixelAbsolute(c11(0),c11(1)) = 255;
            center.pixelAbsolute(c22(0),c22(1)) = 255;

            Vec2 in_p0 = c00 + d0 * l;
            //boundMod(in_p0,in_width,in_height);
            Vec2 in_p1 = c11 + d1 * l;
            //boundMod(in_p1,in_width,in_height);
            Vec2 in_p2 = c22 + d2 * l;
            //boundMod(in_p2,in_width,in_height);

            pix = cov_blending( in.pixelEigenAbsolute(in_p0(0), in_p0(1)), w0,
                                in.pixelEigenAbsolute(in_p1(0), in_p1(1)), w1,
                                in.pixelEigenAbsolute(in_p2(0), in_p2(1)), w2,
                                mean);
            gt.transform_inv(pix);

        });

//        IO::save01_in_u8(out,"out_gaussian.png");
//        gt.transform_inv(out);

        //grid.save("grid.png");
        //center.save("center.png");
        return out;
    }
};

}


#endif // PAVAGE_BRUIT_H
