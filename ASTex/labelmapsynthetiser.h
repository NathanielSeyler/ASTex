#ifndef LABELMAPSYNTHETISER_H
#define LABELMAPSYNTHETISER_H

#include "image_gray.h"
#include "image_rgb.h"
#include "bglam.h"
#include <random>
#include <itkImageRandomNonRepeatingIteratorWithIndex.h>

namespace ASTex {

template<typename IMG>
class LabelMapSynthetiser
{
private:
    struct my_less
    {
        bool operator() (ImageRGBu8::PixelType a, ImageRGBu8::PixelType b) {
            unsigned char ra = a.GetRed();
            unsigned char rb = b.GetRed();
            unsigned char ga = a.GetGreen();
            unsigned char gb = b.GetGreen();
            unsigned char ba = a.GetBlue();
            unsigned char bb = b.GetBlue();
            if(ra != rb)
                return ra < rb;
            else if(ga != gb)
                return ga < gb;
            else
                return ba < bb;
        }

        bool operator()(ImageGrayu8::PixelType a,ImageGrayu8::PixelType b) {return a<b;}

        bool operator()(std::vector<int> ca,std::vector<int> cb)
        {
            int compteur = 0;
            for(int valuea : ca)
            {
                int valueb = cb[compteur++];
                if(valuea != valueb)
                    return valuea < valueb;
            }
            return false;
        }
    };
    std::vector<typename IMG::PixelType> colors;
    std::set<std::vector<int>,my_less> configs;
    int nb_gray_level;
    Bglam bglam_in;
    Bglam bglam_out;

    void computeConfigs()
    {
        bglam_in.getImage().for_all_pixels([&](const ImageGrayu8::PixelType &p,int &x,int &y){
            auto config = bglam_in.getWindow(x,y);
            config.emplace_back(p);
            configs.emplace(config);
        });
        std::cout << configs.size() << " configs" << std::endl;
    }

    ImageRGBu8 computeErrorMap()
    {
        ImageRGBu8 error_map;
        int nb_errors = 0;
        error_map.initItk(bglam_out.getImage().width(),bglam_out.getImage().height(),true);
        bglam_out.getImage().for_all_pixels([&](const ImageGrayu8::PixelType &p,int &x,int &y){
            auto config = bglam_out.getWindow(x,y);
            config.emplace_back(p);
            auto it = configs.find(config);
            if(it==configs.end())
            {
                nb_errors++;
                error_map.pixelAbsolute(x,y) = itkRGBPixel(255,0,0);
            }
        });
        std::cout << nb_errors << " pixels hors config" << std::endl;
        std::cout << float(nb_errors)/(error_map.width()*error_map.height()) * 100 << "% pixels hors config" << std::endl;
        return error_map;
    }

    ImageGrayu8::PixelType aura2DSampling(const ImageGrayu8::PixelType &p, const int &x,const int &y,const float &d)
    {
        std::vector<ImageGrayu8::PixelType> C;
        for(int j = 0; j< nb_gray_level ; j++)
        {
            if(j!=p)
            {
                ImageGrayu8::PixelType pj = ImageGrayu8::itkPixel(j);
                bglam_out.updatePixel(pj,x,y);
                float dj = bglam_in.distance(bglam_out);
                if(dj < d)
                    C.emplace_back(pj);
            }
        }
        bglam_out.updatePixel(p,x,y);
        if(C.empty())
            return p;
        else
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(0,C.size()-1);
            return C[dis(gen)];
        }
    }

    void compute(const float &epsilon)
    {
        float d;
        int nb_modif = -1;
        while( (d = bglam_in.distance(bglam_out)) > epsilon && nb_modif != 0)
        {
            std::cout << d << "\t";
            nb_modif = 0;
            bglam_out.getImage().for_all_random_pixels([&] (ImageGrayu8::PixelType p,const int &x,const int &y){
                ImageGrayu8::PixelType np = aura2DSampling(p,x,y,d);
                if(np!=p)
                {
                    nb_modif++;
                    bglam_out.updatePixel(np,x,y);
                }
                /*if(np!=p)
                    nb_modif++;
                bglam_out.updatePixel(np,x,y);*/
                d = bglam_in.distance(bglam_out);
            });
            std::cout << nb_modif << std::endl;
        }
        std::cout << d << "\t";
        std::cout << nb_modif << std::endl;
    }

public:
    LabelMapSynthetiser() : colors(0),nb_gray_level(0) {}

    LabelMapSynthetiser(const IMG &i,const int &r) : colors(0),
        bglam_in(transfoImg(i),r,nb_gray_level),
        bglam_out(bglam_in.getRandImage(bglam_in.getImage().width()/8,
                                        bglam_in.getImage().height()/8),r,nb_gray_level)
    {
        float epsilon = 0;
        computeConfigs();
        compute(epsilon);
        transfoImgBack(bglam_out.getImage()).save("out.png");
        computeErrorMap().save("error_map.png");

        bglam_out.resolutionUp();
        compute(epsilon);
        transfoImgBack(bglam_out.getImage()).save("out_x2.png");
        computeErrorMap().save("error_x2_map.png");

        bglam_out.resolutionUp();
        compute(epsilon);
        transfoImgBack(bglam_out.getImage()).save("out_x4.png");
        computeErrorMap().save("error_x4_map.png");

        bglam_out.resolutionUp();
        compute(epsilon);
        transfoImgBack(bglam_out.getImage()).save("out_x8.png");
        computeErrorMap().save("error_x8_map.png");
    }

    LabelMapSynthetiser(const IMG &i) : LabelMapSynthetiser(i,1) {}

    ImageGrayu8 transfoImg(const IMG &i)
    {
        std::set<typename IMG::PixelType,my_less> tab;

        i.for_all_pixels([&](const typename IMG::PixelType &p) {
                tab.emplace(p);
        });

        nb_gray_level = tab.size();
        colors.clear();
        for(auto p : tab)
            colors.emplace_back(p);

        ImageGrayu8 img;
        img.initItk(i.width(),i.height());

        img.for_all_pixels([&](ImageGrayu8::PixelType &p,int x,int y){
           p = std::distance(tab.begin(),tab.find(i.pixelAbsolute(x,y)));
        });
        return img;
    }

    IMG transfoImgBack(const ImageGrayu8 &i)
    {
        IMG out;
        out.initItk(i.width(),i.height());
        out.for_all_pixels([&](typename IMG::PixelType &p,int &x,int &y){
            p = colors[i.pixelAbsolute(x,y)];
        });
        return out;
    }
};

}

#endif // LABELMAPSYNTHETISER_H
