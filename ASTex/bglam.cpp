#include "bglam.h"
#include <random>
//#include <itpp/signal/fastica.h>

namespace ASTex {


Bglam::Bglam(){}

//Bglam::Bglam(const int &r) : img(nullptr),bglams(0),ring(r),nb_gray_level(0) {}

Bglam::Bglam(const ImageGrayu8 &i) : img(i),nb_matrix(8), bglams(0)
{
    nbGray();
    radius = {{1,1}};
    for(int i = 0; i < 8;i++)
        structuring_element.push_back(true);
    compute();
}

Bglam::Bglam(const ImageGrayu8 &i, const Size &r) : img(i),radius(r),structuring_element(0), bglams(0)
{
    nbGray();
    nb_matrix = (1+2*radius[0]) * (1+2*radius[1]) -1;
    for(unsigned int i = 0; i < nb_matrix;i++)
        structuring_element.push_back(true);
    compute();
}

Bglam::Bglam(const ImageGrayu8 &i, const Size &r, const int &n) : img(i), nb_gray_level(n), radius(r), bglams(0)
{
    nb_matrix = (1+2*radius[0]) * (1+2*radius[1]) -1;
    for(unsigned int i = 0; i < nb_matrix;i++)
        structuring_element.push_back(true);
    compute();
}

Bglam::Bglam(const ImageGrayu8 &i, const Size &r, const int &n, const std::vector<bool> &se) : img(i), nb_gray_level(n), radius(r),structuring_element(se), bglams(0)
{
    nb_matrix = 0;
    for(auto e : structuring_element)
        if(e)
            nb_matrix++;
    compute();
}

Bglam::~Bglam(){}

std::vector<itk::VariableSizeMatrix<int>> Bglam::getBglams() {return bglams;}

const std::vector<itk::VariableSizeMatrix<int>> Bglam::getBglams() const {return bglams;}

void Bglam::setBglams(std::vector<itk::VariableSizeMatrix<int>> b) {bglams =b;}

ImageGrayu8 Bglam::getImage() {return img;}

ImageGrayu8 Bglam::getRandImage(const int &w,const int&h)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0,nb_gray_level-1);
    ImageGrayu8 r;
    r.initItk(w,h);
    r.for_all_pixels([&] (ImageGrayu8::PixelType &p) {
       p = dis(gen);
    });
    return r;
}

void Bglam::resolutionUp()
{
    ImageGrayu8 i;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> disGraylevel(0,nb_gray_level-1);
    std::uniform_int_distribution<> disNeighbors(4,5);
    i.initItk(2*img.width(),2*img.height());
    img.for_all_pixels([&](const ImageGrayu8::PixelType &p,const int &x,const int &y){
        switch (disNeighbors(gen)) {
        case 0:
            i.pixelAbsolute(2*x,2*y) = disGraylevel(gen);
            i.pixelAbsolute(2*x+1,2*y) = p;
            i.pixelAbsolute(2*x,2*y+1) = p;
            i.pixelAbsolute(2*x+1,2*y+1) = p;
            break;
        case 1:
            i.pixelAbsolute(2*x,2*y) = p;
            i.pixelAbsolute(2*x+1,2*y) = disGraylevel(gen);
            i.pixelAbsolute(2*x,2*y+1) = p;
            i.pixelAbsolute(2*x+1,2*y+1) = p;
            break;
        case 2:
            i.pixelAbsolute(2*x,2*y) = p;
            i.pixelAbsolute(2*x+1,2*y) = p;
            i.pixelAbsolute(2*x,2*y+1) = disGraylevel(gen);
            i.pixelAbsolute(2*x+1,2*y+1) = p;
            break;
        case 3:
            i.pixelAbsolute(2*x,2*y) = p;
            i.pixelAbsolute(2*x+1,2*y) = p;
            i.pixelAbsolute(2*x,2*y+1) = p;
            i.pixelAbsolute(2*x+1,2*y+1) = disGraylevel(gen);
            break;
        default:
            i.pixelAbsolute(2*x,2*y) = p;
            i.pixelAbsolute(2*x+1,2*y) = p;
            i.pixelAbsolute(2*x,2*y+1) = p;
            i.pixelAbsolute(2*x+1,2*y+1) = p;
            break;
        }
    });
    img = i;
    compute();
}

void Bglam::nbGray()
{
    unsigned char max = 0;
    img.for_all_pixels([&](const ImageGrayu8::PixelType &p){
        if(p>max)
            max = p;
    });
    nb_gray_level = max + 1 ;
}

std::vector<int> Bglam::getWindow(const int &x, const int &y)
{
    int rx(radius[0]);
    int ry(radius[1]);
    int tWindow_w = 1 + 2 * rx;
    int tWindow_h = 1 + 2 * ry;
    std::vector<int> window;

    int compteur = 0;
    for(int i= 0;i<tWindow_h;i++)
    {
        for(int j=0;j<tWindow_w;j++)
        {
            if(!(x-rx+j==x && y-ry+i==y))
            {
                if(structuring_element[compteur++])
                {
                    if( x-rx+j >=0 && x-rx+j < img.width()
                            && y-ry+i >=0 && y-ry+i < img.height())
                    {
                        ImageGrayu8::PixelType g = img.pixelAbsolute(x-rx+j,y-ry+i);
                        window.emplace_back(g);
                    }
                    else
                        window.emplace_back(-1);
                }
            }
        }
    }
    return window;
}

void Bglam::compute()
{
    bglams.clear();
    std::cout << nb_matrix << " matrices " << nb_gray_level << "x" << nb_gray_level << std::endl;
    for(unsigned int i=0;i<nb_matrix;i++)
    {
        itk::VariableSizeMatrix<int> bglam(nb_gray_level,nb_gray_level);
        bglam.Fill(0);
        bglams.emplace_back(bglam);
    }
    img.for_all_pixels([&] (const ImageGrayu8::PixelType &p,const int &x,const int &y){
        /*ImageGrayu8::PixelType g;
        int compteur = 0;
        for(int i= 0;i<window;i++)
        {
            for(int j=0;j<window;j++)
            {
                if(!(x-ring+j==x &&y-ring+i==y))
                {
                    if( x-ring+j >=0 && x-ring+j < img.width()
                            && y-ring+i >=0 && y-ring+i < img.height())
                    {
                        g = img.pixelAbsolute(x-ring+j,y-ring+i);
                        bglams[compteur][p][g]++;
                    }
                    compteur++;
                }
            }
        }*/
        int g;
        auto w = getWindow(x,y);
        for(unsigned int i=0;i<nb_matrix;i++)
        {
            g = w[i];
            if(g >= 0)
                bglams[i][p][g]++;
        }
    });
    /*img.for_all_pixels_structuring_element(radius,[&](const ImageGrayu8::PixelType &p,
                                           const ConstShapedNeighborhoodIterator n)
    {
        for(int i =0 ; i <n.GetActiveIndexListSize();i++)
        {
            int g = n.Get();
            if( g>=0)
                bglams[i][p][g]++;

        }
    });*/
    /*std::mutex mutex;
    img.parallel_for_all_pixels([&] (const ImageGrayu8::PixelType &p,int x, int y){
        int g;
        auto w = getWindow(x,y);
        for(int i=0;i<nb_matrix;i++)
        {
            g = w[i];
            if(g >= 0)
            {
                mutex.lock();
                bglams[i][p][g]++;
                mutex.unlock();
            }
        }
    });*/
}

double Bglam::distance(const Bglam &a)
{
    double distance = 0;
    unsigned int cols = bglams[0].Cols();
    itk::VariableSizeMatrix<float> m(cols,cols);
    for(unsigned int i =0;i<nb_matrix;i++)
    {
        m = normalize(bglams[i]) - normalize(a.getBglams()[i]);
        for(unsigned int j =0 ; j< cols;j++)
        {
            for(unsigned int k =0;k<cols;k++)
                distance+= std::abs(m[j][k]);
        }
    }
    return distance / nb_matrix;
}

itk::VariableSizeMatrix<float> Bglam::normalize(const itk::VariableSizeMatrix<int> &m)
{
    int n = norm(m);
    itk::VariableSizeMatrix<float> r(m.Rows(),m.Cols());
    for(unsigned int i =0;i<r.Rows();i++)
        for(unsigned int j=0;j<r.Cols();j++)
            r[i][j] = float(m[i][j])/float(n);
    //std::cout << norm(r) << std::endl;
    return r;
}

std::ostream& operator<<(std::ostream& out, const Bglam & b){
    out << std::endl << "-------------------------------------" << std::endl << std::endl ;
    for(auto bglam : b.bglams)
    {
        out << bglam << std::endl;
        out << std::endl << "-------------------------------------" << std::endl << std::endl ;
    }
    return out;
}

void Bglam::updatePixel(const ImageGrayu8::PixelType &p, const int &x, const int &y)
{
    auto old = img.pixelAbsolute(x,y);
    img.pixelAbsolute(x,y) = p;
    auto window = getWindow(x,y);
    updateBglams(old,p,window);

}

void Bglam::updateBglams(const ImageGrayu8::PixelType &old, const ImageGrayu8::PixelType &p, std::vector<int> &window)
{
    int g;
    int g2;
    for(unsigned int i =0;i<nb_matrix;i++)
    {
        g = window[i];
        g2 = window[nb_matrix-1-i];
        if(g>=0)
        {
            bglams[i][p][g]++;
            bglams[i][old][g]--;
        }
        if(g2>=0)
        {
            bglams[i][g2][p]++;
            bglams[i][g2][old]--;
        }
    }

}

/*Bglam Bglam::clone()
{
   Bglam b(ring);
   b.setBglams(bglams);
   return b;
}*/

}
