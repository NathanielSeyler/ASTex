#include "bglam.h"

namespace ASTex {


Bglam::Bglam(){}

Bglam::Bglam(const ImageGrayu8 &i) : bglams(0), ring(1)
{
    transfoImg(i);
    compute();
}

Bglam::Bglam(const ImageGrayu8 &i,const int &r) : bglams(0), ring(r)
{
    transfoImg(i);
    compute();
}

Bglam::Bglam(const ImageRGBu8 &i) : bglams(0) , ring(1)
{
    transfoImg(i);
    compute();
}

Bglam::Bglam(const ImageRGBu8 &i,const int &r) : bglams(0) , ring(r)
{
    transfoImg(i);
    compute();
}

Bglam::~Bglam(){}

std::vector<itk::VariableSizeMatrix<int>> Bglam::getBglams() {return bglams;}

const std::vector<itk::VariableSizeMatrix<int>> Bglam::getBglams() const {return bglams;}

void Bglam::transfoImg(const ImageGrayu8 &i)
{
    std::set<unsigned char> tab;

    i.for_all_pixels([&](const ImageGrayu8::PixelType &p) {
            tab.emplace(p);
    });

    nb_gray_level = tab.size();

    img.initItk(i.width(),i.height());

    img.for_all_pixels([&](ImageGrayu8::PixelType &p,int x,int y){
       p = std::distance(tab.begin(),tab.find(i.pixelAbsolute(x,y)));
    });

}

void Bglam::transfoImg(const ImageRGBu8 &i)
{
    std::set<ImageRGBu8::PixelType,my_less> tab;

    i.for_all_pixels([&](const ImageRGBu8::PixelType &p) {
        tab.emplace(p);
    });

    nb_gray_level = tab.size();

    img.initItk(i.width(),i.height());

    img.for_all_pixels([&](ImageGrayu8::PixelType &p,int &x,int &y){
       p = std::distance(tab.begin(),tab.find(i.pixelAbsolute(x,y)));
    });

}


void Bglam::compute()
{
    std::cout << nb_gray_level << "x" << nb_gray_level << std::endl;
    int window = 1 + 2*ring;
    int nb_matrix = window * window - 1 ;
    for(int i=0;i<nb_matrix;i++)
    {
        itk::VariableSizeMatrix<int> bglam(nb_gray_level,nb_gray_level);
        bglam.Fill(0);
        bglams.emplace_back(bglam);
    }
    img.for_all_pixels([&] (ImageGrayu8::PixelType &p,int &x,int &y){
        ImageGrayu8::PixelType g;
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
        }
    });
}

double Bglam::distance(const Bglam &a)
{
    double distance = 0;
    unsigned int nb_matrix = (1+2*ring) * (1+2*ring) - 1;
    unsigned int cols = bglams[0].Cols();
    itk::VariableSizeMatrix<float> m;
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

int Bglam::norm(const itk::VariableSizeMatrix<int> &m)
{
    int r =0;
    for(unsigned int i = 0; i < m.Rows();i++)
        for(unsigned int j=0; j < m.Cols();j++)
            r+=m[i][j];
    return r;
}

itk::VariableSizeMatrix<float> Bglam::normalize(const itk::VariableSizeMatrix<int> &m)
{
    int n = norm(m);
    itk::VariableSizeMatrix<float> r(m.Rows(),m.Cols());
    for(unsigned int i =0;i<r.Rows();i++)
        for(unsigned int j=0;j<r.Cols();j++)
            r[i][j] = float(m[i][j])/float(n);
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

}
