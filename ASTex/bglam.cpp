#include "bglam.h"
#include <random>

namespace ASTex {


Bglam::Bglam(){}

Bglam::Bglam(const ImageGrayu8 &i) : img(i), bglams(0), ring(1)
{
    nbGray();
    compute();
}

Bglam::Bglam(const ImageGrayu8 &i,const int &r) : img(i), bglams(0), ring(r)
{
    nbGray();
    compute();
}

Bglam::Bglam(const ImageGrayu8 &i,const int &r,const int &n) : img(i), bglams(0), ring(r), nb_gray_level(n)
{
    compute();
}

Bglam::~Bglam(){}

std::vector<itk::VariableSizeMatrix<int>> Bglam::getBglams() {return bglams;}

const std::vector<itk::VariableSizeMatrix<int>> Bglam::getBglams() const {return bglams;}

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
    int tWindow = 1 + 2*ring;
    std::vector<int> window;

    for(int i= 0;i<tWindow;i++)
    {
        for(int j=0;j<tWindow;j++)
        {
            if(!(x-ring+j==x && y-ring+i==y))
            {
                if( x-ring+j >=0 && x-ring+j < img.width()
                        && y-ring+i >=0 && y-ring+i < img.height())
                {
                    ImageGrayu8::PixelType g = img.pixelAbsolute(x-ring+j,y-ring+i);
                    window.emplace_back(g);
                }
                else
                    window.emplace_back(-1);
            }
        }
    }
    return window;
}

void Bglam::compute()
{
    bglams.clear();
    int window = 1 + 2*ring;
    int nb_matrix = window * window - 1 ;
    std::cout << nb_matrix << " matrices " << nb_gray_level << "x" << nb_gray_level << std::endl;
    for(int i=0;i<nb_matrix;i++)
    {
        itk::VariableSizeMatrix<int> bglam(nb_gray_level,nb_gray_level);
        bglam.Fill(0);
        bglams.emplace_back(bglam);
    }
    img.for_all_pixels([&] (ImageGrayu8::PixelType &p,int &x,int &y){
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
        for(int i=0;i<nb_matrix;i++)
        {
            g = w[i];
            if(g >= 0)
                bglams[i][p][g]++;
        }
    });
}

double Bglam::distance(const Bglam &a)
{
    double distance = 0;
    unsigned int nb_matrix = (1+2*ring) * (1+2*ring) - 1;
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

void Bglam::updatePixel(const ImageGrayu8::PixelType &p, const int &x, const int &y)
{
    ImageGrayu8::PixelType old = img.pixelAbsolute(x,y);
    img.pixelAbsolute(x,y) = p;
    /*ImageGrayu8::PixelType g;

    int window = 1 + 2*ring;
    int compteur = 0;
    for(int i= 0;i<window;i++)
    {
        for(int j=0;j<window;j++)
        {
            if(!(x-ring+j==x && y-ring+i==y))
            {
                if( x-ring+j >=0 && x-ring+j < img.width()
                        && y-ring+i >=0 && y-ring+i < img.height())
                {
                    g = img.pixelAbsolute(x-ring+j,y-ring+i);
                    bglams[compteur][p][g]++;
                    bglams[compteur][old][g]--;
                }
                compteur++;
            }
        }
    }*/
    int twindow = 1 + 2*ring;
    int nb_neighbors = twindow * twindow - 1 ;
    auto window = getWindow(x,y);
    int g;
    int g2;
    for(int i =0;i<nb_neighbors;i++)
    {
        g = window[i];
        g2 = window[nb_neighbors-1-i];
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

}
