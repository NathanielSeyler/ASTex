#include <iostream>
#include <ASTex/image_rgb.h>
#include <ASTex/histogram.h>
#include <itkMatrix.h>
#include <itkVariableSizeMatrix.h>
#include <vector>

using namespace ASTex;
//typedef RGB<unsigned char> RGBType;



int main()
{
    ImageRGBu8 X;
    X.load("../../Label_maps_PNG/flowered_wall_1k/flowered_wall_1k.png");
    auto start_chrono = std::chrono::system_clock::now();
    /*HistogramRGBu8 h(X);
    unsigned int nb_color = h.binsNumber();
    h.saveFullHistogram("histo.txt");
    ImageRGBu8::PixelType p_mean = h.meanPixelType();
    std::cout << nb << std::endl;*/

    // substract the mean color from each RGB color channel
    int n = X.width()*X.height();
    itk::VariableSizeMatrix<float> D(3,n);
    int r = 0;
    int g = 0;
    int b = 0;
    X.for_all_pixels([&] (ImageRGBu8::PixelType &p){
        r += p.GetRed();
        g += p.GetGreen();
        b += p.GetBlue();
    });
    float r_mean = (float)r / n;
    float g_mean = (float)g / n;
    float b_mean = (float)b / n;
    int i = 0;
    X.for_all_pixels([&](ImageRGBu8::PixelType &p){
        D[0][i] = p.GetRed()-r_mean;
        D[1][i] = p.GetGreen()-g_mean;
        D[2][i++] = p.GetBlue()-b_mean;

    });

    //calculate the covariance matrix
    itk::Matrix<float,3,3> C;
    C = D * D.GetTranspose();

    //perform SVD on C

    //calculate T and T's inverse
    /*T = S.GetInverse() * U;
    inverse_T = U * S;*/

    //calculate Y the img with decorrelated color channels
    //Y = T * D;

    //calculate bglams on Y's color channels
    ImageGrayf channel0,channel1,channel2;
    channel0.initItk(X.width(),X.height());
    channel1.initItk(X.width(),X.height());
    channel2.initItk(X.width(),X.height());
    for(int x= 0;x<X.height();x++)
        for(int y=0;y<X.width();y++)
        {
            channel0.pixelAbsolute(y,x) = D[0][X.width() * x + y];
            channel1.pixelAbsolute(y,x) = D[1][X.width() * x + y];
            channel2.pixelAbsolute(y,x) = D[2][X.width() * x + y];
        }

    /*channel0.save(TEMPO_PATH+"tmp0.png");
    channel1.save(TEMPO_PATH+"tmp1.png");
    channel2.save(TEMPO_PATH+"tmp2.png");*/
    //synthetise

    //transfo back

    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "bglam timing: " << elapsed_seconds.count() << " s." << std::endl;
    //itk::VariableSizeMatrix<unsigned int> m(nb_color,nb_color);
    //m.Fill(0);
    //itk::Matrix<double,nb,nb> m() ;
    /*for(int i =0;i<img.height();i++)
    {
        for(int j=0;j<img.width();j++)
        {
            auto p = img.pixelAbsolute(i, j);
*/
            /*unsigned char r = p.GetRed();
            unsigned char g = p.GetGreen();
            unsigned char b = p.GetBlue();*/
            //p = ImageRGBu8::itkPixel(127,127,127);
            //std::cout << p  << std::endl /*<< (int)r << std::endl << (int)g << std::endl << (int)b << std::endl*/;
   /*    }
    }*/
    return EXIT_SUCCESS;
}
