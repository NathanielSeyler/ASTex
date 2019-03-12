#include <iostream>
#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/histogram.h>
#include <ASTex/bglam.h>
#include <ASTex/labelmapsynthetiser.h>
#include <itkMatrix.h>
#include <itkVariableSizeMatrix.h>
#include <set>

using namespace ASTex;
//typedef RGB<unsigned char> RGBType;

int main()
{
    ImageRGBu8 label_map3;
    label_map3.load("../../Label_maps_PNG/PlasterDamaged0232_S/PlasterDamaged0232_S_NH_167.135529.png");

    auto start_chrono = std::chrono::system_clock::now();
    LabelMapSynthetiser<ImageRGBu8> s(label_map3);
    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe timing: " << elapsed_seconds.count() << " s." << std::endl;

    return EXIT_SUCCESS;

    ImageGrayu8 img;
    img.initItk(5,5,true);
    img.pixelAbsolute(1,0) = 255;
    img.pixelAbsolute(2,0) = 255;
    img.pixelAbsolute(3,0) = 255;
    img.pixelAbsolute(4,0) = 255;
    img.pixelAbsolute(3,1) = 255;
    img.pixelAbsolute(4,1) = 255;
    img.pixelAbsolute(1,2) = 255;
    img.pixelAbsolute(2,2) = 255;
    img.pixelAbsolute(4,2) = 255;
    img.pixelAbsolute(2,3) = 255;
    img.pixelAbsolute(0,4) = 255;
    img.pixelAbsolute(4,4) = 255;

    ImageGrayu8 img3;
    img3.initItk(5,5,true);
    img3.pixelAbsolute(1,0) = 125;
    img3.pixelAbsolute(2,0) = 125;
    img3.pixelAbsolute(3,0) = 125;
    img3.pixelAbsolute(4,0) = 125;
    img3.pixelAbsolute(3,1) = 125;
    img3.pixelAbsolute(4,1) = 125;
    img3.pixelAbsolute(1,2) = 125;
    img3.pixelAbsolute(2,2) = 125;
    img3.pixelAbsolute(4,2) = 125;
    img3.pixelAbsolute(2,3) = 125;
    img3.pixelAbsolute(0,4) = 125;
    img3.pixelAbsolute(4,4) = 125;

    start_chrono = std::chrono::system_clock::now();

    LabelMapSynthetiser<ImageGrayu8> synthe;
    Bglam bglams(synthe.transfoImg(img),2);
    ImageGrayu8 img2 = bglams.getRandImage(5,5);
    img2 = synthe.transfoImgBack(img2);
    img2.save("img_random.png");
    Bglam bglams2(synthe.transfoImg(img2),2);
    Bglam bglams3(synthe.transfoImg(img3),2);
    std::cout << "distance img1 et img2 : " << bglams.distance(bglams2) << std::endl;
    std::cout << "distance img1 et img3 : " << bglams.distance(bglams3) << std::endl;

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "bglam timing: " << elapsed_seconds.count() << " s." << std::endl;

    //std::cout << "bglams of img1\n" << bglams << std::endl;
    //std::cout << "bglams of img2\n" << bglams2 << std::endl;
    //std::cout << "bglams of img3\n" << bglams3 << std::endl;

    ImageRGBu8 label_map1;
    //label_map1.load("../../Label_maps_PNG/flowered_wall_1k/flowered_wall_1k_cl_1.uhd.png");
    label_map1.load("../../Label_maps_PNG/PlasterDamaged0232_S/PlasterDamaged0232_S_NH_167.135529.png");
    //label_map1.load("../../Label_maps_PNG/Moss0101_S/Moss0101_S_NH_298.431677.png");

    ImageRGBu8 label_map2;
    label_map2.load("../../Label_maps_PNG/PlasterDamaged0232_S/noise.png");

    LabelMapSynthetiser<ImageRGBu8> synthe2;
    start_chrono = std::chrono::system_clock::now();

    Bglam blab(synthe2.transfoImg(label_map1),2);
    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "bglam timing: 5x5 " << elapsed_seconds.count() << " s." << std::endl;

    start_chrono = std::chrono::system_clock::now();
    Bglam blab2(synthe2.transfoImg(label_map2),2);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "bglam timing: 5x5 " << elapsed_seconds.count() << " s." << std::endl;

    start_chrono = std::chrono::system_clock::now();
    std::cout << blab.distance(blab2) << std::endl;

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "distance timing: " << elapsed_seconds.count() << " s." << std::endl;


    return EXIT_SUCCESS;

    ImageRGBu8 X;
    X.load("../../Label_maps_PNG/flowered_wall_1k/flowered_wall_1k.png");
    start_chrono = std::chrono::system_clock::now();
    /*HistogramRGBu8 h(X);
    unsigned int nb_color = h.binsNumber();
    h.saveFullHistogram("histo.txt");
    ImageRGBu8::PixelType p_mean = h.meanPixelType();
    std::cout << nb_color << std::endl;*/

    //return EXIT_SUCCESS;
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
    float r_mean = float(r) / float(n);
    float g_mean = float(g) / float(n);
    float b_mean = float(b) / float(n);
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
    vnl_svd<float>(C.GetVnlMatrix());

    //calculate T and T's inverse
    /*T = S.GetInverse() * U;
    inverse_T = U * S;*/

    //calculate Y the img with decorrelated color channels
    //Y = T * D;

    //calculate bglams on Y's color channels
    ImageGrayu8 channel0,channel1,channel2;
    channel0.initItk(X.width(),X.height());
    channel1.initItk(X.width(),X.height());
    channel2.initItk(X.width(),X.height());
    for(int i= 0;i<X.height();i++)
        for(int j=0;j<X.width();j++)
        {
            channel0.pixelAbsolute(j,i) = D[0][X.width() * i + j];
            channel1.pixelAbsolute(j,i) = D[1][X.width() * i + j];
            channel2.pixelAbsolute(j,i) = D[2][X.width() * i + j];
        }
    Bglam b0(channel0);
    Bglam b1(channel1);
    Bglam b2(channel2);
    //std::cout << b0.distance(b0) << std::endl;
    std::cout << b0.distance(b1) << std::endl;
    std::cout << b0.distance(b2) << std::endl;

    //std::cout << b1.distance(b0) << std::endl;
    //std::cout << b1.distance(b1) << std::endl;
    std::cout << b1.distance(b2) << std::endl;

    //std::cout << b2.distance(b0) << std::endl;
    //std::cout << b2.distance(b1) << std::endl;
    //std::cout << b2.distance(b2) << std::endl;

    channel0.save("tmp0.png");
    channel1.save("tmp1.png");
    channel2.save("tmp2.png");
    //synthetise

    //transfo back

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "bglam timing: " << elapsed_seconds.count() << " s." << std::endl;

    return EXIT_SUCCESS;
}
