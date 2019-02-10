#ifndef BGLAM_H
#define BGLAM_H

#include "image_gray.h"
#include "image_rgb.h"

namespace ASTex {

class Bglam
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
            if(ra < rb)
                return true;
            else if(ra == rb && ga < gb)
                return true;
            else if(ra == rb && ga == gb && ba < bb)
                return true;
            else
                return false;
        }
    };

    ImageGrayu8 img;
    std::vector<itk::VariableSizeMatrix<int>> bglams;
    int ring;
    int nb_gray_level;
    void transfoImg(const ImageGrayu8 &);
    void transfoImg(const ImageRGBu8 &);
    void compute();
    int norm(const itk::VariableSizeMatrix<int> &);
    itk::VariableSizeMatrix<float> normalize(const itk::VariableSizeMatrix<int>&);
public:
    Bglam();
    Bglam(const ImageGrayu8 &);
    Bglam(const ImageRGBu8 &);
    ~Bglam();
    std::vector<itk::VariableSizeMatrix<int>> getBglams();
    const std::vector<itk::VariableSizeMatrix<int>> getBglams() const;
    double distance(const Bglam &);
    friend std::ostream& operator<<(std::ostream&, const Bglam &);
};

}

#endif // BGLAM_H
