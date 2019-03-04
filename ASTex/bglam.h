#ifndef BGLAM_H
#define BGLAM_H

#include "image_gray.h"

namespace ASTex {

class Bglam
{
private:
    ImageGrayu8 img;
    std::vector<itk::VariableSizeMatrix<int>> bglams;
    int ring;
    int nb_gray_level;
    void nbGray();
    void compute();
    int norm(const itk::VariableSizeMatrix<int> &);
    itk::VariableSizeMatrix<float> normalize(const itk::VariableSizeMatrix<int>&);
public:
    Bglam();
    Bglam(const ImageGrayu8 &);
    Bglam(const ImageGrayu8 &,const int &);
    Bglam(const ImageGrayu8 &,const int &,const int &);
    ~Bglam();
    std::vector<itk::VariableSizeMatrix<int>> getBglams();
    const std::vector<itk::VariableSizeMatrix<int>> getBglams() const;
    ImageGrayu8 getImage();
    ImageGrayu8 getRandImage();
    std::vector<int> getWindow(const int &,const int &);
    double distance(const Bglam &);
    friend std::ostream& operator<<(std::ostream&, const Bglam &);
    void updatePixel(ImageGrayu8::PixelType &,const int &,const int &);
};

}

#endif // BGLAM_H
