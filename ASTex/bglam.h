#ifndef BGLAM_H
#define BGLAM_H

#include "image_gray.h"

namespace ASTex {

class Bglam
{
private:
    ImageGrayu8 img;
    int nb_gray_level;
    Size radius;
    std::vector<bool> structuring_element;
    unsigned int nb_matrix;
    std::vector<itk::VariableSizeMatrix<int>> bglams;
    //std::vector<Offset> offsets;
    void nbGray();
    void compute();
    template<typename T>
    T norm(const itk::VariableSizeMatrix<T> &m)
    {
        T r(0);
        for(unsigned int i = 0; i < m.Rows();i++)
            for(unsigned int j=0; j < m.Cols();j++)
                r+=m[i][j];
        return r;
    }
    itk::VariableSizeMatrix<float> normalize(const itk::VariableSizeMatrix<int>&);
public:
    Bglam();
    //Bglam(const int &);
    Bglam(const ImageGrayu8 &);
    Bglam(const ImageGrayu8 &,const Size &);
    Bglam(const ImageGrayu8 &,const Size &,const int &);
    Bglam(const ImageGrayu8 &,const Size &,const int &,const std::vector<bool>&);
    ~Bglam();
    std::vector<itk::VariableSizeMatrix<int>> getBglams();
    const std::vector<itk::VariableSizeMatrix<int>> getBglams() const;
    void setBglams(std::vector<itk::VariableSizeMatrix<int>>);
    ImageGrayu8 getImage();
    ImageGrayu8 getRandImage(const int &w, const int &h);
    void resolutionUp();
    std::vector<int> getWindow(const int &,const int &);
    double distance(const Bglam &);
    friend std::ostream& operator<<(std::ostream&, const Bglam &);
    void updatePixel(const ImageGrayu8::PixelType &,const int &,const int &);
    void updateBglams(const ImageGrayu8::PixelType &,const ImageGrayu8::PixelType &,std::vector<int> &);
    Bglam clone();
};

}

#endif // BGLAM_H
