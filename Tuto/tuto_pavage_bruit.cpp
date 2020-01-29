#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/pavage_bruit.h>

#include <ASTex/easy_io.h>

using namespace ASTex;

int main()
{
    using IMG = ImageRGBd;
    IMG in;
    IO::loadu8_in_01(in,"../../pavage_bruit/textures/n03.png");
    int width = 2 * in.width();
    int height = 2 * in.height();
    //l compris entre 0 exclus et in.width * 0.5 exclus
    int l = in.width()*0.3 ;

    IMG out = PavageBruit<IMG>::tile(in,width,height,l);
    IO::save01_in_u8(out,"out.png");

    std::cout << "l (pixels) : " << l << std::endl;
    std::cout << "input : " << in.width() << " x " << in.height() << std::endl;
    std::cout << "output : " << out.width() << " x " << out.height() << std::endl;

    return 0;
}
