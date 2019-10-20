#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/histogram.h>
#include <cmath>

using namespace ASTex;

float square_to_triangle_gridX(int r,int x,int y)
{
   return float(x) / r - float(y) / ( r * std::sqrt(3));
}

float square_to_triangle_gridY(int r,int y)
{
   return 2 * float(y) / ( r * std::sqrt(3));
}

int triangle_to_square_gridX(int r ,float u,float v){
    return int(r * u + r * v / 2);
}

int triangle_to_square_gridY(int r,float v){
    return int( r * v * std::sqrt(3) / 2);
}

ImageRGBu8 pavage_bruit(ImageRGBu8 in,int r)
{
    ImageRGBu8 out(in.width(),in.height());

    /*out.parallel_for_all_pixels([&](ImageRGBu8::PixelType &p,int x,int y)
    {
       int l = (2 * y) / (r * std::sqrt(3));
       //int l = y / r;
       if(l%2 == 0)
       {
           float pgx = (int)x/r * r ;
           float pgy = l * r * std::sqrt(3)/2;
           float pdx = pgx + r;
           float pdy = pgy;
           float pmx = pgx + r/2;
           float pmy = (l+1) * r * std::sqrt(3)/2;

           float cp1 = (pgy - pmy) * (x - pmx) - (pgx - pmx) * (y - pmy);
           if(cp1 > 0)
               p = ImageRGBu8::itkPixel(0,0,255);
           else
           {
               float cp2 = (pdy - pmy) * (x - pmx) - (pdx - pmx) * (y - pmy);
               if(cp2 >= 0)
                   p = ImageRGBu8::itkPixel(255,0,0);
               else
                   p = ImageRGBu8::itkPixel(0,0,255);
           }
       }
       else
       {
           float pgx = (int)x/r * r;
           float pgy = (l+1) * r * std::sqrt(3)/2;
           float pdx = pgx + r;
           float pdy = pgy;
           float pmx = pgx + r/2;
           float pmy = l * r * std::sqrt(3)/2;

           float cp1 = (pgy - pmy) * (x - pmx) - (pgx - pmx) * (y - pmy);
           if(cp1 < 0)
               p = ImageRGBu8::itkPixel(255,0,0);
           else
           {
               float cp2 = (pdy - pmy) * (x - pmx) - (pdx - pmx) * (y - pmy);
               if(cp2 <= 0)
                   p = ImageRGBu8::itkPixel(0,0,255);
               else
                   p = ImageRGBu8::itkPixel(255,0,0);
           }
       }
    });*/

    out.parallel_for_all_pixels([&](ImageRGBu8::PixelType &p,int x,int y)
    {
       float tx = square_to_triangle_gridX(r,x,y);
       float ty = square_to_triangle_gridY(r,y);

       float p0x = std::floor(tx);
       float p0y = std::floor(ty);

       float p1x = p0x + 1;
       float p1y = p0y;

       float p2x = p0x;
       float p2y = p0y + 1;

       float cp = (p2y - p1y) * (tx - p1x) - (p2x - p1x) * (ty - p1y);

       if(cp >= 0)
       {
           p0x = p1x;
           p0y = p2y;
       }

       float t0x = tx - p0x;
       float t0y = ty - p0y;

       float t1x = tx - p1x;
       float t1y = ty - p1y;

       float t2x = tx - p2x;
       float t2y = ty - p2y;

       std::srand(triangle_to_square_gridX(r,p0x,p0y));
       int r0x = std::rand() % in.width();
       std::srand(triangle_to_square_gridY(r,p0y));
       int r0y = std::rand() % in.height();
       float c0x = square_to_triangle_gridX(r,r0x,r0y);
       float c0y = square_to_triangle_gridY(r,r0y);

       std::srand(triangle_to_square_gridX(r,p1x,p1y));
       int r1x = std::rand() % in.width();
       std::srand(triangle_to_square_gridY(r,p1y));
       int r1y = std::rand() % in.height();
       float c1x = square_to_triangle_gridX(r,r1x,r1y);
       float c1y = square_to_triangle_gridY(r,r1y);

       std::srand(triangle_to_square_gridX(r,p2x,p2y));
       int r2x = std::rand() % in.width();
       std::srand(triangle_to_square_gridY(r,p2y));
       int r2y = std::rand() % in.height();
       float c2x = square_to_triangle_gridX(r,r2x,r2y);
       float c2y = square_to_triangle_gridY(r,r2y);

       int in_p0x = (triangle_to_square_gridX(r,c0x+t0x,c0y+t0y) + in.width()) % in.width();
       int in_p0y = (triangle_to_square_gridY(r,c0y+t0y) + in.height()) % in.height();

       int in_p1x = (triangle_to_square_gridX(r,c1x+t1x,c1y+t1y) + in.width()) % in.width();
       int in_p1y = (triangle_to_square_gridY(r,c1y+t1y) + in.height()) % in.height();

       int in_p2x = (triangle_to_square_gridX(r,c2x+t2x,c2y+t2y) + in.width()) % in.width();
       int in_p2y = (triangle_to_square_gridY(r,c2y+t2y) + in.height()) % in.height();

       int color_r =
               (in.pixelAbsolute(in_p0x,in_p0y).GetRed() +
               in.pixelAbsolute(in_p1x,in_p1y).GetRed() +
               in.pixelAbsolute(in_p2x,in_p2y).GetRed())/3;
       int color_g =
               (in.pixelAbsolute(in_p0x,in_p0y).GetGreen() +
               in.pixelAbsolute(in_p1x,in_p1y).GetGreen() +
               in.pixelAbsolute(in_p2x,in_p2y).GetGreen())/3;
       int color_b =
               (in.pixelAbsolute(in_p0x,in_p0y).GetBlue() +
               in.pixelAbsolute(in_p1x,in_p1y).GetBlue() +
               in.pixelAbsolute(in_p2x,in_p2y).GetBlue())/3;

       p = ImageRGBu8::itkPixel(color_r,color_g,color_b);

    });

    return out;
}

ImageGrayu8 pavage_bruit(ImageGrayu8 in,int r)
{
    ImageGrayu8 out(in.width(),in.height());

    out.parallel_for_all_pixels([&](ImageGrayu8::PixelType &p,int x,int y)
    {
       float tx = square_to_triangle_gridX(r,x,y);
       float ty = square_to_triangle_gridY(r,y);

       float p0x = std::floor(tx);
       float p0y = std::floor(ty);

       float p1x = p0x + 1;
       float p1y = p0y;

       float p2x = p0x;
       float p2y = p0y + 1;

       float cp = (p2y - p1y) * (tx - p1x) - (p2x - p1x) * (ty - p1y);

       if(cp >= 0)
       {
           p0x = p1x;
           p0y = p2y;
       }

       float t0x = tx - p0x;
       float t0y = ty - p0y;

       float t1x = tx - p1x;
       float t1y = ty - p1y;

       float t2x = tx - p2x;
       float t2y = ty - p2y;

       std::srand(triangle_to_square_gridX(r,p0x,p0y));
       int r0x = std::rand() % in.width();
       std::srand(triangle_to_square_gridY(r,p0y));
       int r0y = std::rand() % in.height();
       float c0x = square_to_triangle_gridX(r,r0x,r0y);
       float c0y = square_to_triangle_gridY(r,r0y);

       std::srand(triangle_to_square_gridX(r,p1x,p1y));
       int r1x = std::rand() % in.width();
       std::srand(triangle_to_square_gridY(r,p1y));
       int r1y = std::rand() % in.height();
       float c1x = square_to_triangle_gridX(r,r1x,r1y);
       float c1y = square_to_triangle_gridY(r,r1y);

       std::srand(triangle_to_square_gridX(r,p2x,p2y));
       int r2x = std::rand() % in.width();
       std::srand(triangle_to_square_gridY(r,p2y));
       int r2y = std::rand() % in.height();
       float c2x = square_to_triangle_gridX(r,r2x,r2y);
       float c2y = square_to_triangle_gridY(r,r2y);

       int in_p0x = (triangle_to_square_gridX(r,c0x+t0x,c0y+t0y) + in.width()) % in.width();
       int in_p0y = (triangle_to_square_gridY(r,c0y+t0y) + in.height()) % in.height();

       int in_p1x = (triangle_to_square_gridX(r,c1x+t1x,c1y+t1y) + in.width()) % in.width();
       int in_p1y = (triangle_to_square_gridY(r,c1y+t1y) + in.height()) % in.height();

       int in_p2x = (triangle_to_square_gridX(r,c2x+t2x,c2y+t2y) + in.width()) % in.width();
       int in_p2y = (triangle_to_square_gridY(r,c2y+t2y) + in.height()) % in.height();

       int color =
               (in.pixelAbsolute(in_p0x,in_p0y) +
               in.pixelAbsolute(in_p1x,in_p1y) +
               in.pixelAbsolute(in_p2x,in_p2y))/3;

       p = ImageGrayu8::itkPixel(color);

    });

    return out;
}



int main()
{
    //int r = 500;

    ImageGrayu8 in;
    in.load("../../textures/tex.jpeg");
    //HistogramGrayu8 h(in);
    //h.saveFullHistogram("h_in");
    ImageGrayu8 out = pavage_bruit(in,in.width()/2);
    out.save("out.png");
    //HistogramGrayu8 ho(out);
    //ho.saveFullHistogram("h_out");

    ImageRGBu8 in2;
    in2.load("../../textures/noise.png");
    //HistogramRGBu8 h2(in2);
    //h2.saveFullHistogram("h2_in");
    ImageRGBu8 out2 = pavage_bruit(in2,in2.width()/2);
    out2.save("out2.png");
    //HistogramRGBu8 ho2(out2);
    //ho2.saveFullHistogram("h2_out");

    ImageRGBu8 in3;
    in3.load("../../textures/Moss0101_S.png");
    //HistogramRGBu8 h3(in3);
    //h3.saveFullHistogram("h3_in");
    ImageRGBu8 out3 = pavage_bruit(in3,in3.width()/2);
    out3.save("out3.png");
    //HistogramRGBu8 ho3(out3);
    //ho3.saveFullHistogram("h3_out");

    return 0;
}
