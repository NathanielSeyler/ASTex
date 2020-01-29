#ifndef SIMPLEX_H
#define SIMPLEX_H

//#include <cmath>
#include <Eigen/Eigen>
#include <vector>
#include <ASTex/image_gray.h>

namespace ASTex {

template<int N>
class Simplex{

public:
    using Vec = Eigen::Matrix<double,N,1>;

    Simplex(){}

    static inline Vec skew(Vec &p){
        assert(N > 0);
        double f((std::sqrt(N + 1) - 1) / N);
        double s(0);

        for(int i =0; i < N; i++)
            s += p(i);
        s *= f;

        Vec p_prime;
        for(int i =0; i < N ; i++)
            p_prime(i) = p(i) + s;

        return p_prime;
    }

    static inline Vec unskew(Vec &p_prime){
        assert(N > 0);
        double g( (N + 1. - std::sqrt(N + 1.)) / (N * (N + 1.)));
        double s(0);
        for(int i = 0; i< N ;i++)
            s+= p_prime(i);
        s *= g;

        Vec p;
        for(int i =0 ; i < N ;i++)
            p(i) = p_prime(i) - s;

        return p;
    }

    static inline std::vector<Vec> vertices(Vec &p,ImageGrayu8 &im ,int x,int y){
        assert(N > 0);
        Vec p_prime = skew(p);

        Vec p0;
        for(int i = 0;i <N ; i++)
            p0[i] = std::floor(p_prime[i]);

        Vec pi = p_prime - p0;
        std::vector<Vec> v;
        v.push_back(p0);

        Vec p_next = p0;
        std::multimap<double,unsigned int> map;
        for(int i = 0;i <N ;i++)
            map.insert(std::pair<double,unsigned int>(pi(i),i));

        for(auto it =map.rbegin();it != map.rend();it++)
        {
            p_next(it->second)++;
            v.push_back(p_next);
        }

        if(pi(0)>pi(1))
            im.pixelAbsolute(x,y)= 0;
        else
            im.pixelAbsolute(x,y)= 255;

        return v;
    }

    static inline double contribution(Vec &d)
    {
        double dd = d.squaredNorm();

        double w = 0.5 - dd;
        if(w>0)
            w = 8 * (w * w * w * w);
        else
            w = 0;

        return w;
    }

};

using Simplex2D = Simplex<2>;
using Vec2 = Simplex2D::Vec;

}
#endif
