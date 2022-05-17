//
// Created by liuja on 30/03/2022.
//

#ifndef LINE2D_CPP_ZBUFFER_H
#define LINE2D_CPP_ZBUFFER_H
#include <vector>
#include <limits>

const double postInf = std::numeric_limits<double>::infinity();
const double negInf = -std::numeric_limits<double>::infinity();

class ZBuffer:public std::vector<std::vector<double>> {
public:
    ZBuffer(const int width, const int height);
};


#endif //LINE2D_CPP_ZBUFFER_H
