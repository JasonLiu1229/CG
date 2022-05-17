//
// Created by liuja on 30/03/2022.
//

#include "ZBuffer.h"


ZBuffer::ZBuffer(const int width, const int height) {
    std::vector<double> test;
    test.resize(height, postInf);
    this->resize(width, test);
    for (int i = 0; i < width; ++i) {
        this->push_back(std::vector<double>(height, postInf));
    }
}
