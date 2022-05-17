//
// Created by liuja on 09/03/2022.
//

#ifndef COLOR_CPP_FACE_H
#define COLOR_CPP_FACE_H
#include <vector>

class Face {
public:
    Face();

    virtual ~Face();

    std::vector<int> point_indexes;

    bool operator == (const Face &face);
};


#endif //COLOR_CPP_FACE_H
