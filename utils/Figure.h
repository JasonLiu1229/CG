//
// Created by liuja on 09/03/2022.
//

#ifndef COLOR_CPP_FIGURE_H
#define COLOR_CPP_FIGURE_H

#include <vector>
#include <list>

#include "Color.h"
#include "vector3d.h"
#include "Face.h"

class Figure {
public:
    Figure();

    virtual ~Figure();

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    double scale;

    double rotationX;
    double rotationY;
    double rotationZ;

    std::vector<double> center;

    std::string name;

    Color ambientReflection;
    Color diffuseReflection;
    Color specularReflection;
    double reflectionCoefficient;
};

typedef std::list<Figure> Figures3D;

#endif //COLOR_CPP_FIGURE_H
