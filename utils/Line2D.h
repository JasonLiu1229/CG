//
// Created by liuja on 02/03/2022.
//

#ifndef INC_2D_L_SYSTEMEN_LINE2D_H
#define INC_2D_L_SYSTEMEN_LINE2D_H

#include "Point2D.h"
#include "Color.h"
#include <iostream>
#include <list>

class Line2D : public std::error_code {
public:
    Point2D p1;
    Point2D p2;
    Color rgb;

    double z1;
    double z2;

    Line2D();
    virtual ~Line2D();
};

using Lines2D = std::list<Line2D>;

#endif //INC_2D_L_SYSTEMEN_LINE2D_H
