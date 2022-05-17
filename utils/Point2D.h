//
// Created by liuja on 02/03/2022.
//

#ifndef INC_2D_L_SYSTEMEN_POINT2D_H
#define INC_2D_L_SYSTEMEN_POINT2D_H
#include <vector>

class Point2D {
    double x;
    double y;
    double z;

public:
    double getZ() const;

    void setZ(double z);

    virtual ~Point2D();

    Point2D();

    void setCoords(double x, double y);

    double getX() const;

    void setX(double x);

    double getY() const;

    void setY(double y);
};

typedef std::vector<Point2D> Points2D;
#endif //INC_2D_L_SYSTEMEN_POINT2D_H
