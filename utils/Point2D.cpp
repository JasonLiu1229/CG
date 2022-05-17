//
// Created by liuja on 02/03/2022.
//

#include "Point2D.h"

double Point2D::getX() const {
    return x;
}

void Point2D::setX(double x) {
    Point2D::x = x;
}

double Point2D::getY() const {
    return y;
}

void Point2D::setY(double y) {
    Point2D::y = y;
}

void Point2D::setCoords(double x, double y) {
    Point2D::x = x;
    Point2D::y = y;
}

Point2D::Point2D() {}

Point2D::~Point2D() {}

double Point2D::getZ() const {
    return z;
}

void Point2D::setZ(double z) {
    Point2D::z = z;
}
