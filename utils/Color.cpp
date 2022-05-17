//
// Created by liuja on 02/03/2022.
//

#include "Color.h"

double doubleModulo(double value, double moduloVal){
    while (value > moduloVal) {
        value -= moduloVal;
    }
    return value;
}

double Color::getRed() const {
    return red;
}

void Color::setRed(double red) {
    Color::red = red;
}

double Color::getGreen() const {
    return green;
}

void Color::setGreen(double green) {
    Color::green = green;
}

double Color::getBlue() const {
    return blue;
}

void Color::setBlue(double blue) {
    Color::blue = blue;
}

void Color::setRGB(double red, double green, double blue) {
    this->red = red;
    this->green = green;
    this->blue = blue;
}

Color::Color() {}

Color::~Color() {}

Color::Color(double red, double green, double blue) : red(red), green(green), blue(blue) {}

Color Color::operator*(const Color &color) {
    Color newColor;
    newColor.setRed((this->red * color.red));
    newColor.setBlue((this->blue * color.blue));
    newColor.setGreen((this->green * color.green));
    return newColor;
}

Color Color::operator+(const Color &color) {
    Color newColor;
    newColor.setRed((this->red + color.red));
    newColor.setBlue((this->blue + color.blue));
    newColor.setGreen((this->green + color.green));
    return newColor;
}

void Color::operator+=(const Color &color) {
    *(this) = *(this) + color;
}
