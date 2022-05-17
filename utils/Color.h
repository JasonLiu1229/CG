//
// Created by liuja on 02/03/2022.
//

#ifndef INC_2D_L_SYSTEMEN_COLOR_H
#define INC_2D_L_SYSTEMEN_COLOR_H


class Color {
    double red = 0;
    double green = 0;
    double blue = 0;

public:
    double getRed() const;

    void setRed(double red);

    double getGreen() const;

    void setGreen(double green);

    double getBlue() const;

    void setBlue(double blue);

    void setRGB(double red, double green, double blue);

    Color operator*(const Color &color);

    Color operator+(const Color &color);

    void operator+=(const Color &color);

    Color();

    Color(double red, double green, double blue);

    virtual ~Color();
};


#endif //INC_2D_L_SYSTEMEN_COLOR_H
