//==== Custom Libs ====//
#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "Line2D.h"
#include "vector3d.h"
#include "Figure.h"
#include "ZBuffer.h"
#include "Light.h"
#include "InfLight.h"
#include "PointLight.h"

//==== Build In Libs ====//
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

//==== Const Variables ====//
    //= math =//
const double pi = 3.14159265358979323846; // pi representation
const std::string pointS = "point";
const std::string lineS = "line";
const std::string figureS = "Figure";

    //= chars =//
const char plus = '+';
const char minus = '-';
const char startBracket = '(';
const char endBracket = ')';
const char ampersand = '&';
const char power = '^';
const char backslash = '\\';
const char slash = '/';
const char collon = '|';

//==== Extra Functions ====//
double doubleModuloE(double value, double moduloVal){
    while (value > moduloVal) {
        value -= moduloVal;
    }
    return value;
}

std::vector<double> matrixMult(std::vector<double> mat1, Matrix mat2){
    std::vector<double> newMat = {0,0,0,0};

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            newMat[i] += mat1[j] * mat2(j + 1,i + 1);
        }
    }
    return newMat;
}

Vector3D vecTo3D(std::vector<double> vect){
    Vector3D newVec;
    newVec.x = vect[0];
    newVec.y = vect[1];
    newVec.z = vect[2];
    return newVec;
}

void toPolar(const Vector3D &point, double &theta, double &phi, double &r){
    r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    theta = std::atan2(point.y, point.x);
    phi = std::acos(point.z/r);
}

Matrix scaleFigure(const double scale){
    Matrix scaleMat;

    for (int i = 0; i < 3; ++i) {
        scaleMat(i+1, i+1) = scale;
    }
    scaleMat(4,4) = 1;

    return scaleMat;
}

Matrix rotXMat(const double angle){
    Matrix rotateX;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j){
                rotateX(i + 1,j + 1)  = cos(angle);
            }
            else if (i == 1 and j == 2){
                rotateX(j + 1, i + 1) = -sin(angle);
            }
            else if (i == 2 and j == 1){
                rotateX(j + 1, i + 1) = sin(angle);
            }
        }
    }

    rotateX(1,1) = 1;
    rotateX(4, 4) = 1;

    return rotateX;
}

Matrix rotYMat(const double angle){
    Matrix rotateY;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j){
                rotateY(j + 1,i + 1)  = cos(angle);
            }
            else if (i == 2 and j == 0){
                rotateY(j + 1, i + 1) = -sin(angle);
            }
            else if (i == 0 and j == 2){
                rotateY(j + 1, i + 1) = sin(angle);
            }
        }
    }

    rotateY(2,2) = 1;
    rotateY(4,4) = 1;

    return rotateY;
}

Matrix rotZMat(const double angle){
    Matrix rotateZ;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j){
                rotateZ(j + 1,i + 1)  = cos(angle);
            }
            else if (i == 1 and j == 0){
                rotateZ(j + 1, i + 1) = sin(angle);
            }
            else if (i == 0 and j == 1){
                rotateZ(j + 1, i + 1) = -sin(angle);
            }
        }
    }

    rotateZ(3,3) = 1;
    rotateZ(4,4) = 1;

    return rotateZ;
}

Matrix translate(const Vector3D &vector){
    Matrix trans;

    for (int i = 0; i < 4; ++i) {
        trans(i+1, i+1) = 1;
    }

    trans(4,1) = vector.x;
    trans(4, 2) = vector.y;
    trans(4, 3) = vector.z;

    return trans;
}

Matrix eyeTransformatieMatrix(const double &theta, const double &phi, const double &r){
    Matrix V;

    // row 1
    V(1,1) = -sin(theta);
    V(1,2) = -cos(theta) * cos(phi);
    V(1,3) = cos(theta) * sin(phi);
    V(1,4) = 0;

    // row 2
    V(2, 1) = cos(theta);
    V(2, 2) = -sin(theta) * cos(phi);
    V(2, 3) = sin(theta) * sin(phi);
    V(2, 4) = 0;

    // row 3
    V(3, 1) = 0;
    V(3, 2) = sin(phi);
    V(3, 3) = cos(phi);
    V(3, 4) = 0;

    // row 4
    V(4, 1) = 0;
    V(4, 2) = 0;
    V(4, 3) = -r;
    V(4, 4) = 1;

    return V;
}

void applyTransformationFig(Figure &fig, const Matrix &m){
    for (auto &point : fig.points) {
        point *= m;
    }
}

void applyTransformationFigs(Figures3D &figs, const Matrix &m){
    for(auto &fig : figs){
        applyTransformationFig(fig, m);
    }
}

std::pair<double, double> determine_max(const Lines2D &lines){
    double x1 = lines.begin()->p1.getX();
    double x2 = lines.begin()->p2.getX();

    double xmax = 0;

    if(x1 >= x2){
        xmax = x1;
    }
    else{
        xmax = x2;
    }

    double y1 = lines.begin()->p1.getY();
    double y2 = lines.begin()->p2.getY();

    double ymax = 0;

    if(y1 > y2){
        ymax = y1;
    }
    else{
        ymax = y2;
    }

    for(auto it : lines){
        x1 = it.p1.getX();
        x2 = it.p2.getX();

        y1 = it.p1.getY();
        y2 = it.p2.getY();

        if(x1 >= x2 and x1 >= xmax){
            xmax = x1;
        }
        else if(x2 > xmax and x2 > x1){
            xmax = x2;
        }

        if(y1 >= y2 and y1 >= ymax){
            ymax = y1;
        }
        else if(y2 > ymax and y2 > y1){
            ymax = y2;
        }
    }
    std::pair<double, double> max;
    max.first = xmax;
    max.second = ymax;

    return max;
}

std::pair<double, double> determine_min(const Lines2D &lines){
    double x1 = lines.begin()->p1.getX();
    double x2 = lines.begin()->p2.getX();

    double xmin;

    if(x1 <= x2){
        xmin = x1;
    }
    else{
        xmin = x2;
    }

    double y1 = lines.begin()->p1.getY();
    double y2 = lines.begin()->p2.getY();

    double ymin;

    if(y1 <= y2){
        ymin = y1;
    }
    else{
        ymin = y2;
    }
    for (auto it:lines) {
        x1 = it.p1.getX();
        x2 = it.p2.getX();

        y1 = it.p1.getY();
        y2 = it.p2.getY();

        if(x1 < xmin){
            xmin = x1;
        }
        if(x2 < xmin){
            xmin = x2;
        }

        if(y1 < ymin){
            ymin = y1;
        }
        if(y2 < ymin)
            ymin = y2;
    }
    std::pair<double, double> min;

    min.first = xmin;
    min.second = ymin;

    return min;
}

double degreesToRadians(double degrees){
    return degrees *(pi/180);
}

void calculateNextCoords(double& curX, double& curY, const double& curAngle){
    const double curAngleRad = degreesToRadians(curAngle);
    curX = curX + cos(curAngleRad);
    curY = curY + sin(curAngleRad);
}

std::pair<std::string, int> extractString(const std::string& source, int startInt, char secondDelimiter){
    std::size_t startP = startInt;
    startP++;
    std::size_t endP = 0;
    int count = 0;
    for (int i = startP; i < source.size(); i++){
        if(source[i] == secondDelimiter and count == 0){
            endP = i;
            break;
        }
        else if(source[i] == startBracket){
            count++;
        }
        else if(source[i] == secondDelimiter){
            count--;
        }
    }
    return {source.substr(startP, endP - startP), endP};
}
//==== Create Figure ====//
Figure createCube(){
    const int PointsT[3][8] = {{1,-1, 1, -1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, -1, -1, 1}, {-1, -1, 1, 1, -1, -1, 1, 1}};
    const int FacesT[4][6] = {{1, 5, 2, 6, 7, 1}, {5, 2, 6, 1, 3, 6}, {3, 8, 4, 7, 8, 2}, {7, 3, 8, 4, 4, 5}};

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    Face face;
    face.point_indexes = {0,0,0,0};
    Vector3D point;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 3; ++j) {
            switch (j){
                case 0:
                    point.x = PointsT[j][i];
                    break;
                case 1:
                    point.y = PointsT[j][i];
                    break;
                case 2:
                    int test = PointsT[j][i];
                    point.z = PointsT[j][i];
                    break;
            }
        }
        points.push_back(point);
    }

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            face.point_indexes[j] = FacesT[j][i] - 1;
        }
        faces.push_back(face);
    }

    Figure figureO;
    figureO.faces = faces;
    figureO.points = points;
    return figureO;
}

Figure createTetraHedron(){
    const int PointsT[3][4] = {{1,-1,1,-1}, {-1, 1, 1, -1}, {-1, -1, 1, 1}};
    const int FacesT[3][4] = {{1,2,1,1}, {2,4,4,3}, {3,3,2,4}};

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    Face face;
    face.point_indexes = {0,0,0};
    Vector3D point;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            switch (j){
                case 0:
                    point.x = PointsT[j][i];
                    break;
                case 1:
                    point.y = PointsT[j][i];
                    break;
                case 2:
                    point.z = PointsT[j][i];
                    break;
            }
        }
        points.push_back(point);
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            face.point_indexes[j] = FacesT[j][i] - 1;
        }
        faces.push_back(face);
    }

    Figure figureO;
    figureO.faces = faces;
    figureO.points = points;
    return figureO;
}

Figure createOctahedron(){
    const int PointsT[3][6] = {{1, 0, -1, 0, 0, 0}, {0, 1, 0, -1, 0, 0}, {0, 0, 0, 0, -1, 1}};
    const int FacesT[3][8] = {{1,2,3,4,2,3,4,1}, {2, 3, 4, 1, 1, 2, 3, 4}, {6, 6, 6, 6, 5, 5, 5, 5}};

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    Face face;
    face.point_indexes = {0,0,0};
    Vector3D point;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            switch (j){
                case 0:
                    point.x = PointsT[j][i];
                    break;
                case 1:
                    point.y = PointsT[j][i];
                    break;
                case 2:
                    point.z = PointsT[j][i];
                    break;
            }
        }
        points.push_back(point);
    }

    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 3; ++j) {
            face.point_indexes[j] = FacesT[j][i] - 1;
        }
        faces.push_back(face);
    }

    Figure figureO;
    figureO.faces = faces;
    figureO.points = points;
    return figureO;
}

Figure createIcosahedron(){
    const int FacesT[3][20] = {{1,1,1,1,1,2,3,3,4,4,5,5,6,6,2,12,12,12,12,12}, {2,3,4,5,6,7,7,8,8,9,9,10,10,11,11,8,9,10,11,7}, {3,4,5,6,2,3,8,4,9,5,10,6,11,2,7,7,8,9,10,11}};

    double PointsT[3][12] = {};

    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 3; ++j) {
            if(i == 0){
                if(j == 0 or j == 1){
                    PointsT[j][i] = 0;
                }
                else{
                    PointsT[j][i] = sqrt(5)/2;
                }
            }
            else if(i == 2 or i == 3 or i == 4 or i == 5 or i == 1){
                if (j == 0){
                    PointsT[j][i] = cos((i - 1) * ((2 * pi)/5));
                }
                else if (j == 1){
                    PointsT[j][i] = sin((i - 1) * ((2 * pi)/5));
                }
                else{
                    PointsT[j][i] = 0.5;
                }
            }
            else if (i == 6 or i == 7 or i == 8 or i == 9 or i == 10){
                if (j == 0){
                    PointsT[j][i] = cos((pi/5) + (i - 6) * ((2*pi)/5));
                }
                else if (j == 1){
                    PointsT[j][i] = sin((pi/5) + (i - 6) * ((2*pi)/5));
                }
                else{
                    PointsT[j][i] = -0.5;
                }
            }
            else if (i == 11){
                if (j == 0 or j == 1){
                    PointsT[j][i] = 0;
                }
                else{
                    PointsT[j][i] = -sqrt(5)/2;
                }
            }
        }
    }

    Face face;
    face.point_indexes = {0,0,0};
    Vector3D point;

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 3; ++j) {
            switch (j){
                case 0:
                    point.x = PointsT[j][i];
                    break;
                case 1:
                    point.y = PointsT[j][i];
                    break;
                case 2:
                    point.z = PointsT[j][i];
                    break;
            }
        }
        points.push_back(point);
    }

    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 3; ++j) {
            face.point_indexes[j] = FacesT[j][i] - 1;
        }
        faces.push_back(face);
    }

    Figure figureO;
    figureO.faces = faces;
    figureO.points = points;
    return figureO;
}

Figure createDodecahedron(){
    Figure ico = createIcosahedron();

    Face face;
    face.point_indexes = {0,0,0,0,0};
    Vector3D point;

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    int FacesT[5][12] = {{1, 1, 2,  3,  4,  5,  20, 20, 19, 18, 17, 16},
                         {2, 6, 8,  10, 12, 14, 19, 15, 13, 11, 9,  7},
                         {3, 7, 9,  11, 13, 15, 18, 14, 12, 10, 8,  6},
                         {4, 8, 10, 12, 14, 6,  17, 13, 11, 9,  7,  15},
                         {5, 2, 3,  4,  5,  1,  16, 19, 18, 17, 16, 20}};


    for (int i = 0; i < 20; ++i) {
        Vector3D temp;
        for (int j = 0; j < 3; ++j){
            int pointI = ico.faces.at(i).point_indexes.at(j);
            temp += ico.points[pointI];
        }
        temp /= 3;
        points.push_back(temp);
    }

    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 5; ++j) {
            face.point_indexes[j] = FacesT[j][i] - 1;
        }
        faces.push_back(face);
    }

    Figure figureO;
    figureO.faces = faces;
    figureO.points = points;
    return figureO;
}

Figure createSphere(const int n){
    Figure icosahedron = createIcosahedron();

    std::vector<Vector3D> points = icosahedron.points;
    Vector3D point;
    Face face1;
    face1.point_indexes = {0,0,0};
    Face face2;
    face2.point_indexes = {0,0,0};
    Face face3;
    face3.point_indexes = {0,0,0};
    Face face4;
    face4.point_indexes = {0,0,0};

    Vector3D A;
    Vector3D B;
    Vector3D C;
    Vector3D D;
    Vector3D E;
    Vector3D F;

    int Dloc = 0;
    int Eloc = 0;
    int Floc = 0;

    for (int i = 0; i < n; ++i) {
        std::vector<Face> faces;
        for (auto faceIt:icosahedron.faces) {
            A = icosahedron.points[faceIt.point_indexes[0]];
            B = icosahedron.points[faceIt.point_indexes[1]];
            C = icosahedron.points[faceIt.point_indexes[2]];

            D = (A + B) / 2;
            E = (A + C) / 2;
            F = (C + B) / 2;

            icosahedron.points.push_back(D);
            Dloc = icosahedron.points.size()-1;
            icosahedron.points.push_back(E);
            Eloc = icosahedron.points.size()-1;
            icosahedron.points.push_back(F);
            Floc = icosahedron.points.size()-1;

            face1.point_indexes[0] = faceIt.point_indexes[0];
            face1.point_indexes[1] = Dloc;
            face1.point_indexes[2] = Eloc;

            face2.point_indexes[0] = faceIt.point_indexes[1];
            face2.point_indexes[1] = Floc;
            face2.point_indexes[2] = Dloc;

            face3.point_indexes[0] = faceIt.point_indexes[2];
            face3.point_indexes[1] = Eloc;
            face3.point_indexes[2] = Floc;

            face4.point_indexes[0] = Dloc;
            face4.point_indexes[1] = Floc;
            face4.point_indexes[2] = Eloc;

            faces.push_back(face1);
            faces.push_back(face2);
            faces.push_back(face3);
            faces.push_back(face4);
        }
        icosahedron.faces = faces;
    }

    // rescale points
    for (auto &point: icosahedron.points) {
        point.normalise();
    }
    return icosahedron;
}

Figure createCone(const int n, const int h){
    Face face;
    face.point_indexes.resize(n, 0);
    Face face2;
    face2.point_indexes.resize(3, 0);
    Vector3D point;

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    for (int i = 0; i < n; ++i) {
        point.x = cos((2 * i * pi) / n);
        point.y = sin((2 * i * pi) / n);
        point.z = 0;
        points.push_back(point);
    }
    point.x = 0;
    point.y = 0;
    point.z = h;
    points.push_back(point);

    for (int i = 0; i < n; ++i) {
        face2.point_indexes[0] = i;
        face2.point_indexes[1] = (i+1) % n;
        face2.point_indexes[2] = n;
        faces.push_back(face2);
    }

    for (int i = 0; i < n; ++i) {
        face.point_indexes[i] = n - 1 - i;
    }
    faces.push_back(face);

    Figure figureO;
    figureO.faces = faces;
    figureO.points = points;
    return figureO;
}

Figure createCylinder(const int n, const double h){
    Face face;
    face.point_indexes.resize(n, 0);
    Face face2;
    face2.point_indexes.resize(4, 0);
    Vector3D point;

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    for (int i = 0; i < n; ++i) {
        point.x = cos((2 * i * pi) / n);
        point.y = sin((2 * i * pi) / n);
        point.z = 0;
        points.push_back(point);
    }

    for (int i = 0; i < n; ++i) {
        point.x = cos((2 * i * pi) / n);
        point.y = sin((2 * i * pi) / n);
        point.z = h;
        points.push_back(point);
    }

    for (int i = 0; i < n; ++i) {
        face2.point_indexes[0] = i;
        face2.point_indexes[1] = (i + 1) % n;
        face2.point_indexes[2] = ((i + 1) % n) + n;
        face2.point_indexes[3] = (n + i);
        faces.push_back(face2);
    }

    for (int i = 0; i < n; ++i) {
        face.point_indexes[i] = n - 1 - i;
    }
    faces.push_back(face);

    for (int i = 0; i < n; ++i) {
        face.point_indexes[i] = n - 1 - i + n;
    }
    faces.push_back(face);

    Figure figureO;
    figureO.faces = faces;
    figureO.points = points;
    return figureO;
}

Figure createTorus(const double r, const double R, const int n, const int m){
    Figure figureO;
    Vector3D point;
    std::pair<int, int> loc;
    Face face;
    face.point_indexes = {0,0,0,0};
    std::map<std::pair<int, int>, int> pointLocations;

    double u = 0;
    double v = 0;

    double xUV = 0;
    double yUV = 0;
    double zUV = 0;

    for (int i = 0; i < n; ++i) {
        u = (2 * pi * i)/n;
        loc.first = i;
        for (int j = 0; j < m; ++j) {
            v = (2 * pi * j)/m;
            xUV = (R + r* cos(v))* cos(u);
            yUV = (R + r* cos(v))* sin(u);
            zUV = r* sin(v);

            point.x = xUV;
            point.y = yUV;
            point.z = zUV;

            figureO.points.push_back(point);
            loc.second = j;
            pointLocations[loc] = figureO.points.size() - 1;
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            loc.first = i;
            loc.second = j;
            face.point_indexes[0] = pointLocations[loc];

            loc.first = (i + 1) % n;
            face.point_indexes[1] = pointLocations[loc];

            loc.second = (j + 1) % m;
            face.point_indexes[2] = pointLocations[loc];

            loc.first = i;
            face.point_indexes[3] = pointLocations[loc];

            figureO.faces.push_back(face);
        }
    }
    return figureO;
}

int checkIfPointExist(std::vector<Vector3D> points, Vector3D point){
    for (int i = 0; i < points.size(); ++i) {
        if (point == points[i]){
            return i;
        }
    }
    return points.size();
}

bool faceInit(const Face &face, int index){
    for (const auto indexF : face.point_indexes){
        if (indexF == index){
            return true;
        }
    }
    return false;
}

Figure createBuckyBall(){
    Figure buckyBall = createIcosahedron();
    Face temp3face;
    temp3face.point_indexes = {0,0,0};

    Face temp6face;
    temp6face.point_indexes = {0,0,0,0,0,0};

    Face temp5face;
    temp5face.point_indexes = {0,0,0,0,0};

    Vector3D A;
    Vector3D B;
    Vector3D C;
    Vector3D D;
    Vector3D E;
    Vector3D F;
    Vector3D G;
    Vector3D H;
    Vector3D I;

    Figure tempFig;

    std::vector<Face> faces;
    std::vector<Face> faces3;

    std::map<int, std::vector<int>> connections;
    std::vector<Vector3D> points;
    for (auto &face : buckyBall.faces) {
        A = buckyBall.points[face.point_indexes[0]];
        B = buckyBall.points[face.point_indexes[1]];
        C = buckyBall.points[face.point_indexes[2]];

        connections[face.point_indexes[0]];
        connections[face.point_indexes[1]];
        connections[face.point_indexes[2]];

        int locD = 0;
        int locE = 0;
        int locF = 0;
        int locG = 0;
        int locH = 0;
        int locI = 0;

        // D en E --> AB
        D = (2*A)/3 + B/3;
        E = (2*B)/3 + A/3;

        locD = checkIfPointExist(points, D);
        if (locD == points.size()){
            points.push_back(D);
            locD = points.size()-1;
        }
        locE = checkIfPointExist(points, E);
        if (locE == points.size()){
            points.push_back(E);
            locE = points.size()-1;
        }

        // F en G --> BC
        F = (2*B)/3 + C/3;
        G = (2*C)/3 + B/3;
        locF = checkIfPointExist(points, F);
        if (locF == points.size()){
            points.push_back(F);
            locF = points.size()-1;
        }

        locG = checkIfPointExist(points, G);
        if (locG == points.size()){
            points.push_back(G);
            locG = points.size()-1;
        }

        // H en I --> CA
        H = (2*C)/3 + A/3;
        I = (2*A)/3 + C/3;
        locH = checkIfPointExist(points, H);
        if (locH == points.size()){
            points.push_back(H);
            locH = points.size()-1;
        }

        locI = checkIfPointExist(points, I);
        if (locI == points.size()){
            points.push_back(I);
            locI = points.size()-1;
        }

        // face CHG
        temp3face.point_indexes = {face.point_indexes[2], locH, locG};
        faces3.push_back(temp3face);
        connections[face.point_indexes[2]].push_back(locH);
        connections[face.point_indexes[2]].push_back(locG);

        // face ADI
        temp3face.point_indexes = {face.point_indexes[0], locD, locI};
        faces3.push_back(temp3face);
        connections[face.point_indexes[0]].push_back(locD);
        connections[face.point_indexes[0]].push_back(locI);

        // face BFE
        temp3face.point_indexes = {face.point_indexes[1], locF, locE};
        faces3.push_back(temp3face);
        connections[face.point_indexes[1]].push_back(locF);
        connections[face.point_indexes[1]].push_back(locE);

        // face IDEFGH
        temp6face.point_indexes = {locI, locD, locE, locF, locG, locH};
        faces.push_back(temp6face);
    }

    // create pentagon
    for (const auto& pentagon10 : connections){
        std::vector<Vector3D> pentagon;
        std::vector<int> pentagonIndex;
        for (const auto& point : pentagon10.second) {
            bool exist = false;
            for (const auto& pointV : pentagon){
                if (pointV == points[point]){
                    exist = true;
                }
            }
            if (!exist){
                pentagon.push_back(points[point]);
                pentagonIndex.push_back(point);
            }
        }

        Face newFace;
        Vector3D firstPoint = pentagon[0];
        Vector3D nextPoint;
        newFace.point_indexes.push_back(pentagonIndex[0]);
        //newFace.point_indexes.push_back(pentagonIndex[1]);
        int nextIndex = 0;
        Vector3D smallest;
        int smallestIndex = 0;
        for (int i = 0; i < 5; ++i) {
            double nearest = postInf;
            for (int j = 1; j < pentagon.size(); ++j) {
                nextPoint = pentagon[j];
                nextIndex = pentagonIndex[j];
                if (nearest > (firstPoint - nextPoint).length() and !faceInit(newFace, nextIndex)){
                    nearest = (firstPoint - nextPoint).length();
                    smallestIndex = nextIndex;
                    smallest = nextPoint;
                }
            }
            if (!faceInit(newFace, smallestIndex)){
                newFace.point_indexes.push_back(smallestIndex);
                firstPoint = smallest;
            }
        }
        faces.push_back(newFace);
    }

    buckyBall.faces = faces;
    buckyBall.points = points;

    for (auto &point : buckyBall.points) {
        point.normalise();
    }

    return buckyBall;
}

//==== Generating Image Functions ====//
img::EasyImage draw2DLines(const Lines2D &lines, const int size, const std::vector<double>& bgClr){
    // Determine min and max of x and y
    std::pair<double,double> max = determine_max(lines);
    std::pair<double,double> min = determine_min(lines);

    double Xmax = max.first;
    double Xmin = min.first;

    double Ymax = max.second;
    double Ymin = min.second;

    double xrange = Xmax - Xmin;
    double yrange = Ymax - Ymin;

    // calculate image size
    double imageX = size * (xrange / std::max(xrange, yrange)); // width
    double imageY = size * (yrange / std::max(xrange, yrange)); // height

    // determine scale factor
    double d = 0.95 * (imageX/xrange);

    // determine Shift amount
    double DCx = d*((Xmin + Xmax)/2);
    double DCy = d*((Ymin+Ymax)/2);

    double dx = (imageX/2)-DCx;
    double dy = (imageY/2)-DCy;

    // Color
    img::Color color;
    img::Color bgColor;

    bgColor.red = bgClr[0] * 255;
    bgColor.green = bgClr[1] * 255;
    bgColor.blue = bgClr[2] * 255;

    // generate image
    img::EasyImage image(imageX, imageY, bgColor);


    for(const auto& it: lines){
        color.red = it.rgb.getRed() * 255;
        color.green = it.rgb.getGreen() * 255;
        color.blue = it.rgb.getBlue() * 255;

        image.draw_line(lround(it.p1.getX()*d + dx), lround(it.p1.getY()*d + dy), lround(it.p2.getX()*d + dx),
                        lround(it.p2.getY()*d + dy), color);
    }
    return image;
}

Line2D calculateLine(double& curX, double& curY, double& curAngle, const std::vector<double>& bgClr){
    Line2D line;
    line.rgb.setRGB(bgClr[0], bgClr[1], bgClr[2]);
    line.p1.setCoords(curX, curY);
    calculateNextCoords(curX, curY, curAngle);
    line.p2.setCoords(curX, curY);
    return line;
}

Lines2D recursionLinesFinder(const std::string& source, const LParser::LSystem2D &l_system, const std::vector<double>& bgClr, double currentAngleDegrees, double curX, double curY){
    Lines2D lines;
    Lines2D tempLines;
    Line2D line;

    double angle = l_system.get_angle();

    double tempX;
    double tempY;
    double tempAngle;

    const std::set<char>& alphabet = l_system.get_alphabet();

    // generate lines
    for (int i = 0; i < source.size(); ++i) {
        if(source[i] == plus){
            currentAngleDegrees += angle;
        }
        else if(source[i] == minus){
            currentAngleDegrees -= angle;
        }
        else if(source[i] == startBracket) {
            std::pair<std::string, int> subStringPair = extractString(source, i, endBracket);
            std::string subString = subStringPair.first;

            tempX = curX;
            tempY = curY;
            tempAngle = currentAngleDegrees;

            tempLines = recursionLinesFinder(subString, l_system, bgClr, tempAngle, tempX, tempY);
            lines.merge(tempLines);

            i = subStringPair.second;
        }
        else{
            line = calculateLine(curX, curY, currentAngleDegrees, bgClr);
            if (l_system.draw(source[i])) {
                lines.push_back(line);
            }
        }
    }

    return lines;
}

Lines2D drawLsystem(const LParser::LSystem2D &l_system , const std::vector<double>& bgClr) {
    Lines2D lines;
    // l_system config
    const std::set<char> &alphabet = l_system.get_alphabet();

    std::string init = l_system.get_initiator();

    double angleDegrees = l_system.get_angle();

    unsigned int iterations = l_system.get_nr_iterations();

    double currentAngleDegrees = l_system.get_starting_angle();

    double curX = 0;

    double curY = 0;

    Line2D line;

    std::string newRep;

    // generate full replacement
    for (int i = 0; i < iterations; ++i) {
        newRep = "";
        for (char c: init) {
            if (c == plus or c == minus or c == startBracket or c == endBracket) {
                newRep += c;
            } else {
                if (alphabet.find(c) != alphabet.end()) {
                    newRep += l_system.get_replacement(c);
                }
            }
        }
        init = newRep;
    }

    if (iterations == 0){
        newRep = init;
    }

    Lines2D tempLines;
    double tempX;
    double tempY;
    double tempAngle;

    // generate lines list
    for (int i = 0; i < newRep.size(); ++i) {
        if (newRep[i] == plus) {
            currentAngleDegrees += angleDegrees;
        } else if (newRep[i] == minus) {
            currentAngleDegrees -= angleDegrees;
        } else if (newRep[i] == startBracket) {
            std::pair<std::string, int> subStringPair = extractString(newRep, i, endBracket);
            std::string subString = subStringPair.first;

            tempX = curX;
            tempY = curY;
            tempAngle = currentAngleDegrees;

            tempLines = recursionLinesFinder(subString, l_system, bgClr, tempAngle, tempX, tempY);
            lines.merge(tempLines);

            i = subStringPair.second;
        } else {
            if (alphabet.find(newRep[i]) != alphabet.end()) {
                line = calculateLine(curX, curY, currentAngleDegrees, bgClr);
                if (l_system.draw(newRep[i])) {
                    lines.push_back(line);
                }
            }
        }
    }

    return lines;
}

Point2D doProjectionP(const Vector3D &point, const double d){
    Point2D p;
    double x = (point.x * d) / (-1*point.z);
    double y = (point.y * d) / (-1*point.z);
    p.setX(x);
    p.setY(y);
    p.setZ(point.z);
    return p;
}

Lines2D doProjectionL(const Figures3D &figures3D, double scale = 1){
    Lines2D lines;
    Point2D tempP;
    Line2D line2D;

    Vector3D newPoint;

    for (const auto& it : figures3D) {
        for (const auto& face : it.faces) {
            Points2D connectedPoints;
            for (const auto& points : face.point_indexes) {
                newPoint = it.points[points];
                tempP = doProjectionP(newPoint, scale);
                connectedPoints.push_back(tempP);
            }
            for (int i = 1; i < connectedPoints.size(); ++i) {
                line2D.p1 = connectedPoints[i-1];
                line2D.z1 = connectedPoints[i-1].getZ();
                line2D.p2 = connectedPoints[i];
                line2D.z2 = connectedPoints[i].getZ();
                line2D.rgb = it.ambientReflection;
                lines.push_back(line2D);
            }
            line2D.p1 = connectedPoints[connectedPoints.size() - 1];
            line2D.z1 = connectedPoints[connectedPoints.size() - 1].getZ();
            line2D.p2 = connectedPoints[0];
            line2D.z2 = connectedPoints[0].getZ();
            lines.push_back(line2D);
        }
    }

    return lines;
}

void applyMatrix(Matrix &m, Figure figureO, const Matrix& eye){
    m *= rotXMat(degreesToRadians(figureO.rotationX));
    m *= rotYMat(degreesToRadians(figureO.rotationY));
    m *= rotZMat(degreesToRadians(figureO.rotationZ));
    m *= translate(Vector3D::point(figureO.center[0], figureO.center[1], figureO.center[2]));
    m *= eye;
}

void recursionLinesFinder3D(const std::string& source, const LParser::LSystem3D &tl_system, Vector3D currentPos, Vector3D H, Vector3D L, Vector3D U, std::vector<Vector3D> &points, std::vector<Face> &faces, int currentPoint){
    Face face;
    face.point_indexes = {0,0};

    double angleDegrees = tl_system.get_angle();

    double angleRadians = degreesToRadians(angleDegrees);

    const std::set<char> &alphabet = tl_system.get_alphabet();

    Vector3D Hnew = H;
    Vector3D Lnew = L;
    Vector3D Unew = U;

    int currentPointIndex = currentPoint;

    for (int i = 0; i < source.size(); ++i) {
        if (source[i] == plus) {
            Hnew = H* cos(angleRadians) + L* sin(angleRadians);
            Lnew = (-H)* sin(angleRadians) + L* cos(angleRadians);

            H = Hnew;
            L = Lnew;
        } else if (source[i] == minus) {
            Hnew = H* cos(-angleRadians) + L* sin(-angleRadians);
            Lnew = (-H)* sin(-angleRadians) + L* cos(-angleRadians);

            H = Hnew;
            L = Lnew;
        } else if (source[i] == backslash){
            Lnew = L* cos(angleRadians) - U* sin(angleRadians);
            Unew = L* sin(angleRadians) + U* cos(angleRadians);

            L = Lnew;
            U = Unew;
        } else if (source[i] == slash) {
            Lnew = L* cos(-angleRadians) - U* sin(-angleRadians);
            Unew = L* sin(-angleRadians) + U* cos(-angleRadians);

            L = Lnew;
            U = Unew;
        } else if(source[i] == ampersand) {
            Hnew = H* cos(-angleRadians) + U* sin(-angleRadians);
            Unew = (-H)* sin(-angleRadians) + U* cos(-angleRadians);

            H = Hnew;
            U = Unew;
        } else if(source[i] == power) {
            Hnew = H* cos(angleRadians) + U* sin(angleRadians);
            Unew = (-H)* sin(angleRadians) + U* cos(angleRadians);

            H = Hnew;
            U = Unew;
        }
        else if (source[i] == startBracket) {
            std::pair<std::string, int> subStringPair = extractString(source, i, endBracket);
            std::string subString = subStringPair.first;

            recursionLinesFinder3D(subString, tl_system, currentPos, H, L, U, points, faces, currentPointIndex);

            i = subStringPair.second;
        }
        else {
            if (alphabet.find(source[i]) != alphabet.end()) {
                currentPos += H;
                points.push_back(currentPos);
                if (tl_system.draw(source[i])) {
                    face.point_indexes[0] = currentPointIndex;
                    face.point_indexes[1] = points.size() - 1;
                    currentPointIndex = points.size() - 1;
                    faces.push_back(face);
                }
            }
        }
    }
}

std::pair<std::vector<Vector3D>, std::vector<Face>> draw3DLsystem(const LParser::LSystem3D &tl_system , const std::vector<double>& bgClr){
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Face face;
    face.point_indexes = {0,0};
    // tl_system config
    const std::set<char> &alphabet = tl_system.get_alphabet();

    std::string init = tl_system.get_initiator();

    double angleDegrees = tl_system.get_angle();

    double angleRadians = degreesToRadians(angleDegrees);

    unsigned int iterations = tl_system.get_nr_iterations();

    std::string newRep;

    Vector3D H = Vector3D::point(1,0,0);
    Vector3D L = Vector3D::point(0,1,0);
    Vector3D U = Vector3D::point(0,0,1);

    Vector3D Hnew = Vector3D::point(1,0,0);
    Vector3D Lnew = Vector3D::point(0,1,0);
    Vector3D Unew = Vector3D::point(0,0,1);

    Vector3D startLoc = Vector3D::point(0,0,0);
    points.push_back(startLoc);

    // generate full replacement
    for (int i = 0; i < iterations; ++i) {
        newRep = "";
        for (char c: init) {
            if (c == plus or c == minus or c == startBracket or c == endBracket or c == backslash or c == slash or c == collon or c == power or c == ampersand) {
                newRep += c;
            } else {
                if (alphabet.find(c) != alphabet.end()) {
                    newRep += tl_system.get_replacement(c);
                }
            }
        }
        init = newRep;
    }

    if (iterations == 0){
        newRep = init;
    }

    std::pair<std::vector<Vector3D>, std::vector<Face>> tempPF;
    std::pair<std::vector<Vector3D>, std::vector<Face>> returnPF;
    int currentPointIndex = 0;

    // generate lines list
    for (int i = 0; i < newRep.size(); ++i) {
        if (newRep[i] == plus) {
            Hnew = H* cos(angleRadians) + L* sin(angleRadians);
            Lnew = (-H)* sin(angleRadians) + L* cos(angleRadians);

            H = Hnew;
            L = Lnew;
        } else if (newRep[i] == minus) {
            Hnew = H* cos(-angleRadians) + L* sin(-angleRadians);
            Lnew = (-H)* sin(-angleRadians) + L* cos(-angleRadians);

            H = Hnew;
            L = Lnew;
        } else if (newRep[i] == backslash){
            Lnew = L* cos(angleRadians) - U* sin(angleRadians);
            Unew = L* sin(angleRadians) + U* cos(angleRadians);

            L = Lnew;
            U = Unew;
        } else if (newRep[i] == slash) {
            Lnew = L* cos(-angleRadians) - U* sin(-angleRadians);
            Unew = L* sin(-angleRadians) + U* cos(-angleRadians);

            L = Lnew;
            U = Unew;
        } else if(newRep[i] == ampersand) {
            Hnew = H* cos(-angleRadians) + U* sin(-angleRadians);
            Unew = (-H)* sin(-angleRadians) + U* cos(-angleRadians);

            H = Hnew;
            U = Unew;
        } else if(newRep[i] == power) {
            Hnew = H* cos(angleRadians) + U* sin(angleRadians);
            Unew = (-H)* sin(angleRadians) + U* cos(angleRadians);

            H = Hnew;
            U = Unew;
        }
        else if (newRep[i] == collon){
            H = -H;
            L = -L;
        }
        else if (newRep[i] == startBracket) {
            std::pair<std::string, int> subStringPair = extractString(newRep, i, endBracket);
            std::string subString = subStringPair.first;

            recursionLinesFinder3D(subString, tl_system, startLoc, H, L, U, points, faces, currentPointIndex);

            i = subStringPair.second;
        }
        else {
            if (alphabet.find(newRep[i]) != alphabet.end()) {
                startLoc += H;
                points.push_back(startLoc);
                if (tl_system.draw(newRep[i])) {
                    face.point_indexes[0] = currentPointIndex;
                    face.point_indexes[1] = points.size() - 1;
                    currentPointIndex = points.size() - 1;
                    faces.push_back(face);
                }
            }
        }
    }
    returnPF.first = points;
    returnPF.second = faces;
    return returnPF;
}

double inverseZ(const double i, const double a, const double z1, const double z2, const double min, const double max){
    double invZ = 0;
    invZ = (i/a)/z1 + (1 - (i/a))/z2;
    return invZ;
}

void draw_zbuf_line(ZBuffer &zbuffer, img::EasyImage &image, unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, img::Color color){
    assert(x0 <= image.get_width() && y0 <= image.get_height());
    assert(x1 <= image.get_width() && y1 <= image.get_height());

    double zInv;

    if (x0 == x1)
    {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            if (y1 > y0){
                std::swap(z1, z0);
            }
            zInv = inverseZ(i, std::max(y0, y1) - std::min(y0, y1), z0, z1, std::min(y0, y1), std::max(y0, y1));
            if (zInv < zbuffer[x0][i]){
                zbuffer[x0][i] = zInv;
                image(x0, i) = color;
            }
        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            if (x1 > x0){
                std::swap(z0, z1);
            }
            zInv = inverseZ(i, std::max(x0, x1) - std::min(x0, x1), z0, z1, std::min(x0, x1), std::max(x0, x1));
            if (zInv < zbuffer[i][y0]){
                zbuffer[i][y0] = zInv;
                image(i, y0) = color;
            }
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                zInv = inverseZ(i, x1 - x0, z1, z0, std::min(x0, x1), std::max(x0, x1));
                if (zInv < zbuffer[x0+i][(unsigned int) round(y0 + m * i)]){
                    zbuffer[x0+i][(unsigned int) round(y0 + m * i)] = zInv;
                    image(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                }

            }
        }
        else if (m > 1.0)
        {
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                zInv = inverseZ(i, y1 - y0, z1, z0, 0, (y1 - y0));
                if (zInv < zbuffer[(unsigned int) round(x0 + (i / m))][y0 + i]){
                    zbuffer[(unsigned int) round(x0 + (i / m))][y0 + i] = zInv;
                    image((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                }

            }
        }
        else if (m < -1.0)
        {
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                zInv = inverseZ(i, y0 - y1, z1, z0, 0, (y0 - y1));
                if(zInv < zbuffer[(unsigned int) round(x0 - (i / m))][y0-i]){
                    zbuffer[(unsigned int) round(x0 - (i / m))][y0-i] = zInv;
                    image((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                }
            }
        }
    }
}

img::EasyImage drawZbufLines(Lines2D &lines, const int size, const std::vector<double>& bgColorV) {
    // Determine min and max of x and y
    std::pair<double,double> max = determine_max(lines);
    std::pair<double,double> min = determine_min(lines);

    double Xmax = max.first;
    double Xmin = min.first;

    double Ymax = max.second;
    double Ymin = min.second;

    double xrange = Xmax - Xmin;
    double yrange = Ymax - Ymin;

    // calculate image size
    double imageX = size * (xrange / std::max(xrange, yrange)); // width
    double imageY = size * (yrange / std::max(xrange, yrange)); // height

    // determine scale factor
    double d = 0.95 * (imageX/xrange);

    // determine Shift amount
    double DCx = d*((Xmin + Xmax)/2);
    double DCy = d*((Ymin+Ymax)/2);

    double dx = (imageX/2)-DCx;
    double dy = (imageY/2)-DCy;

    // Color
    img::Color color;
    img::Color bgColor;

    bgColor.red = bgColorV[0] * 255;
    bgColor.green = bgColorV[1] * 255;
    bgColor.blue = bgColorV[2] * 255;

    // generate image
    img::EasyImage image(imageX, imageY, bgColor);

    ZBuffer zbuf(lround(image.get_width()), lround(image.get_height()));

    for(const auto& it: lines){
        color.red = it.rgb.getRed() * 255;
        color.green = it.rgb.getGreen() * 255;
        color.blue = it.rgb.getBlue() * 255;

        draw_zbuf_line(zbuf, image, lround(it.p1.getX()*d + dx), lround(it.p1.getY()*d + dy), it.p1.getZ(), lround(it.p2.getX()*d + dx), lround(it.p2.getY()*d + dy), it.p2.getZ(), color);
    }
    return image;
}

std::vector<Face> triangulate(const Face &face){
    Face newFace;
    newFace.point_indexes = {face.point_indexes[0],0,0};
    std::vector<Face> faces;

    for (int i = 1; i < face.point_indexes.size()-1; ++i) {
        for (int j = 0; j < 2; ++j) {
            newFace.point_indexes[j + 1] = face.point_indexes[j + i];
        }
        faces.push_back(newFace);
    }
    return faces;
}

void draw_zbuf_triag(ZBuffer &zbuffer, img::EasyImage &image, Vector3D const &A, Vector3D const &B, Vector3D const &C, double d, double dx, double dy, Color ambientReflection, Color diffuseReflection, Color specularReflection, double reflectionCoeff, Lights3D &lights){
    // determine color
    img::Color color;
    Color color1;
    Color color2;
    if (lights.empty()){
        color.red = ambientReflection.getRed() * 255;
        color.green = ambientReflection.getGreen() * 255;
        color.blue = ambientReflection.getBlue() * 255;
    }
    else {
        for (auto light : lights) {
            color2 += (light->ambientLight * ambientReflection);
        }
    }

    // new points
    Point2D newA;
    Point2D newB;
    Point2D newC;

    newA.setX((((A.x * d)/-A.z) + dx));
    newA.setY((((A.y * d)/-A.z) + dy));

    newB.setX((((B.x * d)/-B.z) + dx));
    newB.setY((((B.y * d)/-B.z) + dy));

    newC.setX((((C.x * d)/-C.z) + dx));
    newC.setY((((C.y * d)/-C.z) + dy));

    // y min and y max
    double ymin = newA.getY();
    double ymax = newA.getY();

    // min
    if (ymin > newB.getY()){
        ymin = newB.getY();
    }

    if (ymin > newC.getY()){
        ymin = newC.getY();
    }

    // max
    if (ymax < newB.getY()){
        ymax = newB.getY();
    }

    if (ymax < newC.getY()){
        ymax = newC.getY();
    }

    // zwaartepunt waardes
    double xG = (newA.getX() + newB.getX() + newC.getX())/3;
    double yG = (newA.getY() + newB.getY() + newC.getY())/3;

    double inverZG = 1/(3*A.z) + 1/(3*B.z) + 1/(3*C.z);

    // determine line
    ymin = lround(ymin);
    ymax = lround(ymax);

    // uvw
    Vector3D u;
    Vector3D v;
    Vector3D w;

    u = B - A;
    v = C - A;

    w.x = u.y * v.z - u.z * v.y;
    w.y = u.z * v.x - u.x * v.z;
    w.z = u.x * v.y - u.y * v.x;

    for (double y = ymin + 0.5; y <= ymax - 0.5; ++y){
        double xLab = postInf, xLac = postInf, xLbc = postInf;
        double xRab = negInf, xRac = negInf, xRbc = negInf;

        for (int j = 0; j < 3; ++j) {
            switch (j) {
                case 0: // AB
                    if ((y - newA.getY()) * (y - newB.getY()) <= 0 and newA.getY() != newB.getY()){
                        xLab = xRab = newB.getX() + (newA.getX() - newB.getX())*((y - newB.getY()) / (newA.getY() - newB.getY()));
                    }
                    break;
                case 1: // AC
                    if ((y - newA.getY()) * (y - newC.getY()) <= 0 and newA.getY() != newC.getY()){
                        xLac = xRac = newC.getX() + (newA.getX() - newC.getX())*((y - newC.getY()) / (newA.getY() - newC.getY()));
                    }
                    break;
                case 2: // BC
                    if ((y - newC.getY()) * (y - newB.getY()) <= 0 and newC.getY() != newB.getY()){
                        xLbc = xRbc = newC.getX() + (newB.getX() - newC.getX())*((y - newC.getY()) / (newB.getY() - newC.getY()));
                    }
                    break;
            }
        }
        double xL;
        double xR;

        if (xLab == postInf and xRab == negInf){
            xL = lround(std::min(xLbc, xLac) + 0.5);
            xR = lround(std::max(xRbc, xRac) - 0.5);
        }
        else if (xLac == postInf and xRac == negInf) {
            xL = lround(std::min(xLbc, xLab) + 0.5);
            xR = lround(std::max(xRbc, xRab) - 0.5);
        }
        else if (xLbc == postInf and xRbc == negInf) {
            xL = lround(std::min(xLab, xLac) + 0.5);
            xR = lround(std::max(xRab, xRac) - 0.5);
        }
        else {
            xL = lround(std::min(std::min(xLbc, xLab), xLac) + 0.5);
            xR = lround(std::max(std::max(xRbc, xRab), xRac) - 0.5);
        }

        // calculate inverse z
         // dzdx en dzdy
        Vector3D nW = w;
        nW.normalise();

        double k;
        k = w.x * A.x + w.y * A.y + w.z * A.z;

        double dzdx;
        double dzdy;

        dzdx = w.x/(-d*k);
        dzdy = w.y/(-d*k);

        for (int x = xL; x <= xR; ++x) {
            double inverZ;
            inverZ = 1.0001*inverZG + (x - xG)*dzdx + (y - yG) * dzdy;

            if (zbuffer[x][y] > inverZ){
                color1 = color2;
                // point light
                double newX = (x-dx) / (-inverZ*d);
                double newY = (y-dy) / (-inverZ*d);

                Vector3D newPoint = Vector3D::vector(newX, newY, 1/inverZ);
                for (auto &light : lights) {
                    // scalar product of normalise and lightDir
                    if (light->getLightType() == normalL){
                        continue;
                    }
                    else if (light->getLightType() == pointL){
                        auto pointLight = (PointLight *) light;

                        Vector3D ld = Vector3D::vector(pointLight->location - newPoint);
                        ld.normalise();
                        double alpha = Vector3D::dot(nW, ld);

                        double cosSpotAngle = cos(degreesToRadians(pointLight->spotAngle));

                        Vector3D r = 2 * nW * alpha - ld;
                        r.normalise();
                        double beta = Vector3D::dot(r, -newPoint);

                        if (alpha >= cosSpotAngle){
                            Color temp = diffuseReflection * pointLight->diffuseLight;
                            temp.setRGB(temp.getRed() * (1-((1-alpha)/(1-cosSpotAngle))), temp.getGreen() * (1-((1-alpha)/(1-cosSpotAngle))), temp.getBlue() * (1-((1-alpha)/(1-cosSpotAngle))));
                            color1 += temp;
                        }
                        if (beta > 0){
                            Color temp = specularReflection * pointLight->specularLight;
                            beta = pow(beta, reflectionCoeff);
                            temp.setRGB(temp.getRed() * beta, temp.getGreen() * beta, temp.getBlue() * beta);
                            color1 += temp;
                        }
                    }
                    else {
                        // inf light
                        auto infLight = (InfLight *) light;
                        Vector3D ld = infLight->direction;
                        ld.normalise();
                        ld = -ld;
                        double alpha = Vector3D::dot(nW, ld);

                        Vector3D r = 2 * nW * alpha - ld;
                        r.normalise();

                        double beta = Vector3D::dot(r, -newPoint);

                        if (alpha > 0){
                            Color temp = diffuseReflection * infLight->diffuseLight;
                            temp.setRGB(temp.getRed() * alpha, temp.getGreen() * alpha, temp.getBlue() * alpha);
                            color1 += temp;
                        }
                        if (beta > 0){
                            Color temp = specularReflection * infLight->specularLight;
                            beta = pow(beta, reflectionCoeff);
                            temp.setRGB(temp.getRed() * beta, temp.getGreen() * beta, temp.getBlue() * beta);
                            color1 += temp;
                        }
                    }
                }

                if (color1.getRed() > 1){
                    color1.setRed(1);
                }
                if (color1.getGreen() > 1){
                    color1.setGreen(1);
                }
                if (color1.getBlue() > 1){
                    color1.setBlue(1);
                }

                // convert Color to image color
                color.red = color1.getRed() * 255;
                color.blue = color1.getBlue() * 255;
                color.green = color1.getGreen() * 255;

                zbuffer[x][y] = inverZ;
                image(x, y) = color;
            }
        }
    }
}

void convertFaceToTriFace(Figures3D &figures){
    for (auto &figure : figures){
        std::vector<Face> totFaces;
        for (const auto &face : figure.faces){
            std::vector<Face> newFaces = triangulate(face);
            totFaces.insert(totFaces.end(), newFaces.begin(), newFaces.end());
        }
        figure.faces = totFaces;
    }
}

img::EasyImage drawZbuffTriag(Figures3D &figures, const int size, const std::vector<double>& bgColorV, Lines2D &lSystemLines, Lights3D lights3D = {}){
    Lines2D lines;
    Lines2D newLines;

    // doprojection with d = 1
    lines = doProjectionL(figures);

    // find d
        // Determine min and max of x and y
        std::pair<double,double> max = determine_max(lines);
        std::pair<double,double> min = determine_min(lines);

        double Xmax = max.first;
        double Xmin = min.first;

        double Ymax = max.second;
        double Ymin = min.second;

        double xrange = Xmax - Xmin;
        double yrange = Ymax - Ymin;

        // calculate image size
        double imageX = size * (xrange / std::max(xrange, yrange)); // width
        double imageY = size * (yrange / std::max(xrange, yrange)); // height

        // determine scale factor
        double d = 0.95 * (imageX/xrange);

        // determine Shift amount
        double DCx = d*((Xmin + Xmax)/2);
        double DCy = d*((Ymin+Ymax)/2);

        double dx = (imageX/2)-DCx;
        double dy = (imageY/2)-DCy;

    // doprojection with d = newD
    convertFaceToTriFace(figures);
    newLines = doProjectionL(figures, d);
    newLines.merge(lSystemLines);

    // Color
    img::Color bgColor;

    bgColor.red = bgColorV[0] * 255;
    bgColor.green = bgColorV[1] * 255;
    bgColor.blue = bgColorV[2] * 255;

    // generate img
    img::EasyImage image(imageX, imageY, bgColor);

    ZBuffer zbuf(lround(image.get_width()), lround(image.get_height()));

    Vector3D A;
    Vector3D B;
    Vector3D C;

    for (const auto &figure : figures){
        for (const auto &face : figure.faces){
            A = figure.points[face.point_indexes[0]];
            B = figure.points[face.point_indexes[1]];
            C = figure.points[face.point_indexes[2]];

            draw_zbuf_triag(zbuf, image, A, B, C, d, dx, dy, figure.ambientReflection, figure.diffuseReflection, figure.specularReflection, figure.reflectionCoefficient, lights3D);
        }
    }
    return image;
}

//==== Fractal ====//
void generateFractal(Figure &figure, Figures3D& fractal, const int nr_iterations, const double scale){
    Matrix newScale = scaleFigure(1/scale);

    Figures3D newFractal;

    fractal.push_back(figure);

    for (int i = 0; i < nr_iterations; ++i) {
        for (auto &figure : fractal) {
            for (int j = 0; j < figure.points.size(); ++j) {
                Figure tempfig = figure;

                applyTransformationFig(tempfig, newScale);

                Matrix translation = translate(figure.points[j] - tempfig.points[j]);

                applyTransformationFig(tempfig, translation);

                newFractal.push_back(tempfig);
            }
        }
        fractal = newFractal;
        newFractal = {};
    }
}

Figures3D parseFigures(const ini::Configuration &configuration, const Matrix &eyePointMat, Lines2D &lines, Lines2D& lines3D){
    int figureAm = configuration["General"]["nrFigures"];
    Figures3D figures3D;

    std::vector<Vector3D> pointsV;
    std::vector<double> pointsD;
    std::vector<Face> linesT;

    std::string iString;
    std::string tempString2;

    Figure figureO;

    std::string type = configuration["General"]["type"].as_string_or_die();

    for (int i = 0; i < figureAm; ++i) {
        std::string tempString1;
        tempString1 = figureS + std::to_string(i);
        std::string type2 = configuration[tempString1]["type"];
        figureO.name = type2;
        if (type2 == "LineDrawing") {
            int pointAm = configuration[tempString1]["nrPoints"];
            int lineAm = configuration[tempString1]["nrLines"];

            for (int j = 0; j < pointAm; j++) {
                iString = std::to_string(j);
                tempString2 = pointS + iString;
                pointsD = configuration[tempString1][tempString2].as_double_tuple_or_die();
                pointsV.push_back(vecTo3D(pointsD));
            }

            for (int j = 0; j < lineAm; j++) {
                iString = std::to_string(j);
                tempString2 = lineS + iString;
                Face lineT;
                lineT.point_indexes = configuration[tempString1][tempString2].as_int_tuple_or_die();
                linesT.push_back(lineT);
            }

            figureO.faces = linesT;
            figureO.points = pointsV;

            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);

            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if(type2 == "Cube"){
            figureO = createCube();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if (type2 == "3DLSystem"){
            LParser::LSystem3D tl_system;

            std::ifstream input_stream(configuration[tempString1]["inputfile"]);
            input_stream >> tl_system;
            input_stream.close();

            std::pair<std::vector<Vector3D>, std::vector<Face>> PF;
            PF = draw3DLsystem(tl_system, configuration[tempString1]["color"]);
            figureO.points = PF.first;
            figureO.faces = PF.second;
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if(type2 == "Tetrahedron"){
            figureO = createTetraHedron();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if (type2 == "Icosahedron"){
            figureO = createIcosahedron();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if (type2 == "Octahedron"){
            figureO = createOctahedron();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if (type2 == "Dodecahedron"){
            figureO = createDodecahedron();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if (type2 == "Cone"){
            figureO = createCone(configuration[tempString1]["n"], configuration[tempString1]["height"]);
            figureO.name = type2;
            Color color;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if(type2 == "Cylinder"){
            figureO = createCylinder(configuration[tempString1]["n"], configuration[tempString1]["height"]);
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if (type2 == "Torus"){
            figureO = createTorus(configuration[tempString1]["r"], configuration[tempString1]["R"], configuration[tempString1]["m"], configuration[tempString1]["n"]);
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if (type2 == "Sphere"){
            figureO = createSphere(configuration[tempString1]["n"]);
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
        else if (type2 == "FractalTetrahedron"){
            figureO = createTetraHedron();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);

            Figures3D newListOfFigures;
            generateFractal(figureO, newListOfFigures, configuration[tempString1]["nrIterations"].as_int_or_die(), configuration[tempString1]["fractalScale"].as_double_or_die());

            figures3D.insert(figures3D.end(), newListOfFigures.begin(), newListOfFigures.end());
        }
        else if (type2 == "FractalCube"){
            figureO = createCube();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);

            Figures3D newListOfFigures;
            generateFractal(figureO, newListOfFigures, configuration[tempString1]["nrIterations"].as_int_or_die(), configuration[tempString1]["fractalScale"].as_double_or_die());

            figures3D.insert(figures3D.end(), newListOfFigures.begin(), newListOfFigures.end());
        }
        else if (type2 == "FractalIcosahedron"){
            figureO = createIcosahedron();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);

            Figures3D newListOfFigures;
            generateFractal(figureO, newListOfFigures, configuration[tempString1]["nrIterations"].as_int_or_die(), configuration[tempString1]["fractalScale"].as_double_or_die());

            figures3D.insert(figures3D.end(), newListOfFigures.begin(), newListOfFigures.end());
        }
        else if (type2 == "FractalDodecahedron"){
            figureO = createDodecahedron();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);

            Figures3D newListOfFigures;
            generateFractal(figureO, newListOfFigures, configuration[tempString1]["nrIterations"].as_int_or_die(), configuration[tempString1]["fractalScale"].as_double_or_die());

            figures3D.insert(figures3D.end(), newListOfFigures.begin(), newListOfFigures.end());
        }
        else if (type2 == "FractalOctahedron"){
            figureO = createOctahedron();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);

            Figures3D newListOfFigures;
            generateFractal(figureO, newListOfFigures, configuration[tempString1]["nrIterations"].as_int_or_die(), configuration[tempString1]["fractalScale"].as_double_or_die());

            figures3D.insert(figures3D.end(), newListOfFigures.begin(), newListOfFigures.end());
        }
        else if (type2 == "FractalBuckyBall"){
            figureO = createBuckyBall();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);

            Figures3D newListOfFigures;
            generateFractal(figureO, newListOfFigures, configuration[tempString1]["nrIterations"].as_int_or_die(), configuration[tempString1]["fractalScale"].as_double_or_die());

            figures3D.insert(figures3D.end(), newListOfFigures.begin(), newListOfFigures.end());
        }
        else if(type2 == "BuckyBall") {
            figureO = createBuckyBall();
            figureO.name = type2;
            if (type == "LightedZBuffering"){
                std::vector<double> ambientV = configuration[tempString1]["ambientReflection"].as_double_tuple_or_die();
                figureO.ambientReflection.setRGB(ambientV[0], ambientV[1], ambientV[2]);
            }
            else {
                Color color;
                std::vector<double> colVec = configuration[tempString1]["color"].as_double_tuple_or_die();
                color.setRGB(colVec[0], colVec[1],
                             colVec[2]);

                figureO.ambientReflection = color;
            }
            std::vector<double> diffuseV = configuration[tempString1]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
            std::vector<double> specularV = configuration[tempString1]["specularReflection"].as_double_tuple_or_default({0,0,0});
            figureO.diffuseReflection.setRGB(diffuseV[0],diffuseV[1],diffuseV[2]);
            figureO.specularReflection.setRGB(specularV[0],specularV[1],specularV[2]);
            figureO.reflectionCoefficient = configuration[tempString1]["reflectionCoefficient"].as_double_or_default(0);


            figureO.rotationX = configuration[tempString1]["rotateX"].as_double_or_die();
            figureO.rotationY = configuration[tempString1]["rotateY"].as_double_or_die();
            figureO.rotationZ = configuration[tempString1]["rotateZ"].as_double_or_die();

            figureO.scale = configuration[tempString1]["scale"].as_double_or_die();

            figureO.center = configuration[tempString1]["center"].as_double_tuple_or_die();

            Matrix m = scaleFigure(figureO.scale);

            applyMatrix(m, figureO, eyePointMat);
            applyTransformationFig(figureO, m);
            figures3D.push_back(figureO);
        }
    }

    return figures3D;
}

void setLightSet(Light* light, const ini::Configuration &configuration, int i){
    std::vector<double> ambientV = configuration["Light" + std::to_string(
            i)]["ambientLight"].as_double_tuple_or_default({1, 1, 1});
    std::vector<double> diffuseV = configuration["Light" + std::to_string(
            i)]["diffuseLight"].as_double_tuple_or_default({0, 0, 0});
    std::vector<double> specularV = configuration["Light" + std::to_string(
            i)]["specularLight"].as_double_tuple_or_default({0, 0, 0});

    Color ambi(ambientV[0], ambientV[1], ambientV[2]);
    Color diff(diffuseV[0], diffuseV[1], diffuseV[2]);
    Color spec(specularV[0], specularV[1], specularV[2]);

    light->ambientLight = ambi;
    light->diffuseLight = diff;
    light->specularLight = spec;
}

//==== Generate Image Control Functions ====//
img::EasyImage generate_image(const ini::Configuration &configuration)
{
    // configuration data
        // type
    std::string type = configuration["General"]["type"].as_string_or_die();

    img::EasyImage image(0,0);

    // View frustrum
    bool viewEnable = configuration["General"]["clipping"].as_bool_or_default(false);
    int hfov;
    double aspectRatio;
    double dNear;
    double dFar;
    std::vector<double> viewDir;

    // eye
    std::vector<double> eyePointVec = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eyePoint3D = Vector3D::vector(eyePointVec[0], eyePointVec[1], eyePointVec[2]);

    double theta;
    double phi;
    double r;

    toPolar(eyePoint3D, theta, phi, r);

    if (viewEnable){
        hfov = configuration["General"]["hfov"].as_int_or_die();
        aspectRatio = configuration["General"]["aspectRatio"].as_double_or_die();
        dNear = configuration["General"]["dNear"].as_double_or_die();
        dFar = configuration["General"]["dFar"].as_double_or_die();
        viewDir = configuration["General"]["viewDirection"].as_double_tuple_or_die();

        double theta2;
        double phi2;
        double r2;

        Vector3D viewDir3D = Vector3D::vector(viewDir[0], viewDir[1], viewDir[2]);
        toPolar(-viewDir3D, theta2, phi2, r2);
        theta = theta2;
        phi = phi2;
    }

    Matrix eyePointMat = eyeTransformatieMatrix(theta, phi, r);

    if(type == "2DLSystem"){
        LParser::LSystem2D l_system;

        std::ifstream input_stream(configuration["2DLSystem"]["inputfile"]);
        input_stream >> l_system;
        input_stream.close();

        Lines2D lines = drawLsystem(l_system, configuration["2DLSystem"]["color"]);
        image = draw2DLines(lines, configuration["General"]["size"].as_int_or_die(), configuration["General"]["backgroundcolor"].as_double_tuple_or_die());
    }
    else if (type == "Wireframe"){
        int figureAm = configuration["General"]["nrFigures"];
        int size = configuration["General"]["size"].as_int_or_die();
        std::vector<double> bgc = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

        Lines2D lines;
        Lines2D lines3D;
        Figures3D figures3D = parseFigures(configuration, eyePointMat, lines, lines3D);

        Lines2D Flines = doProjectionL(figures3D);
        lines.merge(Flines);
        lines.merge(lines3D);
        image = draw2DLines(lines, size, bgc);
    }
    else if (type == "ZBufferedWireframe"){
        int size = configuration["General"]["size"].as_int_or_die();
        std::vector<double> bgc = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

        Lines2D lines;
        Lines2D lines3D;

        Figures3D figures3D = parseFigures(configuration, eyePointMat, lines, lines3D);
        Lines2D Flines = doProjectionL(figures3D);
        lines.merge(Flines);
        lines.merge(lines3D);
        image = drawZbufLines(lines, size, bgc);
    }
    else if (type == "ZBuffering"){
        int figureAm = configuration["General"]["nrFigures"];
        int size = configuration["General"]["size"].as_int_or_die();
        std::vector<double> bgc = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

        Figures3D figures3D;
        Lines2D lines;
        Lines2D lines3D;

        figures3D = parseFigures(configuration, eyePointMat, lines, lines3D);
        image = drawZbuffTriag(figures3D, size, bgc, lines3D);
    }
    else if (type == "LightedZBuffering"){
        int figureAm = configuration["General"]["nrFigures"];
        int size = configuration["General"]["size"].as_int_or_die();
        std::vector<double> bgc = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int amountLights = configuration["General"]["nrLights"].as_int_or_die();

        Lights3D lights3D;
        for (int i = 0; i < amountLights; ++i) {
            bool inf_val;
            bool inf = configuration["Light" + std::to_string(i)]["infinity"].as_bool_if_exists(inf_val);
            if (inf){
                if (inf_val){
                    auto* infLight = new InfLight();
                    std::vector<double> dirV = configuration["Light" + std::to_string(i)]["direction"].as_double_tuple_or_die();
                    Vector3D direction = Vector3D::vector(dirV[0],dirV[1],dirV[2]);
                    direction *= eyePointMat;
                    direction.normalise();
                    infLight->direction = direction;
                    setLightSet(infLight, configuration, i);
                    infLight->setLightType(infL);
                    lights3D.push_back(infLight);
                }
                else {
                    auto* pointLight = new PointLight();
                    std::vector<double> locV = configuration["Light" + std::to_string(i)]["location"].as_double_tuple_or_die();
                    Vector3D location = Vector3D::point(locV[0], locV[1], locV[2]);
                    location *= eyePointMat;
                    pointLight->location = location;
                    pointLight->spotAngle = configuration["Light" + std::to_string(i)]["spotAngle"].as_double_or_default(90);
                    setLightSet(pointLight, configuration, i);
                    pointLight->setLightType(pointL);
                    lights3D.push_back(pointLight);
                }
            }
            else {
                auto* light = new Light();
                setLightSet(light, configuration, i);
                lights3D.push_back(light);
            }
        }

        Figures3D figures3D;
        Lines2D lines;
        Lines2D lines3D;

        figures3D = parseFigures(configuration, eyePointMat, lines, lines3D);
        image = drawZbuffTriag(figures3D, size, bgc, lines3D, lights3D);
    }
    return image;
}

//==== Main ====//
int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
