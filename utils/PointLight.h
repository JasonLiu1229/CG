//
// Created by liuja on 11/05/2022.
//

#ifndef ZBUFFER_CPP_POINTLIGHT_H
#define ZBUFFER_CPP_POINTLIGHT_H
#include "Light.h"
#include "vector3d.h"

class PointLight : public Light{
public:
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;

    PointLight();

    PointLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight);
};


#endif //ZBUFFER_CPP_POINTLIGHT_H
