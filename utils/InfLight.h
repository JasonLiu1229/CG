//
// Created by liuja on 11/05/2022.
//

#ifndef ZBUFFER_CPP_INFLIGHT_H
#define ZBUFFER_CPP_INFLIGHT_H
#include "Light.h"
#include "vector3d.h"

class InfLight : public Light{
public:
    //de richting waarin het
    //licht schijnt
    Vector3D direction;

    InfLight();

    InfLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight);
};


#endif //ZBUFFER_CPP_INFLIGHT_H
