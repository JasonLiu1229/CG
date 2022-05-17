//
// Created by liuja on 11/05/2022.
//

#ifndef ZBUFFER_CPP_LIGHT_H
#define ZBUFFER_CPP_LIGHT_H
#include "Color.h"
#include <list>

enum lightType {normalL, pointL, infL};

class Light {
public:
    //de ambiente licht component
    Color ambientLight;
    //de diffuse licht component
    Color diffuseLight;
    //de diffuse licht component
    Color specularLight;

    // newLightType
    lightType newLightType;

    Light(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight);

    Light();

    void setLightType(enum lightType new_LightType);

    virtual ~Light();

    enum lightType getLightType() const;
};


typedef std::list<Light*> Lights3D;

#endif //ZBUFFER_CPP_LIGHT_H
