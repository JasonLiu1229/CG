//
// Created by liuja on 11/05/2022.
//

#include "Light.h"

Light::Light(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight) : ambientLight(
        ambientLight), diffuseLight(diffuseLight), specularLight(specularLight), newLightType(normalL) {}

Light::Light() {}

Light::~Light() {

}

lightType Light::getLightType() const {
    return newLightType;
}

void Light::setLightType(enum lightType new_LightType) {
    Light::newLightType = new_LightType;
}
