//
// Created by liuja on 11/05/2022.
//

#include "PointLight.h"

PointLight::PointLight() {}

PointLight::PointLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight) : Light(
        ambientLight, diffuseLight, specularLight){
}
