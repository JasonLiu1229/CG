//
// Created by liuja on 11/05/2022.
//

#include "InfLight.h"

InfLight::InfLight() {}

InfLight::InfLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight) : Light(
        ambientLight, diffuseLight, specularLight) {
}
