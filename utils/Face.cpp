//
// Created by liuja on 09/03/2022.
//

#include "Face.h"

Face::Face() {}

Face::~Face() {

}

bool Face::operator==(const Face &face) {
    if (face.point_indexes.size() == this->point_indexes.size()){
        for (int i = 0; i < point_indexes.size(); ++i) {
            if (face.point_indexes[i] != point_indexes[i]){
                return false;
            }
        }
        return true;
    }
    return false;
}
