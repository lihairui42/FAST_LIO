#pragma once

#define INS_TOOLS_H
#ifdef INS_TOOLS_H

#include <iostream>
#include <Eigen/Core>

Eigen::Matrix3d Askew(Eigen::Vector3d a)
{
    Eigen::Matrix3d     b;

    b << 0   , -a.z(),  a.y(),
         a.z(),     0, -a.x(),
        -a.y(),  a.x(),     0;

    return b;
}

#endif



