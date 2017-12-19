//
// Created by neo on 11/28/17.
//
#pragma once

#include <cstdlib>
#include <vector>
#include "Vector.cpp"

constexpr double NUMROWS{20};
constexpr double NUMCOLUMNS{20};
constexpr double CLOTHHEIGHT{280};
constexpr double CLOTHWIDTH{400};
constexpr double CLOTHMASS{400};
constexpr double NUMSTRUCTURALSPRINGS{NUMCOLUMNS * (NUMROWS + 1) + NUMROWS * (NUMCOLUMNS + 1) + NUMCOLUMNS * NUMROWS * 2};
constexpr double NUMFACES{(NUMCOLUMNS * NUMROWS) * 2};
constexpr double NUMVERTICES{(NUMCOLUMNS + 1) * (NUMROWS + 1)};
constexpr double MASSPERFACE{CLOTHMASS / NUMFACES};
constexpr double CSTEP{CLOTHWIDTH / NUMCOLUMNS};
constexpr double RSTEP{CLOTHHEIGHT / NUMROWS};

constexpr double GRAVITY{-32.174};
constexpr double SPRINGTENSIONCONSTANT{500};
constexpr double SPRINGSHEARCONSTANT{300};
constexpr double SPRINGDAMPINGCONSTANT{2};
constexpr double YOFFSET{120};
constexpr double DRAGCOEFFICIENT{0.01};
constexpr double WINDFACTOR{100};
constexpr double FLAGPOLEHEIGHT{200};
constexpr double FLAGPOLERADIUS{10};
constexpr double COLLISIONTOLERANCE{0.05};
constexpr double KRESTITUTION{0.25};
constexpr double FRICTIONFACTOR{0.5};
constexpr double WindForceFactor{WINDFACTOR};

constexpr double NOCOLLISION{0};
constexpr double COLLISION{1};
constexpr double PENETRATING{-1};

constexpr double rho{0.0023769};
constexpr double tol{0.0001};
