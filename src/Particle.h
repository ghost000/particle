//
// Created by neo on 11/28/17.
//
#pragma once

#include "ConstantValues.h"
#include "ofMain.h"

struct Particle : public ofNode {
    double Mass;
    double InvMass;
    Vector Position;
    Vector Velocity;
    Vector Acceleration;
    Vector Forces;
    bool Locked;
    ofColor color;
};

struct ParticleRef {
    int r;
    int c;
};

struct Spring {
    ParticleRef p1;
    ParticleRef p2;
    double k;
    double d;
    double L;
};

struct Collision {
    ParticleRef p1;
    Vector n;
};