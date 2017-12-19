//
// Created by neo on 11/28/17.
//
#pragma once

#include <cmath>
#include "ofMain.h"

class Vector : public ofVec3f {
  public:
    double x;
    double y;
    double z;

    Vector();

    Vector(const double &, const double &, const double &);

    ~Vector();

    double Magnitude();

    void Normalize();

    void Reverse();

    Vector &operator+=(const Vector &);

    Vector &operator-=(const Vector &);

    Vector &operator*=(const double &);

    Vector &operator/=(const double &);

    Vector operator-();
};

extern inline Vector operator+(const Vector &, const Vector &);

extern inline Vector operator-(const Vector &, const Vector &);

extern inline Vector operator^(const Vector &, const Vector &);

extern inline double operator*(const Vector &, const Vector &);

extern inline Vector operator*(const double &, const Vector &);

extern inline Vector operator*(const Vector &, const double &);

extern inline Vector operator/(const Vector &, const double &);