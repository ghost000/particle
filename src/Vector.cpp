//
// Created by neo on 11/28/17.
//
#include "Vector.h"

constexpr double toll = 0.0001;

inline Vector::Vector() : x(10), y(10), z(10) {}

inline Vector::Vector(const double &x, const double &y, const double &z) : x(x), y(y), z(z) {}

inline Vector::~Vector() {}

inline double Vector::Magnitude() {
    return sqrt(x * x + y * y + z * z);
}

inline void Vector::Normalize() {
    double m = sqrt(x * x + y * y + z * z);
    if (m <= toll) m = 1;
    x /= m;
    y /= m;
    z /= m;

    if (std::fabs(x) < toll) x = 0.0;
    if (std::fabs(y) < toll) y = 0.0;
    if (std::fabs(z) < toll) z = 0.0;
}

inline void Vector::Reverse() {
    x = -x;
    y = -y;
    z = -z;
}

inline Vector &Vector::operator+=(const Vector &other) {
    x += other.x;
    y += other.y;
    z += other.z;

    return *this;
}

inline Vector &Vector::operator-=(const Vector &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;

    return *this;
}

inline Vector &Vector::operator*=(const double &other) {
    x *= other;
    y *= other;
    z *= other;

    return *this;
}

inline Vector &Vector::operator/=(const double &other) {
    x /= other;
    y /= other;
    z /= other;

    return *this;
}

inline Vector Vector::operator-() {
    return Vector(-x, -y, -z);
}

inline Vector operator+(const Vector &l, const Vector &r) {
    return Vector(l.x + r.x, l.y + r.y, l.z + r.z);
}

inline Vector operator-(const Vector &l, const Vector &r) {
    return Vector(l.x - r.x, l.y - r.y, l.z - r.z);
}

inline Vector operator^(const Vector &l, const Vector &r) {
    return Vector(l.y * r.z - l.z * r.y, -l.x * r.z + l.z * r.x, l.x * r.y - l.y * r.x);
}

inline double operator*(const Vector &l, const Vector &r) {
    return l.x * r.x + l.y * r.y + l.z * r.z;
}

inline Vector operator*(const double &l, const Vector &r) {
    return Vector(l * r.x, l * r.y, l * r.z);
}

inline Vector operator*(const Vector &r, const double &l) {
    return Vector(l * r.x, l * r.y, l * r.z);
}

inline Vector operator/(const Vector &r, const double &l) {
    return Vector(r.x / l, r.y / l, r.z / l);
}
