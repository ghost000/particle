//
// Created by neo on 11/28/17.
//
#pragma once

#include "Particle.h"
#include "ofMain.h"

class Physics : public ofNode {
public:
    Physics();

    ~Physics();

    void Initialize();

    void CalcForces(std::vector <std::vector<Particle>> &);

    void StepSimulation(const double &);

    void ReleaseLockedParticles();

    void SetWindVector(const double &, const double &, const double &);

    void SetWindForceFactor(const double &);

    void UpdateClothGeometry();

    void CopyParticles(const std::vector <std::vector<Particle>> &, std::vector <std::vector<Particle>> &);

    double CheckForCollisions(const std::vector <std::vector<Particle>> &);

    void ResolveCollisions(std::vector <std::vector<Particle>> &);

    int GetParticleSize1();

    int GetParticleSize2(const size_t &);

    int tb_Rnd(const int &min, const int &max);

    std::vector <std::vector<Particle>> Particles;

    void customDraw();

    double GetWind();

    void NextFag();

    ofLight light;

private:
    std::vector<double> ClothFaces;
    std::vector<double> ClothVertices;
    std::vector <Spring> Springs;
    std::vector <Collision> Collisions;
    Vector WindVector;
    double wind;
    int flag;
    std::vector <std::vector<int>> flags;

};
