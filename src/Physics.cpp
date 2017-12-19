//
// Created by neo on 11/28/17.
//
#include "Physics.h"

inline Physics::Physics() : Particles(NUMROWS + 1, std::vector<Particle>(NUMCOLUMNS + 1)),
    Springs(NUMSTRUCTURALSPRINGS), Collisions(NUMVERTICES), wind(WindForceFactor), flag(0),
    flags{
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0},
    {0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0},
    {0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0},
    {0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0},
    {0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0},
    {0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0},
    {0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0},
    {0,0,0,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0},
    {0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,1,0,1,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}} {
    light.setAmbientColor(ofColor(0, 0, 0));
}

inline Physics::~Physics() = default;

inline void Physics::Initialize() {

    double f;
    Vector L;
    int count;

    ofColor color;
    wind = WindForceFactor;

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {

            if ((r == 0) && (c == 0))
                f = 1;
            else if ((r == NUMROWS) && (c == 0))
                f = 2;
            else if ((r == 0) && (c == NUMCOLUMNS))
                f = 2;
            else if ((r == NUMROWS) && (c == NUMCOLUMNS))
                f = 1;
            else if (((r == 0) || (r == NUMROWS)) && ((c != 0) && (c != NUMCOLUMNS)))
                f = 3;
            else
                f = 6;
            if( flag == 0 ) {
                color.r = ofRandom(255);
                color.g = ofRandom(255);
                color.b = ofRandom(255);
                color.a = 255;
            }
            if( flag == 1 ) {
                if(r < NUMROWS / 2) {

                    color.r = 255;
                    color.g = 255;
                    color.b = 255;
                } else {

                    color.r = 255;
                    color.g = 0;
                    color.b = 0;
                }
                color.a = 255;
            }
            if( flag == 2 ) {
                if(flags.at(r).at(c) == 0 ) {

                    color.r = 255;
                    color.g = 255;
                    color.b = 255;
                } else {

                    color.r = 0;
                    color.g = 0;
                    color.b = 0;
                }
                color.a = 255;
            }
            if( flag == 3 ) {
                if(flags.at(r).at(c) == 1 ) {

                    color.r = 255;
                    color.g = 255;
                    color.b = 255;
                } else {

                    color.r = 0;
                    color.g = 0;
                    color.b = 0;
                }
                color.a = 255;
            }

            Particles.at(r).at(c).Mass = (f * MASSPERFACE) / 3;
            Particles.at(r).at(c).InvMass = 1 / Particles.at(r).at(c).Mass;

            Particles.at(r).at(c).Position.x = c * CSTEP;
            Particles.at(r).at(c).Position.y = (CLOTHHEIGHT - (r * RSTEP)) + YOFFSET;
            Particles.at(r).at(c).Position.z = 0.0;

            Particles.at(r).at(c).Velocity.x = 0.0;
            Particles.at(r).at(c).Velocity.y = 0.0;
            Particles.at(r).at(c).Velocity.z = 0.0;

            Particles.at(r).at(c).Acceleration.x = 0.0;
            Particles.at(r).at(c).Acceleration.y = 0.0;
            Particles.at(r).at(c).Acceleration.z = 0.0;

            Particles.at(r).at(c).Forces.x = 0.0;
            Particles.at(r).at(c).Forces.y = 0.0;
            Particles.at(r).at(c).Forces.z = 0.0;
            Particles.at(r).at(c).Locked = (c == 0) && (r == 0 || r == NUMROWS);
            Particles.at(r).at(c).color = color;
        }
    }

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            ClothVertices.emplace_back(Particles.at(r).at(c).Position.x);
            ClothVertices.emplace_back(Particles.at(r).at(c).Position.y);
            ClothVertices.emplace_back(Particles.at(r).at(c).Position.z);
        }
    }
    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            ClothVertices.emplace_back(Particles.at(r).at(c).Position.x);
            ClothVertices.emplace_back(Particles.at(r).at(c).Position.y);
            ClothVertices.emplace_back(Particles.at(r).at(c).Position.z);
        }
    }


    for (auto r = 0; r < NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            if (c == 0) {
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + c); // vertex 1
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + (c + 1)); // vertex 2
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + c); // vertex 3
            } else if (c == NUMCOLUMNS) {
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + c); // vertex 1
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + c); // vertex 2
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + (c - 1)); // vertex 3
            } else {
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + c); // vertex 1
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + c); // vertex 2
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + (c - 1)); // vertex 3

                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + c); // vertex 1
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + (c + 1)); // vertex 2
                ClothFaces.emplace_back(((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + c); // vertex 3
            }
        }
    }

    for (auto r = 0; r < NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            if (c == 0) {
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + c); // vertex 3
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + (c + 1)); // vertex 2
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + c); // vertex 1
            } else if (c == NUMCOLUMNS) {
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + (c - 1)); // vertex 3
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + c); // vertex 2
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + c); // vertex 1
            } else {
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + (c - 1)); // vertex 3
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + c); // vertex 2
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + c); // vertex 1

                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + (NUMCOLUMNS + 1) + c); // vertex 3
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + (c + 1)); // vertex 2
                ClothFaces.emplace_back(NUMVERTICES + ((NUMCOLUMNS + 1) * r) + c); // vertex 1
            }
        }
    }

    count = 0;

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            if (c < NUMCOLUMNS) {
                Springs.at(count).p1.r = r;
                Springs.at(count).p1.c = c;
                Springs.at(count).p2.r = r;
                Springs.at(count).p2.c = c + 1;
                Springs.at(count).k = SPRINGTENSIONCONSTANT;
                Springs.at(count).d = SPRINGDAMPINGCONSTANT;
                L = Particles.at(r).at(c).Position - Particles.at(r).at(c + 1).Position;
                Springs.at(count).L = L.Magnitude();
                count++;
            }
            if (r < NUMROWS) {
                Springs.at(count).p1.r = r;
                Springs.at(count).p1.c = c;
                Springs.at(count).p2.r = r + 1;
                Springs.at(count).p2.c = c;
                Springs.at(count).k = SPRINGTENSIONCONSTANT;
                Springs.at(count).d = SPRINGDAMPINGCONSTANT;
                L = Particles.at(r).at(c).Position - Particles.at(r + 1).at(c).Position;
                Springs.at(count).L = L.Magnitude();
                count++;
            }
            if (r < NUMROWS && c < NUMCOLUMNS) {
                Springs.at(count).p1.r = r;
                Springs.at(count).p1.c = c;
                Springs.at(count).p2.r = r + 1;
                Springs.at(count).p2.c = c + 1;
                Springs.at(count).k = SPRINGSHEARCONSTANT;
                Springs.at(count).d = SPRINGDAMPINGCONSTANT;
                L = Particles.at(r).at(c).Position - Particles.at(r + 1).at(c + 1).Position;
                Springs.at(count).L = L.Magnitude();
                count++;
            }
            if (c > 0 && r < NUMROWS) {
                Springs.at(count).p1.r = r;
                Springs.at(count).p1.c = c;
                Springs.at(count).p2.r = r + 1;
                Springs.at(count).p2.c = c - 1;
                Springs.at(count).k = SPRINGSHEARCONSTANT;
                Springs.at(count).d = SPRINGDAMPINGCONSTANT;
                L = Particles.at(r).at(c).Position - Particles.at(r + 1).at(c - 1).Position;
                Springs.at(count).L = L.Magnitude();
                count++;
            }
        }
    }

    WindVector.x = 10.0;
    WindVector.y = 0.0;
    WindVector.z = 1.0;
}

inline int Physics::tb_Rnd(const int &min, const int &max) {
    int number= {(std::abs(std::rand()) % (max - min + 1)) + min};

    if (number > max) {
        number = max;
    }

    if (number < min) {
        number = min;
    }

    return number;
}

inline void Physics::CalcForces(std::vector <std::vector<Particle>> &particles) {
    Vector dragVector;
    Vector f1, f2, d, v;
    double L;

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            particles.at(r).at(c).Forces.x = 0;
            particles.at(r).at(c).Forces.y = 0;
            particles.at(r).at(c).Forces.z = 0;
        }
    }

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {

            if (!particles.at(r).at(c).Locked) {
                particles.at(r).at(c).Forces.y += GRAVITY * particles.at(r).at(c).Mass;

                dragVector = -particles.at(r).at(c).Velocity;
                dragVector.Normalize();
                particles.at(r).at(c).Forces += dragVector * (particles.at(r).at(c).Velocity.Magnitude() *
                                                particles.at(r).at(c).Velocity.Magnitude()) *
                                                DRAGCOEFFICIENT;

                SetWindVector(tb_Rnd(0, 10), 0, tb_Rnd(0, 1));
                WindVector.Normalize();
                particles.at(r).at(c).Forces += WindVector * tb_Rnd(0, wind);
            }
        }
    }

    for (auto i = 0; i < NUMSTRUCTURALSPRINGS; i++) {
        auto r1 = Springs.at(i).p1.r;
        auto c1 = Springs.at(i).p1.c;
        auto r2 = Springs.at(i).p2.r;
        auto c2 = Springs.at(i).p2.c;

        d = particles.at(r1).at(c1).Position - particles.at(r2).at(c2).Position;
        v = particles.at(r1).at(c1).Velocity - particles.at(r2).at(c2).Velocity;
        L = Springs.at(i).L;

        f1 = -(Springs.at(i).k * (d.Magnitude() - L) + Springs.at(i).d * ((v * d) / d.Magnitude())) *
             (d / d.Magnitude());
        f2 = -f1;

        if (!particles.at(r1).at(c1).Locked)
            particles.at(r1).at(c1).Forces += f1;

        if (!particles.at(r2).at(c2).Locked)
            particles.at(r2).at(c2).Forces += f2;
    }
}

inline void Physics::StepSimulation(const double &dt) {
    Vector Ae;
    std::vector <std::vector<Particle>> p(NUMROWS + 1, std::vector<Particle>(NUMCOLUMNS + 1));
    double dtime = dt;
    bool tryAgain = true;
    double check = 0;
    bool didPen = false;


    while (tryAgain && (dtime > tol)) {
        tryAgain = false;
        CopyParticles(Particles, p);

        // calculate all of the forces
        CalcForces(p);

        for (auto r = 0; r <= NUMROWS; r++) {
            for (auto c = 0; c <= NUMCOLUMNS; c++) {
                Ae = p.at(r).at(c).Forces * p.at(r).at(c).InvMass;
                p.at(r).at(c).Acceleration = Ae;
                p.at(r).at(c).Velocity += Ae * dtime;
                p.at(r).at(c).Position += p.at(r).at(c).Velocity * dtime;
            }
        }

        // check for collisions
        check = CheckForCollisions(p);

        if (check == PENETRATING) {
            dtime = dtime / 2;
            tryAgain = true;
            didPen = true;
        } else if (check == COLLISION) {
            ResolveCollisions(p);
            didPen = false;
        }
    }

    CopyParticles(p, Particles);
    // update cloth object's geometry
    UpdateClothGeometry();

}

inline void Physics::ReleaseLockedParticles() {

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            Particles.at(r).at(c).Locked = false;
        }
    }
}

inline void Physics::SetWindVector(const double &x, const double &y, const double &z) {
    WindVector.x = x;
    WindVector.y = y;
    WindVector.z = z;
}

inline void Physics::SetWindForceFactor(const double &f) {
    wind += f;
}

inline void Physics::UpdateClothGeometry() {
    // fill the vertex array
    auto vertices = ClothVertices.begin();
    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            // setup vertices
            *vertices = Particles.at(r).at(c).Position.x;
            vertices++;
            *vertices = Particles.at(r).at(c).Position.y;
            vertices++;
            *vertices = Particles.at(r).at(c).Position.z;
            vertices++;
        }
    }
    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            // setup vertices
            *vertices = Particles.at(r).at(c).Position.x;
            vertices++;
            *vertices = Particles.at(r).at(c).Position.y;
            vertices++;
            *vertices = Particles.at(r).at(c).Position.z;
            vertices++;
        }
    }
}

inline void
Physics::CopyParticles(const std::vector <std::vector<Particle>> &src, std::vector <std::vector<Particle>> &dst) {

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            dst.at(r).at(c).Mass = src.at(r).at(c).Mass;
            dst.at(r).at(c).InvMass = src.at(r).at(c).InvMass;
            dst.at(r).at(c).Position = src.at(r).at(c).Position;
            dst.at(r).at(c).Velocity = src.at(r).at(c).Velocity;
            dst.at(r).at(c).Forces = src.at(r).at(c).Forces;
            dst.at(r).at(c).Locked = src.at(r).at(c).Locked;
        }
    }

}

inline double Physics::CheckForCollisions(const std::vector <std::vector<Particle>> &p) {
    int count{0};
    auto state{NOCOLLISION};
    double d;
    Vector n;
    double Vn;

    for (auto i = 0; i < NUMVERTICES; i++) {
        Collisions.at(i).p1.r = -1;
        Collisions.at(i).p1.c = -1;
        Collisions.at(i).n.x = 0;
        Collisions.at(i).n.y = 0;
        Collisions.at(i).n.z = 0;
    }

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            if (!p.at(r).at(c).Locked) {
                if ((p.at(r).at(c).Position.y <= COLLISIONTOLERANCE) && (p.at(r).at(c).Velocity.y < tol)) {
                    state = COLLISION;
                    Collisions.at(count).p1.r = r;
                    Collisions.at(count).p1.c = c;
                    Collisions.at(count).n.x = 0.0f;
                    Collisions.at(count).n.y = 1.0f;
                    Collisions.at(count).n.z = 0.0f;
                    Collisions.at(count).n.Normalize();
                    count++;
                }
            }
        }
    }

    for (auto r = 0; r <= NUMROWS; r++) {
        for (auto c = 0; c <= NUMCOLUMNS; c++) {
            if (!p.at(r).at(c).Locked) {
                d = sqrt(p.at(r).at(c).Position.x * p.at(r).at(c).Position.x +
                         p.at(r).at(c).Position.z * p.at(r).at(c).Position.x);
                n.x = p.at(r).at(c).Position.x;
                n.y = 0.0f;
                n.z = p.at(r).at(c).Position.z;
                Vn = (p.at(r).at(c).Velocity * n);

                if ((d <= (FLAGPOLERADIUS + COLLISIONTOLERANCE)) &&
                        (p.at(r).at(c).Position.y < FLAGPOLEHEIGHT) &&
                        (p.at(r).at(c).Position.y > 0.0f) &&
                        (Vn < tol)
                   ) {
                    state = COLLISION;
                    Collisions.at(count).p1.r = r;
                    Collisions.at(count).p1.c = c;
                    Collisions.at(count).n = n;
                    Collisions.at(count).n.Normalize();
                    count++;
                }
            }
        }
    }

    return state;
}

inline void Physics::ResolveCollisions(std::vector <std::vector<Particle>> &p) {

    for (auto i = 0; i < NUMVERTICES; i++) {
        if (Collisions.at(i).p1.r != -1) {
            auto r = Collisions.at(i).p1.r;
            auto c = Collisions.at(i).p1.c;
            Vector Vn = (Collisions.at(i).n * p.at(r).at(c).Velocity) * Collisions.at(i).n;
            Vector Vt = p.at(r).at(c).Velocity - Vn;
            p.at(r).at(c).Velocity = (-(KRESTITUTION + 1) * Vn) + (FRICTIONFACTOR * Vt);
        }
    }
}

inline int Physics::GetParticleSize1() {
    return Particles.size();
}

inline int Physics::GetParticleSize2(const size_t &count) {
    return Particles.at(count).size();
}

inline void Physics::customDraw() {
    // We run the update ourselves manually. ofNode does
    //  not do this for us.
    //update();



    //--
    // Draw particles

    // We use the position of the first
    //  particle as the position of the
    //  light.
    ofPushStyle();
    light.enable();
    light.setPosition(Particles.at(0).at(0).Position);

    for (auto i = 0; i < GetParticleSize1(); i++) {
        for (auto j = 0; j < GetParticleSize2(i); j++) {
            ofPushStyle();
            ofSetColor(Particles.at(i).at(j).color);

            ofDrawSphere(Particles.at(i).at(j).Position.x,
                         Particles.at(i).at(j).Position.y,
                         Particles.at(i).at(j).Position.z, 4.0);

            ofPopStyle();
        }
    }

    light.disable();
    ofDisableLighting();

    //
    //--



    // Render light as white sphere
    ofSetColor(255, 255, 255);
    ofDrawSphere(light.getPosition(), 2.0);
    ofSetDrawBitmapMode(OF_BITMAPMODE_MODEL);
    ofDrawBitmapString(" light", 
					   Particles.at(0).at(0).Position.x,
                       Particles.at(0).at(0).Position.y,
                       Particles.at(0).at(0).Position.z);
    ofPopStyle();
}

inline double Physics::GetWind() {
    return wind;
}

inline void Physics::NextFag() {
    flag++;
    if(flag > 3)
        flag = 0;
}
