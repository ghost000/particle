#pragma once

#include "ofMain.h"

class orthoCamera : public ofCamera {
  public:
    orthoCamera();

    void begin(ofRectangle rect = ofGetWindowRect());

    float scale;
};
