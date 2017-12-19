#pragma once
#include "ofMain.h"

// Custom objects for this example
#include "Physics.cpp"
#include "Grid.h"
#include "OrthoCamera.h"

#define N_CAMERAS 4

class ofApp : public ofBaseApp {

  public:
    void setup();

    void update();

    void draw();

    void setupViewports();

    void drawScene(int iCameraDraw);

    void updateMouseRay();

    void keyPressed(int key);

    void keyReleased(int key);

    void mouseMoved(int x, int y);

    void mouseDragged(int x, int y, int button);

    void mousePressed(int x, int y, int button);

    void mouseReleased(int x, int y, int button);

    void mouseEntered(int x, int y);

    void mouseExited(int x, int y);

    void windowResized(int w, int h);

    //cameras (all these inherit from ofCamera)
    ofEasyCam camEasyCam;
    orthoCamera camFront;
    orthoCamera camTop;
    orthoCamera camLeft;

    //cameras have parent?
    bool bCamParent;

    //camera pointers
    ofCamera *cameras[N_CAMERAS];
    int iMainCamera;

    //viewports
    ofRectangle viewMain;
    ofRectangle viewGrid[N_CAMERAS];

    Physics pysics;
    grid nodeGrid;

    //ray drawn under mouse cursor [start,end]
    ofVec3f ray[2];
};
