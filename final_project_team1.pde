import gab.opencv.*;
import processing.video.*;
import java.util.ArrayList;

import picking.*;

Picker picker;

// final boolean MARKER_TRACKER_DEBUG = false;
boolean MARKER_TRACKER_DEBUG = true;

final boolean BALL_DEBUG = false;

final boolean USE_SAMPLE_IMAGE = true;

// We've found that some Windows build-in cameras (e.g. Microsoft Surface)
// cannot work with processing.video.Capture.*.
// Instead we use DirectShow Library to launch these cameras.
final boolean USE_DIRECTSHOW = false;


// final double kMarkerSize = 0.036; // [m]
final double kMarkerSize = 0.024; // [m]

Capture cap;
DCapture dcap;
OpenCV opencv;

// Variables 
// **************************************************************
float fov = 45; // for camera capture

// Marker codes to draw snowmans example
// final int[] ExistenceList = {0x1228, 0x0690};
// int towards = 0x1228; // the target marker that the ball flies towards
int towardscnt = 0;   // if ball reached, +1 to change the target

int[] ExistenceList;
int[] ExistenceState;
int point = 0;
int towards = 0x005A; //for throwing ball example
boolean calibration_boolean = false;

final float GA = 9.80665;

PVector snowmanLookVector;
PVector ballPos;
float ballAngle = 45;
float ballspeed = 0;

final int ballTotalFrame = 10;
final float snowmanSize = 0.020;
int frameCnt = 0;

// HashMap<Integer, PMatrix3D> markerPoseMap;

MarkerTracker markerTracker;
PImage img;

KeyState keyState;

int moleKeyDebug = -1;
boolean moleHitDebug = false;

float globalRotateAngle = 0;

GameState gameState;


// Circle button setup to start the game
// from https://processing.org/examples/button.html
int circleX, circleY;  // Position of circle button
int circleSize = 93;   // Diameter of circle
color circleColor, baseColor;
color circleHighlight;
color currentColor;
boolean circleOver = false;

void selectCamera() {
  String[] cameras = Capture.list();

  if (cameras == null) {
    println("Failed to retrieve the list of available cameras, will try the default");
    cap = new Capture(this, 640, 480);
  } else if (cameras.length == 0) {
    println("There are no cameras available for capture.");
    exit();
  } else {
    println("Available cameras:");
    printArray(cameras);

    // The camera can be initialized directly using an element
    // from the array returned by list():
    // cap = new Capture(this, cameras[5]);

    // Or, the settings can be defined based on the text in the list
    // cap = new Capture(this, 1280, 720, "USB2.0 HD UVC WebCam", 30);
    cap = new Capture(this, 1280, 720, 10);
  }
}

void setupCamera() {
  if (!USE_SAMPLE_IMAGE) {
    selectCamera();
    opencv = new OpenCV(this, cap.width, cap.height);
  }
}

void setupButton(){
  circleColor = color(0);
  circleHighlight = color(204);
  baseColor = color(102);
  currentColor = baseColor;
  circleX = opencv.width/2+circleSize/2;
  circleY = opencv.height/2;
  ellipseMode(CENTER);
}


void settings() {
  if (USE_SAMPLE_IMAGE) {
    // Here we introduced a new test image in Lecture 6 (20/05/27)
    size(1280, 720, P3D);
    opencv = new OpenCV(this, "./marker_test2.jpg");
    // size(1000, 730, P3D);
    // opencv = new OpenCV(this, "./marker_test.jpg");
  } else {
    if (USE_DIRECTSHOW) {
      dcap = new DCapture();
      size(dcap.width, dcap.height, P3D);
      opencv = new OpenCV(this, dcap.width, dcap.height);
    } else {
      selectCamera();
      size(cap.width, cap.height, P3D);
      opencv = new OpenCV(this, cap.width, cap.height);
    }
  }
}

void setup() {
  background(0);
  smooth();
  // frameRate(10);

  markerTracker = new MarkerTracker(kMarkerSize);

  // if (!USE_DIRECTSHOW){
  //   cap.start();
  // }
  if (!USE_DIRECTSHOW && !USE_SAMPLE_IMAGE) {
    setupCamera();
    cap.start();
  }

  img = createImage(opencv.width, opencv.height, RGB);

  // picker = new Picker(this);

  // draw circle button as start button
  setupButton();

  // Align the camera coordinate system with the world coordinate system
  // (cf. drawSnowman.pde)
  PMatrix3D cameraMat = ((PGraphicsOpenGL)g).camera;
  cameraMat.reset();

  keyState = new KeyState();

  ballPos = new PVector();  // ball position
  // markerPoseMap = new HashMap<Integer, PMatrix3D>();  // hashmap (code, pose)

  calibration_boolean = false;
}

void draw() {

  // if we use video camera
  if (!USE_SAMPLE_IMAGE) {
    if (USE_DIRECTSHOW) {
      img = dcap.updateImage();
      opencv.loadImage(img);
    } else {
      if (cap.width <= 0 || cap.height <= 0) {
        println("Incorrect capture data. continue");
        return;
      }
      opencv.loadImage(cap);
    }
  }

  // setup marker detection
  ArrayList<Marker> markers = new ArrayList<Marker>();
  markers.clear();
  // markers = calibration(markers);
  calibration(markers);


  // calibration case
  if (calibration_boolean == false) {
    keyState.getKeyEvent();

    // startButton
    update(mouseX, mouseY);
    if (circleOver) {
      fill(circleHighlight);
    } else {
      fill(circleColor);
    }
    // TODO can someone help to make an ellipse appear to indicate button presence
    strokeWeight(3);
    stroke(255, 255, 255);
    ellipse(circleX, circleY, circleSize, circleSize);

    println("In calibration case");

    println("Number of Markers detected in this calibration:" + markers.size());

    // if(mousePressed == true && circleOver){
    if (key == ENTER && markers.size() > 0) {
      gameState = new GameState(markers);
      println("finish calibration");
      calibration_boolean = true;
      MARKER_TRACKER_DEBUG = false;
    }
    System.gc();
  }
  //game start case
  else {
    println("point:" + point);

    // use perspective camera
    perspective(radians(fov), float(width)/float(height), 0.01, 1000.0);

    // setup light
    // (cf. drawSnowman.pde)
    ambientLight(180, 180, 180);
    directionalLight(180, 150, 120, 0, 1, 0);
    lights();

    println("mouseX : " + mouseX + " , mouseY : " + mouseY);
    gameState.updateHoleState(markers);

    gameState.drawGame();


    // TODO @Daphne modifed from here to generate a function to call molePopUp
    if (key != TAB){
      moleKeyDebug = int(key) % markers.size();
    }
    gameState.updateMoleExistence();
    gameState.molePopUp(moleKeyDebug);

    // int id = picker.get(mouseX, mouseY);
    // println(id);
    // if(contains(ExistenceList, id)){
    //   println("mouse over hole " + id);
    // }
  
    // Your Code for Homework 6 (20/06/03) - End
    // **********************************************

    noLights();
    keyState.getKeyEvent();

    System.gc();
  }
}

void update(int x, int y) {
  if ( overCircle(circleX, circleY, circleSize) ) {
    circleOver = true;
  } else {
    circleOver = false;
  }
  // println(circleOver);
}

boolean overCircle(int x, int y, int diameter) {
  float disX = x - mouseX;
  float disY = y - mouseY;
  if (sqrt(sq(disX) + sq(disY)) < diameter/2 ) {
    return true;
  } else {
    return false;
  }
}

boolean overObject(PMatrix3D thisObject) {
  // Someone implement this
  PVector relativeVector = new PVector();
  relativeVector.x = thisObject.m03 - mouseX;
  relativeVector.y = thisObject.m13 - mouseY;
  relativeVector.z = thisObject.m23 ;
  float relativeLen = relativeVector.mag();

  return true;
}

void captureEvent(Capture c) {
  PGraphics3D g;
  if (!USE_DIRECTSHOW && c.available())
      c.read();
}

