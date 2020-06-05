// sample code for keyevents
// currently we manipulate binary thresholds of the marker tracker

class KeyState {
  HashMap<Integer, Boolean> key;

  KeyState() {
    key = new HashMap<Integer, Boolean>();

    key.put(RIGHT, false);
    key.put(LEFT,  false);
    key.put(UP,    false);
    key.put(DOWN,  false);
  }

  void putState(int code, boolean state) {
    key.put(code, state);
  }

  boolean getState(int code) {
    return key.get(code);
  }

  void getKeyEvent() {
    if (getState(LEFT)) {
      markerTracker.thresh -= 1;
    }

    if (getState(RIGHT)) {
      markerTracker.thresh += 1;
    }

    if (getState(UP)) {
      markerTracker.bw_thresh += 1;
    }

    if (getState(DOWN)) {
      markerTracker.bw_thresh -= 1;
    }
  }
}

void keyPressed() {
  keyState.putState(keyCode, true);
}

void keyReleased() {
  keyState.putState(keyCode, false);
}