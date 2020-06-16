import java.util.Random;

class GameState {
  int holeNum = 0;
  int threshold = 5;
  int[] holeState; 
  int[] holeExistence;

  // list of int state {0,1,2}, 
  // 0 means Mole is hiding (popdown), 
  // 1 means Mole is appears (popup), 
  // 2 means Mole is Hit
  int[] moleState; 
  int[] moleExistence;

  int point;

  int score = 0;

  int gameDuration= 20;

  HashMap<Integer, PMatrix3D> markerPoseMap;
  HashMap<Integer, PMatrix3D> frameMarkerPoseMap;

  GameState (ArrayList<Marker> markers) {
    markerPoseMap = new HashMap<Integer, PMatrix3D>();
    frameMarkerPoseMap = new HashMap<Integer, PMatrix3D>();
    for (int i = 0; i < markers.size(); i++) {
      Marker m = markers.get(i);
      markerPoseMap.put(m.code, m.pose);
      // frameMarkerPoseMap.put(m.code, m.pose);
    }

    holeNum = markerPoseMap.size();
    holeExistence = new int[holeNum];
    for (int j = 0; j < holeExistence.length; j++){
      Marker m = markers.get(j);
      holeExistence[j] = m.code; 
    }

    holeState = new int[holeNum];

    moleExistence = new int[holeNum];
    moleState = new int[holeNum];
  }

  void resetState(){
    markerPoseMap.clear();
    holeNum = 0;
  }


  void updateHoleState(ArrayList<Marker> markers){
    // don't clear markerPoseMap to display the current/last marker pose
    // markerPoseMap.clear();

    // to detect hit using marker loss, clear frameMarkerPose instead;
    frameMarkerPoseMap.clear();
    for (int i = 0; i < markers.size(); i++) {
      Marker m = markers.get(i);
      markerPoseMap.put(m.code, m.pose);
      frameMarkerPoseMap.put(m.code, m.pose);
    }
    // for each marker, update loss interval if the marker detection lost.
    for (int i = 0; i < holeExistence.length; i++) {
      if (holeHitOnMarkerLoss(i) || holeHitOnClick(i)){
        holeState[i] += 1;
      }
      else {
        //moleState[i] = 1; @Fariz it makes molestate always be 1 and can't popdown the mole
        holeState[i] = 0;
      }

      if (holeState[i] > threshold && moleExistence[i] == 1) {
        moleState[i] = 2;
        point += 1;
      }
    }
    debugDisplay(holeExistence, holeState);
  }

  void updateMoleExistence(){
    // for each mole, check if any mole appear in the mole boolean state
    for (int i = 0; i < moleState.length; i++) {
      // When moleState is Appear and holeState is still normal (not loss), update moleExistence
      if (moleState[i] == 1 && holeState[i] <= threshold) {
        moleExistence[i] = 1;
      }
      else {
        moleExistence[i] = 0;
      }
    }
    // println(Arrays.toString(array));
    printArray(moleExistence);
  }

  void drawGame(){
    // The hole will be always animated, so loop for detected marker to generate hole
    for (int i = 0; i < holeExistence.length; i++) {
      // current marker position to draw
      PMatrix3D pose_this = markerPoseMap.get(holeExistence[i]);
      // pose current marker to another marker direction
      PMatrix3D pose_look = markerPoseMap.get(holeExistence[(i+1)%2]);

      if (pose_this == null || pose_look == null)
        continue;

      float angle = rotateToMarker(pose_this, pose_look, holeExistence[i]);

      pushMatrix();
        applyMatrix(pose_this);
        // rotateZ(angle);
        rotateZ(angle);
        drawHole(snowmanSize);
      popMatrix();
    }
  }

  void molePopUp(int moleIndex, long passTime){
    // PMatrix3D pose_this = markerPoseMap.get(ExistenceList[moleIndex]);
    // PMatrix3D pose_look = markerPoseMap.get(ExistenceList[(moleIndex+1)%2]);
    PMatrix3D pose_this = markerPoseMap.get(holeExistence[moleIndex]);
    PMatrix3D pose_look = markerPoseMap.get(holeExistence[(moleIndex+1)%2]);
    
    if (pose_this != null && pose_look != null ){
      float angle = rotateToMarker(pose_this, pose_look, holeExistence[moleIndex]);
      pushMatrix();
        // apply matrix (cf. drawPrimitives.pde)
        applyMatrix(pose_this);
        rotateZ(angle);
        // globalRotateAngle += 0.005;
        // rotateZ(globalRotateAngle);

        // draw snowman
        // println("moleKey : " + moleIndex);

        if (moleHitDebug || moleState[moleIndex] == 2){
            angle -= 40;
            rotateZ(angle);
            // println("mole hit : "+moleHitDebug);
            drawMole(calcMoleSize(moleIndex, passTime),1);
        } else {
            drawMole(calcMoleSize(moleIndex, passTime),0);
        }

        noFill();
        strokeWeight(3);
        stroke(255, 0, 0);
        line(0, 0, 0, 0.02, 0, 0); // draw x-axis
        stroke(0, 255, 0);
        line(0, 0, 0, 0, 0.02, 0); // draw y-axis
        stroke(0, 0, 255);
        line(0, 0, 0, 0, 0, 0.02); // draw z-axis
      popMatrix();
    }
  }

  void molePopDown(int moleIndex){
    moleState[moleIndex] = 0;
  }

  /**
  * Param int @i : index of checked marker still detected
  **/
  boolean holeHitOnMarkerLoss(int i){
    // return markerPoseMap.get(holeExistence[i]) == null;
    return frameMarkerPoseMap.get(holeExistence[i]) == null;
  }

  boolean holeHitOnClick(int i){
    if (markerPoseMap.get(holeExistence[i]) == null){
        // return true same as holeHitByMarkerLoss because cannot get the PMatrix3D
        return true;
    } else {
        // TODO check if x and y collide in 2D Projection of marker i PMatrix3D
        PMatrix3D holePose = markerPoseMap.get(holeExistence[i]);
        
        PVector camera = new PVector(-width/2, -height/2,-(height/2)/tan(radians(fov)));
        PVector result = new PVector();
        holePose.mult(camera, result);
        println(result.x, result.y);
        float x_proj = screenX(result.x, result.y, result.z);
        float y_proj = screenY(result.x, result.y, result.z);
        // println(x_proj, y_proj);
        // by default return false
    }
    return false;
  }

  /**
  * TODO holeHitOnKeyboard Tab are for debugging only until mouse click is available
  **/
  boolean holeHitOnKeyboard(int i){
    if (moleHitDebug){
      return true;
    }
    else {
      return false;
    }
  }

  /**
  * update mole at corresponding hole index
  * @param holePos : 
  **/
  void updateMoleState(int holeIndex, int state){
    moleState[holeIndex] = state;
  }

  int getMoleState(int holeIndex){
    return moleState[holeIndex];
  }

  int getMoleExistence(int holeIndex){
    return moleExistence[holeIndex];
  }

  int getNumberofHole(){
    return holeExistence.length;
  }

  int getScore(){
      return score;
  }

  void addScore(){
      score += 10;
  }

  void minusScore(){
      score -= 5;
  }

  boolean timeup(long passTime, long dueTime){
      if(passTime >= dueTime) return true;
      return false;
  }

  int randAppearDuration(){
      return rand.nextInt(2000)+1000;
  }

  int randHideDuration(){
      return rand.nextInt(3000)+2000;
  }

  float calcMoleSize(int moleIndex, long passTime){
    float moleSize = 0;
    int frame = int(passTime/(moleAppearDuration[moleIndex]/40));

    switch(getMoleState(moleIndex)){
      case 0:
        moleSize = -0.00005*frame*frame+0.002*frame;
        break;
      case 1:
        moleSize = -0.00005*frame*frame+0.002*frame;
        break;
      case 2:
        frame += 20;
        moleSize = -0.00005*frame*frame+0.002*frame;
        break;
      default:
        moleSize = 0;
        break;
    }
    if(moleSize<0) moleSize = 0;
    return moleSize;
  }
}

