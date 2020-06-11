
class GameState {
    int holeNum = 0;
    int threshold = 30;
    int[] HoleState; 
    int[] HoleExistence;

    boolean[] MoleState; 
    int[] MoleExistence;

    int point;

    HashMap<Integer, PMatrix3D> markerPoseMap;

    GameState (ArrayList<Marker> markers) {
        markerPoseMap = new HashMap<Integer, PMatrix3D>();
        for (int i = 0; i < markers.size(); i++) {
            Marker m = markers.get(i);
            markerPoseMap.put(m.code, m.pose);
        }

        holeNum = markerPoseMap.size();
        HoleState = new int[holeNum];
        HoleExistence = new int[holeNum];

        MoleState = new boolean[holeNum];
        MoleExistence = new int[holeNum];
    }

    void resetState(){
        markerPoseMap.clear();
        holeNum = 0;
    }


    void updateHoleExistence(ArrayList<Marker> markers){
        markerPoseMap.clear();
        for (int i = 0; i < markers.size(); i++) {
            Marker m = markers.get(i);
            markerPoseMap.put(m.code, m.pose);
        }
        // for each marker, update loss interval if the marker detection lost.
        for (int i = 0; i < HoleExistence.length; i++) {
            if (markerPoseMap.get(HoleExistence[i]) == null){
                HoleState[i] += 1;
            }
            else {
                HoleState[i] = 0;
            }
        }
    }

    void updateMoleExistence(){
        // for each mole, check if any mole appear in the mole boolean state
        for (int i = 0; i < MoleState.length; i++) {
            if (MoleState[i] && HoleState[i] == 0) {
                MoleExistence[i] += 1;
            }
            else {
                MoleExistence[i] = 0;
            }
        }
    }

    boolean mouseOverHole(){
        return true;
    }

    /**
    * update mole at corresponding hole index
    * @param holePos : 
    **/
    void updateMoleState(int holeIndex){
        
    }

    boolean getMoleState(int holeIndex){
        return true;
    }

    boolean getMoleExistence(int holeIndex){
        return true;
    }

}

