//all method used for first calibration
import java.util.Arrays;
import java.util.Iterator;

ArrayList<Integer> calibration(ArrayList<Integer> ExistenceList) {
    ArrayList<Marker> markers_calibrate = new ArrayList<Marker>();
    ortho();
    pushMatrix();
    translate(-width/2, -height/2,-(height/2)/tan(radians(fov)));
    markerTracker.findMarker(markers_calibrate);
    popMatrix();
//   println("marker_calibrate.size in loop:" + markers_calibrate.size());
//   println("ExistenceList length in loop:" + ExistenceList.length);
  for (int i = 0; i < markers_calibrate.size(); i++) {
      Marker m = markers_calibrate.get(i);
      ExistenceList.add(m.code);
    //   println("ExistenceList.length:" + ExistenceList.length);
  }
  return ExistenceList;
}

ArrayList<Integer> initializeState(ArrayList<Integer> ExistenceState, ArrayList<Integer> ExistenceList){
     //int in the array have values 0,1,..,5 ( == get point when marker detection loss 5 iteration). each int correspond to each matrix in ExistenceList with same index.
    for (int i = 0; i < ExistenceList.size(); i++) {
      ExistenceState.add(0);
    }
  return ExistenceState;
}

int[] convertIntegersArray(ArrayList<Integer> integers)
{
    int[] ret = new int[integers.size()];
    Iterator<Integer> iterator = integers.iterator();
    for (int i = 0; i < ret.length; i++)
    {
        ret[i] = iterator.next().intValue();
    }
    return ret;
}

void debug_display(int[] ExistenceList, int[] ExistenceState){
    println("ExistenceList:" + Arrays.toString(ExistenceList));
    println("ExistenceState:" + Arrays.toString(ExistenceState));
}