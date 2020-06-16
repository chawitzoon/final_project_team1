//all method used for first calibration
import java.util.Arrays;
import java.util.Iterator;

void calibration(ArrayList<Marker> markers){
  markers.clear();
  ortho();
  pushMatrix();
    translate(-width/2, -height/2,-(height/2)/tan(radians(fov)));
    markerTracker.findMarker(markers);
  popMatrix();
  println("marker.size in loop:" + markers.size());

}

void debugDisplay(int[] ExistenceList, int[] ExistenceState){
    println("ExistenceList:" + Arrays.toString(ExistenceList));
    println("ExistenceState:" + Arrays.toString(ExistenceState));
}