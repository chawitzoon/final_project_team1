
float rotateToMarker(PMatrix3D thisMarker, PMatrix3D lookAtMarker, int markernumber) {
  PVector relativeVector = new PVector();
  relativeVector.x = lookAtMarker.m03 - thisMarker.m03;
  relativeVector.y = lookAtMarker.m13 - thisMarker.m13;
  relativeVector.z = lookAtMarker.m23 - thisMarker.m23;
  float relativeLen = relativeVector.mag();

  relativeVector.normalize();

  float[] defaultLook = {1, 0, 0, 0};
  snowmanLookVector = new PVector();
  snowmanLookVector.x = thisMarker.m00 * defaultLook[0];
  snowmanLookVector.y = thisMarker.m10 * defaultLook[0];
  snowmanLookVector.z = thisMarker.m20 * defaultLook[0];

  snowmanLookVector.normalize();

  float angle = PVector.angleBetween(relativeVector, snowmanLookVector);
  if (relativeVector.x * snowmanLookVector.y - relativeVector.y * snowmanLookVector.x < 0)
    angle *= -1;

  return angle;
}