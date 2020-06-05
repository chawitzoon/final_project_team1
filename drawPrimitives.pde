void drawCylinder(float topRadius, float bottomRadius, float tall, int sides) {
  float angle = 0;
  float angleIncrement = TWO_PI / sides;

  beginShape(QUAD_STRIP);
  for (int i = 0; i < sides + 1; ++i) {
    vertex(topRadius*cos(angle), 0, topRadius*sin(angle));
    vertex(bottomRadius*cos(angle), tall, bottomRadius*sin(angle));
    angle += angleIncrement;
  }
  endShape();

  // If it is not a cone, draw the circular top cap
  if (topRadius != 0) {
    angle = 0;
    beginShape(TRIANGLE_FAN);

    // Center point
    vertex(0, 0, 0);
    for (int i = 0; i < sides + 1; i++) {
      vertex(topRadius * cos(angle), 0, topRadius * sin(angle));
      angle += angleIncrement;
    }
    endShape();
  }

  // If it is not a cone, draw the circular bottom cap
  if (bottomRadius != 0) {
    angle = 0;
    beginShape(TRIANGLE_FAN);

    // Center point
    vertex(0, tall, 0);
    for (int i = 0; i < sides + 1; i++) {
      vertex(bottomRadius * cos(angle), tall, bottomRadius * sin(angle));
      angle += angleIncrement;
    }
    endShape();
  }
}

void drawSnowman(float s) {
    pushMatrix();
      translate(0, 0, -s/2);
      noFill();
      strokeWeight(2);
      stroke(0,255,255);
      box(s);
    popMatrix();

    // draw snowman
    pushMatrix();
      translate(0, 0, -s * 0.4);
      noStroke();
      fill(255);
      sphere(s * 0.4);
      translate(0, 0, -s * 0.7);
      sphere(s * 0.3);
      translate(0, 0, -s * 0.5);
      sphere(s * 0.2);
    popMatrix();

    // draw eyes
    pushMatrix();
      translate(s * 0.1, s * 0.08, -s * 1.7);
      fill(0);
      sphere(s * 0.05);

      translate(0,  -s * 0.16, 0);
      fill(0);
      sphere(s * 0.05);
    popMatrix();

    // draw nose
    pushMatrix();
      translate(s * 0.2,  0, -s * 1.6);
      rotateZ(radians(-90));
      fill(255, 120, 0);
      drawCylinder(s * 0.05, 0, s * 0.15, 32);
    popMatrix();
}