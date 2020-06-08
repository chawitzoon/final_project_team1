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

void drawHole(float r) {
  r *= 0.7;
  pushMatrix();
    fill(0);
    noStroke();
    ellipse(0,0,2*r,2*r);
  popMatrix();
  pushMatrix();
    float theta = 0;
    float s = r/12;
    int n = 12*4;
    for (int i = 0; i < n; i++) {
      theta += 2*PI/n;
      translate(r*cos(theta), r*sin(theta), 0);
      fill(100);
      noStroke();
      sphere(s);//stone size
      translate(-r*cos(theta), -r*sin(theta), 0);
    }
  popMatrix();
}

void drawWisker(float s,float l) {
  fill(250);
  box(s*l, s*0.01, s*0.01);
}

void drawMole(float s, int state) {
    pushMatrix();
      translate(0, 0, -s/2);
      noFill();
      strokeWeight(2);
      stroke(0,255,255);
      box(s);
    popMatrix();

    // body
    pushMatrix();
      rotateX(-PI/2);
      noStroke();
      fill(#8b4513);
      drawCylinder(s*0.5,s*0.5,s,100);
    popMatrix();

    // head
    pushMatrix();
      translate(0, 0, -s);
      noStroke();
      fill(#8b4513);
      sphere(s*0.5);
    popMatrix();

    // eyebrows
    pushMatrix();
      translate(s*1,0,0);
      fill(255);
      beginShape();
      curveVertex(84,  91);
      curveVertex(84,  91);
      curveVertex(68,  19);
      curveVertex(21,  17);
      curveVertex(32, 100);
      curveVertex(32, 100);
      endShape();
    popMatrix();

    // draw eyes red green blue
    if (state == 0) {
      pushMatrix();
        translate(s*0.28, s*0.15, -s*1.05);
        fill(30);
        sphere(s * 0.2);
        translate(s*0.091, s*0.01, -s*0.035);
        fill(255);
        sphere(s * 0.105);
      popMatrix();
      pushMatrix();
        translate(s*0.28, -s*0.15, -s*1.05);
        fill(30);
        sphere(s * 0.2);
        translate(s*0.06, -s*0.07, -s*0.035);
        fill(255);
        sphere(s * 0.105);
      popMatrix();
    }else{
      pushMatrix();
        translate(s*0.43, s*0.2, -s*1.09);
        rotateZ(PI/2);
        rotateZ(PI/6);
        rotateY(PI/6);
        fill(30);
        box(s*0.2, s*0.05, s*0.05);
      popMatrix();
      pushMatrix();
        translate(s*0.44, s*0.2, -s*1.01);
        rotateZ(PI/2);
        rotateZ(PI/8);
        rotateY(-PI/6);
        fill(30);
        box(s*0.2, s*0.05, s*0.05);
      popMatrix();
      pushMatrix();
        translate(s*0.43, -s*0.2, -s*1.09);
        rotateZ(PI/2);
        rotateZ(-PI/6);
        rotateY(-PI/6);
        fill(30);
        box(s*0.2, s*0.05, s*0.05);
      popMatrix();
      pushMatrix();
        translate(s*0.44, -s*0.2, -s*1.01);
        rotateZ(PI/2);
        rotateZ(-PI/8);
        rotateY(PI/6);
        fill(30);
        box(s*0.2, s*0.05, s*0.05);
      popMatrix();
    }
    

    // draw nose
    pushMatrix();
      translate(s*0.28,0,-s*0.85);
      rotateZ(radians(-90));
      fill(255, 120, 0);
      drawCylinder(s * 0.1, s*0.03, s * 0.25, 32);
    popMatrix();

    // draw whisker
    pushMatrix();
      translate(s*0.4, s*0.25, -s*0.85);
      rotateZ(PI/2);
      drawWisker(s,0.5);
    popMatrix();
    pushMatrix();
      translate(s*0.4, s*0.25, -s*0.86);
      rotateZ(PI/2);
      rotateY(PI/12);
      drawWisker(s,0.5);
    popMatrix();
    pushMatrix();
      translate(s*0.4, s*0.25, -s*0.84);
      rotateZ(PI/2);
      rotateY(-PI/12);
      drawWisker(s,0.5);
    popMatrix();

    pushMatrix();
      translate(s*0.4, -s*0.25, -s*0.85);
      rotateZ(PI/2);
      drawWisker(s,0.5);
    popMatrix();
    pushMatrix();
      translate(s*0.4, -s*0.25, -s*0.84);
      rotateZ(PI/2);
      rotateY(PI/12);
      drawWisker(s,0.5);
    popMatrix();
    pushMatrix();
      translate(s*0.4, -s*0.25, -s*0.86);
      rotateZ(PI/2);
      rotateY(-PI/12);
      drawWisker(s,0.5);
    popMatrix();

}