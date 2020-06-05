import org.opencv.core.Point;
import org.opencv.core.Size;
import org.opencv.core.Mat;

import org.opencv.core.CvType;

PVector pointToPVector(Point p) {
    return OpenCV.pointToPVector(p);
}

PVector[] pointToPVectorArray(Point[] p) {
    PVector[] v = new PVector[p.length];
    for (int i = 0; i < p.length; i++)
        v[i] = OpenCV.pointToPVector(p[i]);
    return v;
}

Point pVectorToPoint(PVector v) {
    return new Point(v.x, v.y);
}

Point[] pVectorToPointArray(PVector[] v) {
    Point[] p = new Point[v.length];
    for (int i = 0; i < v.length; i++)
        p[i] = pVectorToPoint(v[i]);
    return p;
}

Mat pMatrixToMat(PMatrix pm) {
    Mat mat;
    float[] fm_all = new float[16];
    float[] fm;

    pm.get(fm_all);

    if (pm instanceof PMatrix2D) {
        mat = new Mat(new Size(3, 3), CvType.CV_32FC1);
        fm = new float[] {
                    fm_all[0], fm_all[1], fm_all[2],
                    fm_all[4], fm_all[5], fm_all[6],
                    fm_all[8], fm_all[9], fm_all[10], 
                };
    } else {
        mat = new Mat(new Size(4, 4), CvType.CV_32FC1);
        fm = fm_all;
    }

    mat.put(0, 0, fm);
    return mat;
}

PMatrix2D fArrayToPMatrix2D(float[][] fm) {
    return fArrayToPMatrix2D(fm[0], fm[1]);
}

PMatrix2D fArrayToPMatrix2D(float[] a, float[] b) {
    return new PMatrix2D(
        a[0], a[1], 0, 
        b[0], b[1], 0
    );
}

PMatrix3D fArrayToPMatrix3D(float[][] fm) {
    return fArrayToPMatrix3D(fm[0], fm[1], fm[2]);
}

PMatrix3D fArrayToPMatrix3D(float[] a, float[] b, float[] c) {
    return new PMatrix3D(
        a[0], a[1], a[2], 0, 
        b[0], b[1], b[2], 0, 
        c[0], c[1], c[2], 0,
        0, 0, 0, 1);
}

PMatrix2D matToPMatrix2D(Mat m) {
    float[] a = new float[2];
    float[] b = new float[2];
    m.get(0, 0, a);
    m.get(1, 0, b);

    return fArrayToPMatrix2D(a, b);
}

PMatrix3D matToPMatrix3D(Mat m) {
    float[] a = new float[3];
    float[] b = new float[3];
    float[] c = new float[3];
    m.get(0, 0, a);
    m.get(1, 0, b);
    m.get(2, 0, c);

    return fArrayToPMatrix3D(a, b, c);
}
