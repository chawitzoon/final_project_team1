import gab.opencv.*;
import org.opencv.imgproc.Imgproc;

import org.opencv.core.Core;

import org.opencv.core.Mat;
import org.opencv.core.MatOfPoint;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.CvType;

import org.opencv.core.Point;
import org.opencv.core.Size;
import org.opencv.core.Rect;

final int QX = 0;
final int QY = 1;
final int QZ = 2;
final int QW = 3;

final boolean PRINT_OPTIMIZATION = false;

// /**
//  * implementation based on proposal by Horn, idea:
//  * first determine the largest entry of unit quaternion,
//  * then use this to get other entries
//  * adapted from dwarfutil.cpp
//  */
float[] matrixToQuaternion(float[][] m, float[] q) {
	// get entry of q with largest absolute value
	// note: we compute here 4 * q[..]^2 - 1
	float[] tmp = new float[4];
	tmp[QW] =  m[0][0] + m[1][1] + m[2][2];
	tmp[QX] =  m[0][0] - m[1][1] - m[2][2];
	tmp[QY] = -m[0][0] + m[1][1] - m[2][2];
	tmp[QZ] = -m[0][0] - m[1][1] + m[2][2];

	int max = QW;
	if (tmp[QX] > tmp[max]) max = QX;
	if (tmp[QY] > tmp[max]) max = QY;
	if (tmp[QZ] > tmp[max]) max = QZ;

	// depending on largest entry compute the other values
	// note: these formulae can be derived very simply from the
	//       matrix representation computed in quaternionToMatrix
	switch(max) {
	case QW:
		q[QW] = sqrt(tmp[QW]+1) * 0.5;
		q[QX] = (m[2][1] - m[1][2]) / (4 * q[QW]);
		q[QY] = (m[0][2] - m[2][0]) / (4 * q[QW]);
		q[QZ] = (m[1][0] - m[0][1]) / (4 * q[QW]);
		break;

	case QX:
		q[QX] = sqrt(tmp[QX]+1) * 0.5f;
		q[QW] = (m[2][1] - m[1][2]) / (4 * q[QX]);
		q[QY] = (m[1][0] + m[0][1]) / (4 * q[QX]);
		q[QZ] = (m[0][2] + m[2][0]) / (4 * q[QX]);
		break;

	case QY:
		q[QY] = sqrt(tmp[QY]+1) * 0.5;
		q[QW] = (m[0][2] - m[2][0]) / (4 * q[QY]);
		q[QX] = (m[1][0] + m[0][1]) / (4 * q[QY]);
		q[QZ] = (m[2][1] + m[1][2]) / (4 * q[QY]);
		break;

	case QZ:
		q[QZ] = sqrt(tmp[QZ]+1) * 0.5;
		q[QW] = (m[1][0] - m[0][1]) / (4 * q[QZ]);
		q[QX] = (m[0][2] + m[2][0]) / (4 * q[QZ]);
		q[QY] = (m[2][1] + m[1][2]) / (4 * q[QZ]);
		break;
	}

	normalizeQuaternion(q);

	return q;
}

/**
 * computes where a point should be on the screen, given a pose
 * @param p2D output 2d vector
 * @param p3d input 3D vector
 * @param pRotation rotation as quaternione
 * @param pTranslation 3-element translation
 * @param f focal length
 */
 PVector projectPoint(PVector p3D, float[] pRot, PVector pTrans, float f) {
	PVector point = new PVector();
	PVector point3D = new PVector();
	point3D.set(p3D);

	// rotate
	point = rotateQuaternion(pRot, point3D);

	// translate
	point = PVector.add(point, pTrans);

	// project
	// TODO: check for division by zero
    point.x = f * point.x / -point.z;
	point.y = f * point.y / -point.z;
	point.z = 0;

	return point;
}

/**
 * rotate a vector around a quaternion
 * implementation based on precomputation, correct mult
 * matrix can be found in Horn, Closed-form solution of
 * absolute orientation using unit quaternions. (1987)
 * adapted from dwarfutil.cpp
 */
PVector rotateQuaternion(float[] q, PVector p) {
	// precomputation of some values
	float xy = q[QX]*q[QY];
	float xz = q[QX]*q[QZ];
	float yz = q[QY]*q[QZ];
	float ww = q[QW]*q[QW];
	float wx = q[QW]*q[QX];
	float wy = q[QW]*q[QY];
	float wz = q[QW]*q[QZ];

    PVector r = new PVector();
	r.x = p.x * ( 2.0*(q[QX]*q[QX] + ww) - 1.0 ) + p.y * 2.0 * (xy - wz) + p.z * 2.0 * (wy + xz);
	r.y = p.x * 2.0 * (xy + wz) + p.y * ( 2.0*(q[QY]*q[QY] + ww) - 1.0 ) + p.z * 2.0 * (yz - wx);
	r.z = p.x * 2.0 * (xz - wy) + p.y * 2.0 * (wx + yz) + p.z * ( 2.0*(q[QZ]*q[QZ] + ww) - 1.0 );

	return r;
}

/**
 * factors a non-unit-length of a quaternion into the translation
 */
void normalizePose(float[] pRot, PVector pTrans) {
	// compute length of quaternion
	float fQuatLenSq = 0.0f;

	for (int i = 0; i < 4; i++)
		fQuatLenSq += pRot[i] * pRot[i];
	float fQuatLen = sqrt(fQuatLenSq);

	// normalize quaternion
	for (int i = 0; i < 4; i++)
		pRot[i] /= fQuatLen;

	// scale translation
    PVector.div(pTrans, fQuatLenSq);
}

/** 
 * normalizes a quaternion (makes it a unit quaternion)
 */
float[] normalizeQuaternion(float[] q) {
	float norm = 0.0;

    for (int i = 0; i < 4; i++)
		norm += q[i] * q[i];

	norm = sqrt(1.0 / norm);

	for (int i = 0; i < 4; i++)
		q[i] *= norm;

	return q;
}

/**
 * computes the orientation and translation of a square using homography
 * @param pRot : result as quaternion
 * @param pTrans : result position
 * @param p2D : four input coordinates. the origin is assumed to be at the camera's center of projection
 * @param fMarkerSize : side-length of marker. Origin is at marker center.
 * @param f : focal length
 */
void getInitialPose(float[] pRot, PVector pTrans, PVector[] p2D, float fMarkerSize, float f) {
	// compute homography
	float[][] hom = new float[3][3];
	calcHomography(hom, p2D);

	// compute rotation matrix by multiplying with inverse of camera matrix 
    // and inverse marker scaling:
	// R = C^-1 H S^-1
	float[][] fRotMat = new float[3][3];
	float[] fScaleLeft  = { 1.0/f, 1.0/f, -1.0 };
	float[] fScaleRight = { 1.0/fMarkerSize, 1.0/fMarkerSize, 1.0 };

	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
			fRotMat[r][c] = hom[r][c] * fScaleLeft[r] * fScaleRight[c];

	// check sign of z-axis translation, multiply matrix with -1 if necessary
	if (fRotMat[2][2] > 0.0f)
    	for (int r = 0; r < 3; r++)
	    	for (int c = 0; c < 3; c++)
		    	fRotMat[r][c] *= -1;

    PVector[] cols = new PVector[3];
    for (int i = 0; i < 3; i++)
        cols[i] = new PVector(fRotMat[0][i], fRotMat[1][i], fRotMat[2][i]);

	// compute length of the first two colums
	float fXLen = cols[0].mag();
	float fYLen = cols[1].mag();

	// copy & normalize translation
	float fTransScale = 2.0 / (fXLen + fYLen);
	pTrans.x = fRotMat[0][2] * fTransScale;
	pTrans.y = fRotMat[1][2] * fTransScale;
	pTrans.z = fRotMat[2][2] * fTransScale;

	// normalize first two colums
	cols[0].normalize();
	cols[1].normalize();

	// compute  third row as vector product
	cols[2] = (cols[0].cross(cols[1])).normalize();

	// recompute y vector from x and z
	cols[1] = (cols[0].cross(cols[2])).mult(-1);

	for (int i = 0; i < 3; i++) {
		fRotMat[0][i] = cols[i].x;
		fRotMat[1][i] = cols[i].y;
		fRotMat[2][i] = cols[i].z;
	}

	// compute rotation quaternion from matrix
	matrixToQuaternion(fRotMat, pRot);
}

/**
 * computes the Jacobian for optimizing the pose
 * @param pResult 2x7 matrix
 * @param pParam rotation parameters: 4 * quaternion rotation + 3 * translation
 * @param p3D the 3d input vector
 * @param f focal length
 */
void computeJacobian(float[] pResult, float[] pRot, PVector pTrans, PVector p3D, float f ) {
	// TODO: check for division by zero
	// maple-generated code
	float t4  = pRot[0]*p3D.x+pRot[1]*p3D.y+pRot[2]*p3D.z;
	float t10 = pRot[3]*p3D.x+pRot[1]*p3D.z-pRot[2]*p3D.y;
	float t15 = pRot[3]*p3D.y-pRot[0]*p3D.z+pRot[2]*p3D.x;
	float t20 = pRot[3]*p3D.z+pRot[0]*p3D.y-pRot[1]*p3D.x;
	float t22 = -t4*pRot[2]+t10*pRot[1]-t15*pRot[0]-t20*pRot[3]-pTrans.z;
	float t23 = 1.0/t22;
	float t24 = 2.0*f*t4*t23;
	float t30 = f*(t4*pRot[0]+t10*pRot[3]-t15*pRot[2]+t20*pRot[1]+pTrans.x);
	float t31 = t22*t22;
	float t32 = 1.0/t31;
	float t33 = -2.0*t32*t15;
	float t38 = 2.0*t32*t10;
	float t43 = -2.0*t32*t4;
	float t47 = 2.0*f*t10*t23;
	float t48 = -2.0*t32*t20;
	float t51 = f*t23;
	float t60 = f*(t4*pRot[1]+t10*pRot[2]+t15*pRot[3]-t20*pRot[0]+pTrans.y);

    pResult[0] = t24-t30*t33;
	pResult[1] = 2.0*f*t20*t23-t30*t38;
	pResult[2] = -2.0*f*t15*t23-t30*t43;
	pResult[3] = t47-t30*t48;
	pResult[4] = t51;
	pResult[5] = 0.0;
	pResult[6] = t30*t32;
	pResult[7+0] = -2.0*f*t20*t23-t60*t33;
	pResult[7+1] = t24-t60*t38;
	pResult[7+2] = t47-t60*t43;
	pResult[7+3] = 2.0*f*t15*t23-t60*t48;
	pResult[7+4] = 0.0;
	pResult[7+5] = t51;
	pResult[7+6] = t60*t32;
}

/**
 * computes the reprojection error for each point
 * @param pError output: two entries (x,y) for each input point = (measured - reprojected)
 * @param p3D 3D coordinates of the points
 * @param p2D measured 2D coordinates
 * @param nPoints number of input points
 * @param pRot rotation
 * @param pTrans translation
 * @param f focal length
 * @returns absolute squared error
 */
float computeReprojectionError(PVector[] pError, PVector[] p3D, PVector[] p2D, int nPoints,
                               float[] pRot, PVector pTrans, float f) {
	float fAbsErrSq = 0.0;

	for (int i = 0; i < nPoints; i++) {
		// println("pTrans:", i, pTrans);
		// reproject
		PVector projected = projectPoint(p3D[i], pRot, pTrans, f);
		// println("proj:", i, projected);

		// compute deviation
        pError[i].x = p2D[i].x - projected.x;
        pError[i].y = p2D[i].y - projected.y;

		// println("err:", pError[i]);
		// update absolute error
		fAbsErrSq += pError[i].x * pError[i].x + pError[i].y * pError[i].y;
	}

    return fAbsErrSq;
}

/**
 * optimize a pose with levenberg-marquardt
 * @param pRotation rotation as quaternion, both used as output and initial value
 * @param pTranslation 3-element translation, both used as output and initial value
 * @param nPoints number of correspondences
 * @param p2D pointer to camera coordinates
 * @param p3D pointer to object coordinates
 * @param f focal length
 */
void optimizePose(float[] pRot_, PVector pTrans_, int nPoints, 
                  PVector[] p2D, PVector[] p3D, float f) {
	float[] pRot = new float[4];
	for (int i = 0; i < 4; i++)
		pRot[i] = pRot_[i];

	PVector pTrans = new PVector();
	pTrans.set(pTrans_);

	// compute initial error
	PVector[] measurementDiffPrev = new PVector[nPoints];
    PVector[] measurementDiffNew  = new PVector[nPoints];

	for (int i = 0; i < nPoints; i++) {
		measurementDiffPrev[i] = new PVector();
		measurementDiffNew[i]  = new PVector();
	}

    float[] jacobian = new float[nPoints * 2 * 7];
    float[] jacobiSquare = new float[7 * 7];
    float fLambda = 1.0;

    Mat matJacobian = new Mat(nPoints * 2, 7, CvType.CV_32F);
    Mat matJacobiSquare = new Mat(7, 7, CvType.CV_32F);
    Mat matMeasurementDiffPrev = new Mat(nPoints * 2, 1, CvType.CV_32F);
    Mat matParamDiff = new Mat(7, 1, CvType.CV_32F);
    Mat matMDiff2 = new Mat(7, 1, CvType.CV_32F);

	float fPreviousErr = computeReprojectionError(measurementDiffPrev, p3D, p2D, nPoints, pRot, pTrans, f);
	float[] fMeasurementDiffPrev = new float[2 * nPoints];

	for (int i = 0; i < nPoints; i++) {
		fMeasurementDiffPrev[2 * i]     = measurementDiffPrev[i].x;
		fMeasurementDiffPrev[2 * i + 1] = measurementDiffPrev[i].y;
	}
	matMeasurementDiffPrev.put(0, 0, fMeasurementDiffPrev);
	// println("matMeasurementDiffPrev---");
	// println(matMeasurementDiffPrev.dump());

    if (PRINT_OPTIMIZATION)
    	println("initial error: ", fPreviousErr);

	// iterate (levenberg-marquardt)
	int nMaxIterations = 3;
	for (int iIteration = 0; iIteration < nMaxIterations; iIteration++) {
		// println("--- itr:", iIteration, "-------------------------------------");

		// create jacobian
		for (int i = 0; i < nPoints; i++) {
            float[] tmpJ = new float[2 * 7];
			computeJacobian(tmpJ, pRot, pTrans, p3D[i], f);

            for (int j = 0; j < 14; j++) {
                jacobian[i * 2 * 7 + j] = tmpJ[j];
			}

			// println("jacobian",i);
			// println(tmpJ[0],tmpJ[1],tmpJ[2],tmpJ[3],tmpJ[4],tmpJ[5],tmpJ[6]);
        }
        matJacobian.put(0, 0, jacobian);

        // multiply both sides with J^T
        Core.gemm(matJacobian.t(), matJacobian, 1, 
                  Mat.zeros(7, 7, CvType.CV_32F), 0, matJacobiSquare);
        Core.gemm(matJacobian.t(), matMeasurementDiffPrev, 1, 
                  Mat.zeros(7, 1, CvType.CV_32F), 0, matMDiff2);
        // Core.gemm(matJacobiSquare, Mat.eye(7, 7, CvType.CV_32F), 1,
        //           Mat.eye(7, 7, CvType.CV_32F), fLambda, matJacobiSquare);

		for (int i = 0; i < 7; i++) {
			double tmp = matJacobiSquare.get(i, i)[0];
			matJacobiSquare.put(i, i, tmp + fLambda);
		}

		// println("MatJacobiSquare ---" );
		// println(matJacobiSquare.dump());

		// println("matMDiff2 ---" );
		// println(matMDiff2.dump());

        // do least squares
        // performance improvement: use cholesky decomp
        Core.solve(matJacobiSquare, matMDiff2, matParamDiff, Core.DECOMP_CHOLESKY);

		// println("matParamDiff ---");
		// println(matParamDiff.dump());

        float[] fParamDiff = new float[7];
        matParamDiff.get(0, 0, fParamDiff);

        // update parameters
		float[] pRotNew = new float[4];
 		for (int i = 0; i < 4; i++)
			pRotNew[i] = pRot[i] + fParamDiff[i];

        PVector pTransNew = PVector.add(pTrans, new PVector(fParamDiff[4], fParamDiff[5], fParamDiff[6]));

        // factor the quaternion length into the translation
        normalizePose(pRotNew, pTransNew);

		// println("pRotNew 2:", pRotNew[0], pRotNew[1], pRotNew[2], pRotNew[3]);
		// println("pTransNew 2:", pTransNew);

    	// compute new error
		float fErr = computeReprojectionError(measurementDiffNew, p3D, p2D, nPoints, pRotNew, pTransNew, f);

		// for (int i = 0; i < 4; i++)
		// 	println("mdiffNew[",i,"]:",measurementDiffNew[i]);

		if (fErr >= fPreviousErr)
			fLambda *= 10.0;
		else {
			fLambda /= 10.0;

			// update parameters
			for (int i = 0; i < 4; i++)
				pRot[i] = pRotNew[i];
            pTrans.set(pTransNew);

			// copy measurement error
			for ( int i = 0; i < nPoints; i++)
				measurementDiffPrev[i].set(measurementDiffNew[i]);

			fPreviousErr = fErr;
		}

        if (PRINT_OPTIMIZATION)
            println(" it", iIteration, ": fErr=", fErr, "lambda=", fLambda);
    }
}

/**
 * @param mat : result as 4x4 matrix in row-major format
 * @param p2D : coordinates of the four corners in counter-clock-wise order. 
 *              the origin is assumed to be at the camera's center of projection
 * @param markerSize : side-length of marker. Origin is at marker center.
 */
PMatrix3D estimateSquarePose(PVector[] p2D_, float markerSize) {
	PVector[] p2D = new PVector[4];
	for (int i = 0; i < 4; i++) {
		p2D[i] = new PVector();
		p2D[i].set(p2D_[i]);
	}

	// approximate focal length for logitech quickcam 4000 at 320*240 resolution
	// float fFocalLength = 900.0f;

	// Adjust fFocalLength according to your camera setup
	float fFocalLength = (height/2)/tan(radians(fov)) * 2.5;

	// compute initial pose
	float[] rot = new float[4];
    PVector trans = new PVector();
	getInitialPose(rot, trans, p2D, markerSize, fFocalLength);
    // println(rot[0], rot[1], rot[2], rot[3]);

	// corner 3D coordinates
	float fCp = markerSize / 2;
	PVector[] p3D = new PVector[4];

    // counter-clock-wise
    p3D[0] = new PVector(-fCp,  fCp, 0);
    p3D[1] = new PVector(-fCp, -fCp, 0);
    p3D[2] = new PVector( fCp, -fCp, 0);
    p3D[3] = new PVector( fCp,  fCp, 0);

	PVector[] points = new PVector[4];
	for (int i = 0; i < 4; i++) {
		points[i] = new PVector();
		points[i].set(p2D[i]);
	}

	// refine pose using nonlinear optimization
	optimizePose(rot, trans, 4, points, p3D, fFocalLength);

    // convert quaternion to matrix
	// convert quaternion to matrix
	float X = -rot[0];
	float Y = -rot[1];
	float Z = -rot[2];
	float W =  rot[3];

	float xx = X * X;
	float xy = X * Y;
	float xz = X * Z;
	float xw = X * W;
	float yy = Y * Y;
	float yz = Y * Z;
	float yw = Y * W;
	float zz = Z * Z;
	float zw = Z * W;

	PMatrix3D result = new PMatrix3D();
	result.m00 = 1 - 2 * ( yy + zz );
	result.m01 =     2 * ( xy + zw );
	result.m02 =     2 * ( xz - yw );
	result.m10 =     2 * ( xy - zw );
	result.m11 = 1 - 2 * ( xx + zz );
	result.m12 =     2 * ( yz + xw );
	result.m20 =     2 * ( xz + yw );
	result.m21 =     2 * ( yz - xw );
	result.m22 = 1 - 2 * ( xx + yy );

    result.m03 = trans.x;
	result.m13 = trans.y;
	result.m23 = trans.z;
	result.m30 = result.m31 = result.m32 = 0;
	result.m33 = 1;

    // right -> left hand conversion
    result.m13 = -result.m13;
    result.m01 = -result.m01;
    result.m10 = -result.m10;
    result.m12 = -result.m12;
    result.m21 = -result.m21;

    if (PRINT_OPTIMIZATION) {
    	println("rot:", rot[0], rot[1], rot[2], rot[3]);
    	println("tra: ", trans);
    }

	return result;
}

// Returns Matrix in Row-major format
void calcHomography(float[][] pResult, PVector[] pQuad) {
	// homography computation from Ela Harker & O'Leary, simplified for squares

	// subtract mean from points
	PVector[] c = new PVector[4];
    for (int i = 0; i < 4; i++)
        c[i] = new PVector();

	PVector mean = new PVector();
	for (int i = 0; i < 4; i++)
        mean = PVector.add(mean, pQuad[i]);
    mean = PVector.div(mean, 4);

	for (int i = 0; i < 4; i++)
        c[i] = PVector.sub(pQuad[i], mean);

	// build simplified matrix A
	float[] fMatA = new float[12];
    fMatA[0]  =  c[0].x - c[1].x - c[2].x + c[3].x;
	fMatA[1]  = -c[0].x - c[1].x + c[2].x + c[3].x;
	fMatA[2]  = -2 * (c[0].x + c[2].x);
	fMatA[3]  = -fMatA[0];
    fMatA[4]  = -fMatA[1];
	fMatA[5]  = -2 * (c[1].x + c[3].x);

	fMatA[6]  =  c[0].y - c[1].y - c[2].y + c[3].y;
	fMatA[7]  = -c[0].y - c[1].y + c[2].y + c[3].y;
	fMatA[8]  = -2 * (c[0].y + c[2].y);
	fMatA[9]  = -fMatA[6];
	fMatA[10] = -fMatA[7];
	fMatA[11] = -2 * (c[1].y + c[3].y);

	Mat matA = new Mat(4, 3, CvType.CV_32F);
	matA.put(0, 0, fMatA);

	// compute SVD
	// Idea: replace the whole thing with an analytical solution of the line at infinity
	Mat matW = new Mat(3, 1, CvType.CV_32F);
	Mat matU = new Mat(3, 3, CvType.CV_32F);
	Mat matV = new Mat(3, 3, CvType.CV_32F);
	Core.SVDecomp(matA, matW, matU, matV, Core.SVD_MODIFY_A);

	// copy bottom line of homography
	matV.get(2, 0, pResult[2]);

	// compute entries 1,1 and 1,2, multiply by 2 to compensate scaling
	pResult[0][0] = (( c[0].x + c[1].x + c[2].x + c[3].x) * pResult[2][0] +
		             (-c[0].x + c[1].x - c[2].x + c[3].x) * pResult[2][1] +
                     (-c[0].x - c[1].x + c[2].x + c[3].x) * pResult[2][2]) / 2;
	pResult[0][1] = ((-c[0].x + c[1].x - c[2].x + c[3].x) * pResult[2][0] +
                     ( c[0].x + c[1].x + c[2].x + c[3].x) * pResult[2][1] +
                     ( c[0].x - c[1].x - c[2].x + c[3].x) * pResult[2][2] ) / 2;

	// compute entries 2,1 and 2,2, multiply by 2 to compensate scaling
	pResult[1][0] = (( c[0].y + c[1].y + c[2].y + c[3].y) * pResult[2][0] +
                     (-c[0].y + c[1].y - c[2].y + c[3].y) * pResult[2][1] +
                     (-c[0].y - c[1].y + c[2].y + c[3].y) * pResult[2][2]) / 2;
	pResult[1][1] = ((-c[0].y + c[1].y - c[2].y + c[3].y) * pResult[2][0] +
                     ( c[0].y + c[1].y + c[2].y + c[3].y) * pResult[2][1] +
                     ( c[0].y - c[1].y - c[2].y + c[3].y) * pResult[2][2]) / 2;

	// compute entries 1,3 and 2,3
	pResult[0][2] = (( c[0].x + c[1].x - c[2].x - c[3].x) * pResult[2][0] +
                     (-c[0].x + c[1].x + c[2].x - c[3].x) * pResult[2][1]) / -4;
	pResult[1][2] = (( c[0].y + c[1].y - c[2].y - c[3].y) * pResult[2][0] +
                     (-c[0].y + c[1].y + c[2].y - c[3].y) * pResult[2][1]) / -4;

	// now multiply last row with factor 2 to compensate scaling
	pResult[2][0] *= 2;
	pResult[2][1] *= 2;

	// multiply with shift to compensate mean subtraction
	for (int i = 0; i < 3; i++) {
		pResult[0][i] = pResult[0][i] + pResult[2][i] * mean.x;
		pResult[1][i] = pResult[1][i] + pResult[2][i] * mean.y;
	}
}

