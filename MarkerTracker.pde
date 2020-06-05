import gab.opencv.*;
import org.opencv.imgproc.Imgproc;

import org.opencv.core.Core;

import org.opencv.core.Mat;
import org.opencv.core.MatOfPoint;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.MatOfFloat;

import org.opencv.core.CvType;

import org.opencv.core.Point;
import org.opencv.core.Size;
import org.opencv.core.Rect;


class Marker {
    int code;
    PMatrix3D pose;
    // float[] pose;

    Marker() {
        pose = new PMatrix3D();
    }

    void print_matrix() {
        pose.print();
    }
}

class MarkerTracker {
    int thresh;     // Threshold: gray to mono
    int bw_thresh;  // Threshold for gray marker to ID image
    double kMarkerSizeLength;

    Mat image_bgr;
	Mat image_gray;
	Mat image_gray_filtered;

    MarkerTracker(double _kMarkerSizeLength) {
        thresh = 80;
        bw_thresh = 100;
        kMarkerSizeLength = _kMarkerSizeLength;
        init();
    }

    MarkerTracker(double _kMarkerSizeLength, int _thresh, int _bw_thresh) {
        thresh = _thresh;
        bw_thresh = _bw_thresh;
        kMarkerSizeLength = _kMarkerSizeLength;
        init();
    }

    void init() {
        println("Startup");
    }

    void cleanup() {
        println("Finished");
    }

    int subpixSampleSafe(Mat pSrc, PVector p) {
        int x = (int)(floor(p.x));
        int y = (int)(floor(p.y));

        if (x < 0 || x >= pSrc.cols() - 1 || y < 0 || y >= pSrc.rows() - 1)
            return 127;

        int dx = (int)(256 * (p.x - floor(p.x)));
        int dy = (int)(256 * (p.y - floor(p.y)));

        int i   = (int)(pSrc.get(y,   x  )[0]);
        int ix  = (int)(pSrc.get(y,   x+1)[0]);
        int iy  = (int)(pSrc.get(y+1, x  )[0]);
        int ixy = (int)(pSrc.get(y+1, x+1)[0]);

        int a = i  + ((dx * (ix  - i )) >> 8);
        int b = iy + ((dx * (ixy - iy)) >> 8);

        return a + ((dy * (b - a)) >> 8);
    }

    // A function to check if a contour is a good candidate
    boolean checkContourCondition(MatOfPoint2f contour_approx, Mat image_bgr, int kNumOfCorners) {
        // Get a bounding rectangle of the approximated contour
        Rect bounding_rectangle = Imgproc.boundingRect(new MatOfPoint(contour_approx.toArray()));
        double marker_size = bounding_rectangle.area();

        // Filter bad contours
        int kImageSize = image_bgr.rows() * image_bgr.cols();
        double kMarkerSizeMin = kImageSize * 0.01;
        double kMarkerSizeMax = kImageSize * 0.99;
        boolean is_contour_valid = (marker_size > kMarkerSizeMin) 
            && (marker_size < kMarkerSizeMax)
            && contour_approx.size().height == kNumOfCorners
            && Imgproc.isContourConvex(new MatOfPoint(contour_approx.toArray()));

        return is_contour_valid;
    }

	void findMarker(ArrayList<Marker> markers) {
        boolean isFirstStripe = true;
        int markerCnt = 0;
        // boolean isFirstMarker = true;

        image_bgr = OpenCV.imitate(opencv.getColor());
        opencv.getColor().copyTo(image_bgr);
        int image_height = image_bgr.rows();
        int image_width  = image_bgr.cols();

        image_gray = OpenCV.imitate(opencv.getGray());
        opencv.getGray().copyTo(image_gray);

        image_gray_filtered = OpenCV.imitate(opencv.getGray());
        Imgproc.threshold(image_gray, image_gray_filtered, thresh, 255, Imgproc.THRESH_BINARY);

        ArrayList<MatOfPoint> contours = new ArrayList<MatOfPoint>();
        Mat hierarchy = new Mat();

        Imgproc.findContours(image_gray_filtered, contours, hierarchy, Imgproc.RETR_LIST, Imgproc.CHAIN_APPROX_SIMPLE);

        // show image
        PImage dst = createImage(image_width, image_height, ARGB);
        opencv.toPImage(image_bgr, dst);

        loadPixels();
        dst.loadPixels();
        for (int i = 0; i < image_width * image_height; i++)
            pixels[i] = dst.pixels[i];
        updatePixels();

        for (MatOfPoint contour : contours) {
            MatOfPoint2f contour_approx = new MatOfPoint2f();
            double kEpsilon = 0.05 * Imgproc.arcLength(new MatOfPoint2f(contour.toArray()), true);

            Imgproc.approxPolyDP(new MatOfPoint2f(contour.toArray()), contour_approx, kEpsilon, true);

            Rect kRectangle = Imgproc.boundingRect(new MatOfPoint(contour_approx.toArray()));
            double kMarkerSize = kRectangle.area();

            int kNumOfCorners = 4;
            boolean is_valid = checkContourCondition(contour_approx, image_bgr, kNumOfCorners);

            if (!is_valid)
                continue;

            if (MARKER_TRACKER_DEBUG) {
                // Draw lines
                noFill();
                strokeWeight(4);
                stroke(random(255), random(255), random(255));

                beginShape();
                Point[] p = contour_approx.toArray();
                for (int i = 0; i < p.length; i++)
                    vertex((float)p[i].x, (float)p[i].y);
                endShape(CLOSE);
            }

            MatOfFloat[] line_parameters = new MatOfFloat[4];

            for (int i = 0; i < kNumOfCorners; i++) {
                int kNumOfEdgePoints = 7;

                if (MARKER_TRACKER_DEBUG) {
                    // Draw corner points
                    float kCircleSize = 10;
                    fill(0, 255, 0);
                    noStroke();
                    Point[] p = contour_approx.toArray();
                    circle((float)p[i].x, (float)p[i].y, kCircleSize);
                }

                Point[] approx_points = contour_approx.toArray();
                PVector pa = pointToPVector(approx_points[(i+1) % kNumOfCorners]);
                PVector pb = pointToPVector(approx_points[i]);
                PVector kEdgeDirectionVec = PVector.div(PVector.sub(pa, pb), kNumOfEdgePoints);
                float kEdgeDirectionVecNorm = kEdgeDirectionVec.mag();

                // Added in Exercise 2 - Start
                // ************************************
                int stripe_length = (int)(0.8 * kEdgeDirectionVecNorm);
                if (stripe_length < 5)
                    stripe_length = 5;
                stripe_length |= 1;

                int kStop = stripe_length >> 1;
                int kStart = -kStop;

                int kStripeWidth = 3;
                Size kStripeSize = new Size(kStripeWidth, stripe_length);

                PVector kStripeVecX = PVector.div(kEdgeDirectionVec, kEdgeDirectionVecNorm);
                PVector kStripeVecY = new PVector(-kStripeVecX.y, kStripeVecX.x);

                Mat stripe_image = new Mat(kStripeSize, CvType.CV_8UC1);
                // Added in Exercise 2 - End
                // ************************************

                // Array for edge point centers
                PVector[] edge_points = new PVector[kNumOfEdgePoints - 1];

                for (int j = 1; j < kNumOfEdgePoints; j++) {
                    PVector edge_point = PVector.add(pb, PVector.mult(kEdgeDirectionVec, j));

                    if (MARKER_TRACKER_DEBUG) {
                        // Draw line delimeters
                        fill(0, 0, 255);
                        noStroke();
                        circle(edge_point.x, edge_point.y, 2);
                    }

                    // Added in Exercise 2 - Start
                    // **********************************************

                    // Construct a stripe image
                    for (int m = -1; m <= 1; m++) {
                        for (int n = kStart; n <= kStop; n++) {
                            PVector subpixel = PVector.add(
                                PVector.add(edge_point, PVector.mult(kStripeVecX, m)), 
                                PVector.mult(kStripeVecY, n)
                            );

                            if (MARKER_TRACKER_DEBUG) {
                                noStroke();
                                if (isFirstStripe) {
                                    fill(255, 0, 255);
                                    circle(subpixel.x, subpixel.y, 2);
                                } else {
                                    fill(0, 255, 255);
                                    circle(subpixel.x, subpixel.y, 2);
                                }
                            }

                            // Fetch subpixel value
                            int kSubpixelValue = subpixSampleSafe(image_gray, subpixel);
                            int kStripeX = m + 1;
                            int kStripeY = n + (stripe_length >> 1);
                            stripe_image.put(kStripeY, kStripeX, kSubpixelValue);
                        }
                    }

                    // use sobel operator on stripe
				    // ( -1 , -2, -1 )
	    			// (  0 ,  0,  0 )
		    		// (  1 ,  2,  1 )

                    double[] sobelValues = new double[stripe_length - 2];

                    for (int n = 1; n < stripe_length - 1; n++) {
                        byte[] p = new byte[3];

                        stripe_image.get(n - 1, 0, p);
                        double r1 = -(p[0] & 0xFF) - 2.0 * (p[1] & 0xFF) - (p[2] & 0xFF);

                        stripe_image.get(n + 1, 0, p);
                        double r3 =  (p[0] & 0xFF) + 2.0 * (p[1] & 0xFF) + (p[2] & 0xFF);

                        sobelValues[n - 1] = -(r1 + r3);
                    }

                    double maxVal = -1;
                    int maxIndex = 0;
                    for (int n = 0; n < stripe_length - 2; n++) {
                        if (sobelValues[n] > maxVal) {
                            maxVal = sobelValues[n];
                            maxIndex = n;
                        }
                    }

                    double y0, y1, y2;
    				y0 = (maxIndex <= 0) ? 0 : sobelValues[maxIndex - 1];
	    			y1 = sobelValues[maxIndex];
		    		y2 = (maxIndex >= stripe_length - 3) ? 0 : sobelValues[maxIndex + 1];

                    // formula for calculating the x-coordinate of the vertex of a parabola,
                    // given 3 points with equal distances
				    // (xv means the x value of the vertex, d the distance between the points):
				    // xv = x1 + (d / 2) * (y2 - y0)/(2*y1 - y0 - y2)
                    double pos = (y2 - y0) / (4 * y1 - 2 * y0 - 2 * y2);

                    // exact point with subpixel accuracy
                    int maxIndexShift = maxIndex - (stripe_length >> 1);
                    PVector edgeCenter = PVector.add(edge_point, PVector.mult(kStripeVecY, maxIndexShift + (float)pos));

                    if (MARKER_TRACKER_DEBUG) {
                        fill(0, 0, 255);
                        noStroke();
                        circle(edgeCenter.x, edgeCenter.y, 4);
                    }
                    edge_points[j - 1] = new PVector(edgeCenter.x, edgeCenter.y);

                    if (isFirstStripe) {
                        if (MARKER_TRACKER_DEBUG) {
                            // TODO: move stripe_image to another window
                            PImage dst_stripe = createImage(100, 300, ARGB);
                            Mat iplTmp = new Mat(new Size(100, 300), CvType.CV_8UC1);
                            Imgproc.resize(stripe_image, iplTmp, new Size(100, 300), 0.0, 0.0, Imgproc.INTER_NEAREST);
                            opencv.toPImage(iplTmp, dst_stripe);
                            image(dst_stripe, 0, 0);
                        }
                        isFirstStripe = false;
                    }
                } // --- end of loop over edge points of one edge

                // Added in Homework 4 (2020/5/28), Part 1 - Start
                // **************************************************************
                // Derive line parameters from subpixel-precise edge points

                MatOfPoint2f mat = new MatOfPoint2f();
                mat.fromArray(pVectorToPointArray(edge_points));

                line_parameters[i] = new MatOfFloat(new float[4]);
                Imgproc.fitLine(mat, line_parameters[i], Imgproc.CV_DIST_L2, 0, 0.01, 0.01);

                if (MARKER_TRACKER_DEBUG) {
                    int length = 50;
                    float[] tmpLine = line_parameters[i].toArray();
                    PVector p1 = new PVector(tmpLine[2] - length * tmpLine[0], tmpLine[3] - length * tmpLine[1]);
                    PVector p2 = new PVector(tmpLine[2] + length * tmpLine[0], tmpLine[3] + length * tmpLine[1]);

                    strokeWeight(2);
                    stroke(0, 255, 255);
                    line(p1.x, p1.y, p2.x, p2.y);
                }
                // Added in Homework 4 (2020/5/28), Part 1 - End
                // **************************************************************
            } // end of loop over the 4 edges

            // Added in Homework 4 (2020/5/28) Part 2 - Start
            // ****************************************************************
            // so far we stored the exact line parameters and show the lines in the image
            // now we have to calculate the exact corners

            Point[] corners = new Point[kNumOfCorners];
            for (int i = 0; i < kNumOfCorners; i++) {
                int j = (i + 1) % 4;
                float[] l1 = line_parameters[i].toArray();
                float[] l2 = line_parameters[j].toArray();

                float x0 = l1[2]; float y0 = l1[3];
                float x1 = l2[2]; float y1 = l2[3];
                float u0 = l1[0]; float v0 = l1[1];
                float u1 = l2[0]; float v1 = l2[1];

    			// (x|y) = p + s * vec
	    		// s = Ds / D (see cramer's rule)
		    	// (x|y) = p + (Ds / D) * vec
    			// (x|y) = (p * D / D) + (Ds * vec / D)
	    		// (x|y) = (p * D + Ds * vec) / D
	    		// (x|y) = a / c;
		    	float a =  x1 * u0 * v1 - y1 * u0 * u1 - x0 * u1 * v0 + y0 * u0 * u1;
		    	float b = -x0 * v0 * v1 + y0 * u0 * v1 + x1 * v0 * v1 - y1 * v0 * u1;
			    float c =  v1 * u0 - v0 * u1;

                if (abs(c) < 0.001) { // lines parallel ?
                    println("lines parallel");
                    continue;
                }

                a /= c;
                b /= c;
                corners[i] = new Point(a, b);

                if (MARKER_TRACKER_DEBUG) {
                    fill(255/(i+1), 255/(i+1), 0);
                    noStroke();
                    circle((float)corners[i].x, (float)corners[i].y, 10);
                }
            } // finished the calculation of the exact corners

            Point[] targetCorners = new Point[kNumOfCorners];
            int scale = 10;
            targetCorners[0] = new Point(-0.5,             -0.5);
            targetCorners[1] = new Point(-0.5 + 6 * scale, -0.5);
            targetCorners[2] = new Point(-0.5 + 6 * scale, -0.5 + 6 * scale);
            targetCorners[3] = new Point(-0.5,             -0.5 + 6 * scale);

            MatOfPoint2f cornersMat = new MatOfPoint2f();
            MatOfPoint2f targetCornersMat = new MatOfPoint2f();
            cornersMat.fromArray(corners);
            targetCornersMat.fromArray(targetCorners);

    		// create and calculate the matrix of perspective transform
            Mat projMat = Imgproc.getPerspectiveTransform(cornersMat, targetCornersMat);

            // change the perspective in the marker image using the previously calculated matrix
            Mat iplMarker = new Mat(new Size(6, 6), CvType.CV_8UC1);
            Mat iplMarkerScaled = new Mat(new Size(6 * scale, 6 * scale), CvType.CV_8UC1);
            Imgproc.warpPerspective(image_gray, iplMarkerScaled, projMat, new Size(6 * scale, 6 * scale));

            int bw_thresh = 100;
            Imgproc.threshold(iplMarkerScaled, iplMarkerScaled, bw_thresh, 255, Imgproc.THRESH_BINARY | Imgproc.THRESH_OTSU);

            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    Mat subMat = iplMarkerScaled.submat(i*scale, (i+1)*scale-1, j*scale, (j+1)*scale-1);
                    int blackPixNum = (int)(Core.sumElems(subMat).val[0] / 255);

                    if (blackPixNum > scale * scale / 2)
                        iplMarker.put(i, j, 255);
                    else
                        iplMarker.put(i, j, 0);
                }
            }

            // now, we have a Black/White image of a supposed Marker
            // check if border is black
            int code = 0;
            for (int i = 0; i < 6; i++) {
                int pixel1 = (int)(iplMarker.get(0, i)[0]);
                int pixel2 = (int)(iplMarker.get(5, i)[0]);
                int pixel3 = (int)(iplMarker.get(i, 0)[0]);
                int pixel4 = (int)(iplMarker.get(i, 5)[0]);
                if ((pixel1 > 0) || (pixel2 > 0) || (pixel3 > 0) || (pixel4 > 0)) {
                    code = -1;
                    break;
                }
            }

            if (code < 0)
                continue;

            // copy the BW values into cP
            int[][] cP = new int[4][4];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    cP[i][j] = (int)(iplMarker.get(i+1, j+1)[0]);
                    cP[i][j] = (cP[i][j] == 0) ? 1 : 0; // if black then 1 else 0
                }
            }

            // save the ID of the marker
            int[] codes = new int[4];
            codes[0] = codes[1] = codes[2] = codes[3] = 0;
            for (int i = 0; i < 16; i++) {
    			int row = i >> 2;
	    		int col = i % 4;

    			codes[0] <<= 1;
	    		codes[0] |= cP[row][col]; // 0 deg

    			codes[1] <<= 1;
	    		codes[1] |= cP[3 - col][row]; // 90 deg

    			codes[2] <<= 1;
	    		codes[2] |= cP[3 - row][3 - col]; // 180 deg

    			codes[3] <<= 1;
	    		codes[3] |= cP[col][3 - row]; // 270 deg
		    }

    	    if ((codes[0] == 0) || (codes[0] == 0xffff))
	    		continue;

            // account for symmetry
            code = codes[0];
            int angle = 0;
            for (int i = 1; i < 4; i++) {
                if (codes[i] < code) {
                    code = codes[i];
                    angle = i;
                }
            }

            if (MARKER_TRACKER_DEBUG) {
                Mat iplTmp = new Mat();
                int dispSize = 100;
                Imgproc.resize(iplMarker, iplTmp, new Size(dispSize, dispSize), 
                               0.0, 0.0, Imgproc.INTER_NEAREST);

                PImage dst_marker = createImage(dispSize, dispSize, ARGB);
                opencv.toPImage(iplTmp, dst_marker);

                int dispX = markerCnt % 3;
                int dispY = markerCnt / 3;
                image(dst_marker, dispX * dispSize, dispY * dispSize);

                strokeWeight(1);
                stroke(255, 0, 0);
                for (int i = 0; i < 6; i++) {
                    line((dispX + i / 6.0) * dispSize,      dispY  * dispSize,
                         (dispX + i / 6.0) * dispSize, (dispY + 1) * dispSize);
                    line( dispX      * dispSize, (dispY + i / 6.0) * dispSize,
                         (dispX + 1) * dispSize, (dispY + i / 6.0) * dispSize);
                }
                markerCnt += 1;
            }

            // Added in Homework 5 (2020/6/3) - Start
            // **************************************************************
            // correct the order of the corners
            if (angle != 0) {
                Point[] corrected_corners = new Point[4];
                for (int i = 0; i < 4; i++)
                    corrected_corners[i] = new Point();

                for (int i = 0; i < 4; i++)
                    corrected_corners[(i + angle) % 4] = corners[i].clone();
                for (int i = 0; i < 4; i++)
                    corners[i] = corrected_corners[i].clone();
            }

            // transfer screen coords to camera coords
            for (int i = 0; i < 4; i++) {
                // here you have to use your own camera resolution (x) * 0.5
                corners[i].x -= image_width / 2;
                // here you have to use your own camera resolution (y) * 0.5
                corners[i].y = -corners[i].y + image_height / 2;
            }

            Marker marker = new Marker();
            marker.code = code;

            marker.pose = estimateSquarePose(pointToPVectorArray(corners), (float)kMarkerSizeLength);
		    markers.add(marker);

            if (MARKER_TRACKER_DEBUG) {
                println("Found: " + hex(code));
                // marker.print_matrix();
            }
            // Added in Homework 5 (2020/6/3) - End
            // **************************************************************
        } // end of loop over contour candidates

        // If debugging, draw thresholds to be used for marker tracking
        if (MARKER_TRACKER_DEBUG) {
            fill(255, 0, 0);
            textSize(20);
            text("thresh : "    + thresh,    width-200, 50);
            text("bw_thresh : " + bw_thresh, width-200, 80);
        }
    }
}