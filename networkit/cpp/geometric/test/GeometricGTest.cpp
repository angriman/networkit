/*
 * GeometricGTest.cpp
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#include <cmath>
#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/geometric/HyperbolicSpace.hpp>
#include <networkit/viz/Point.hpp>

#include <algorithm>
#include <cmath>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/geometric/HyperbolicSpace.hpp>
#include <networkit/geometric/Point2DWithIndex.hpp>

namespace NetworKit {

class GeometricGTest : public testing::Test {};

/**
 * test conversion of polar coordinates to cartesian coordinates and back
 */
TEST_F(GeometricGTest, testConversion) {
    Point2DWithIndex<double> a(6, 7);
    double epsilon = 10E-6;
    double angle, radius;
    HyperbolicSpace::cartesianToPolar(a, angle, radius);
    Point2DWithIndex<double> back = HyperbolicSpace::polarToCartesian(angle, radius);
    EXPECT_NEAR(a[0], back[0], epsilon);
    EXPECT_NEAR(a[1], back[1], epsilon);
    count n = 1000;
    vector<double> angles(n);
    vector<double> radii(n);
    HyperbolicSpace::fillPoints(angles, radii, 1, 1);
    for (index i = 0; i < n; i++) {
        Point2DWithIndex<double> point = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
        double phi, r;
        HyperbolicSpace::cartesianToPolar(point, phi, r);
        EXPECT_GE(phi, 0) << "Point (" << point[0] << "," << point[1]
                          << ") was not converted correctly";
        EXPECT_GE(r, 0);
        EXPECT_LE(std::abs(phi - angles[i]), epsilon);
        EXPECT_LE(std::abs(r - radii[i]), epsilon);
    }
}

/**
 * Test numeric stability of Euclidean Circle conversion
 */
TEST_F(GeometricGTest, testEuclideanCircleConsistency) {
    const count n = 100000;
    vector<double> angles(n);
    vector<double> radii(n);
    double stretch = 2;
    double alpha = 3;
    double epsilon = 10E-6;
    double R = stretch * HyperbolicSpace::hyperbolicAreaToRadius(n);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);

    for (index i = 0; i < n; i++) {
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
    }

    for (index i = 0; i < n; i++) {
        Point2DWithIndex<double> cartesianPoint =
            HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
        double r_e, euRadius;
        HyperbolicSpace::getEuclideanCircle(radii[i], R, r_e, euRadius);
        double mirrorangle = fmod(angles[i] + PI, 2 * PI);
        double mirrorradiusInside = std::abs(r_e - euRadius) - epsilon;
        Point2DWithIndex<double> counterPointInside =
            HyperbolicSpace::polarToCartesian(mirrorangle, mirrorradiusInside);
        EXPECT_LE(HyperbolicSpace::poincareMetric(cartesianPoint, counterPointInside), R)
            << "(" << cartesianPoint.getX() << ", " << cartesianPoint.getY() << ")"
            << " and (" << counterPointInside.getX() << ", " << counterPointInside.getY() << ")"
            << " are " << HyperbolicSpace::poincareMetric(cartesianPoint, counterPointInside)
            << " apart from each other, which is more than " << R << ".";

        double mirrorradiusOutside = std::abs(r_e - euRadius) + epsilon;
        Point2DWithIndex<double> counterPointOutside =
            HyperbolicSpace::polarToCartesian(mirrorangle, mirrorradiusOutside);
        EXPECT_GE(HyperbolicSpace::poincareMetric(cartesianPoint, counterPointOutside), R);
    }
}

TEST_F(GeometricGTest, testHyperbolicTargetRadius) {
    const count n = 5000;
    const double T = 10;
    const double k = 5;
    const double alpha = 1;
    double R = HyperbolicSpace::getTargetRadius(n, (n * k) / 2, alpha, T);
    EXPECT_NEAR(R, 167.08503, 1e-4);
}

TEST_F(GeometricGTest, testNativeDistanceOfNearbyPointsFarFromOrigin) {
    // Reference values were computed with 120 significant digits.
    EXPECT_NEAR(HyperbolicSpace::nativeDistance(0, 30, 1e-12, 29.5), 2.971469618934, 1e-9);
    EXPECT_NEAR(HyperbolicSpace::nativeDistance(1e-12, 28, 0, 28), 0.708233522592, 1e-9);
}

TEST_F(GeometricGTest, testExpectedAverageDegree) {
    // Reference values were computed with a much finer quadrature.
    EXPECT_NEAR(HyperbolicSpace::getExpectedAverageDegree(10000, 20, 1, 0.5), 1.814882, 0.01);
    EXPECT_NEAR(HyperbolicSpace::getExpectedAverageDegree(10000, 22, 10. / 11, 1.1), 7.630244,
                0.01);
    EXPECT_NEAR(HyperbolicSpace::getExpectedAverageDegree(1000000, 30, 0.75, 0.9), 11.519722, 0.02);
    EXPECT_NEAR(HyperbolicSpace::getExpectedAverageDegree(2000, 60, 0.2, 5), 21.884188, 0.03);
}

TEST_F(GeometricGTest, testCalibrateTargetRadius) {
    const double n = 5000;
    const double k = 10;
    const double alpha = 1;
    for (const double T : {0.1, 0.5, 0.9, 0.99, 1.01, 1.1, 2., 8.}) {
        const double initialR = HyperbolicSpace::getTargetRadius(n, n * k / 2, alpha, T);
        const double R =
            HyperbolicSpace::calibrateTargetRadius(n, k, alpha / std::max(T, 1.), T, initialR);
        EXPECT_NEAR(HyperbolicSpace::getExpectedAverageDegree(n, R, alpha / std::max(T, 1.), T), k,
                    0.02)
            << "T = " << T;
    }

    // Graphs of two nodes cannot reach this degree at any radius.
    EXPECT_EQ(HyperbolicSpace::calibrateTargetRadius(2, 1.5, alpha, 0.5, 3), 3);
}

} /* namespace NetworKit */
