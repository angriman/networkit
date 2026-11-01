/*
 * HyperbolicSpace.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/geometric/HyperbolicSpace.hpp>

using std::abs;
using std::max;

namespace NetworKit {

namespace {

/**
 * @return Nodes and weights of the Gauss-Legendre quadrature rule on [lower, upper]
 */
std::pair<vector<double>, vector<double>> gaussLegendre(count size, double lower, double upper) {
    vector<double> nodes(size), weights(size);
    for (index i = 0; i < size; ++i) {
        double x = std::cos(PI * (i + 0.75) / (size + 0.5));
        double derivative = 0;
        for (index iteration = 0; iteration < 100; ++iteration) {
            double previous = 1;
            double current = x;
            for (index degree = 2; degree <= size; ++degree) {
                const double next =
                    ((2 * degree - 1) * x * current - (degree - 1) * previous) / degree;
                previous = current;
                current = next;
            }
            derivative = size * (x * current - previous) / (x * x - 1);
            const double step = current / derivative;
            x -= step;
            if (std::abs(step) < 1e-14) {
                break;
            }
        }
        nodes[i] = lower + (x + 1) / 2 * (upper - lower);
        weights[i] = (upper - lower) / ((1 - x * x) * derivative * derivative);
    }
    return {std::move(nodes), std::move(weights)};
}

} // namespace

double HyperbolicSpace::getExpectedAverageDegree(double n, double R, double alpha, double T) {
    assert(n > 1);
    assert(R > 0);
    assert(alpha > 0);
    assert(T > 0);

    // Most of the degree comes from small angles and radii, so integrate over the logarithm of the
    // angle and over the radius itself rather than over its much more skewed distribution.
    constexpr count radialPoints = 24;
    constexpr count angularPoints = 64;
    const auto [radii, radialWeights] = gaussLegendre(radialPoints, 0, R);
    const auto [logAngles, logAngleWeights] =
        gaussLegendre(angularPoints, std::log(1e-7 * std::min(1., 50. / n)), std::log(PI));

    const double logTwo = std::log(2.);
    const double normalization = std::pow(1 - std::exp(-alpha * R), 2);
    vector<double> radialMass(radialPoints), logSinh(radialPoints);
    for (index i = 0; i < radialPoints; ++i) {
        const double density = alpha * std::exp(alpha * (radii[i] - R))
                               * (1 - std::exp(-2 * alpha * radii[i])) / normalization;
        radialMass[i] = radialWeights[i] * density;
        logSinh[i] = radii[i] + std::log1p(-std::exp(-2 * radii[i])) - logTwo;
    }

    vector<double> angularMass(angularPoints), logTwoSinSquared(angularPoints);
    for (index i = 0; i < angularPoints; ++i) {
        const double angle = std::exp(logAngles[i]);
        angularMass[i] = logAngleWeights[i] * angle;
        logTwoSinSquared[i] = logTwo + 2 * std::log(std::sin(angle / 2));
    }

    double expected = 0;
    for (index i = 0; i < radialPoints; ++i) {
        for (index j = 0; j <= i; ++j) {
            const double radialDifference = std::abs(radii[i] - radii[j]);
            const double coshDifference = std::cosh(radialDifference);
            double probability = 0;
            for (index a = 0; a < angularPoints; ++a) {
                // cosh(d) = cosh(r1 - r2) + 2 sinh(r1) sinh(r2) sin^2(angle / 2). Large values are
                // added in log space because the products would overflow at high temperatures.
                const double logAngularTerm = logSinh[i] + logSinh[j] + logTwoSinSquared[a];
                double distance;
                if (std::max(logAngularTerm, radialDifference) < 30) {
                    distance = std::acosh(coshDifference + std::exp(logAngularTerm));
                } else {
                    const double logDoubled = logAngularTerm + logTwo;
                    const double larger = std::max(logDoubled, radialDifference);
                    const double smaller = std::min(logDoubled, radialDifference);
                    distance = larger + std::log1p(std::exp(smaller - larger));
                }
                probability += angularMass[a] / (std::exp((distance - R) / (2 * T)) + 1);
            }
            expected += (i == j ? 1 : 2) * radialMass[i] * radialMass[j] * probability;
        }
    }
    return (n - 1) * expected / PI;
}

double HyperbolicSpace::calibrateTargetRadius(double n, double k, double alpha, double T,
                                              double initialR) {
    const auto logDegreeRatio = [n, k, alpha, T](double R) -> double {
        return std::log(getExpectedAverageDegree(n, R, alpha, T) / k);
    };

    // The degree falls with the radius, so first enclose the target.
    double smallR = std::max(initialR, 1.);
    double largeR = smallR;
    double smallRatio = logDegreeRatio(smallR);
    double largeRatio = smallRatio;
    while (smallRatio < 0) {
        if (smallR < 1e-3) {
            return initialR;
        }
        largeR = smallR;
        largeRatio = smallRatio;
        smallR /= 2;
        smallRatio = logDegreeRatio(smallR);
    }
    while (largeRatio > 0) {
        smallR = largeR;
        smallRatio = largeRatio;
        largeR *= 2;
        largeRatio = logDegreeRatio(largeR);
    }

    // The logarithm of the degree is almost linear in the radius, so false position is fast.
    // Damping the bound that stays in place keeps it from stalling on the remaining curvature.
    constexpr double tolerance = 1e-3;
    double R = smallR;
    for (index iteration = 0; iteration < 100; ++iteration) {
        R = smallR + (largeR - smallR) * smallRatio / (smallRatio - largeRatio);
        const double ratio = logDegreeRatio(R);
        if (std::abs(ratio) < tolerance) {
            break;
        }
        if (ratio > 0) {
            smallR = R;
            smallRatio = ratio;
            largeRatio /= 2;
        } else {
            largeR = R;
            largeRatio = ratio;
            smallRatio /= 2;
        }
    }
    return R;
}

double HyperbolicSpace::nativeDistance(double firstangle, double firstR, double secondangle,
                                       double secondR) {
    assert(firstR >= 0);
    assert(secondR >= 0);
    assert(firstangle >= 0);
    assert(firstangle < 2 * PI);
    assert(secondangle >= 0);
    assert(secondangle < 2 * PI);
    double result;
    if (firstangle == secondangle) {
        result = std::abs(firstR - secondR);
    } else {
        // The usual law of cosines subtracts two huge numbers for points far from the origin, which
        // rounds the distance of nearby points to zero.
        double sinHalfDelta = std::sin((firstangle - secondangle) / 2);
        double coshDist =
            std::cosh(firstR - secondR)
            + 2 * std::sinh(firstR) * std::sinh(secondR) * sinHalfDelta * sinHalfDelta;
        result = std::acosh(coshDist);
    }
    assert(result >= 0);
    return result;
}

/**
 * This distance measure is taken from the Poincare disc model.
 */
double HyperbolicSpace::poincareMetric(double phi_a, double r_a, double phi_b, double r_b) {
    assert(r_a < 1);
    assert(r_b < 1);
    return poincareMetric(polarToCartesian(phi_a, r_a), polarToCartesian(phi_b, r_b));
}

double HyperbolicSpace::poincareMetric(Point2DWithIndex<double> a, Point2DWithIndex<double> b) {
    assert(a.length() < 1);
    assert(b.length() < 1);
    double result = std::acosh(
        1 + 2 * a.squaredDistance(b) / ((1 - a.squaredLength()) * (1 - b.squaredLength())));
    assert(result >= 0);
    return result;
}

void HyperbolicSpace::fillPoints(vector<double> &angles, vector<double> &radii, double R,
                                 double alpha) {
    fillPoints(angles, radii, 0, 2 * PI, 0, R, alpha);
}

void HyperbolicSpace::fillPoints(vector<double> &angles, vector<double> &radii, double minPhi,
                                 double maxPhi, double minR, double maxR, double alpha) {
    uint64_t n = radii.size();
    assert(angles.size() == n);

    double mincdf = std::cosh(alpha * minR);
    double maxcdf = std::cosh(alpha * maxR);
    std::uniform_real_distribution<double> phidist{minPhi, maxPhi};
    std::uniform_real_distribution<double> rdist{mincdf, maxcdf};

    for (uint64_t i = 0; i < n; i++) {
        angles[i] = phidist(Aux::Random::getURNG());
        /**
         * for the radial coordinate distribution, I took the probability density from Greedy
         * Forwarding in Dynamic Scale-Free Networks Embedded in Hyperbolic Metric Spaces f (r) =
         * sinh r/(cosh R  1) \int sinh = cosh+const
         */
        double random = rdist(Aux::Random::getURNG());
        radii[i] = (std::acosh(random) / alpha);
        assert(radii[i] <= maxR);
        assert(radii[i] >= minR);
        assert(angles[i] <= maxPhi);
        assert(angles[i] >= minPhi);
        if (radii[i] == maxR)
            radii[i] = std::nextafter(radii[i], 0);
        assert(radii[i] < maxR);
    }
}
void HyperbolicSpace::fillPointsSorted(vector<double> &angles, vector<double> &radii, double R,
                                       double alpha) {
    fillPointsSorted(angles, radii, 0, 2 * PI, 0, R, alpha);
}

void HyperbolicSpace::fillPointsSorted(vector<double> &angles, vector<double> &radii, double minPhi,
                                       double maxPhi, double minR, double maxR, double alpha) {
    uint64_t n = radii.size();
    assert(angles.size() == n);

    double mincdf = std::cosh(alpha * minR);
    double maxcdf = std::cosh(alpha * maxR);
    std::uniform_real_distribution<double> uniformDist{0, 1};
    std::uniform_real_distribution<double> rdist{mincdf, maxcdf};
    double sum = 0;

    for (uint64_t i = 0; i < n; i++) {
        double value = -std::log(uniformDist(Aux::Random::getURNG()));
        sum += value;
        angles[i] = sum;
        /**
         * for the radial coordinate distribution, I took the probability density from Greedy
         * Forwarding in Dynamic Scale-Free Networks Embedded in Hyperbolic Metric Spaces f (r) =
         * sinh r/(cosh R  1) \int sinh = cosh+const
         */
        double random = rdist(Aux::Random::getURNG());
        radii[i] = (std::acosh(random) / alpha);
        assert(radii[i] <= maxR);
        assert(radii[i] >= minR);
        if (radii[i] == maxR)
            radii[i] = std::nextafter(radii[i], 0);
        assert(radii[i] < maxR);
    }

    sum += -std::log(uniformDist(Aux::Random::getURNG()));
    double inverseSum = (maxPhi - minPhi) / sum;

    for (uint64_t i = 0; i < n; i++) {
        angles[i] = minPhi + angles[i] * inverseSum;
        assert(angles[i] <= maxPhi);
        assert(angles[i] >= minPhi);
    }
}

Point2DWithIndex<double> HyperbolicSpace::polarToCartesian(double phi, double r) {
    return Point2DWithIndex<double>(r * std::cos(phi), r * std::sin(phi));
}

std::map<index, Point<float>> HyperbolicSpace::polarToCartesian(const vector<double> &angles,
                                                                const vector<double> &radii) {
    assert(radii.size() == angles.size());
    std::map<index, Point<float>> result;
    for (index i = 0; i < angles.size(); i++) {
        Point2DWithIndex<double> coord = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
        result.emplace(i, Point<float>{static_cast<float>(coord[0]), static_cast<float>(coord[1])});
    }
    return result;
}

void HyperbolicSpace::cartesianToPolar(Point2DWithIndex<double> a, double &phi, double &r) {
    r = a.length();
    if (r == 0)
        phi = 0;
    else if (a[1] >= 0) {
        phi = std::acos(a[0] / r);
    } else {
        phi = -std::acos(a[0] / r);
    }
    if (phi < 0)
        phi += 2 * PI;
}

void HyperbolicSpace::getEuclideanCircle(Point2DWithIndex<double> hyperbolicCenter,
                                         double hyperbolicRadius,
                                         Point2DWithIndex<double> &euclideanCenter,
                                         double &euclideanRadius) {
    double phi_h, r_h;
    HyperbolicSpace::cartesianToPolar(hyperbolicCenter, phi_h, r_h);
    double r_c;
    HyperbolicSpace::getEuclideanCircle(r_h, hyperbolicRadius, r_c, euclideanRadius);
    euclideanCenter = HyperbolicSpace::polarToCartesian(phi_h, r_c);
}

void HyperbolicSpace::getEuclideanCircle(double r_h, double hyperbolicRadius,
                                         double &radialCoordOfEuclideanCenter,
                                         double &euclideanRadius) {
    double a = std::cosh(hyperbolicRadius) - 1;
    double b = 1 - (r_h * r_h);
    radialCoordOfEuclideanCenter = (2 * r_h) / (b * a + 2);
    euclideanRadius = std::sqrt(radialCoordOfEuclideanCenter * radialCoordOfEuclideanCenter
                                - (2 * r_h * r_h - b * a) / (b * a + 2));
}

double HyperbolicSpace::hyperbolicRadiusToEuclidean(double hyperbolicRadius) {
    double ch = std::cosh(hyperbolicRadius);
    return std::sqrt((ch - 1) / (ch + 1));
}

double HyperbolicSpace::EuclideanRadiusToHyperbolic(double euclideanRadius) {
    double eusq = euclideanRadius * euclideanRadius;
    double result = std::acosh(1 + 2 * eusq / ((1 - eusq)));
    return result;
}

double HyperbolicSpace::maxRinSlice(double minPhi, double maxPhi, double phi_c, double r_c,
                                    double euRadius) {
    double maxCos =
        std::max(std::cos(std::abs(minPhi - phi_c)), std::cos(std::abs(maxPhi - phi_c)));

    if (minPhi < phi_c && phi_c < maxPhi)
        maxCos = 1;
    // applying law of cosines here
    double p = r_c * maxCos;
    double maxR = p + std::sqrt(p * p - r_c * r_c + euRadius * euRadius);
    return maxR;
}

double HyperbolicSpace::hyperbolicSpaceInEuclideanCircle(double r_c, double d_c, double r_max) {
    double result = 0;
    assert(r_c >= 0);
    assert(d_c >= 0);
    assert(r_c <= r_max);
    double min = r_c - d_c;
    double max = std::min(r_c + d_c, r_max);

    if (d_c > r_c) {
        // the query circle overlaps the origin

        if (d_c - r_c < r_max) {
            // remaining query circle is smaller than the disk representation
            // adding small circle around origin
            result += 2 * PI * (std::cosh(EuclideanRadiusToHyperbolic(d_c - r_c)) - 1);
        } else {
            // adding small circle around origin
            result += 2 * PI * (std::cosh(EuclideanRadiusToHyperbolic(r_max)) - 1);
        }
        assert(result <= 2 * PI * (std::cosh(EuclideanRadiusToHyperbolic(r_max)) - 1));
        // correcting integral start to exclude circle
        min = std::nextafter(d_c - r_c, std::numeric_limits<double>::max());
    }

    /**
     * Now, the integral.
     * It is 4\int_{min}^{max} \text{acos}(\frac{r_c^2-d_c^2+r^2}{2r_c\cdot r})  \cdot
     * \frac{1}{1-r^2} \cdot (\sinh (\text{acosh}( 1 + 2\frac{r^2}{1 - r^2})))\,dr This solution was
     * computed by WolframAlpha
     */

    if (max < min)
        return result;

    auto realpart = [](double r, double d, double c) {
        double result = std::acos((c * c - d * d + r * r) / (2 * c * r)) / (r * r - 1);
        return result;
    };

    auto firstlogpart = [](double r, double d, double c) {
        double s = (c * c - d * d);
        double rsqs = r * r + s;
        double real = -2 * s * std::sqrt(4 * c * c * r * r - rsqs * rsqs);
        double imag = -4 * c * c * r * r + 2 * s * r * r + 2 * s * s;
        return std::atan2(imag, real) / 2;
    };

    auto secondlogpart = [](double r, double d, double c) {
        double s = (c * c - d * d);
        double rsqs = r * r + s;
        double real = std::sqrt(4 * c * c * r * r - rsqs * rsqs);
        double imag = 2 * c * c * (r * r + 1) - (s + 1) * rsqs;
        imag = imag / std::sqrt((s + 1) * (s + 1) - (4 * c * c));
        return (s - 1) * std::atan2(imag, real) / (2 * std::sqrt((s + 1) * (s + 1) - (4 * c * c)));
    };

    double lower =
        -realpart(min, d_c, r_c) - firstlogpart(min, d_c, r_c) + secondlogpart(min, d_c, r_c);
    double upper =
        -realpart(max, d_c, r_c) - firstlogpart(max, d_c, r_c) + secondlogpart(max, d_c, r_c);
    result += 4 * (upper - lower);
    return result;
}
} // namespace NetworKit
