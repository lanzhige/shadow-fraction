#include "fractioncalculator.h"

#include <cmath>

namespace calculate {
FractionCalculator::FractionCalculator(double width, double height) {
  this->width = width;
  this->height = height;
  double du = 2 / width;
  double dv = 2 / height;

  mask.resize(width * height);

  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      double u = i / width * 2 - 1;
      double v = j / height * 2 - 1;

      std::vector<Point3d> corners = getPixelCorners(u, v, du, dv);
      double p1 = getSurfaceArea(corners[0], corners[1], corners[2]);
      double p2 = getSurfaceArea(corners[1], corners[2], corners[3]);
      mask[i * width + j] = (p1 + p2) / (M_PI * 4);
    }
  }
}

double FractionCalculator::getLat(const Point3d &p) {
  double r = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
  return asin(p.z / r);
}

double FractionCalculator::getLon(const Point3d &p) {
  return atan(p.y / p.x);
}

double FractionCalculator::getD(double lat0, double lon0,
                                double lat1, double lon1) {
  return 2 * asin(
      sqrt(
          (1 - cos(lat1 - lat0)) / 2 
        + (1 - cos(lon1 - lon0)) / 2 * cos(lat1) * cos(lat0)
      )
  );
}

double FractionCalculator::getA(double x, double y, double z) {
  return acos( (cos(x) - cos(y) * cos(z)) / (sin(y) * sin(z)) );
}

double FractionCalculator::getE(double x, double y, double z) {
  return x + y + z - M_PI;
}

Point3d FractionCalculator::p(double u, double v) {
  double r = sqrt(u * u + v * v + 1);
  return {r, r*u, r*v};
}

std::vector<Point3d> FractionCalculator::getPixelCorners(double u, double v,
                                                         double du, double dv) {
  return{
    p(u, v),
    p(u + du, v),
    p(u, v + dv),
    p(u + du, v + dv)
  };
}

double FractionCalculator::getSurfaceArea(const Point3d &x, const Point3d &y,
                                          const Point3d &z) {
  double dz = getD(getLat(x), getLon(x), getLat(y), getLon(y));
  double dy = getD(getLat(x), getLon(x), getLat(z), getLon(z));
  double dx = getD(getLat(y), getLon(y), getLat(z), getLon(z));

  double X = getA(dx, dy, dz);
  double Y = getA(dy, dz, dx);
  double Z = getA(dz, dx, dy);

  return getE(X, Y, Z);
}

}