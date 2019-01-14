#ifndef SHADOWFRACTION_FRACTIONCALCULATOR_H_
#define SHADOWFRACTION_FRACTIONCALCULATOR_H_

#include <vector>

namespace calculate {
struct Point3d {
  double x;
  double y;
  double z;
  Point3d(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  Point3d& operator=(const Point3d& p) {
    if (this != &p) {
      this->x = p.x;
      this->y = p.y;
      this->z = p.z;
    }
    return *this;
  }
};


class FractionCalculator {
#define M_PI 3.14159265358979323846

public:
  FractionCalculator(double width = 512, double height = 512);
private:
  double getLat(const Point3d &p);
  double getLon(const Point3d &p);
  double getD(double lat0, double lon0, double lat1, double lon1);
  double getA(double x, double y, double z);
  double getE(double x, double y, double z);
  Point3d p(double u, double v);
  std::vector<Point3d> getPixelCorners(double u, double v, double du, double dv);
  double getSurfaceArea(const Point3d &x, const Point3d &y, const Point3d &z);

  std::vector<double> mask;
  double width, height;
};
}

#endif // !SHADOWFRACTION_FRACTIONCALCULATOR_H_
