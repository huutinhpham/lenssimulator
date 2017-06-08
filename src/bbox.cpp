#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO Part 2, task 2:
  // Implement ray - bounding box intersection test
  // If the ray intersected the bounding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  Vector3D origin = r.o;
  Vector3D direction = r.d;
  Vector3D inv_d = r.inv_d;
  double tmin, tmax,
         x_tmin, x_tmax,
         y_tmin, y_tmax,
         z_tmin, z_tmax;

  //intersection along x-axis
  if (not r.sign[0]) {
    x_tmin = (min.x - origin.x)*inv_d.x;
    x_tmax = (max.x - origin.x)*inv_d.x;
  } else {
    x_tmin = (max.x - origin.x)*inv_d.x;
    x_tmax = (min.x - origin.x)*inv_d.x;
  }

  //intersection along y-axis
  if (not r.sign[1]){
    y_tmin = (min.y - origin.y)*inv_d.y;
    y_tmax = (max.y - origin.y)*inv_d.y;
  } else {
    y_tmin = (max.y - origin.y)*inv_d.y;
    y_tmax = (min.y - origin.y)*inv_d.y;
  }

  //intersection along z-axis
  if (not r.sign[2]){
    z_tmin = (min.z - origin.z)*inv_d.z;
    z_tmax = (max.z - origin.z)*inv_d.z;
  } else {
    z_tmin = (max.z - origin.z)*inv_d.z;
    z_tmax = (min.z - origin.z)*inv_d.z;
  }

  if (x_tmin > y_tmax || y_tmin > x_tmax) {
    return false;
  }

  tmin = std::fmax(x_tmin, y_tmin);
  tmax = std::fmin(x_tmax, y_tmax);

  if (tmin > z_tmax || tmax < z_tmin) {
    return false;
  }

  tmin = std::fmax(z_tmin, tmin);
  tmax = std::fmin(z_tmax, tmax);
  

  if (tmin > t1 || tmax < t0){
    return false;
  }

  t0 = std::fmax(t0, tmin);
  t1 = std::fmin(t1, tmax);

  return true;
}

void BBox::draw(Color c) const {

  glColor4f(c.r, c.g, c.b, c.a);

	// top
	glBegin(GL_LINE_STRIP);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
	glEnd();

	// bottom
	glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glEnd();

	// side
	glBegin(GL_LINES);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
	glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
	glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
	glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
