#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.

  return false;

}

bool Sphere::intersect(const Ray& r) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double a = dot(r.d, r.d);
  double b = dot(2*(r.o - o), r.d);
  double c = dot((r.o - o), (r.o - o)) - r2;

  if ((pow(b, 2) - 4*a*c) < 0) {
    return false;
  }
  double t0 = (-b - sqrt(pow(b, 2) - 4*a*c))/(2*a);
  double t1 = (-b + sqrt(pow(b, 2) - 4*a*c))/(2*a);
  double min_t = min(t0, t1);
  double max_t = max(t0, t1);
  double t = min_t;
  if(t <= 0) {
    t = max_t;
  }

  if (t > r.max_t || t < r.min_t) {
    return false;
  }

  r.max_t = t;
  
  return true;
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO Part 1m task 4:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double a = dot(r.d, r.d);
  double b = dot(2*(r.o - o), r.d);
  double c = dot((r.o - o), (r.o - o)) - r2;

  if ((pow(b, 2) - 4*a*c) < 0) {
    return false;
  }

  double t0 = (-b - sqrt(pow(b, 2) - 4*a*c))/(2*a);
  double t1 = (-b + sqrt(pow(b, 2) - 4*a*c))/(2*a);
  double min_t = min(t0, t1);
  double max_t = max(t0, t1);
  double t = min_t;
  if(t <= 0) {
    t = max_t;
  }

  if (t > r.max_t || t < r.min_t) {
    return false;
  }

  r.max_t = t;
  i->t = t;
  i->n = ((r.o + t*r.d) - o).unit();
  i->primitive = this;
  i->bsdf = get_bsdf();

  
  return true;

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
