#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return albedo * (1.0 / PI);
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi = sampler.get_sample(pdf);
  return albedo * (1.0 / PI);
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return reflectance;
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Implement MirrorBSDF
  reflect(wo, wi);
  *pdf = 1;
  return reflectance;
}

// Glossy BSDF //

/*
Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0f;
  return reflect(wo, wi, reflectance);
}
*/

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return transmittance;
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Implement RefractionBSDF
  if (refract(wo, wi, ior)) {
    *pdf = 1;
    return transmittance;
  }
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return transmittance;
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Compute Fresnel coefficient and either reflect or refract based on it.
  if (!refract(wo, wi, ior)){
    reflect(wo, wi);
    *pdf = 1;
    return reflectance/cos_theta(wo);
  } else {
    double ni = 1.0;
    double no = ior;
    if (wo.z < 0){
      ni = ior;
      no = 1.0;
    }
    double R = pow((ni - no)/(ni + no), 2);
    if (coin_flip(R)){
      reflect(wo, wi);
      *pdf = R;
      return R*reflectance/cos_theta(wo);
    } else {
      refract(wo, wi, ior);
      *pdf = 1 - R;
      return transmittance*(1 - R)*pow((no/ni), 2)/cos_theta(wo);
    }
  }
  
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO Part 5:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.

  Vector3D n = Vector3D(0, 0, 1);
  *wi = -wo + 2*(dot(wo, n))*n;
  

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO Part 5:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  double ni = 1.0;
  double no = ior;
  Vector3D w_in = -wo;
  double sin_the_i = sin_theta(w_in);
  Vector3D n = Vector3D(0.0, 0.0, 1.0);
  if (wo.z < 0) {
    ni = ior;
    no = 1.0;
    n = Vector3D(0.0 , 0.0, -1.0);
  }

  if (ni*sin_the_i >= no) {
    return false;
  }

  double r = ni/no;
  double c = dot(-n, w_in);

  *wi = r*w_in + (r*c - sqrt(1 - pow(r, 2)*(1 - pow(c, 2))))*n;

  return true;

}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL