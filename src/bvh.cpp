#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // TODO Part 2, task 1:
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox centroid_box, bbox;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    Vector3D c = bb.centroid();
    centroid_box.expand(c);
  }

  // You'll want to adjust this code.
  // Right now we just return a single node containing all primitives.
  BVHNode *node = new BVHNode(bbox);
  node->prims = new vector<Primitive *>(prims);

  if (prims.size() <= max_leaf_size) {
  	return node;
  }

  //finding an axis to split along
  int split_axis = 1;
  double midPoint = centroid_box.min.x + centroid_box.extent.x/2.0;
  if (centroid_box.extent.x < centroid_box.extent.y) {
    split_axis = 2;
  	midPoint = centroid_box.min.y + centroid_box.extent.y/2.0;
  	if (centroid_box.extent.y < centroid_box.extent.z) {
    	split_axis = 3;
  		midPoint = centroid_box.min.z + centroid_box.extent.z/2.0;
  	}
  }

  //looping over all prims to split
  std::vector<Primitive *> l_prims;
  std::vector<Primitive *> r_prims;

  for (Primitive *p : prims) {
	BBox bb = p->get_bbox();
	Vector3D c_prim = bb.centroid();
	if (split_axis == 1) {
	  	if (c_prim.x <= midPoint) {
	  		l_prims.insert(l_prims.end(), p);
	  	} else {
	 		r_prims.insert(r_prims.end(), p);
	  	}
	}
	if (split_axis == 2) {
	  	if (c_prim.y <= midPoint) {
	  		l_prims.insert(l_prims.end(), p);
	  	} else {
	  		r_prims.insert(r_prims.end(), p);
		}
	}
	if (split_axis == 3) {
	  	if (c_prim.z <= midPoint) {
	  		l_prims.insert(l_prims.end(), p);
	  	} else {
	  		r_prims.insert(r_prims.end(), p);
	 	}
	}
  }


  node->l = construct_bvh(l_prims, max_leaf_size);
  node->r = construct_bvh(r_prims, max_leaf_size);

  return node;  

}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {

  // TODO Part 2, task 3:
  // Implement BVH intersection.
  // Currently, we just naively loop over every primitive.
	double t0 = ray.min_t;
  	double t1 = ray.max_t;  
  	if (node->bb.intersect(ray, t0, t1)){
  		if (node->isLeaf()){
	  		for (Primitive *p : *(node->prims)){
				if (p->intersect(ray)) {
					return true;
				}
			}
		} else {
			return intersect(ray, node->l) || intersect(ray, node->r);
		}
  	}
  	return false;	
}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {

  // TODO Part 2, task 3:
  // Implement BVH intersection.
  bool hit = false;
  bool hitR = false;
  bool hitL = false;
  Vector3D o = ray.o;
  Vector3D d = ray.d;
  double t0 = ray.min_t;
  double t1 = ray.max_t;  
  if (node->bb.intersect(ray, t0, t1)){
  	if (node->isLeaf()){
	 	for (Primitive *p : *(node->prims)){
			if (p->intersect(ray, i)) {
				hit = true;
				++total_isects;
			}
		}
	} else {
		Intersection iL = Intersection();
		Intersection iR = Intersection();
		Intersection* piL = &iL;
		Intersection* piR = &iR;
		Ray rayR = Ray(o, d);
		rayR.min_t = ray.min_t;
		rayR.max_t = ray.max_t;
		Ray rayL = Ray(o, d);
		rayL.min_t = ray.min_t;
		rayL.max_t = ray.max_t;
		Ray* p_rayR = &rayR;
		Ray* p_rayL = &rayL;
		hitR = intersect(*p_rayR, piR, node->r);
		hitL = intersect(*p_rayL, piL, node->l);
		if (hitR && hitL){
			hit = true;
			if (piR->t < piL->t){
				i->t = piR->t;
				i->primitive = piR->primitive;
				i->n = piR->n;
				i->bsdf = piR->bsdf;
			} else {
				i->t = piL->t;
				i->primitive = piL->primitive;
				i->n = piL->n;
				i->bsdf = piL->bsdf;
			}
		} else if (hitR) {
			hit = true;
			i->t = piR->t;
			i->primitive = piR->primitive;
			i->n = piR->n;
			i->bsdf = piR->bsdf;
		} else if (hitL){
			hit = true;
			i->t = piL->t;
			i->primitive = piL->primitive;
			i->n = piL->n;
			i->bsdf = piL->bsdf;
		}
	}
  }
  return hit;

}

}  // namespace StaticScene
}  // namespace CGL
