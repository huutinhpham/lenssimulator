#include "lenscamera.h"

#include "image.h"

using namespace std;

namespace CGL {


/****** Helpers ******/
  

// Extract the R, G, or B channel out of an RGBA color stored in a single 32bit integer
static uint32_t red_channel(uint32_t color) {
    return (255 & (color >> (0)));
}

static uint32_t green_channel(uint32_t color) {
    return (255 & (color >> (8)));
}

static uint32_t blue_channel(uint32_t color) {
    return (255 & (color >> (16)));
}

// Convert from millimeters to meters
static const double scale = .001;







/****** LensElement functions ******/


bool LensElement::pass_through(Ray &r, double &prev_ior) const {
  // Part 1 Task 1: Implement this. It takes r and passes it through this lens element.
  Vector3D hit_p;
  if (intersect(r, &hit_p)) {
    if (refract(r, hit_p, prev_ior)){
      prev_ior = ior;
      return true;
    }
  }

  return false;
}

bool LensElement::intersect(const Ray &r, Vector3D *hit_p) const {
  // Part 1 Task 1: Implement this. It intersects r with this spherical lens elemnent 
  // (or aperture diaphragm). You'll want to reuse some sphere intersect code.

  if (radius == 0) {
    double t = (center - r.o.z)/r.d.z;
    *hit_p = r.o + t*r.d;
    double hit_x = hit_p->x;
    double hit_y = hit_p->y;

  if (sqrt(pow(hit_x,2) + pow(hit_y, 2)) > aperture/2) {
      return false;
  }

  return true;

  }else{
    Vector3D v_center = Vector3D(0.0, 0.0, center);
    double a = dot(r.d, r.d);
    double b = dot(2*(r.o - v_center), r.d);
    double c = dot((r.o - v_center), (r.o - v_center)) - pow(radius, 2);

    if ((pow(b, 2) - 4*a*c) < 0) {
      return false;
    }

    double t0 = (-b - sqrt(pow(b, 2) - 4*a*c))/(2*a);
    double t1 = (-b + sqrt(pow(b, 2) - 4*a*c))/(2*a);
    double min_t = min(t0, t1);
    double max_t = max(t0, t1);

    double t = min_t;
    if (radius*r.d.z < 0) {
      t = max_t;
    }

    *hit_p = r.o + t*r.d;
    double hit_x = hit_p->x;
    double hit_y = hit_p->y;

    if (sqrt(pow(hit_x, 2) + pow(hit_y, 2)) > aperture/2) {
      return false;
    }

    return true;
  }
  
}
bool LensElement::refract(Ray& r, const Vector3D& hit_p, const double& prev_ior) const {
  // Part 1 Task 1: Implement this. It refracts the Ray r with this lens element or 
  // does nothing at the aperture element.
  // You'll want to consult your refract function from the previous assignment.

  if (radius == 0) {
    return true;
  }

  double ior_in = prev_ior;
  double ior_out = ior;

  Vector3D v_center = Vector3D(0.0, 0.0, center);
  Vector3D hit_p_normal = (hit_p - v_center).unit();
  if (dot(r.d, hit_p_normal)*r.d.z < 0) {
    hit_p_normal = -hit_p_normal;
  }

  Matrix3x3 o2w;
  make_coord_space(o2w, hit_p_normal);
  Matrix3x3 w2o = o2w.T();

  Vector3D w_in = w2o*(r.d);
  double sin_the_i = sin_theta(w_in);
  Vector3D n = Vector3D(0.0, 0.0, 1.0);

  if (ior_in*sin_the_i >= ior_out) {
    return false;
  }

  double const_r = ior_in/ior_out;
  if (r.d.z > 0) {
    const_r = 1/const_r;
    n.z = -1.0;
  }
  double c = dot(-n, w_in);

  Vector3D w_out = const_r*w_in + (const_r*c - sqrt(1 - pow(const_r, 2)*(1 - pow(c, 2))))*n;
  w_out = o2w*(w_out);

  r.o = hit_p;
  r.d = w_out.unit();

  return true;
}






/****** Lens functions ******/



void Lens::parse_lens_file(std::string filename) {

  ifstream infile(filename);
  string line;
  double z_coord = 0;
  double z_ap;
  vector<LensElement> backwards;
  elts.clear();
  bool first = true;
  while (getline(infile, line)) {
    if (first) {
      cout << "[Lens] Loading lens file " << line << endl;
      first = false;
    }
    if (line[0] == '#')
      continue;
    stringstream ss(line);
    LensElement lens;
    double offset;
    ss >> lens.radius >> offset >> lens.ior >> lens.aperture;
    lens.center = z_coord;
    if (!lens.radius) {
      z_ap = z_coord;
    }
    z_coord += offset;
    backwards.push_back(lens);
  }
  for (int i = backwards.size() - 1; i >= 0; --i) {
    LensElement l = backwards[i];
    l.center = (l.center - z_ap) + l.radius;
    if (i) l.ior = backwards[i-1].ior;
    else l.ior = 1;
    if (!l.ior) l.ior = 1;
    elts.push_back(l);
    if (!l.radius)
      ap_i = elts.size()-1;
    // cout << "Lens element edge first " << (l.center - l.radius) << " " 
    //   << l.radius << " " << l.center << " " << l.ior << " " << l.aperture << endl;
  }
  double c = elts.front().center, r = elts.front().radius, a = elts.front().aperture * .5;
  back_elt = c - (r>0?1:-1) * sqrt(r*r-a*a);
  ap_radius = ap_original = elts[ap_i].aperture;

  // Get infinity and close focus depths, also get focal length.
  set_focus_params();
  // Focus at infinity to start.
  sensor_depth = infinity_focus;
       
}


void Lens::set_focus_params() {

  // Part 1 Task 2: Implement this. 
  // After this function is called, the three variables
  // infinity_focus, near_focus, and focal_length
  // should be set correctly.

  double paraxial_dist = exp(1.0);
  Ray inf_inc_ray = Ray(Vector3D(paraxial_dist, 0, -99999), Vector3D(0, 0, 1));
  while (not trace_backwards(inf_inc_ray)) {
    paraxial_dist = paraxial_dist/2;
    inf_inc_ray = Ray(Vector3D(paraxial_dist, 0, -99999), Vector3D(0, 0, 1));
  }
  double inf_t = (0 - inf_inc_ray.o.x)/inf_inc_ray.d.x;

  infinity_focus = (inf_inc_ray.o + inf_t*inf_inc_ray.d).z;

  double t_p = (paraxial_dist - inf_inc_ray.o.x)/inf_inc_ray.d.x;
  double p_prime = (inf_inc_ray.o + t_p*inf_inc_ray.d).z;

  focal_length = abs(p_prime - infinity_focus);



  double near_inc_ray_z = elts.front().center - elts.front().radius -
                          (1 + log(focal_length))*focal_length;
  Vector3D near_o = Vector3D(0, 0, near_inc_ray_z);
  Vector3D near_d = (Vector3D(elts.back().aperture/10.0, 0, 0) - near_o).unit();
  Ray near_inc_ray = Ray(near_o, near_d);
  while (not trace_backwards(near_inc_ray)){
    near_d.x = near_d.x/2;
    near_inc_ray = Ray(near_o, near_d.unit());
  }

  double near_t = (0.0 - near_inc_ray.o.x)/near_inc_ray.d.x;

  near_focus = (near_inc_ray.o + near_t*near_inc_ray.d).z;


  cout << "[Lens] Infinity focus depth is " << infinity_focus << endl;
  cout << "[Lens] Close focus depth is " << near_focus << endl;
  cout << "[Lens] True focal length is " << focal_length << endl;
}




bool Lens::trace(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the sensor out into the world.
  double prev_ior = 1.0;
  for (LensElement elt: elts) {
    if (elt.pass_through(r, prev_ior)){
      if (trace){
        trace->push_back(r.o);
      }
    } else {
      return false;
    } 
  }

  return true;
}

bool Lens::trace_backwards(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the world backwards through 
  // the lens towards the sensor.
  double prev_ior;
  for (int i = elts.size() - 1; i >= 0; --i) {
    if (i > 0) {
      prev_ior = elts[i - 1].ior;
    } else {
      prev_ior = 1.0;
    }
    if (elts[i].pass_through(r, prev_ior)){
      if (trace) {
        trace->push_back(r.o);
      }
    } else {
      return false;
    }
  }

  return true;
}

float Lens::focus_depth(float d) const {

  // Part 1 Task 2: Implement this. Should find the conjugate of a ray
  // starting from the sensor at depth d.
  double offset = elts.back().aperture/10.0;
  Ray r = Ray(Vector3D(0, 0, d), Vector3D(offset, 0, -1).unit());
  while (not trace(r)) {
    offset = offset/2;
    r = Ray(Vector3D(0, 0, d), Vector3D(offset, 0, -1).unit());
  }

  double t = (0.0 - r.o.x)/r.d.x;

  return (r.o + t*r.d).z;
}

Vector3D Lens::back_lens_sample() const {

  // Part 1 Task 2: Implement this. Should return a point randomly sampled
  // on the back element of the lens (the element closest to the sensor)

  double center = elts[0].center;
  double radius = elts[0].aperture*0.5;
  double diameter = radius*2;
  double z = center - (radius>0?1:-1) * sqrt(radius*radius - diameter*diameter*0.25);
  double sampled_x = random_uniform()*diameter - radius;
  double sampled_y = random_uniform()*diameter - radius;
  while (sqrt(pow(sampled_x, 2) + pow(sampled_y, 2)) > radius) {
    sampled_x = random_uniform()*diameter - radius;
    sampled_y = random_uniform()*diameter - radius;
  }

  return Vector3D(sampled_x, sampled_y, z);

}



/****** LensCamera functions ******/


LensCamera::LensCamera(): pt(NULL) {
  string path = string(__FILE__).substr(0,string(__FILE__).find_last_of('/')+1) + "../lenses/";
  static const vector<string> lens_files = {"dgauss.50mm.dat", "wide.22mm.dat", "telephoto.250mm.dat", "fisheye.10mm.dat"};
  for (string lens_file : lens_files)
    lenses.emplace_back(path + lens_file);

  mount_lens(0);
}


Ray LensCamera::generate_ray(double x, double y) const {

  Ray r = Ray(Vector3D(),Vector3D() );
  if (lens_ind >= 0) {

    // Part 1 Task 2: Implement this. It generates a ray from sensor pixel (x,y)
    // pointing toward the back element of the lens (use back_lens_sample) and traces
    // it through the Lens (using your "trace" function)

    double film_d = sqrt(24*24 + 36*36);
    double screen_d = sqrt(screenW*screenW + screenH*screenH);
    double film_w = film_d * screenW / screen_d;
    double film_h = film_d * screenH / screen_d;
    Vector3D sensor_point= Vector3D(-(x-0.5)*film_w, -(y-0.5)*film_h, lenses[lens_ind].sensor_depth);

    Vector3D sample = lenses[lens_ind].back_lens_sample();
    r.o = sensor_point;
    r.d = (sample - sensor_point).unit();

    if (lenses[lens_ind].trace(r)) {
      r.o = pos + c2w * r.o * scale;
      r.d = (c2w * r.d).unit();
    } else {
      r.o = pos + c2w * r.o * scale;
      r.d = (c2w * Vector3D(0, 0, 1)).unit();
    }

    r.inv_d.x = 1/r.d.x;
    r.inv_d.y = 1/r.d.y;
    r.inv_d.z = 1/r.d.z;

    r.sign[0] = (r.inv_d.x < 0);
    r.sign[1] = (r.inv_d.y < 0);
    r.sign[2] = (r.inv_d.z < 0);


    /***** end of your code ******/


    // This code converts the ray you traced through the lens into world coordinates.

  } else {

    // Generate ray for a pinhole camera. Same as in the previous assignment.
    x = 2*(x-.5); y = 2*(y-.5);
    r = Ray(pos,(c2w*Vector3D(x*tan(radians(hFov)*.5),y*tan(radians(vFov)*.5),-1)).unit());

  }

  r.min_t = nClip; r.max_t = fClip;
  return r;
}



void LensCamera::move_sensor(float delta) {
  if (lens_ind < 0) return;
  curr_lens().sensor_depth += delta;
  cout << "[LensCamera] Sensor plane moved to " << curr_lens().sensor_depth
       << ", focus now at " << lenses[lens_ind].focus_depth(lenses[lens_ind].sensor_depth) << endl;
}

void LensCamera::stop_down(float ratio) {
  float ap = curr_lens().ap_radius * ratio;
  if (ap > curr_lens().ap_original) ap = curr_lens().ap_original;
  curr_lens().ap_radius = ap;
  cout << "[LensCamera] Aperture is now " << curr_lens().ap_radius << "mm" << endl;
}

void LensCamera::mount_lens(int i) {
  lens_ind = i;
  if (i >= 0) {
    cout << "[LensCamera] Switched to lens #" << (i+1) 
         << " with focal length " << curr_lens().focal_length << "mm" << endl;
  } else {
    cout << "[LensCamera] Switched to pinhole camera" << endl;
  }
}



// A dummy function to demonstrate how to work with the image buffer.
// Calculates the average value of the green color channel in the image.
// You'll have to remember your 2D array indexing in order to take differences
// of neighboring pixels in a more sophisticated metric function.
static double mean_green(const ImageBuffer& ib) {
  double sum = 0;
  for (int i = 0; i < ib.w * ib.h; ++i) {
      sum += green_channel(ib.data[i]);
  }
  double mean = sum / (ib.w * ib.h);
  
  return mean;
}

static double mean_red(const ImageBuffer& ib) {
  double sum = 0;
  for (int i = 0; i < ib.w * ib.h; ++i) {
      sum += red_channel(ib.data[i]);
  }
  double mean = sum / (ib.w * ib.h);
  
  return mean;
}

static double mean_blue(const ImageBuffer& ib) {
  double sum = 0;
  for (int i = 0; i < ib.w * ib.h; ++i) {
      sum += blue_channel(ib.data[i]);
  }
  double mean = sum / (ib.w * ib.h);
  
  return mean;
}

double LensCamera::focus_metric(const ImageBuffer& ib) const {

  // Part 2 Task 1: Implement this. Design a metric to judge how "in-focus"
  // the image patch stored in the provided ImageBuffer is.
  pt->bump_settings();

  double green_avg = mean_green(ib);
  double red_avg = mean_red(ib);
  double blue_avg = mean_blue(ib);

  double green_var;
  double red_var;
  double blue_var;

  for (int i = 0; i < ib.w * ib.h; i++) {
    green_var += pow(((double)green_channel(ib.data[i]) - green_avg), 2);
    blue_var += pow(((double)blue_channel(ib.data[i]) - blue_avg), 2);
    red_var += pow(((double)red_channel(ib.data[i]) - red_avg), 2);
  }

  double variance = (green_var + red_var + blue_var)/(3*ib.w*ib.h);

  return variance; //  A meaningless standin
}


void LensCamera::autofocus() {


  // Part 2 Task 2: Implement this. Design a global search using your 
  // focus metric to set the sensor to be at the depth where the 
  // render cell is most "in focus". Provided code shows how to 
  // move the sensor, request a render of the cell, and evaluate the focus metric.

  // This call ensures that your pathtracer is rendering at high enough quality.
  // Increase samples per pixel to 16 and samples per light to 16.
  pt->bump_settings();

  //double step = sqrt(36*36 + 24*24) / sqrt(screenW*screenW + screenH*screenH)*2;
  double step = sqrt(36*36 + 24*24) / sqrt(screenW*screenW + screenH*screenH)*2;

  double curr_depth = curr_lens().infinity_focus;
  double best_metric = -1.0;
  double best_depth = curr_depth;
  std::vector<double> metric;
  std::vector<double> depth;
  while (curr_depth <= curr_lens().near_focus){
    ImageBuffer ib;
    curr_lens().sensor_depth = curr_depth;
    pt->raytrace_cell(ib);
    double curr_metric = focus_metric(ib);
    metric.push_back(curr_metric);
    depth.push_back(curr_depth);
    if (curr_metric > best_metric) {
      best_metric = curr_metric;
      best_depth = curr_depth;
    }
    curr_depth += step;
  }


  // Example code. Nothing to do with your actual implementation except to 
  // demonstrate functionality.



  // ImageBuffer ib;
  // curr_lens().sensor_depth += 1;
  // pt->raytrace_cell(ib);
  // cout << "[LensCamera] The mean green is " << focus_metric(ib) << endl;


  
}





void LensCamera::dump_settings(string filename) {
  ofstream file(filename);
  file << hFov << " " << vFov << " " << ar << " " << nClip << " " << fClip << endl;
  for (int i = 0; i < 3; ++i)
    file << pos[i] << " ";
  for (int i = 0; i < 3; ++i)
    file << targetPos[i] << " ";
  file << endl;
  file << phi << " " << theta << " " << r << " " << minR << " " << maxR << endl;
  for (int i = 0; i < 9; ++i)
    file << c2w(i/3, i%3) << " ";
  file << endl;
  file << screenW << " " << screenH << " " << screenDist << endl;

  file << lens_ind << endl;
  for (Lens &lens : lenses) {
    file << lens.sensor_depth << " ";
  }
  file << endl;

  cout << "[LensCamera] Dumped settings to " << filename << endl;
}

void LensCamera::load_settings(string filename) {
  ifstream file(filename);

  file >> hFov >> vFov >> ar >> nClip >> fClip;
  for (int i = 0; i < 3; ++i)
    file >> pos[i];
  for (int i = 0; i < 3; ++i)
    file >> targetPos[i];
  file >> phi >> theta >> r >> minR >> maxR;
  for (int i = 0; i < 9; ++i)
    file >> c2w(i/3, i%3);
  file >> screenW >> screenH >> screenDist;

  file >> lens_ind;
  for (Lens &lens : lenses) {
    file >> lens.sensor_depth;
  }

  cout << "[LensCamera] Loaded settings from " << filename << endl;
}


} // namespace CGL

