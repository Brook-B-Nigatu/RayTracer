#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <cmath>

#include "iostream"
#define PI 3.1415296

double square(double x){
    return x * x;
}
class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
    Vector direction;
    Vector origin;
    Ray(Vector direction, Vector origin){
        this->direction = direction;
        this->origin = origin;
    }
};
 
class Sphere {
public:
    Vector center, albedo;
    double radius;
    Sphere(Vector center, Vector albedo, double radius){
        this->center = center;
        this->albedo = albedo;
        this->radius = radius;
    }
    bool intersect(Ray &ray, Vector& P, Vector& N, double& t){
        Vector temp = ray.origin - this->center;
        double delta = square(dot(ray.direction, temp)) - (temp.norm2() - square(this->radius));
        if (delta < 0){
            return false;
        }
        double sqrt_delta = sqrt(delta);
        t = dot(ray.direction, temp * -1) - sqrt_delta;
        if (t < 0){
            t += 2 * sqrt_delta;
        }
        if (t < 0 ){
            return false;
        }
        P = ray.origin + t * ray.direction;
        N = P - this->center;
        N.normalize();
        return true;
    }
    Vector computeColor(double lightIntensity, const Vector& lightSource, const Vector& P, const Vector& N){
        Vector tmp = lightSource - P; 
        return lightIntensity / (4 * PI * tmp.norm2() * PI) * std::max(0., dot(N, (1 / tmp.norm()) * tmp)) * this->albedo;

    }




};

class Scene {
public:
    Vector lightSource;
    double lightIntensity;
    std::vector<Sphere> objects;
    double epsilon = 0.001;
    Scene(double lightIntensity, Vector lightSource){
        this->lightSource = lightSource;
        this->lightIntensity = lightIntensity;
    }
    void addObject(Sphere S){
        objects.push_back(S);
    }
    bool intersect(Ray& ray, Vector& P, Vector& N, size_t& sphere_id){
        Vector Ptemp, Ntemp;
        double t = INT64_MAX, t_temp;
        bool intersected = false; 
        for (size_t i = 0; i < size(objects); ++i){
            if (objects[i].intersect(ray, Ptemp, Ntemp, t_temp)){
                intersected = true;
                if (t_temp < t){
                    t = t_temp;
                    P = Ptemp;
                    N = Ntemp;
                    sphere_id = i;
                }
            }
           
        }
        return intersected;
    }
    Vector getColor(Ray& ray, int recDepth){
        if (recDepth == 0){
            return Vector(0., 0., 0.);
        }
        Vector N, P;
        size_t sphere_id;
        if (this->intersect(ray, P, N, sphere_id)){
            P = P + epsilon * N;
            Vector P2, N2;
            size_t si;
            Vector rayDir = this->lightSource - P;
            rayDir.normalize();
            Ray rayPL(rayDir, P);
            if (this->intersect(rayPL, P2, N2, si) && (this->lightSource - P).norm2() > (P2 - P).norm2()){
                return Vector(0., 0., 0.);
            }
            return objects[sphere_id].computeColor(this->lightIntensity, this->lightSource, P, N);
        }
        return Vector(0., 0., 0.);
    }

};
 

 
int main() {
    int W = 512;
    int H = 512;
    Vector camera_pos = Vector(0., 0., 55.); 
    double angle = PI / 3.;
    Scene scene = Scene(1E10, Vector(-10., 20., 40.));
    scene.addObject(Sphere(Vector(0., 0., 0.), Vector(0.3, 0.2, 0.1), 10.));
    scene.addObject(Sphere(Vector(0., 1000., 0.), Vector(1., 0., 0.), 940.));
    scene.addObject(Sphere(Vector(0., 0., 1000.), Vector(1., 1., 1.), 940.));
    scene.addObject(Sphere(Vector(0., 0., -1000.), Vector(0., 1., 0.), 940.));
    scene.addObject(Sphere(Vector(0., -1000., 0.), Vector(0., 0., 1.), 990.));
    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector ray_dir = Vector(j - W / 2 + 0.5, H / 2 - i - 0.5, - W / (2 * tan(angle / 2)));
            ray_dir.normalize();
            Ray ray = Ray(ray_dir, camera_pos);
            Vector P, N;
            Vector col = scene.getColor(ray, 5);
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(col[0], 1 / 2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(col[1], 1 / 2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(col[2], 1 / 2.2));
            
            
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    
 
    return 0;
}
