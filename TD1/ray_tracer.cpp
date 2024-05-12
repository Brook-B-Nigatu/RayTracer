#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../lib/stb_image.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "../lib/stb_image_write.h"
#include <cmath>

#include <iostream>
#define PI 3.1415296
#define EPS 0.001

#include <random>
static std::default_random_engine engine (10); // random seed = 10
static std::uniform_real_distribution<double> uniform (0, 1);

#include <omp.h>
#include <string>
#include <stdio.h>
#include <algorithm>


class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};




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
Vector operator*(const Vector& a, const Vector& b){
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
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
Vector randomVect(const Vector& N){
        double r1 = uniform(engine);
        double r2 = uniform(engine);
        double x = cos(2 * PI * r1) * sqrt(1 - r2);
        double y = sin(2 * PI * r1) * sqrt(1 - r2);
        double z = sqrt(r2);

        Vector T1(N[1], -N[0], 0.);
        if (N[1] < N[0] && N[1] < N[2]){
            T1 = Vector(N[2], 0., -N[0]);
        }
        else if (N[0] < N[1] && N[0] < N[2]){
            T1 = Vector(0., N[2], -N[1]);
        }
        T1.normalize();
        Vector T2 = cross(N, T1);

        return x * T1 + y * T2 + z * N;
    }
Vector boxMuller(double stdev){
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = sqrt(-2 * log(r1)) * cos(2 * PI * r2) * stdev;
    double y = sqrt(-2 * log(r1)) * sin(2 * PI * r2) * stdev;
    return Vector(x, y, 0.);
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
Ray reflect(const Ray& incident, Vector P, const Vector& N){
    P = P + EPS * N;
    Vector reflectedDir = incident.direction - 2 * dot(incident.direction, N) * N;
    reflectedDir.normalize();
    Ray reflected(reflectedDir, P);
    return reflected;
}

struct Intersection{
    bool intersect;
    Vector P;  // point of intersection
    Vector N;  // normal to the surface at P
    void *object;
    double t; 
    Intersection(){
        object = nullptr;
        intersect = false;
        t = std::numeric_limits<double>::max();
    }
}; 


class BoundingBox {
    public:
        Vector B_min;
        Vector B_max;

        BoundingBox(){
            B_min = Vector(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
            B_max = Vector(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
        }
        BoundingBox(Vector B_min, Vector B_max) : B_min{B_min}, B_max{B_max}{}

        int getLongestAxis(){
            Vector diff = B_max - B_min;
            if (diff[0] <= diff[1] && diff[0] <= diff[2]){
                return 0;
            }
            else if (diff[1] <= diff[0] && diff[1] <= diff[2]){
                return 1;
            }
            return 2;
        }
        
        bool intersect(const Ray& ray){
            bool originInBox = true;
            
            for (int i = 0; i < 3; ++i){
                if (ray.origin[i] < B_min[i] || ray.origin[i] > B_max[i]){
                    originInBox = false;
                    break;
                }
            }
            if (originInBox){
                return true;
            }
            if (ray.direction[0] == 0 || ray.direction[1] == 0 || ray.direction[2] == 0){
                return false;
            }

            double ts[6];
            for (int i = 0; i < 6; ++i){
                int j = i / 2;
                if (i % 2 == 0){
                    ts[i] = (B_min[j] - ray.origin[j]) / ray.direction[j];
                }
                else{
                    ts[i] = (B_max[j] - ray.origin[j]) / ray.direction[j];
                }
            }

            double mint1, maxt0;

            if (ts[0] < ts[1]){
                mint1 = ts[1];
                maxt0 = ts[0];
            }
            else{
                maxt0 = ts[1];
                mint1 = ts[0];
            }
            if (ts[2] < ts[3]){
                mint1 = std::min(mint1, ts[3]);
                maxt0 = std::max(maxt0, ts[2]);
            }
            else{
                maxt0 = std::max(maxt0, ts[3]);
                mint1 = std::min(mint1, ts[2]);
            }
            if (ts[4] < ts[5]){
                mint1 = std::min(mint1, ts[5]);
                maxt0 = std::max(maxt0, ts[4]);
            }
            else{
                maxt0 = std::max(maxt0, ts[5]);
                mint1 = std::min(mint1, ts[4]);
            }
            
            return mint1 > maxt0;

        }


};

class BVH {
    public:
        BoundingBox bbox;
        BVH *left;
        BVH *rightss;
};

class Geometry {
    public:
        double refractionIndex; 
        Vector albedo;
        bool isMirror, isTransparent, invertNormal;
        virtual Intersection intersect(const Ray& ray) = 0; 
        
};


class TriangleMesh : Geometry{
public:
    std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    ~TriangleMesh() {}
	TriangleMesh() {
        this->albedo = Vector(0.3, 0.2, 0.25);
        this->isMirror = false;
    };

    void translateMesh(Vector translationVector){
        for (int i = 0; i < vertices.size(); ++i){
            vertices[i] = vertices[i] + translationVector;
        }
    }

    void scaleMesh(double scale){
        for (int i = 0; i < vertices.size(); ++i){
            vertices[i] = vertices[i] * scale;
        }
    }


    void computeBoundingBox(){
        Vector &B_max = bbox.B_max;
        Vector &B_min = bbox.B_min;
        for(const Vector &vertex : vertices){

            B_max[0] = std::max(B_max[0], vertex[0]);
            B_max[1] = std::max(B_max[1], vertex[1]);
            B_max[2] = std::max(B_max[2], vertex[2]);

            B_min[0] = std::min(B_min[0], vertex[0]);
            B_min[1] = std::min(B_min[1], vertex[1]);
            B_min[2] = std::min(B_min[2], vertex[2]);
        }
    }
	
    Intersection intersect(const Ray& ray) override {
        Intersection info;
        if (!bbox.intersect(ray)){
            return info;
        }
        double alpha, beta, gamma, t = std::numeric_limits<double>::max(), alphaTemp, betaTemp, gammaTemp, tTemp;
        size_t intersectionIndex;
        
        for (size_t i = 0; i < indices.size(); ++i){
            if (mollerTrumbore(ray, indices[i], alphaTemp, betaTemp, gammaTemp, tTemp) && tTemp < t){
                alpha = alphaTemp;
                beta = betaTemp;
                gamma = gammaTemp;
                t = tTemp;
                intersectionIndex = i;
                info.intersect = true;
                info.object = (void *) this;
            }
        }
        if (!info.intersect){
            return info;
        }
        info.t = t; 
        TriangleIndices &ti = indices[intersectionIndex];
        info.P = alpha * vertices[ti.vtxi] + beta * vertices[ti.vtxj] + gamma * vertices[ti.vtxk];
        info.N = alpha * normals[ti.ni] + beta * normals[ti.nj] + gamma * normals[ti.nk];
        info.N.normalize();
        
        return info;

    }

    bool mollerTrumbore(const Ray &ray, const TriangleIndices& inds, double &alpha, double &beta, double &gamma, double &t){
        Vector A = vertices[inds.vtxi], B = vertices[inds.vtxj], C = vertices[inds.vtxk]; 
        
        Vector e1 = B - A, e2 = C - A;
        
        Vector N = cross(e1, e2), O = ray.origin, u = ray.direction;
        
        double udN = dot(u, N);
        if (udN == 0)
            return false;
        Vector OAcru = cross(A - O, u);
        
        beta = dot(e2, OAcru) / udN, gamma = - dot(e1, OAcru) / udN, t = dot(A - O, N) / udN;

        alpha = 1 - beta - gamma;

        return t > 0 && alpha <= 1 && alpha >= 0 && beta <= 1 && beta >= 0 && gamma <= 1 && gamma >= 0;

    }
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}


	
	
};

class Sphere : Geometry{
public:
    Vector center;
    double radius;

    

    Sphere(Vector center, double radius, Vector albedo, bool isMirror) : center{center}, radius{radius}{
        this->albedo = albedo;
        this->isMirror = isMirror;
    }

    Intersection intersect(const Ray &ray) override {
        Vector P, N;
        double t;
        Intersection info;
        Vector temp = ray.origin - this->center;
        double delta = square(dot(ray.direction, temp)) - (temp.norm2() - square(this->radius));
        if (delta < 0){
            return info;
        }
        double sqrt_delta = sqrt(delta);
        t = dot(ray.direction, temp * -1) - sqrt_delta;
        if (t < 0){
            t += 2 * sqrt_delta;
        }
        if (t < 0 ){
            return info;
        }
        P = ray.origin + t * ray.direction;
        N = getNormal(P);
        info.P = P;
        info.N = N;
        info.t = t;
        info.object = (void *)this;
        info.intersect = true;
        return info;
    }
    Vector getNormal(const Vector& P){
        Vector N = P - this->center;
        N.normalize();
        return N;
    }

};





class Scene {
public:
    Vector lightSource;
    double lightIntensity;
    std::vector<Geometry*> objects;
    int recursionDepth;
    Scene(double lightIntensity, Vector lightSource, int recursionDepth = 5){
        this->lightSource = lightSource;
        this->lightIntensity = lightIntensity;
        this->recursionDepth = recursionDepth;
    }
    void addObject(Geometry* S){
        objects.push_back(S);
    }
    Intersection intersect(const Ray& ray){
        Intersection info;
        for (size_t i = 0; i < objects.size(); ++i){
            Intersection infoTemp = objects[i]->intersect(ray);
            if (infoTemp.intersect && infoTemp.t < info.t){
                info = infoTemp;
            }
        }
        return info;
    }
    Vector directLighting(const Intersection& info){
        Vector P = info.P, N = info.N;
        P = P + EPS * N;

        Vector rayDir = lightSource - P;
        rayDir.normalize();

        Ray PtoL(rayDir, P);
        Intersection info2 = intersect(PtoL);

        if (info2.intersect && (lightSource - P).norm2() > (info2.P - P).norm2()){
            return Vector(0., 0., 0.);
        }
        else{
            Vector tmp = lightSource - P; 
            return lightIntensity / (4 * PI * tmp.norm2() * PI) * std::max(0., dot(N, (1 / tmp.norm()) * tmp)) * ((Geometry *)info.object)->albedo;
        }

        
    }
    Vector getColor(const Ray& ray){
        return getColorRec(ray, recursionDepth);
    }
    Vector getColorRec(const Ray& ray, int recDepth){
        if (recDepth == 0){
            return Vector(0., 0., 0.);
        }
        Intersection info = intersect(ray);
        if (info.intersect){
            if(((Geometry *)info.object)->isMirror){
                return getColorRec(reflect(ray, info.P, info.N), recDepth - 1);
            }
            // if(this->objects[sphere_id]->isTransparent){
            //     double dotProd = dot(ray.direction, N);
            //     Geometry* obj = this->objects[sphere_id];
            //     double ratio = 1 / obj->refractionIndex;
            //     if (dotProd > 0){
            //         N = -1 * N;
            //         ratio = 1 / ratio;
            //     }
            //     if (obj->invertNormal){
            //         ratio = 1 / ratio;
            //     }
            //     double temp = std::max(1 - std::pow(ratio, 2) * (1 - std::pow(dot(ray.direction, N), 2)), 0.);
            //     if (temp < 0){
            //         return getColor(reflect(ray, P, N), recDepth - 1);
            //         //return Vector(0., 0., 0.);
            //     }
            //     Vector tangential = ratio * (ray.direction - dot(ray.direction, N) * N);
            //     Vector normal = -sqrt(temp) * N;
            //     P = P - EPS * N;
            //     Vector transmittedDir = tangential + normal;
            //     transmittedDir.normalize();
            //     Ray transmitted(transmittedDir, P);
            //     return getColor(transmitted, recDepth - 1);
            // }
            
            Vector col = directLighting(info);     
            return col + (((Geometry *)info.object)->albedo * getColorRec(Ray(randomVect(info.N), info.P + EPS * info.N), recDepth - 1));
        }
        return Vector(0., 0., 0.);
    }

};
 

 
int main() {
    int W = 512;
    int H = 512;
    Vector camera_pos = Vector(0., 0., 55.); 
    double angle = PI / 3.;
    Scene scene = Scene(3E10, Vector(-10., 20., 40.), 5);

    // Sphere diffuseSphere = Sphere(Vector(0., 0., 0.), 10., Vector(0., 0., 0.5), false);
    // scene.addObject((Geometry*)&diffuseSphere);

    TriangleMesh cat;
    cat.readOBJ("../CSE306/objs/cat.obj");
    cat.translateMesh(Vector(0., -17., 0.));
    cat.scaleMesh(0.6);
    cat.computeBoundingBox();
    scene.addObject((Geometry *) &cat);

    

    Sphere leftWall = Sphere(Vector(-1000., 0., 0.), 940., Vector(0.5, 0.5, 1.), false);
    scene.addObject((Geometry *) &leftWall);

    Sphere rightWall = Sphere(Vector(1000., 0., 0.), 940., Vector(1., 1., 0.), false);
    scene.addObject((Geometry *) &rightWall);

    Sphere topWall = Sphere(Vector(0., 1000., 0.), 940., Vector(1., 0., 0.), false);
    scene.addObject((Geometry *) &topWall);

    Sphere bottomWall = Sphere(Vector(0., -1000., 0.), 990., Vector(0., 0., 1.), false);
    scene.addObject((Geometry *) &bottomWall);

    Sphere frontWall = Sphere(Vector(0., 0., -1000.), 940., Vector(0., 1., 0.), false);
    scene.addObject((Geometry *) &frontWall);

    Sphere backWall = Sphere(Vector(0., 0., 1000.), 940., Vector(1., 0.5, 0.5), false);
    scene.addObject((Geometry *) &backWall);


    // scene.addObject(Sphere(Vector(0., 0., 0.), 10., 1.5));
    // scene.addObject(Sphere(Vector(20., 0., 0.), 9., 1.5, true));
    // scene.addObject(Sphere(Vector(20., 0., 0.), 10., 1.5));
    
    std::vector<unsigned char> image(W * H * 3, 0);
    
    int ray_count = 5;
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            
            Vector col(0., 0., 0.);
            for (int k = 0; k < ray_count; ++k){
                Vector ray_dir = Vector(j - W / 2 + 0.5, H / 2 - i - 0.5, - W / (2 * tan(angle / 2))) + boxMuller(1);
                ray_dir.normalize();
                Ray ray = Ray(ray_dir, camera_pos);
                col = col + scene.getColor(ray);
            }
            col = col / ray_count;
            image[(i * W + j) * 3 + 0] = (unsigned char) std::min(255., std::pow(col[0], 1 / 2.2));
            image[(i * W + j) * 3 + 1] = (unsigned char) std::min(255., std::pow(col[1], 1 / 2.2));
            image[(i * W + j) * 3 + 2] = (unsigned char) std::min(255., std::pow(col[2], 1 / 2.2));
            
            
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    
 
    return 0;
}
