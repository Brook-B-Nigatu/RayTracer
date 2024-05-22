
#include <random>
static std::default_random_engine engine (10); // random seed = 10
static std::uniform_real_distribution<double> uniform (0, 1);

#include<vector>
#include<string>
#include<iostream>
#include "../Assignment1/VectorClass.h"
// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
	std::vector<Vector> vertices;
};	

class Voronoi {
	public:
		std::vector<Vector> points;
		std::vector<Polygon> computeVoronoi(){
			size_t N = points.size();
			std::vector<Polygon> res(N);
			Polygon square;
			square.vertices.push_back(Vector(0., 0., 0.));
			square.vertices.push_back(Vector(1., 0., 0.));
			square.vertices.push_back(Vector(1., 1., 0.));
			square.vertices.push_back(Vector(0., 1., 0.));
		#pragma omp parallel for schedule(dynamic, 1) 
			for (size_t i = 0; i < N; ++i){
				Polygon cell = square;
				for (size_t j = 0; j < N; ++j){
					if (i == j) continue;
					cell = sutherlandHodgman(cell, points[i], points[j]);
				}
				res[i] = cell;
			}
			return res;
		}

		Polygon sutherlandHodgman(const Polygon& cell, const Vector& p1, const Vector& p2){
			Polygon res;
			size_t n = cell.vertices.size();
			Vector M = (p1 + p2) / 2; 
			for (size_t i = 0; i < n; ++i){
				Vector curr = cell.vertices[i];
				Vector prev = cell.vertices[(n + i - 1) % n];
				double denom = dot(curr - prev, p2 - p1);
				double t;
				if (denom == 0.){
					t = 1e9;
				}
				else{
					t = dot(M - prev, p2 - p1) / denom;
				}	
				Vector P = prev + t * (curr - prev);
				if (closer(curr, p1, p2)){
					if (!closer(prev, p1, p2)){
						res.vertices.push_back(P);
					}
					res.vertices.push_back(curr);
				}
				else if (closer(prev, p1, p2)){
					res.vertices.push_back(P);
				}
			}
			return res;

		}
		bool closer(const Vector &point, const Vector &p1, const Vector &p2){
			Vector p1V = point - p1;
			Vector p2V = point - p2;

			return p1V.norm2() <= p2V.norm2();
		}


};

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
	FILE* f = fopen(filename.c_str(), "w+"); 
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
	for (int i=0; i<polygons.size(); i++) {
		fprintf(f, "<g>\n");
		fprintf(f, "<polygon points = \""); 
		for (int j = 0; j < polygons[i].vertices.size(); j++) {
			fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
		}
		fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		fprintf(f, "</g>\n");
	}
	fprintf(f, "</svg>\n");
	fclose(f);
}


// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
	FILE* f;
	if (frameid == 0) {
		f = fopen(filename.c_str(), "w+");
		fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
		fprintf(f, "<g>\n");
	} else {
		f = fopen(filename.c_str(), "a+");
	}
	fprintf(f, "<g>\n");
	for (int i = 0; i < polygons.size(); i++) {
		fprintf(f, "<polygon points = \""); 
		for (int j = 0; j < polygons[i].vertices.size(); j++) {
			fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
		}
		fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
	}
	fprintf(f, "<animate\n");
	fprintf(f, "	id = \"frame%u\"\n", frameid);
	fprintf(f, "	attributeName = \"display\"\n");
	fprintf(f, "	values = \"");
	for (int j = 0; j < nbframes; j++) {
		if (frameid == j) {
			fprintf(f, "inline");
		} else {
			fprintf(f, "none");
		}
		fprintf(f, ";");
	}
	fprintf(f, "none\"\n	keyTimes = \"");
	for (int j = 0; j < nbframes; j++) {
		fprintf(f, "%2.3f", j / (double)(nbframes));
		fprintf(f, ";");
	}
	fprintf(f, "1\"\n	dur = \"5s\"\n");
	fprintf(f, "	begin = \"0s\"\n");
	fprintf(f, "	repeatCount = \"indefinite\"/>\n");
	fprintf(f, "</g>\n");
	if (frameid == nbframes - 1) {
		fprintf(f, "</g>\n");
		fprintf(f, "</svg>\n");
	}
	fclose(f);
}


int main(){
	size_t N = 1000;
	Voronoi vor;
	for (size_t i = 0; i < N; ++i){
		vor.points.push_back(Vector(uniform(engine), uniform(engine), 0.));
	}

	save_svg(vor.computeVoronoi(), "test.svg");
	return 0;
}