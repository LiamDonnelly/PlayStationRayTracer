#include <stdio.h>
#include <stdlib.h>
#include <scebase.h>
#include <kernel.h>
#include <gnmx.h>
#include <video_out.h>

#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

#include <algorithm>
#include <sstream>
#include <string.h>

// Time precision
#include <chrono>

// Threading

#include <tuple>
#include "allocator.h"
#include <ult.h>

static const size_t kOnionMemorySize = 64 * 1024 * 1024;

static const int num_threads = 10;

using namespace sce;
using namespace sce::Gnmx;

template<typename T>
class Vec3
{
public:
	T x, y, z;
	Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
	Vec3(T xx) : x(xx), y(xx), z(xx) {}
	Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
	Vec3& normalize()
	{
		T nor2 = length2();
		if (nor2 > 0) {
			T invNor = 1 / sqrt(nor2);
			x *= invNor, y *= invNor, z *= invNor;
		}
		return *this;
	}
	Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
	Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
	T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
	Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
	Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
	Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
	Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
	Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
	T length2() const { return x * x + y * y + z * z; }
	T length() const { return sqrt(length2()); }
	friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
	{
		os << "[" << v.x << " " << v.y << " " << v.z << "]";
		return os;
	}
};

typedef Vec3<float> Vec3f;

class Sphere
{
public:
	Vec3f center;                           /// position of the sphere
	float radius, radius2;                  /// sphere radius and radius^2
	Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
	float transparency, reflection;         /// surface transparency and reflectivity
	Sphere(
		const Vec3f &c,
		const float &r,
		const Vec3f &sc,
		const float &refl = 0,
		const float &transp = 0,
		const Vec3f &ec = 0) :
		center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
		transparency(transp), reflection(refl)
	{ /* empty */
	}
	//[comment]
	// Compute a ray-sphere intersection using the geometric solution
	//[/comment]
	bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
	{
		Vec3f l = center - rayorig;
		float tca = l.dot(raydir);
		if (tca < 0) return false;
		float d2 = l.dot(l) - tca * tca;
		if (d2 > radius2) return false;
		float thc = sqrt(radius2 - d2);
		t0 = tca - thc;
		t1 = tca + thc;

		return true;
	}
};

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
class Rectangle {
	public:
	float x,y,z,width,length,height;
	float centerX,centerY,centerZ;
	Rectangle(){}
	Rectangle(int xx, int yy, int zz, float ww, float hh, float ll) : x(xx), y(yy), z(zz), width(ww), height(hh), length(ll) {
		 centerX = x + (width / 2);
		 centerY = y + (height / 2);
		 centerZ = z + (length / 2);
	}
};

class OctTree {
	private:
	const int MAX_LEVELS = 5;
	const int MAX_SPHERES = 10;
	int p_level;
	std::vector<Sphere> p_objects;
	Rectangle p_bounds;
	OctTree* p_nodes[8];

	public:
	OctTree(int level, Rectangle bounds)
	{
		p_level = level;
		p_bounds = bounds;

	}

	void Clear()
	{
		for (int i = 0; i < 8; i++)
		{
			p_nodes[i] = nullptr;
		}
	}

	void Split()
	{
		float nw = p_bounds.width / 2;
		float nh = p_bounds.height / 2;
		float nl = p_bounds.length / 2;
		int x = p_bounds.x;
		int y = p_bounds.y;
		int z = p_bounds.z;

		p_nodes[0] = new OctTree(p_level +1, Rectangle ( x	 , y   , z, nw, nh, nl));
		p_nodes[1] = new OctTree(p_level + 1, Rectangle(x + nw, y, z, nw, nh, nl));
		p_nodes[2] = new OctTree(p_level + 1, Rectangle(x , y, z + nl, nw, nh, nl));
		p_nodes[3] = new OctTree(p_level + 1, Rectangle(x + nw, y, z + nl, nw, nh, nl));
		p_nodes[4] = new OctTree(p_level + 1, Rectangle( x	 , y+nh, z, nw, nh, nl));
		p_nodes[5] = new OctTree(p_level + 1, Rectangle( x+nw, y+nh, z, nw, nh, nl));
		p_nodes[6] = new OctTree(p_level + 1, Rectangle(x, y + nh, z + nl, nw, nh, nl));
		p_nodes[7] = new OctTree(p_level + 1, Rectangle(x + nw, y + nh, z + nl, nw, nh, nl));
	}

	int GetIndex(Sphere object)
	{
		float dimension = object.radius + object.radius;
		Rectangle objRec = Rectangle(object.center.x - object.radius, object.center.y - object.radius, object.center.z - object.radius, dimension, dimension, dimension);
		int index = -1;

		bool fitTopY = objRec.y > p_bounds.centerY && (objRec.y + objRec.height) < (p_bounds.y + p_bounds.length);
		bool fitBotY = objRec.y > p_bounds.y && (objRec.y + objRec.height) < p_bounds.centerY;

		bool fitRightX = objRec.x > p_bounds.centerX && (objRec.x + objRec.width) < (p_bounds.x + p_bounds.width);
		bool fitLeftX = objRec.x > p_bounds.x && (objRec.x + objRec.width) < p_bounds.centerX;

		bool fitBackZ = objRec.z > p_bounds.centerZ && (objRec.z + objRec.length) < (p_bounds.z + p_bounds.length);
		bool fitFrontZ = objRec.z > p_bounds.z && (objRec.z + objRec.length) < p_bounds.centerZ;

		if (fitBotY)
		{
			//either 0,1,2,3
			if (fitRightX)
			{
				//either 1,3
				if (fitBackZ)
				{
					return 3;
				}
				else if (fitFrontZ)
				{
					return 1;
				}

			}
			else if(fitLeftX)
			{
				//either 0,2
				if (fitBackZ)
				{
					return 2;
				}
				else if (fitFrontZ)
				{
					return 0;
				}
			}
		}
		else if(fitTopY)
		{
			//either 4,5,6,7
			if (fitRightX)
			{
				//either 5,7
				if (fitBackZ)
				{
					return 7;
				}
				else if (fitFrontZ)
				{
					return 5;
				}

			}
			else if (fitLeftX)
			{
				//either 4,6
				if (fitBackZ)
				{
					return 6;
				}
				else if (fitFrontZ)
				{
					return 4;
				}
			}
		}
		
		return index;
	}

	void insert(Sphere obj)
	{
		if (p_nodes[0] != nullptr)
		{
			int index = GetIndex(obj);

			if (index != -1)
			{
				p_nodes[index]->insert(obj);
				return;
			}
		}

		p_objects.push_back(obj);
		if (p_objects.size() > MAX_SPHERES && p_level < MAX_LEVELS)
		{
			if (p_nodes[0] == nullptr)
			{
				Split();
			}
		}

		int i = 0;
		while (i < p_objects.size())
		{
			int index = GetIndex(p_objects[i]);
			if (index != -1)
			{
				p_nodes[index]->insert(p_objects[i]);
				p_objects.erase(p_objects.begin() + i);
			}
			else
			{
				i++;
			}
		}
	}

	std::vector<Sphere> Intersectable(OctTree* tree, const Vec3f &rayorig, const Vec3f &raydir)
	{
		std::vector<Sphere> spheres;
		
		if (I(Vec3f(tree->p_bounds.centerX, tree->p_bounds.centerY, tree->p_bounds.centerZ), rayorig, raydir))
		{
			spheres.insert(spheres.end(), p_objects.begin(), p_objects.end());
			if (tree->p_nodes[0] != nullptr)
			{
				for (int i = 0; i < 8; i++)
				{
					std::vector<Sphere> moreSpheres = Intersectable(p_nodes[i], rayorig, raydir);
					spheres.insert(spheres.end(), moreSpheres.begin(), moreSpheres.end());
				}
			}
		}
 
		return spheres;
	}

	bool I(const Vec3f& center, const Vec3f &rayorig, const Vec3f &raydir)
	{

		Vec3f o = center - rayorig;
		float de = o.dot(raydir);
		if (de > 0.0001f)
		{
			float de2 = o.dot(o) - de*de;
			if (de2 >= 0)
			{
				return true;
			}
		}
		return false;
	}
};
//[comment]
// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]
#define MAX_RAY_DEPTH 5

float mix(const float &a, const float &b, const float &mix)
{
	return b * mix + a * (1 - mix);
}

Vec3f trace(
	const Vec3f &rayorig,
	const Vec3f &raydir,
	const std::vector<Sphere> &spheres,
	const int &depth)
{
	//if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
	float tnear = INFINITY;
	const Sphere* sphere = NULL;
	// find intersection of this ray with the sphere in the scene
	for (unsigned i = 0; i < spheres.size(); ++i) {
		float t0 = INFINITY, t1 = INFINITY;
		if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
			if (t0 < 0) t0 = t1;
			if (t0 < tnear) {
				tnear = t0;
				sphere = &spheres[i];
			}
		}
	}
	// if there's no intersection return black or background color
	if (!sphere) return Vec3f(2);
	Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
	Vec3f phit = rayorig + raydir * tnear; // point of intersection
	Vec3f nhit = phit - sphere->center; // normal at the intersection point
	nhit.normalize(); // normalize normal direction
					  // If the normal and the view direction are not opposite to each other
					  // reverse the normal direction. That also means we are inside the sphere so set
					  // the inside bool to true. Finally reverse the sign of IdotN which we want
					  // positive.
	float bias = 1e-4; // add some bias to the point from which we will be tracing
	bool inside = false;
	if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
	if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
		float facingratio = -raydir.dot(nhit);
		// change the mix value to tweak the effect
		float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
		// compute reflection direction (not need to normalize because all vectors
		// are already normalized)
		Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
		refldir.normalize();
		Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
		Vec3f refraction = 0;
		// if the sphere is also transparent compute refraction ray (transmission)
		if (sphere->transparency) {
			float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
			float cosi = -nhit.dot(raydir);
			float k = 1 - eta * eta * (1 - cosi * cosi);
			Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
			refrdir.normalize();
			refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
		}
		// the result is a mix of reflection and refraction (if the sphere is transparent)
		surfaceColor = (
			reflection * fresneleffect +
			refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
	}
	else {
		// it's a diffuse object, no need to raytrace any further
		for (unsigned i = 0; i < spheres.size(); ++i) {
			if (spheres[i].emissionColor.x > 0) {
				// this is a light
				Vec3f transmission = 1;
				Vec3f lightDirection = spheres[i].center - phit;
				lightDirection.normalize();
				for (unsigned j = 0; j < spheres.size(); ++j) {
					if (i != j) {
						float t0, t1;
						if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
							transmission = 0;
							break;
						}
					}
				}
				surfaceColor += sphere->surfaceColor * transmission *
					std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColor;
			}
		}
	}

	return surfaceColor + sphere->emissionColor;
}

const int NUMTHREADS = 6;
static void *ThreadedRender(void* sphere)
{
	std::tuple<std::vector<Sphere>, Vec3f*, int>* a = static_cast<std::tuple<std::vector<Sphere>, Vec3f*, int>*>(sphere);

	unsigned width = 1920, height = 1080;
	int numThreads = NUMTHREADS;
	int endheight = std::get<2>(*a) + (height / numThreads);
	float invWidth = 1 / float(width), invHeight = 1 / float(height);
	float fov = 30, aspectratio = width / float(height);
	float angle = tan(M_PI * 0.5 * fov / 180.);

	int interation = std::get<2>(*a) * width;
	for (unsigned y = std::get<2>(*a); y < endheight; ++y) {
		for (unsigned x = 0; x < width; ++x) {
			float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
			float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
			Vec3f raydir(xx, yy, -1);
			raydir.normalize();
			std::get<1>(*a)[interation] = trace(Vec3f(0), raydir, std::get<0>(*a), 0);
			interation++;
		}
	}

	return SCE_OK;
}
//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.
//[/comment]
void render(const std::vector<Sphere> &spheres)
{



	// Initialize the WB_ONION memory allocator

	LinearAllocator onionAllocator;
	int ret = onionAllocator.initialize(
		kOnionMemorySize, SCE_KERNEL_WB_ONION,
		SCE_KERNEL_PROT_CPU_RW | SCE_KERNEL_PROT_GPU_ALL);

	//if (ret != SCE_OK)
	//	return ret;

	unsigned width = 1920, height = 1080;
	size_t totalSize = sizeof(Vec3f)* width * height;

	void * buffer = onionAllocator.allocate(totalSize, Gnm::kAlignmentOfBufferInBytes);

	Vec3f *image = reinterpret_cast<Vec3f *>(buffer);
	Vec3f *pixel = image;

	int numThreads = NUMTHREADS;
	int rowHeight = height / numThreads;
	// Trace rays
	
	int count = 0;
	ScePthread threads[NUMTHREADS];
	ScePthreadAttr attribute;
	SceKernelSchedParam param;
	SceKernelCpumask mask;
	
	param.sched_priority = SCE_KERNEL_PRIO_FIFO_HIGHEST;

	int RESULT = 0;

	//Init Attribute
	RESULT = scePthreadAttrInit(&attribute);

	//Set scheduling to FIFO
	RESULT = scePthreadAttrSetschedpolicy(&attribute, 1);

	//Set inheritance to explict
	RESULT = scePthreadAttrSetinheritsched(&attribute, SCE_PTHREAD_EXPLICIT_SCHED);

	//Set prioty to highest
	RESULT = scePthreadAttrSetschedparam(&attribute, &param);
	
	
	//mask = 100;
	//RESULT = sceKernelSetFsstParam(255,mask);

	int cpumode = sceKernelGetCpumode();


	for (int i = 0; i < numThreads; i++)
	{ 
		switch (i)
		{
			case(0) :
			{
				mask = 001;
				break;
			}
			case(1):
			{
				mask = 002;
				break;
			}
			case(2):
			{
				mask = 004;
				break;
			}
			case(3):
			{
				mask = 010;
				break;
			}
			case(4):
			{
				mask = 020;
				break;
			}
			case(5):
			{
				mask = 040;
				break;
			}
		}
		//RESULT = scePthreadAttrSetaffinity(&attribute, mask);

		void* pointer = new std::tuple<std::vector<Sphere>, Vec3f*, int>(spheres, pixel, count);
		int value = scePthreadCreate(&threads[i], &attribute, ThreadedRender, pointer, "Tracing");
		count += rowHeight;
	}

	for (int i = 0; i < numThreads; i++)
	{
		scePthreadJoin(threads[i], NULL);
	}

	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::stringstream ss;
	ss << "/app0/spheres"  << ".ppm";
	std::string tempString = ss.str();
	char* filename = (char*)tempString.c_str();

	std::ofstream ofs(filename, std::ios::out | std::ios::binary);
	ofs << "P6\n" << width << " " << height << "\n255\n";
	for (unsigned i = 0; i < width * height; ++i) {
	ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
	(unsigned char)(std::min(float(1), image[i].y) * 255) <<
	(unsigned char)(std::min(float(1), image[i].z) * 255);
	}
	ofs.close();
	//delete[] image;
}



void BasicRender()
{
	std::vector<Sphere> spheres;

	spheres.push_back(Sphere(Vec3f(0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
	spheres.push_back(Sphere(Vec3f(1, 1, -20), 1, Vec3f(1.00, 0.32, 0.36), 1, 0.5)); // Radius++ change here
	spheres.push_back(Sphere(Vec3f(3.0, 2, -15), 2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
	spheres.push_back(Sphere(Vec3f(5.0, 3, -25), 3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));

	OctTree tree =  OctTree(0, Rectangle(-15,-15,0,30,30,30));
	tree.Clear();
	tree.Split();
	for (int i = 0; i < 4; i++)
	{
		tree.insert(spheres[i]);
	}

	render(spheres);


}
int main(int argc, char **argv)
{
	auto start = std::chrono::high_resolution_clock::now();
	BasicRender();
	auto end = std::chrono::high_resolution_clock::now();

	//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	//std::cout << "Total Time: " << duration;
	return 0;
}