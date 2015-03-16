#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#define NOMINMAX
#include <windows.h>
#endif // Win32 platform

#include <GL/gl.h>
#include <GL/glu.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/
#include <GL/glut.h>

#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>
#include <algorithm>
#include "perlin.h"
#include "LightSource.h"
#include "Camera.h"
#include "Material.h"
#include "Ray.h"

// For procedural texturing
Perlin perlin;

// for quadrics, so that we do not need a float4.cpp

const float4x4 float4x4::identity(
	1, 0, 0, 0,
	0,1, 0, 0,
	0, 0,1, 0,
	0, 0, 0, 1);

//float4x4 A = float4x4::identity;

// TO BE CREATED AT PRACTICAL
class DirectionalLightSource : public LightSource
{
private:
	float3 powerDensity;
	float3 lightDir;
public:
	DirectionalLightSource::DirectionalLightSource(float3 pD, float3 d){
		powerDensity = pD;
		lightDir = d;
	}
	float3 getPowerDensityAt(float3 x){
		//powerDensity = float3(0, 0, 3);
		return powerDensity;
	}
	float3 getLightDirAt(float3 x){
		//lightDir = float3(0, 0, 1);
		return lightDir;
	}
	float getDistanceFrom(float3 x){
		//return float3(FLT_MAX, FLT_MAX, FLT_MIN);
		//float3 maxFloat = FLT_MAX;
		return FLT_MAX;

	}

};


// TO BE CREATED AT PRACTICAL
class PointLightSource : public LightSource
{
private:
	float3 center;
	float3 density;

public:
	PointLightSource::PointLightSource(const float3 c, float3 d){
		center = c;
		density = d;

	}
	float3 getPowerDensityAt(float3 x){
		//C/4(pi)* distance ^2
		//	density = float3(1, 1, 1);
		return density *(1 / (4 * 3.1459) * 1 / (getDistanceFrom(x)* getDistanceFrom(x)));
		//1 / 4 * (getPowerDensityAt(x)* getDistanceFrom(x);
	}
	float3 getLightDirAt(float3 x){
		return -(x - center).normalize();
	}
	float getDistanceFrom(float3 x){
		//center = float3(0, 0, 1);
		//return sqrt(pow(x.x-center.x, 2), pow(ce,2));
		return sqrt(pow(x.x - center.x, 2) + pow(x.y - center.y, 2) + pow(x.z - center.z, 2));
	}
};





// Hit record structure. Contains all data that describes a ray-object intersection point.
class Hit
{
public:
	Hit()
	{
		t = -1;
	}
	float t;
	float3 position;
	float3 normal;
	Material* material;
};

// Object abstract base class.
class Intersectable
{
protected:
	Material* material;
public:
	Intersectable(Material* material) :material(material) {}
	virtual Hit intersect(const Ray& ray) = 0;
};

// Object realization.
class Sphere : public Intersectable
{
	float3 center;
	float radius;
public:
	Sphere(const float3& center, float radius, Material* material) :
		Intersectable(material),
		center(center),
		radius(radius)
	{
	}

	
	Hit intersect(const Ray& ray)
	{
		float3 diff = ray.origin - center;
		double a = ray.dir.dot(ray.dir);
		double b = diff.dot(ray.dir) * 2.0;
		double c = diff.dot(diff) - radius * radius;

		double discr = b * b - 4.0 * a * c;
		if (discr < 0)
			return Hit();
		double sqrt_discr = sqrt(discr);
		double t1 = (-b + sqrt_discr) / 2.0 / a;
		double t2 = (-b - sqrt_discr) / 2.0 / a;

		float t = (t1<t2) ? t1 : t2;
		if (t < 0)
			t = (t1<t2) ? t2 : t1;
		if (t < 0)
			return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		h.normal = h.position - center;
		h.normal.normalize();

		return h;

	}
};

// TO BE CREATED AT PRACTICAL
class Plane : public Intersectable 
{
	//Plane is parralel to xAx and Zaz, so the y-intercept is at this value
	float3 yInt;
	float limit;
public:
	Plane(const float3& yInt, float limit, Material* material) :
		Intersectable(material),
		yInt(yInt),
		limit(limit)
	{
	}
	Hit  intersect(const Ray& ray)
	{
		

		float bot = yInt.dot(ray.dir);
		if (bot == 0)
			return Hit();


		float t =  -(ray.origin.dot(yInt) + limit) / bot;

		if (t < 0)
			return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir*t;
		h.normal = (bot > 0) ? -yInt : yInt;
		h.normal.normalize();

		return h;
	}
};

class Ellipsoid :public Intersectable
{
	
	float3 pos;
	
	float3 scale;
	float4x4 A;


public:
	Ellipsoid(Material* material, int i ) :
		Intersectable(material)
		
		//theta((theta *3.14)/180),
		//phi((phi * 3.14)/180)
	{

		if (i == 1){
			A = float4x4::identity;


			A._00 = 2;
			A._11 = 4;
			A._22 = 2;
			A._33 = -0.7;
		}
		else if (i == 2){
			
			A = float4x4::identity;
		

			A._00 = 2;
			A._11 = 4;
			//A._22 = 2;
			A._33 = -1.5;
		
		
		}
		else if (i == 3){
			A = float4x4::identity;

			A._00 = 2;
			A._11 = 4;
			//A._22 = 2;
			A._33 = -2;
		}
		
	}


	Ellipsoid* translation(Ellipsoid* ell, float3 offset){
		float4x4 trans = float4x4::translation(offset);
		float4x4 matrix = trans*A*trans.transpose();
		ell->A = matrix;
		return ell;

	}
	Ellipsoid* rotation(Ellipsoid* ell, float3 axis, float angle){
		float4x4 rot = float4x4::rotation(axis, angle);
		float4x4 matrix = rot*A*rot.invert();
		ell->A = matrix;
		return ell;

	}
	Ellipsoid* scaling(Ellipsoid* ell, float3 factors){
		float4x4 scal = float4x4::scaling(factors);
		float4x4 matrix = scal*A*scal.transpose();
		ell->A = matrix;
		return ell;
	}
	Hit  intersect(const Ray& ray)
	{
	

		 
		// ray in homo coords
		float4 e = float4(ray.origin.x,
			ray.origin.y, ray.origin.z, 1);
		float4 d = float4(ray.dir.x,
			ray.dir.y, ray.dir.z, 0);
		// quadratic coeffs.
		double a = d.dot(A * d);
		double b = e.dot(A * d)
			+ d.dot(A * e);
		double c = e.dot(A * e);


		double discr = b * b - 4.0 * a * c;
		if (discr < 0)
			return Hit();

		double sqrt_discr = sqrt(discr);
		double t1 = (-b + sqrt_discr) / 2.0 / a;
		double t2 = (-b - sqrt_discr) / 2.0 / a;

		float t = (t1<t2) ? t1 : t2;
		if (t < 0)
			t = (t1<t2) ? t2 : t1;
		if (t < 0)
			return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		//h.normal = h.position - center;
		// homo position
		float4 hPos = float4(h.position.x,
			h.position.y, h.position.z, 1);
		// homo normal per quadric normal formula
		float4 hNormal = A * hPos + hPos * A;
		// Cartesian normal
		h.normal = float3(hNormal.x, hNormal.y, hNormal.z).normalize();

		h.normal.normalize();

	//	if (abs(h.position-center) )
		return h;
	}
};


class Cylinder :public Intersectable
{
	float4x4 A;
	float4x4 B;
public:
	Cylinder(Material* material, int i):Intersectable(material) {
		// infinite cylinder hardwired
		if (i == 1){
			A = float4x4::identity;
			A._00 = 0;
			A._11 = 8;
			A._33 = -.5;

			
			// sphere or radius 2 hardwired
			B = float4x4::identity;
			B._33 = -1;
		}
		else if (i == 2){
			A = float4x4::identity;
			A._11= 0;
			A._33 = -.25;

			//A._31 = 2;
			// sphere or radius 2 hardwired
			B = float4x4::identity;
			B._33 = -.5;
		}
	} // add methods 


	Cylinder* translation(Cylinder* cyl, float3 offset){
		float4x4 trans = float4x4::translation(offset);
		float4x4 matrix2 = trans*A*(trans.transpose());
		float4x4 matrix = trans*B*(trans.transpose());
		cyl->A = matrix2;

		cyl->B = matrix;
		return cyl;

	}
	Cylinder* rotation(Cylinder* cyl, float3 axis, float angle){
		float4x4 rot = float4x4::rotation(axis, angle);
		float4x4 matrix = rot*A*rot.invert();
		cyl->A = matrix;
		return cyl;;

	}
	Cylinder* scaling(Cylinder* cyl, float3 factors){
		float4x4 scal = float4x4::scaling(factors);
		float4x4 matrix = scal*A*scal.transpose();
		float4x4 matrix2 = scal*B*scal.transpose();
		cyl->A = matrix;
		cyl->B = matrix2;
		return cyl;
	}
	Hit  intersect(const Ray& ray)
	{
		// ray in homo coords
		float4 e = float4(ray.origin.x,
			ray.origin.y, ray.origin.z, 1);
		float4 d = float4(ray.dir.x,
			ray.dir.y, ray.dir.z, 0);
		// quadratic coeffs.
		double a = d.dot(A * d);
		double b = e.dot(A * d)
			+ d.dot(A * e);
		double c = e.dot(A * e);



		double discr = b * b - 4.0 * a * c;
		if (discr < 0)
			return Hit();

		double sqrt_discr = sqrt(discr);
		double t1 = (-b + sqrt_discr) / 2.0 / a;
		double t2 = (-b - sqrt_discr) / 2.0 / a;

		float4 hit1 = e + d * t1;
		if (hit1.dot(B * hit1) > 0) // if not in B
			t1 = -1;				 // invalidate
		float4 hit2 = e + d * t2;
		if (hit2.dot(B * hit2) > 0) // if not in B
			t2 = -1;

		float t = (t1<t2) ? t1 : t2;
		if (t < 0)
			t = (t1<t2) ? t2 : t1;
		if (t < 0)
			return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		//h.normal = h.position - center;
		// homo position
		float4 hPos = float4(h.position.x,
			h.position.y, h.position.z, 1);
		// homo normal per quadric normal formula
		float4 hNormal = A * hPos + hPos * A;

		// Cartesian normal
		h.normal = float3(hNormal.x, hNormal.y, hNormal.z).normalize();

		h.normal.normalize();

		return h;
	}
};
class Paraboloid :public Intersectable
{

	float3 pos;

	float3 scale;
	float4x4 A;
	float4x4 B;

public:
	Paraboloid(Material* material, int i) :
		Intersectable(material)

	
	{

		if (i == 1){
			A = float4x4::identity;
			//A = float4x4::identity;
			A._00 =10;
			A._11 = 5;
			A._12 = -1; 
			A._23 = -1; 
			A._32 = -1;
			
			A._33 =-4;
			B = float4x4::identity;
			B._33 = -5;
			
		}
		if (i == 2){

			A = float4x4::identity;
			//A = float4x4::identity;
			A._00 = 30;
			A._11 = 10;
			A._12 = -1;
			A._23 = 1;
			A._32 = 1;

			A._33 = -4;
			B = float4x4::identity;
			B._33 = -5;
		}
	}
	
	Paraboloid* translation(Paraboloid* par, float3 offset){
		float4x4 trans = float4x4::translation(offset);
		float4x4 matrix = trans*A*trans.transpose();
		float4x4 matrix2 = trans*B*trans.transpose();
		par->B = matrix2;
		par->A = matrix;
		return par;

	}
	
	Paraboloid* rotation(Paraboloid* ell, float3 axis, float angle){
		float4x4 rot = float4x4::rotation(axis, angle);
		float4x4 matrix = rot*A*rot.invert();
		float4x4 matrix2 = rot*B*rot.invert();
		ell->A = matrix;
		ell->B = matrix2;
		return ell;

	}

	Paraboloid* scaling(Paraboloid* ell, float3 factors){
		float4x4 scal = float4x4::scaling(factors);
		float4x4 matrix = scal*A*scal.transpose();
		float4x4 matrix2 = scal*B*scal.transpose();
		ell->B = matrix2;
		ell->A = matrix;
		return ell;
	}

	Hit  intersect(const Ray& ray)
	{



		// ray in homo coords
		float4 e = float4(ray.origin.x,
			ray.origin.y, ray.origin.z, 1);
		float4 d = float4(ray.dir.x,
			ray.dir.y, ray.dir.z, 0);
		// quadratic coeffs.
		double a = d.dot(A * d);
		double b = e.dot(A * d)
			+ d.dot(A * e);
		double c = e.dot(A * e);


		double discr = b * b - 4.0 * a * c;
		if (discr < 0)
			return Hit();

		double sqrt_discr = sqrt(discr);
		double t1 = (-b + sqrt_discr) / 2.0 / a;
		double t2 = (-b - sqrt_discr) / 2.0 / a;
	
		float4 hit1 = e + d * t1;
		if (hit1.dot(B * hit1) > 0) // if not in B
			t1 = -1;	
	
		float4 hit2 = e + d * t2;
		if (hit2.dot(B * hit2) > 0) // if not in B
	t2 = -1;

		float t = (t1<t2) ? t1 : t2;
		if (t < 0)
			t = (t1<t2) ? t2 : t1;
		if (t < 0)
			return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		//h.normal = h.position - center;
		// homo position
		float4 hPos = float4(h.position.x,
			h.position.y, h.position.z, 1);
		// homo normal per quadric normal formula
		float4 hNormal = A * hPos + hPos * A;
		// Cartesian normal
		h.normal = float3(hNormal.x, hNormal.y, hNormal.z).normalize();

		h.normal.normalize();

		//	if (abs(h.position-center) )
		return h;
	}
};



// Object realization.
class Cone : public Intersectable
{
	float3 center;
	float radius;

	float4x4 A;
	float4x4 B;
	//Perlin perlin;
public:
	Cone(Material* material, int i) :
		Intersectable(material)

	{
		if (i == 1){
		
			
			A = float4x4::identity;

			A._00 = 2;
			A._11 = 2;
		
			A._22 = -1;
			A._33 =0;

		

			
			B = float4x4::identity;
			B._33 = -1;

	
		}
	}
	Cone* translation(Cone* cone, float3 offset){
		float4x4 trans = float4x4::translation(offset);
		float4x4 matrix2 = trans*A*trans.transpose();
	float4x4 matrix = trans*B*trans.transpose();
		cone->A = matrix2;

		cone->B = matrix;
		return cone;

	}
	Cone* rotation(Cone* cone, float3 axis, float angle){
		float4x4 rot = float4x4::rotation(axis, angle);
		float4x4 matrix = rot*A*rot.invert();
		float4x4 matrix2 = rot*B*rot.invert();
		cone->A = matrix;
		cone->B = matrix2;
		return cone;;

	}
	Cone* scaling(Cone* cone, float3 factors){
		float4x4 scal = float4x4::scaling(factors);
		float4x4 matrix = scal*A*scal.transpose();
		float4x4 matrix2 = scal*B*scal.transpose();
		cone->A = matrix;
		cone->B = matrix2;
		return cone;
	}
	Hit  intersect(const Ray& ray)
	{



		// ray in homo coords
		float4 e = float4(ray.origin.x,
			ray.origin.y, ray.origin.z, 1);
		float4 d = float4(ray.dir.x,
			ray.dir.y, ray.dir.z, 0);
		// quadratic coeffs.
		double a = d.dot(A * d);
		double b = e.dot(A * d)
			+ d.dot(A * e);
		double c = e.dot(A * e);


		double discr = b * b - 4.0 * a * c;
		if (discr < 0)
			return Hit();

		double sqrt_discr = sqrt(discr);
		double t1 = (-b + sqrt_discr) / 2.0 / a;
		double t2 = (-b - sqrt_discr) / 2.0 / a;
		
	float4 hit1 = e + d* t1;
		if (hit1.dot(B * hit1) > 0) // if not in B
		t1 =-1;	
		// invalidate
		float4 hit2 = e + d * t2;
		if (hit2.dot(B * hit2) >0) // if not in B
			t2 = -1;
	
		float t = (t1<t2) ? t1 : t2;
		if (t < 0)
			t = (t1<t2) ? t2 : t1;
		if (t < 0)
			return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		//h.normal = h.position - center;
		// homo position
		float4 hPos = float4(h.position.x,
			h.position.y, h.position.z, 1);
		// homo normal per quadric normal formula

		float4 hNormal = A * hPos + hPos * A;
		// Cartesian normal
		h.normal = float3(hNormal.x, hNormal.y, hNormal.z).normalize();

		h.normal.normalize();

	
		return h;
	}
};
// TO BE CREATED AT PRACTICAL
//class Quadric : public Intersectable
// {};

class Scene
{
	Camera camera;
	std::vector<LightSource*> lightSources;
	std::vector<Intersectable*> objects;
	std::vector<Material*> materials;
public:
	Scene()
	{
		
		for (int i = 0; i < 15; i++){
			materials.push_back(new Material());
		}
		perlin.noise(float3(1, 1, 0));
		//Directional(float3 powerDensity, float3 Direction)
		lightSources.push_back(new DirectionalLightSource(float3(2, 3,3), float3(1,1, 1)));

		//PointLightSource(center,density)
		lightSources.push_back(new PointLightSource(float3(0, -1.3, 4), float3(50,100, 50)));

		/***************************Plane*****************************/
		
	//	objects.push_back(new Plane(float3(0, 200, 2), 150, materials[6]));
		//materials[6]->kd = float3(1, 0, 1);
		
		/************************************************************/
		/**/
		Ellipsoid* sphere = new Ellipsoid(materials[0],1);
		//materials[0]->reflective = true;
		//materials[0]->minReflectance = float3(1, 1, 1);
		materials[0]->kd = float3(1, 0, 0);
		materials[0]->ks = float3(0.6, 0.5, -0.3);
		sphere->translation(sphere, float3(0, -.75,0));
		objects.push_back(sphere);
	
		/***************************Snowman*****************************/
		/**
		Ellipsoid* ell = new Ellipsoid(materials[1], 1);
		materials[1]->refractive = true;
		materials.at(1)->minReflectance = float3(0.5, 0.5, 0.5);
		materials.at(1)->refractiveIndex = 2;
		materials[1]->kd = float3(1, 0, 0);
		ell->translation(ell, float3(-1, -1, -1.0));
		//ell->translation(ell, float3(-1, -.80, -1.0));

		//ell->scaling(ell, float3(1.25, 1.25, 1.25));
		objects.push_back(ell);
		**/
		
		Ellipsoid* ell = new Ellipsoid(materials[0], 1);
		materials[0]->kd = float3(1, 1, 1);
		ell->translation(ell, float3(0, -1.0, 0.2));
		//ell->scaling(ell, float3(1.25, 1.25, 1.25));
		objects.push_back(ell);
		
		
		Ellipsoid* ell2 = new Ellipsoid(materials[1], 2);
		materials[1]->kd = float3(1, 1, 1);
		ell2->translation(ell2, float3(0, 0, 0.2));
		//ell->scaling(ell, float3(1.25, 1.25, 1.25));
		objects.push_back(ell2);
		
		Ellipsoid* ell3 = new Ellipsoid(materials[2], 3);
		materials[2]->kd = float3(1, 1, 1);
		ell3->translation(ell3, float3(0, 1.2, 0.5));
		objects.push_back(ell3);

	
		/************************************************************/
		

		/***********************Hat******************************/
	
		Cylinder* cyl1 = new Cylinder(materials[3], 1);
		materials[3]->reflective = true;
		materials[3]->minReflectance = float3(1, 1, 1);
		materials[3]->kd = float3(1, 0, 0);
		materials[3]->ks = float3(0.6, 0.5, -0.3);
		cyl1->translation(cyl1, float3(0, -1.5, 0.5));
		objects.push_back(cyl1);
		
		Cylinder* cyl2 = new Cylinder(materials[4], 2);
		materials[4]->reflective = true;
		materials[4]->minReflectance = float3(1, 1, 1);
		materials[4]->kd = float3(1, 0, 0);
		materials[4]->ks = float3(0.6, 0.5, -0.3);
		cyl2->scaling(cyl2, float3(1, 1, 2));
		cyl2->translation(cyl2, float3(0, -2.4, 1));
			objects.push_back(cyl2); 
		
		
		/*****************************************************/
		
		/************************Cone****************************/
		
		Cone* cone = new Cone(materials[5], 1);
		materials[5]->kd = float3(1, 0.3, 0);
		cone->translation(cone, float3(0, 3, 0));
	cone->rotation(cone, float3(-1, 0, 0),180);
		cone->scaling(cone, float3(4,4,4));
		objects.push_back(cone);
		
	
		/*****************************************************/
		
		/***********************Paraboloid**********************/
		
		Paraboloid* par = new Paraboloid(materials[7], 1);
		materials[7]->refractive = true;
		materials.at(7)->minReflectance = float3(0.5, 0.5, 0.5);
		materials.at(7)->refractiveIndex = 1.33;
		materials[7]->kd = float3(0.3, 0, 0.5);
		par->rotation(par, float3(-1, 0, 0),90);
		par->translation(par, float3(-5, -1, 5));
		par->scaling(par, float3(2,2,2));
		objects.push_back(par);
	
		Paraboloid* par2 = new Paraboloid(materials[8], 2);
		materials[8]->refractive = true;
		materials.at(8)->minReflectance = float3(0.5, 0.5, 0.5);
		materials.at(8)->refractiveIndex = 1.33;
		materials[8]->kd = float3(0, 0.6, 0.9);
		par2->rotation(par2, float3(-1, 0, 0), 90);
		par2->translation(par2, float3(2.5, -5,0));
		par2->scaling(par2, float3(2, 2, 2));
		objects.push_back(par2);

		Paraboloid* par3 = new Paraboloid(materials[9], 2);
		materials[9]->refractive = true;
		materials.at(9)->minReflectance = float3(0.5, 0.5, 0.5);
		materials.at(9)->refractiveIndex = 1.33;
		materials[9]->kd = float3(0, 0.6, 0.9);
		par3->rotation(par3, float3(-1, 0, 0), 90);
		par3->translation(par3, float3(-2.5, -5, 0));
		par3->scaling(par3, float3(2, 2, 2));
		objects.push_back(par3);

		Paraboloid* par4 = new Paraboloid(materials[10], 1);
		materials[10]->refractive = true;
		materials.at(10)->minReflectance = float3(0.5, 0.5, 0.5);
		materials.at(10)->refractiveIndex = 1.33;
		materials[10]->kd = float3(0.3, 0, 0.5);
		par4->rotation(par4, float3(-1, 0, 0), 90);
		par4->translation(par4, float3(5, -1, 5));
		par4->scaling(par4, float3(2, 2, 2));
		objects.push_back(par4);
		
		/*****************************************************/
		
	}
	~Scene()
	{
		for (std::vector<LightSource*>::iterator iLightSource = lightSources.begin(); iLightSource != lightSources.end(); ++iLightSource)
			delete *iLightSource;
		for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
			delete *iMaterial;
		for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
			delete *iObject;
	}

public:
	Camera& getCamera()
	{
		return camera;
	}

	// IMPLEMENTED FOR YOUR CONVENIENCE, CALL THIS WHEN APPROPRIATE
	Hit firstIntersect(const Ray& ray)
	{
		Hit bestHit;
		bestHit.t = FLT_MAX;
		for (int oi = 0; oi < objects.size(); oi++)
		{
			Hit hit = objects[oi]->intersect(ray);
			if (hit.t > 0 && hit.t < bestHit.t)
				bestHit = hit;
		}
		if (bestHit.t == FLT_MAX)
			return Hit();
		return bestHit;
	}

	float3 trace(const Ray& ray, int depth)
	{
		if (depth == 0){
			glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			glEnd();
			glutSwapBuffers();
		}

		
		Hit hit = firstIntersect(ray);
		
	
		if (hit.t < 0)
		return ray.dir * ray.dir;
	
		float3 outRadiance(0, 0, 0);
		if (hit.material->reflective){
			float3 reflectionDir = hit.material->reflect(ray.dir, hit.normal);
			Ray reflectedRay(hit.position + hit.normal, reflectionDir);
			outRadiance += trace(reflectedRay, depth+1)*hit.material->getReflectance(ray.dir, hit.normal);
		}

		if (hit.material->refractive){
			float3 refractionDir = hit.material->refract(ray.dir, hit.normal);
			Ray refractedRay(hit.position + hit.normal, refractionDir);
			outRadiance += trace(refractedRay, depth+1)*(float3(1, 1, 1) - hit.material->getReflectance(ray.dir, hit.normal));
		}
		for (int i = 0; i < lightSources.size(); i++){
			//outRadiance += outRadiance + lightSources[i]->getPowerDensityAt(hit.position);

			Ray shadowRay(hit.position + hit.normal*.01, lightSources[i]->getLightDirAt(hit.position));
			Hit shadowHit = firstIntersect(shadowRay);
			if (shadowHit.t < 0 || shadowHit.t >abs(lightSources[i]->getDistanceFrom(hit.position))){
				outRadiance += hit.material->shade(hit.position, hit.normal, -ray.dir, lightSources[i]->getLightDirAt(hit.position),
					lightSources[i]->getPowerDensityAt(hit.position));
				//outRadiance +=lightSources[i]->getPowerDensityAt(hit.position)*(lightSources[i]->getLightDirAt(ray.origin).dot(hit.normal));
			}
			return outRadiance;
		}

		
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// global application data

// screen resolution
const int screenWidth = 640;
const int screenHeight = 640;
// image to be computed by ray tracing
float3 image[screenWidth*screenHeight];

Scene scene;

bool computeImage()
{
	static unsigned int iPart = 0;

	if (iPart >= 64)
		return false;
	for (int j = iPart; j < screenHeight; j += 64)
	{
		for (int i = 0; i < screenWidth; i++)
		{
			float3 pixelColor = float3(0, 0, 0);
			float2 ndcPixelCentre((2.0 * i - screenWidth) / screenWidth, (2.0 * j - screenHeight) / screenHeight);

			Camera& camera = scene.getCamera();
			Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcPixelCentre));

			image[j*screenWidth + i] = scene.trace(ray,1);
		}
	}
	iPart++;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL starts here. In the ray tracing example, OpenGL just outputs the image computed to the array.

// display callback invoked when window needs to be redrawn
void onDisplay() {
	glClearColor(0.0f, 0.2f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear screen

	if (computeImage())
		glutPostRedisplay();
	glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);

	glutSwapBuffers(); // drawing finished
}

int main(int argc, char **argv) {
	glutInit(&argc, argv);						// initialize GLUT
	glutInitWindowSize(640, 640);				// startup window size 
	glutInitWindowPosition(100, 100);           // where to put window on screen
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    // 8 bit R,G,B,A + double buffer + depth buffer

	glutCreateWindow("Ray caster");				// application window is created and displayed

	glViewport(0, 0, screenWidth, screenHeight);

	glutDisplayFunc(onDisplay);					// register callback

	glutMainLoop();								// launch event handling loop

	return 0;
}

