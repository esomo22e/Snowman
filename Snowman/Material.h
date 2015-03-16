#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>
#include <algorithm>
#include "perlin.h"

// Skeletal Material class. Feel free to add methods e.g. for illumination computation (shading).
class Material
{
public:
	bool reflective;
	bool refractive;
	bool textured;
	float3 minReflectance;		// Fresnel coefficient
	float refractiveIndex;			// index of refraction
	float3 kd;			// diffuse reflection coefficient
	float3 ks;			// specular reflection coefficient
	float shininess;	// specular exponent
	Material()
	{
		reflective = false;
		refractive = false;
		textured = false;
		minReflectance = float3(0.93, 0.85, 0.4);
		refractiveIndex = 1;
		kd = float3(0.5, 0.5, 0.5) + kd * 0.5;
		ks = float3(1, 1, 1);
		shininess = 15;
	}

	float3 reflect(float3 inDir, float3 normal) {
		return inDir - normal * normal.dot(inDir) * 2;
	}
	float3 refract(float3 inDir, float3 normal) {
		float ri = refractiveIndex;
		float cosa = -normal.dot(inDir);
		if (cosa < 0) { cosa = -cosa; normal = -normal; ri = 1 / ri; }
		float disc = 1 - (1 - cosa * cosa) / ri / ri;
		if (disc < 0) return reflect(inDir, normal);
		return inDir * (1.0 / ri) + normal * (cosa / ri - sqrt(disc));
	}
	float3 getReflectance(float3 inDir, float3 normal) {
		float cosa = fabs(normal.dot(inDir));
		return minReflectance + (float3(1, 1, 1) - minReflectance) * pow(1 - cosa, 5);
	}

	virtual float3 shade(
		float3 position,
		float3 normal,
		float3 viewDir,
		float3 lightDir,
		float3 lightPowerDensity)
	{
		float cosTheta = normal.dot(lightDir);
		if (cosTheta < 0) return float3(0, 0, 0);
		float3 halfway = (viewDir + lightDir).normalize();
		float cosDelta = normal.dot(halfway);
		float3 diffuse = kd * lightPowerDensity * cosTheta;
		if (cosDelta < 0) return diffuse;
		return diffuse + lightPowerDensity * ks * pow(cosDelta, shininess);
		// TO BE IMPLEMENTED AT PRACTICAL
	}
};