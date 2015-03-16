#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>
#include <algorithm>
#include "perlin.h"

// Abstract base class for light sources
class LightSource
{
public:
	virtual float3 getPowerDensityAt(float3 x) = 0;
	virtual float3 getLightDirAt(float3 x) = 0;
	virtual float  getDistanceFrom(float3 x) = 0;

};