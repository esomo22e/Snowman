#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>
#include <algorithm>
#include "perlin.h"
// Ray structure.
class Ray
{
public:
	float3 origin;
	float3 dir;
	Ray(float3 o, float3 d)
	{
		origin = o;
		dir = d;
	}
};
