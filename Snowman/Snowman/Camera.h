#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>
#include <algorithm>
#include "perlin.h"

// Skeletal Camera class. Feel free to add custom initialization, set aspect ratio to fit viewport dimensions, or animation.
class Camera
{
	float3 eye;

	float3 lookAt;
	float3 right;
	float3 up;

public:
	float3 getEye()
	{
		return eye;
	}
	Camera()
	{
		eye = float3(0, 0, 3);
		lookAt = float3(0, 0, 2);
		right = float3(1, 0, 0);
		up = float3(0, 1, 0);
	}

	float3 rayDirFromNdc(const float2 ndc) {
		return (lookAt - eye
			+ right * ndc.x
			+ up    * ndc.y
			).normalize();
	}
};
