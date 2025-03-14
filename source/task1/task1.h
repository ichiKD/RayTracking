


#ifndef INCLUDED_TASK1B
#define INCLUDED_TASK1B

#pragma once

#include <cstddef>
#include <array>

#include <framework/math/vector.h>
#include <framework/math/matrix.h>
#include <framework/image.h>


struct Plane;
class Scene;

using Triangle = std::array<int, 3>;

struct HitPoint
{
	float3 position;
	float3 normal;
	float3 k_d;
	float3 k_s;
	float m;
};

struct Camera
{
	float w_s;
	float f;
	float3 eye;
	float3 lookat;
	float3 up;
};

struct Pointlight
{
	float3 position;
	float3 color;
};

bool intersectRayPlane(const float3& p, const float3& d, const float4& plane, float& t);
bool intersectsRayPlane(const float3& p, const float3& d, const Plane* planes, std::size_t num_planes, float t_min, float t_max);
const Plane* findClosestHitPlanes(const float3& p, const float3& d, const Plane* planes, std::size_t num_planes, float& t);

bool intersectsRayTriangle(const float3& p, const float3& d, const Triangle* triangles, std::size_t num_triangles, const float3* vertices, float t_min, float t_max);
const Triangle* findClosestHitTriangles(const float3& p, const float3& d, const Triangle* triangles, std::size_t num_triangles, const float3* vertices, float& t, float& lambda_1, float& lambda_2);

float3 shade(const float3& p, const float3& d, const HitPoint& hit, const Scene& scene, const Pointlight* lights, std::size_t num_lights);

void render(image2D<float3>& framebuffer, int left, int top, int right, int bottom, const Scene& scene, const Camera& camera, const Pointlight* lights, std::size_t num_lights, const float3& background_color, int max_bounces);

#endif  // INCLUDED_TASK1B
