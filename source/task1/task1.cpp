

#include <limits>
#include <iostream>
#include <iomanip>
#include "Scene.h"
#include "task1.h"

constexpr float epsilon = 0.001f;




bool isZero(float a){
    if(a<=epsilon && a>=-epsilon){
        return true;
    }
    else{
        return false;
    }
}

void printfloat3(math::vector<float, 3> a){
    std::cout<<"{"<<std::setprecision(8)<<a.x<<", "<<std::setprecision(8)<<a.y<<", "<<std::setprecision(8)<<a.z<<"}\n";
}


math::vector<float, 3> normalized(math::vector<float, 3> a){
  float magnitude = a.x* a.x + a.y* a.y +  a.z* a.z;
  magnitude = math::sqrt(magnitude);
  return {a.x/magnitude, a.y/magnitude, a.z/magnitude};
}

float dot(math::vector<float, 3> a, math::vector<float, 3> b){
  return a.x* b.x + a.y* b.y +  a.z* b.z;
}
float dot(math::vector<float, 3> a, math::vector<float, 4> b){
  return a.x* b.x + a.y* b.y +  a.z* b.z;
}
float dot(math::vector<float, 4> a, math::vector<float, 4> b){
  return a.x* b.x + a.y* b.y +  a.z* b.z;
}

float determinant(math::vector<float, 3> a, math::vector<float, 3> b, math::vector<float, 3> c){
    float ans = a.x *(b.y*c.z -b.z*c.y);
    ans -= a.y *(b.x*c.z -b.z*c.x);
    ans += a.z *(b.x*c.y -b.y*c.x);
    return ans;
}

math::vector<float, 3> cross(math::vector<float, 3> a, math::vector<float, 3> b){
    math::vector<float, 3> ans;
    ans.x = a.y*b.z - a.z*b.y;
    ans.y = a.z*b.x - a.x*b.z;
    ans.z = a.x*b.y - a.y*b.x;
    return ans;
}




math::vector<float, 3> vec_product(math::vector<float, 3> a , math::vector<float, 3> b){
  return { a.x* b.x, a.y* b.y, a.z* b.z};
}
float3 reflect(const float3& incident, const float3& normal) {
    return incident - 2.0f * dot(incident, normal) * normal;
}



bool intersectRayPlane(const float3 &p, const float3 &d, const float4 &plane, float &t)
{
  // TODO implement the intersection test between a ray and a plane.
  //
  // A plane is defined by the plane normal n and the offset w along the normal.
  // The plane is given as parameter plane where (plane.x, plane.y, plane.z) represents the plane normal
  // and plane.w is the offset w.
  //
  // If there is no intersection (Hint: or one we do not care about), return false.
  // Otherwise, compute and set the parameter t such that p + t * d yields the intersection point and return true.
    // cout<<
    float dot_ray_n = dot(d, {plane.x, plane.y, plane.z});
    if(isZero(dot_ray_n)){
        //The dot product is zero
        return false;

        //Check if the point p is in the plane
        float T = dot(p, plane) - plane.w;
        if(isZero(T)){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        t = plane.w - dot(p, plane);
        t = t/dot(d, plane);
        return true;
    }
}

bool intersectsRayPlane(const float3 &p, const float3 &d, const Plane *planes, 
    std::size_t num_planes, float t_min, float t_max)
{
  // TODO: implement intersection test between a ray and a set of planes.
  // This method only has to detect whether there is an intersection with ANY
  // plane along the given subset of the ray. The ray is given by its start
  // point p and direction d.
  // A plane is defined by the plane normal n and the offset w along the normal.
  // Each plane in planes contains a float4 member p where the plane normal n is
  // (p.x, p.y, p.z) and w is p.w.
  // If an intersection is found that falls on a point on the ray
  // between t_min and t_max, return true. Otherwise, return false.


    for (int i = 0; i < num_planes; ++i) {
        auto plane = planes[i].p;
        float t;
        bool check = intersectRayPlane(p, d, plane, t);
        if(check){
            if( t >= (t_min- epsilon) && t <= (t_max + epsilon) ){
                return true;
            }
        }
    }
    return false;
}

const Plane *findClosestHitPlanes(const float3 &p, const float3 &d, const Plane *planes, std::size_t num_planes, float &t)
{
  // TODO: implement intersection test between a ray and a set of planes.
  // This function should find the CLOSEST intersection with a plane along
  // the ray. The ray is given by its start point p and direction d.
  // A plane is defined by the plane normal n and the offset w along the normal.
  // Each plane in planes contains a float4 member p where the plane normal n is
  // (p.x, p.y, p.z) and w is p.w.
  // If an intersection is found, set t to the ray parameter and
  // return a pointer to the hit plane.
  // If no intersection is found, return nullptr.
    const Plane * ans = nullptr;
    float ans_t;
    for (int i = 0; i < num_planes; ++i) {
        auto plane = planes[i].p;
        float t;
        bool check = intersectRayPlane(p, d, plane, t);
        if(check){
            if(ans == nullptr){
                ans = &planes[i];
                ans_t = t;
            }
            else{
                if(t < ans_t){
                    ans = &planes[i];
                    ans_t = t;
                }
            }
        }
    }
    t = ans_t;
    return ans;
}

bool intersectRayTriangle(const float3 &p, const float3 &d, const float3 &p1,
                          const float3 &p2, const float3 &p3, float &t,
                          float &lambda_1, float &lambda_2)
{
  // TODO implement the intersection test between a ray and a triangle.
  //
  // The triangle is defined by the three vertices p1, p2 and p3
  //
  // If there is no intersection return false.
  // Otherwise, compute and set the parameters lambda_1 and lambda_2
  // to the barycentric coordinates corresponding to the
  // closest point of intersection

    float A = determinant(p1- p2, p1-p3, d);
    if(isZero(A)){
        // ray is parallel to triangle
        return false;
    }
    lambda_1 = determinant(p1- p, p1-p3, d) / A;
    lambda_2 = determinant(p1- p2, p1-p, d) / A;
    float3 n = cross(p2-p1, p3-p1);
    float w = dot(n, p1);
    t = w - dot(n, p);
    t = t / dot(n, d);
    if(lambda_1 >= 0 && lambda_2 >= 0 && (lambda_1 + lambda_2) < 1){
        return true;
    }
    else{
        return false;
    }
}
const Triangle *findClosestHitTriangles(const float3 &p, const float3 &d,
                               const Triangle *triangles,
                               std::size_t num_triangles,
                               const float3 *vertices, float &t,
                               float &lambda_1, float &lambda_2)
{
  // TODO: implement intersection test between a ray and a set of triangles.
  // This function should find the CLOSEST intersection with a triangle along
  // the ray. The ray is given by its start point p and direction d. A triangle
  // is represented as an array of three vertex indices. The position of each
  // vertex can be looked up from the vertex array via the vertex index.
  // triangles points to the first element of an array of num_triangles
  // triangles. If an intersection is found, set t to the ray parameter and
  // lambda_1 and lambda_2 to the barycentric coordinates corresponding to the
  // closest point of intersection, and return a pointer to the hit triangle.
  // If no intersection is found, return nullptr.


    const Triangle * ans = nullptr;
    float ans_t=std::numeric_limits<float>::infinity(), ans_lambda_1, ans_lambda_2;
    for (int i = 0; i < num_triangles; ++i) {
        auto p1 = vertices[triangles[i][0]];
        auto p2 = vertices[triangles[i][1]];
        auto p3 = vertices[triangles[i][2]];
        float t, lambda_1, lambda_2;
        bool check = intersectRayTriangle(p, d, p1, p2, p3, t, lambda_1, lambda_2);
        if(check){
            if(ans == nullptr){
                ans = &triangles[i];
                ans_t = t;
                ans_lambda_1 = lambda_1;
                ans_lambda_2 = lambda_2;
            }
            else{
                if(t < ans_t){
                    ans = &triangles[i];
                    ans_t = t;
                    ans_lambda_1 = lambda_1;
                    ans_lambda_2 = lambda_2;
                }
            }
        }
        else{
            if(ans == nullptr){
                if(t < ans_t){
                    ans_t = t;
                    ans_lambda_1 = lambda_1;
                    ans_lambda_2 = lambda_2;
                }
            }
        }
    }
    t= ans_t;
    lambda_1= ans_lambda_1;
    lambda_2= ans_lambda_2;
    return ans;
}

bool intersectsRayTriangle(const float3 &p, const float3 &d,
                           const Triangle *triangles, std::size_t num_triangles,
                           const float3 *vertices, float t_min, float t_max)
{
  // TODO: implement intersection test between a ray and a set of triangles.
  // This method only has to detect whether there is an intersection with ANY
  // triangle along the given subset of the ray. The ray is given by its start
  // point p and direction d. A triangle is represented as an array of three
  // vertex indices. The position of each vertex can be looked up from the
  // vertex array via the vertex index. triangles points to an array of
  // num_triangles. If an intersection is found that falls on a point on the ray
  // between t_min and t_max, return true. Otherwise, return false.


    for (int i = 0; i < num_triangles; ++i) {
        auto p1 = vertices[triangles[i][0]];
        auto p2 = vertices[triangles[i][1]];
        auto p3 = vertices[triangles[i][2]];
        float t, lambda_1, lambda_2;
        bool check = intersectRayTriangle(p, d, p1, p2, p3, t, lambda_1, lambda_2);
        if(check){
            if( t > (t_min- epsilon) && t < (t_max + epsilon) ){
                return true;
            }
        }
    }
    return false;
}
float3 shade(const float3 &p, const float3 &d, const HitPoint &hit,
             const Scene &scene, const Pointlight *lights,
             std::size_t num_lights)
{
    float3 shade = {0, 0, 0};
    for (int i = 0; i < num_lights; i++) {
        Pointlight light = lights[i];
        float3 n = normalized(light.position - hit.position);
        float3 d_n = normalized(-d);
        float3 h_n = normalized(hit.normal + d_n); 
        float t_max = length(light.position - hit.position);
        if (scene.intersectsRay(hit.position, n, epsilon, t_max)) {
            continue;
        }
        else {
            float cos_angle = std::max(dot(h_n, n), 0.0f);
            float3 diffuse = vec_product(hit.k_d * cos_angle, light.color);
            float3 r = reflect(-n, hit.normal); 
            float cos_alpha = std::max(dot(r, d_n), 0.0f);
            float3 specular = vec_product(hit.k_s * std::pow(cos_alpha, hit.m), light.color);
            shade += diffuse + specular;
        }
    }
    return shade * 255.0f;
}

void render(image2D<float3> &framebuffer, int left, int top, int right,
            int bottom, const Scene &scene, const Camera &camera,
            const Pointlight *lights, std::size_t num_lights,
            const float3 &background_color, int max_bounces)
{
    // Calculate camera basis vectors
    float3 w = normalized(camera.eye - camera.lookat);
    float3 u = normalized(cross(camera.up, w));
    float3 v = normalized(cross(w, u));

    // Calculate image plane dimensions
    float magnificationFactor = 1; 
    float magnificationFactor2 =  ((float)framebuffer.height) / framebuffer.width; 
    float ws = camera.w_s * magnificationFactor;
    float hs = ws * magnificationFactor2;

    // Loop through the specified region of the framebuffer
 /*   for (int y = top; y < bottom; ++y) {
        for (int x = left; x < right; ++x) {
            framebuffer(x, framebuffer.height - 1 - y) = background_color;
        }
    }
*/
    for (float y = top; y < bottom; ++y) {
        for (float x = left; x < right; ++x) {
            // Compute ray direction through the center of pixel (x, y) on the image plane
            float u_offset = (x +1) * ws / framebuffer.width - ws / 2.0f;
            float v_offset = (y +1) * hs / framebuffer.height - hs / 2.0f;
            float3 ray_direction =100.0f* normalized(-camera.f * w + u_offset * u + v_offset * v);

            // std::cout<< "the lookat is: ";
            // printfloat3(normalized(camera.eye - camera.lookat));
            // std::cout<< "the -camera.f*w  is: ";
            // printfloat3(normalized(-camera.f*w ));
            // Trace ray from camera eye 
            float3 ray_origin = camera.eye;
            float T_MIN=900000;
            std::optional<HitPoint> hit_plane;
            float t_plane;
            // std::cout<< "the ray origin is at: ";
            // printfloat3(ray_origin);
            // std::cout<< "the ray direction is at: ";
            // printfloat3(ray_direction);
            // float ws2 = camera.w_s;
            // float hs2 = ws * framebuffer.height /framebuffer.width ;
            // float3 L = camera.lookat - u * (ws2/2)- v * (hs2/2);
            // std::cout<< "the new ray direction is at: ";
            // printfloat3(normalized(L + u * (x*ws2/(2*framebuffer.width)) + v * ((framebuffer.height - y)*hs2/(2*framebuffer.height))));

           /* const Plane* plane_hit = scene.findClosestHitPlane(ray_origin, ray_direction, t_plane);
            if (plane_hit != nullptr && t_plane >= 0.0f) {
                // Calculate hit point position and normal
                float3 hit_position = ray_origin + t_plane * ray_direction;
                float3 hit_normal = normalized({ plane_hit->p.x, plane_hit->p.y, plane_hit->p.z }); // Plane normal

                Material hit_material = scene.materials.get()[plane_hit->material];
                float3 k_d = hit_material.diffuse;     //  diffuse color
                float3 k_s = hit_material.specular;    //  specular color
                float shine = hit_material.shininess;  //  shininess

                hit_plane = HitPoint{ hit_position, hit_normal, k_d, k_s, shine };
                 if (hit_plane.has_value()) {
                if(t_plane< 0.0f){
                    t_plane = 9000 + 10000;
                }
                if(t_plane < T_MIN){
                    T_MIN = t_plane;
                    float3 color = shade(ray_origin, ray_direction, *hit_plane, scene, lights, num_lights);
                    framebuffer(x, framebuffer.height - 1 - y) = color;
                }
            }
            }*/
           




            std::optional<HitPoint> hit_triangle;
            float t_triangle = 9000, lambda_1, lambda_2;
            
            auto tri = scene.findClosestHitTriangle(ray_origin, ray_direction, t_triangle, lambda_1, lambda_2);
            const Triangle* triangle_hit = std::get<0>(tri);
            int triangle_meterial_id = std::get<1>(tri);
            if (triangle_hit != nullptr && t_triangle > 0.0f)
            {   
                float3 hit_position = ray_origin + (t_triangle * ray_direction);
                float3 hit_position_adj={hit_position.x,hit_position.y,hit_position.z+epsilon};
                float3 v1 = scene.positions.get()[(*triangle_hit)[0]];
                float3 v2 = scene.positions.get()[(*triangle_hit)[1]];
                float3 v3 = scene.positions.get()[(*triangle_hit)[2]];
                float3 hit_normal = normalized(cross(v2-v1, v3-v1)); //Triangle Normal

                Material hit_meterial = scene.materials.get()[triangle_meterial_id];
                float3 k_d = hit_meterial.diffuse;     //  diffuse color
                float3 k_s = hit_meterial.specular;    //  specular color
                float shine = hit_meterial.shininess;  //  shininess
                hit_triangle = HitPoint{ hit_position, hit_normal
                , k_d, k_s, shine };
                
                float3 color = shade(ray_origin, ray_direction, *hit_triangle, scene, lights, num_lights);
                framebuffer(x, framebuffer.height - 1 - y) = color;
            }
         


            std::optional<HitPoint> hit_cone = scene.findClosestHit(ray_origin, ray_direction);
          
            if (1) {
                
                float3 color = shade(ray_origin, ray_direction
                , *hit_cone, scene, lights, num_lights);
                       
                        framebuffer(x, framebuffer.height - 1 - y) = color;}
        }
    }
}