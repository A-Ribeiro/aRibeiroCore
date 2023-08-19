#include "Frustum.h"

#include <aRibeiroCore/all_math.h>

namespace aRibeiro {

namespace collision {
    
    void Frustum::computePlanes(const mat4& matrix, bool depth_zero_one) {
        mat4 clipMatrix = transpose(matrix);
        /*
                    0 _11,1 _21,2 _31,3 _41,
                    4 _12,5 _22,6 _32,7 _42,
                    8 _13,9 _23,10_33,11_43,
                    12_14,13_24,14_34,15_44;
        */
        
        vec4 aux = clipMatrix[3] - clipMatrix[0];
        rightPlane.normal = toVec3(aux);
        rightPlane.distance = -aux.w;
        rightPlane.normalize();

        aux = clipMatrix[3] + clipMatrix[0];
        leftPlane.normal = toVec3(aux);
        leftPlane.distance = -aux.w;
        leftPlane.normalize();

        aux = clipMatrix[3] + clipMatrix[1];
        bottomPlane.normal = toVec3(aux);
        bottomPlane.distance = -aux.w;
        bottomPlane.normalize();

        aux = clipMatrix[3] - clipMatrix[1];
        topPlane.normal = toVec3(aux);
        topPlane.distance = -aux.w;
        topPlane.normalize();

        if (depth_zero_one)
            aux = clipMatrix[2];
        else
            aux = clipMatrix[3] + clipMatrix[2];
        nearPlane.normal = toVec3(aux);
        nearPlane.distance = -aux.w;
        nearPlane.normalize();

        aux = clipMatrix[3] - clipMatrix[2];
        farPlane.normal = toVec3(aux);
        farPlane.distance = -aux.w;
        farPlane.normalize();
        
        //
        // Compute vertices
        //
        
        Plane::intersectPlanes(rightPlane, topPlane, nearPlane,
                               &vertices[FrustumVertex_Near_Right_Top]);
        Plane::intersectPlanes(rightPlane, bottomPlane, nearPlane,
                               &vertices[FrustumVertex_Near_Right_Bottom]);
        Plane::intersectPlanes(leftPlane, bottomPlane, nearPlane,
                               &vertices[FrustumVertex_Near_Left_Bottom]);
        Plane::intersectPlanes(leftPlane, topPlane, nearPlane,
                               &vertices[FrustumVertex_Near_Left_Top]);
        
        Plane::intersectPlanes(rightPlane, topPlane, farPlane,
                               &vertices[FrustumVertex_Far_Right_Top]);
        Plane::intersectPlanes(rightPlane, bottomPlane, farPlane,
                               &vertices[FrustumVertex_Far_Right_Bottom]);
        Plane::intersectPlanes(leftPlane, bottomPlane, farPlane,
                               &vertices[FrustumVertex_Far_Left_Bottom]);
        Plane::intersectPlanes(leftPlane, topPlane, farPlane,
                               &vertices[FrustumVertex_Far_Left_Top]);
        
        
        aabb = AABB(vertices[0],vertices[1]);
        for (int i = 2; i < 8 ; i ++ ) {
            aabb.min_box = minimum(aabb.min_box, vertices[i]);
            aabb.max_box = maximum(aabb.max_box, vertices[i]);
        }

        //aabb = AABB::joinAABB(aabb, AABB(vertices[i],vertices[i+1]) );
        
        //pre-calculate the vertex projections over the frustum plane
        for (int i = 0; i < 6; i++)
            projectOnAxis(vertices, 8, (*this)[i].normal , &minProjections[i], &maxProjections[i]);
        

        //
        // OBB calculation...
        //
        
        //vec3 center = nearPlane.normal * ( (nearPlane.distance - farPlane.distance) * 0.5f );
        vec3 right = normalize( vertices[FrustumVertex_Far_Right_Bottom] - vertices[FrustumVertex_Far_Left_Bottom] );
        //vec3 up = normalize( vertices[FrustumVertex_Far_Right_Top] - vertices[FrustumVertex_Far_Right_Bottom] );
        vec3 front = nearPlane.normal;

        const vec4 _0001 = vec4(0, 0, 0, 1);
        mat4 rotationBase = mat4(toVec4(right), toVec4(cross(front,right)), toVec4(front), _0001);

        quat rotationBase_quat = extractQuat(
            rotationBase
        );

        float dimension_min[3];
        float dimension_max[3];

        float center_local[3];
        vec3 dimension;

        vec3 center;// = inv( rotationBase_quat ) * center_local;

        for (int i=0;i<3;i++) {
            projectOnAxis(vertices, 8, toVec3(rotationBase[i]), &dimension_min[i], &dimension_max[i]);
            dimension[i] = dimension_max[i] - dimension_min[i];
            center_local[i] = (dimension_min[i] + dimension_max[i]) * 0.5f;
            center += toVec3(rotationBase[i]) * center_local[i];
        }

        obb = OBB( center, dimension, rotationBase_quat );
    }

    Plane& Frustum::operator[](int idx){
        return ((&rightPlane)[idx]);
    }
    
    const Plane& Frustum::operator[](int idx)const {
        return ((&rightPlane)[idx]);
    }

    Frustum::Frustum() {

    }
    
    Frustum::Frustum(const mat4& projection, bool depth_zero_one) {
        computePlanes(projection, depth_zero_one);
    }

    Frustum::Frustum(const mat4& projection, const mat4& camera, bool depth_zero_one) {
        computePlanes(projection * camera, depth_zero_one);
    }


    bool Frustum::pointInsideFrustum(const vec3 &p, const Frustum &frustum) {
        for (int i = 0; i < 6; i++) {
            if (Plane::pointDistanceToPlane(p, frustum[i]) < 0)
                return false;
        }
        return true;
    }

    bool Frustum::sphereOverlapsFrustum(const vec3 &center, const float &radius, const Frustum &frustum) {
        for (int i = 0; i < 6; i++) {
            if (Plane::pointDistanceToPlane(center, frustum[i]) < -radius)
                return false;
        }
        return true;
    }

    bool Frustum::sphereOverlapsFrustum(const Sphere &s, const Frustum &frustum) {
        return Frustum::sphereOverlapsFrustum(s.center, s.radius, frustum);
    }
    
    bool Frustum::aabbOverlapsFrustum(const AABB &aabb, const Frustum &frustum) {
        
        if (!AABB::aabbOverlapsAABB(aabb, frustum.aabb))
            return false;

        /*
        // Test the box normals (x-, y- and z-axes)
        vec3 boxNormals[] = {
            vec3(1,0,0),
            vec3(0,1,0),
            vec3(0,0,1)
        };
        
        float frustumMin,frustumMax;
        for (int i = 0; i < 3; i++)
        {
            projectOnAxis(frustum.vertices, 8, boxNormals[i], &frustumMin, &frustumMax);
            if (frustumMax < aabb.min_box[i] - EPSILON || frustumMin > aabb.max_box[i] + EPSILON)
                return false; // No intersection possible.
        }
        */
        
        
        vec3 box_Vertices[] = {
            vec3(aabb.min_box.x,aabb.min_box.y,aabb.min_box.z),// 000
            vec3(aabb.min_box.x,aabb.min_box.y,aabb.max_box.z),// 001
            vec3(aabb.min_box.x,aabb.max_box.y,aabb.min_box.z),// 010
            vec3(aabb.min_box.x,aabb.max_box.y,aabb.max_box.z),// 011
            vec3(aabb.max_box.x,aabb.min_box.y,aabb.min_box.z),// 100
            vec3(aabb.max_box.x,aabb.min_box.y,aabb.max_box.z),// 101
            vec3(aabb.max_box.x,aabb.max_box.y,aabb.min_box.z),// 110
            vec3(aabb.max_box.x,aabb.max_box.y,aabb.max_box.z)// 111
        };
        
        float boxMin, boxMax;
        
        for (int i = 0; i < 6; i++) {
            projectOnAxis(box_Vertices, 8, frustum[i].normal , &boxMin, &boxMax);
            //projectOnAxis(frustum.vertices, 8, frustum[i].normal , &frustum.minProjections[i], &frustum.maxProjections[i]);
            if (boxMax < frustum.minProjections[i] - EPSILON || boxMin > frustum.maxProjections[i] + EPSILON)
                return false; // No intersection possible.
        }
        
        return true;
        
        /*
        vec3 center = (aabb.min_box + aabb.max_box) * 0.5;
        vec3 halfbounds = (aabb.max_box - aabb.min_box) * 0.5;

        for (int i = 0; i < 6; i++) {
            
            if (Plane::pointDistanceToPlane(center + halfbounds * vec3(1, 1, 1), frustum[i]) >= 0.0f) continue; // 000
            if (Plane::pointDistanceToPlane(center + halfbounds * vec3(1, 1, -1), frustum[i]) >= 0.0f) continue; // 001
            if (Plane::pointDistanceToPlane(center + halfbounds * vec3(1, -1, 1), frustum[i]) >= 0.0f) continue; // 010
            if (Plane::pointDistanceToPlane(center + halfbounds * vec3(1, -1, -1), frustum[i]) >= 0.0f) continue; // 011
            if (Plane::pointDistanceToPlane(center + halfbounds * vec3(-1, 1, 1), frustum[i]) >= 0.0f) continue; // 100
            if (Plane::pointDistanceToPlane(center + halfbounds * vec3(-1, 1, -1), frustum[i]) >= 0.0f) continue; // 101
            if (Plane::pointDistanceToPlane(center + halfbounds * vec3(-1, -1, 1), frustum[i]) >= 0.0f) continue; // 110
            if (Plane::pointDistanceToPlane(center + halfbounds * vec3(-1, -1, -1), frustum[i]) >= 0.0f) continue; // 111

            return false;
        }
        return true;
        */
    }


    bool Frustum::obbOverlapsFrustum(const OBB &obb, const Frustum &frustum) {

        float frustumMin,frustumMax;
        for (int i = 0; i < 3; i++)
        {
            projectOnAxis(frustum.vertices, 8, obb.right_up_depth_vec[i], &frustumMin, &frustumMax);
            if (frustumMax < obb.center_right_up_depth_proj_min[i] || frustumMin > obb.center_right_up_depth_proj_max[i])
                return false; // No intersection possible.
        }

        float boxMin, boxMax;
        for (int i = 0; i < 6; i++) {
            projectOnAxis(obb.box_vertices, 8, frustum[i].normal , &boxMin, &boxMax);
            if (boxMax < frustum.minProjections[i] || boxMin > frustum.maxProjections[i])
                return false; // No intersection possible.
        }
        
        return true;

    }

}

}
