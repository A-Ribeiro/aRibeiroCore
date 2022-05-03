#include "OBB.h"
#include "../all_math.h"

#include "AABB.h"

namespace aRibeiro {

namespace collision {

    void OBB::computeOptimizationData() {

        right_up_depth_vec[0] = orientation * vec3(1,0,0);
        right_up_depth_vec[1] = orientation * vec3(0,1,0);
        right_up_depth_vec[2] = orientation * vec3(0,0,1);

        for (int i=0;i<3;i++) {
            float center_projected = dot(center, right_up_depth_vec[i]);
            center_right_up_depth_proj_min[i] = center_projected - dimension_2[i];
            center_right_up_depth_proj_max[i] = center_projected + dimension_2[i];
        }

        vec3 right = right_up_depth_vec[0] * dimension_2.x;
        vec3 top = right_up_depth_vec[1] * dimension_2.y;
        vec3 depth = right_up_depth_vec[2] * dimension_2.z;

        box_vertices[0] = center - right - top - depth; //000
        box_vertices[1] = center - right - top + depth; //001
        box_vertices[2] = center - right + top - depth; //010
        box_vertices[3] = center - right + top + depth; //011
        box_vertices[4] = center + right - top - depth; //100
        box_vertices[5] = center + right - top + depth; //101
        box_vertices[6] = center + right + top - depth; //110
        box_vertices[7] = center + right + top + depth; //111

       /*

        box_vertices[0] = center + orientation * vec3(-dimension_2.x,-dimension_2.y,-dimension_2.z);// 000
        box_vertices[1] = center + orientation * vec3(-dimension_2.x,-dimension_2.y, dimension_2.z);// 001
        box_vertices[2] = center + orientation * vec3(-dimension_2.x, dimension_2.y,-dimension_2.z);// 010
        box_vertices[3] = center + orientation * vec3(-dimension_2.x, dimension_2.y, dimension_2.z);// 011
        box_vertices[4] = center + orientation * vec3( dimension_2.x,-dimension_2.y,-dimension_2.z);// 100
        box_vertices[5] = center + orientation * vec3( dimension_2.x,-dimension_2.y, dimension_2.z);// 101
        box_vertices[6] = center + orientation * vec3( dimension_2.x, dimension_2.y,-dimension_2.z);// 110
        box_vertices[7] = center + orientation * vec3( dimension_2.x, dimension_2.y, dimension_2.z);// 111
*/
    }

    OBB::OBB() {
        computeOptimizationData();
    }
    OBB::OBB(const vec3 &_center, const vec3 &_dimension, const quat &_orientation) {
        setOBB(_center, _dimension, _orientation);
    }
    OBB::OBB(const mat4 &_transform) {
        setOBB(_transform);
    }

    void OBB::setOBB(const vec3 &_center, const vec3 &_dimension, const quat &_orientation) {
        
        center = _center;
        dimension_2 = _dimension * 0.5f;
        orientation = _orientation;

        computeOptimizationData();
    }

    void OBB::setOBB(const mat4 &_transform) {
        
        center = toVec3(_transform[3]);
        dimension_2 = vec3(
            length(_transform[0]),
            length(_transform[1]),
            length(_transform[2])
        ) * 0.5f;
        orientation = extractQuat( _transform );

        computeOptimizationData();
    }

    //--------------------------------------------------------------------------
    // Static methods
    //--------------------------------------------------------------------------

    bool OBB::pointInsideOBB(const vec2& ptn,const OBB& b) {
        for (int i=0;i<2;i++) {
            float ptn_projected = dot(ptn, vec2( b.right_up_depth_vec[i].x, b.right_up_depth_vec[i].y ) );
            if ( ( ptn_projected < b.center_right_up_depth_proj_min[i] ) || ( ptn_projected > b.center_right_up_depth_proj_max[i] ) )
                return false;
        }
        return true;
    }

    bool OBB::pointInsideOBB(const vec3& ptn,const OBB& b) {
        for (int i=0;i<3;i++) {
            float ptn_projected = dot(ptn, b.right_up_depth_vec[i]);
            if ( ( ptn_projected < b.center_right_up_depth_proj_min[i] ) || ( ptn_projected > b.center_right_up_depth_proj_max[i] ) )
                return false;
        }
        return true;
    }

    bool OBB::obbOverlapsOBB(const OBB& a,const OBB& b) {

        float min, max;

        //test vertices from 'b' into 'a' orientation
        for (int i = 0; i < 3; i++)
        {
            projectOnAxis(b.box_vertices, 8, a.right_up_depth_vec[i], &min, &max);
            if ((max < a.center_right_up_depth_proj_min[i]) || (min > a.center_right_up_depth_proj_max[i]) )
                return false; // No intersection possible.
        }

        //test vertices from 'a' into 'b' orientation
        for (int i = 0; i < 3; i++)
        {
            projectOnAxis(a.box_vertices, 8, b.right_up_depth_vec[i], &min, &max);
            if ((max < b.center_right_up_depth_proj_min[i]) || (min > b.center_right_up_depth_proj_max[i]) )
                return false; // No intersection possible.
        }

        return true;
    }

    OBB OBB::fromAABB(const AABB& a) {
        vec3 center = (a.min_box + a.max_box) * 0.5f;
        vec3 dimension = (a.max_box - a.min_box);
        return OBB(center, dimension, quat());
    }

    OBB OBB::fromTriangle(const vec3& a, const vec3& b, const vec3& c) {
        vec3 ab = b-a;
        vec3 ab_unit = normalize(ab);

        vec3 normal = cross(ab, c-a);
        vec3 normal_unit = normalize( normal );

        const vec4 _0001 = vec4(0, 0, 0, 1);
        mat4 rotationBase = mat4(toVec4(ab_unit), toVec4(cross(normal_unit,ab_unit)), toVec4(normal_unit), _0001);

        quat rotationBase_quat = extractQuat(
            rotationBase
        );

        float dimension_min[3];
        float dimension_max[3];

        vec3 vertices[3] = {a,b,c};
        
        float center_local[3];
        vec3 dimension;

        vec3 center;// = inv( rotationBase_quat ) * center_local;

        for (int i=0;i<3;i++) {
            projectOnAxis(vertices, 3, toVec3(rotationBase[i]), &dimension_min[i], &dimension_max[i]);
            dimension[i] = dimension_max[i] - dimension_min[i];
            center_local[i] = (dimension_min[i] + dimension_max[i]) * 0.5f;
            center += toVec3(rotationBase[i]) * center_local[i];
        }

        //vec3 center = inv( rotationBase_quat ) * center_local;

        return OBB( center, dimension, rotationBase_quat );
    }

    OBB OBB::fromTriangle(const Triangle& triangle) {
        return OBB::fromTriangle(triangle.a,triangle.b,triangle.c);
    }

    OBB OBB::fromSphere(const vec3& pos, const float &radius) {
        return OBB::fromAABB( AABB::fromSphere(pos, radius) );
    }
    
    OBB OBB::fromSphere(const Sphere& sphere) {
        return OBB::fromSphere(sphere.center, sphere.radius);
    }

    OBB OBB::fromLineSegment(const vec3& a, const vec3& b) {
        vec3 center = (a+b) * 0.5f;
        vec3 dimension;
        dimension.z = distance(a,b);
        
        vec3 ab = normalize( b-a );
        vec3 up = vec3(0,1,0);

        // 0º case construction...
        if (angleBetween(ab,up) <= EPSILON){
            const quat fixedRotation = quatFromEuler(DEG2RAD(0.5f), DEG2RAD(0.5f), DEG2RAD(0.5f));
            up = fixedRotation * vec3(0,1,0);
        }

        quat orientation = quatLookAtRotationLH( ab, up );

        return OBB( center, dimension, orientation );
    }

    OBB OBB::fromLineSegment(const LineSegment& ls) {
        return OBB::fromLineSegment(ls.a,ls.b);
    }

    OBB OBB::fromFrustum(const Frustum& frustum) {
        return frustum.obb;
    }
    
    bool OBB::raycastOBB(const Ray &r, const OBB& a, float *outTmin, vec3 *outNormal) {
        Ray r_transformed = Ray( 
            r.origin - a.center,
            inv(a.orientation) * r.dir
        );
        AABB aabb = AABB( -a.dimension_2, a.dimension_2 );

        float t;
        vec3 normal;
        if (AABB::raycastAABB( r_transformed, aabb, &t, &normal ) ){
            //vec3 pt = r_transformed.dir * t; //r_transformed.origin +
            //pt = a.orientation * pt;// + a.center - r.origin;
            //float t_2 = dot(pt, r.dir);

            normal = a.orientation * normal;
            
            *outTmin = t;
            *outNormal = normal;

            return true;
        }

        return false;
    }

    bool OBB::segmentIntersectsOBB(const vec3& p0, const vec3& p1, const OBB &obb) {

        vec3 p0_p1 = p1 - p0;
        float lgth = length(p0_p1);
        vec3 p0_p1_dir = p0_p1;
        if (lgth <= EPSILON){
            lgth = 0.0f;
            p0_p1_dir = vec3(1,0,0);
        } else {
            p0_p1_dir = p0_p1 / lgth;
        }

        vec3 p0_transformed = p0 - obb.center;
        vec3 p0_p1_dir_transformed = inv(obb.orientation) * p0_p1_dir;
        vec3 p1_transformed = p0_transformed + p0_p1_dir_transformed * lgth;

        AABB aabb = AABB( -obb.dimension_2, obb.dimension_2 );

        return AABB::segmentIntersectsAABB( p0_transformed, p1_transformed, aabb );
    }

    bool OBB::segmentIntersectsOBB(const LineSegment& ls, const OBB &obb) {
        return OBB::segmentIntersectsOBB(ls.a, ls.b, obb);
    }

    vec3 OBB::closestPointToOBB(const vec3 &p, const OBB &obb) {

        vec3 d = p - obb.center;
        // Start result at center of box; make steps from there
        vec3 q = obb.center;
        // For each OBB axis...
        for (int i = 0; i < 3; i++) {
            // ...project d onto that axis to get the distance
            // along the axis of d from the box center
            float dist = dot(d, obb.right_up_depth_vec[i]);
            // If distance farther than the box extents, clamp to the box
            dist = clamp(dist, obb.center_right_up_depth_proj_min[i], obb.center_right_up_depth_proj_max[i]);

            //if (dist > b.e[i]) dist = b.e[i];
            //if (dist < -b.e[i]) dist = -b.e[i];
            // Step that distance along the axis to get world coordinate
            q += dist * obb.right_up_depth_vec[i];
        }

        return q;
    }

    bool OBB::aabbOverlapsOBB(const AABB& _a,const OBB& b) {
        //OBB a_OBB = OBB::fromAABB(_a);
        //return OBB::OBBOverlapsOBB(a_OBB, b);

        float min, max;

        /*
        vec3 a_normals[3] = {
            vec3(1,0,0),
            vec3(0,1,0),
            vec3(0,0,1)
        };*/

        vec3 b_min, b_max;
        b_min = b_max = b.box_vertices[0];
        for (int i = 1; i < 8; i++){
            b_min = minimum(b_min, b.box_vertices[i]);
            b_max = maximum(b_max, b.box_vertices[i]);
        }

        //test vertices from 'b' into 'a' orientation
        for (int i = 0; i < 3; i++)
        {
            //projectOnAxis(b.box_vertices, 8, a_normals[i], &min, &max);
            //if ((max < _a.min_box[i]) || (min > _a.max_box[i]) )
            if ((b_max[i] < _a.min_box[i]) || (b_min[i] > _a.max_box[i]) )
                return false; // No intersection possible.
        }

        vec3 a_box_vertices[] = {
            vec3(_a.min_box.x,_a.min_box.y,_a.min_box.z),// 000
            vec3(_a.min_box.x,_a.min_box.y,_a.max_box.z),// 001
            vec3(_a.min_box.x,_a.max_box.y,_a.min_box.z),// 010
            vec3(_a.min_box.x,_a.max_box.y,_a.max_box.z),// 011
            vec3(_a.max_box.x,_a.min_box.y,_a.min_box.z),// 100
            vec3(_a.max_box.x,_a.min_box.y,_a.max_box.z),// 101
            vec3(_a.max_box.x,_a.max_box.y,_a.min_box.z),// 110
            vec3(_a.max_box.x,_a.max_box.y,_a.max_box.z),// 111
        };

        //test vertices from 'a' into 'b' orientation
        for (int i = 0; i < 3; i++)
        {
            projectOnAxis(a_box_vertices, 8, b.right_up_depth_vec[i], &min, &max);
            if ((max < b.center_right_up_depth_proj_min[i]) || (min > b.center_right_up_depth_proj_max[i]) )
                return false; // No intersection possible.
        }

        return true;
    }

    bool OBB::sphereOverlapsOBB(const vec3 &center, const float &radius, const OBB &obb, vec3 *penetration) {
        vec3 closestPointInTriangle = closestPointToOBB(center, obb);

        vec3 SphereToOBB = closestPointInTriangle - center;
        float sqrLength_SphereToOBB = dot(SphereToOBB, SphereToOBB);

        float Max_Radius_sqr = radius * radius;

        if (sqrLength_SphereToOBB > 0.00002f && sqrLength_SphereToOBB < Max_Radius_sqr) {

            float Length_SphereToTriangle = sqrtf(sqrLength_SphereToOBB);
            vec3 SphereToTriangleNorm = SphereToOBB * (1.0f / Length_SphereToTriangle); //normalize(SphereToTriangle);
            //Vector3 triangleNormal = Vectormath::Aos::normalize( Vectormath::Aos::cross( p2-p1 , p3-p1 ) );

            //EPSILON - to avoid process the same triangle again...
            const float EPSILON = 0.002f;
            *penetration = SphereToTriangleNorm * (radius - Length_SphereToTriangle + EPSILON);
            return true;
        }
        return false;
    }

    bool OBB::sphereOverlapsOBB(const Sphere &sphere, const OBB &obb, vec3 *penetration) {
        return OBB::sphereOverlapsOBB(sphere.center, sphere.radius, obb, penetration);
    }

    bool OBB::planeIntersectsOBB(const Plane &plane, const OBB &obb) {
        float min, max;
        projectOnAxis( obb.box_vertices, 8, plane.normal, &min, &max );
        if ( plane.distance < min || plane.distance > max )
            return false;
        return true;
    }

    bool OBB::triangleIntersectsOBB(const vec3 &v0, const vec3 &v1, const vec3 &v2, const OBB &obb) {
        const float EPSILON = 1e-6f;

        float triangleMin, triangleMax;
        float boxMin, boxMax;

        vec3 triangle_Vertices[] = {
            v0,v1,v2
        };

        for (int i = 0; i < 3; i++)
        {
            projectOnAxis(triangle_Vertices, 3, obb.right_up_depth_vec[i], &triangleMin, &triangleMax);
            if ((triangleMax < obb.center_right_up_depth_proj_min[i] - EPSILON) || (triangleMin > obb.center_right_up_depth_proj_max[i] + EPSILON) )
                return false; // No intersection possible.
        }

        // Test the triangle normal
        vec3 triangle_Normal = normalize(cross(v1 - v0, v2 - v0));
        float triangleOffset = dot(triangle_Normal, v0);

        projectOnAxis(obb.box_vertices, 8, triangle_Normal, &boxMin, &boxMax);

        if ( (boxMax < triangleOffset - EPSILON) || (boxMin > triangleOffset + EPSILON) )
            return false; // No intersection possible.

        //
        // TODO: test the nine edge cross-products or single edge case
        //
        // Test the single edge
        vec3 triangleEdges[] = {
            v0 - v1,
            v1 - v2,
            v2 - v0
        };
        
        for (int i = 0; i < 3; i++)
        {
            // The box normals are the same as it's edge tangents
            vec3 axis = normalize( cross(triangleEdges[i], triangle_Normal) );
            projectOnAxis(obb.box_vertices, 8, axis, &boxMin, &boxMax);
            projectOnAxis(triangle_Vertices, 3, axis, &triangleMin, &triangleMax);
            if ( (boxMax < triangleMin - EPSILON) || (boxMin > triangleMax + EPSILON) )
                return false; // No intersection possible
        }
        
        /*
        // Test the nine edge cross-products
        vec3 triangleEdges[] = {
            v0 - v1,
            v1 - v2,
            v2 - v0
        };

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                // The box normals are the same as it's edge tangents
                vec3 axis = cross(triangleEdges[i], boxNormals[j]);
                projectOnAxis(box_Vertices, 8, axis, &boxMin, &boxMax);
                projectOnAxis(triangle_Vertices, 3, axis, &triangleMin, &triangleMax);
                if (boxMax < triangleMin - EPSILON || boxMin > triangleMax + EPSILON)
                    return false; // No intersection possible
            }
        */
        
        // No separating axis found.
        return true;
    }

    bool OBB::triangleIntersectsOBB(const Triangle &t, const OBB &obb) {
        return OBB::triangleIntersectsOBB(t.a, t.b, t.c, obb);
    }

    //
    // Cloned methods from other collision classes
    //

    bool OBB::frustumOverlapsOBB(const Frustum &f, const OBB &obb) {
        return Frustum::obbOverlapsFrustum(obb, f);
    }

}

}
