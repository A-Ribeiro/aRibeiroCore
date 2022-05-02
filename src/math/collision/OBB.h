#ifndef OBB_hpp__
#define OBB_hpp__

#include <aRibeiroCore/common.h>

//I need a complete definition to instantiate a vec3 as class field
#include <aRibeiroCore/vec3.h>
#include <aRibeiroCore/vec2.h>
#include <aRibeiroCore/quat.h>
#include <aRibeiroCore/mat4.h>

namespace aRibeiro {

//class vec2;

namespace collision {

class Sphere;
class LineSegment;
class Plane;
class Ray;
class Triangle;
class Frustum;
class AABB;

/// \brief Oriented Bounding Box (OBB)
///
/// Stores 3D points to represent an Oriented Bounding Box. <br/>
/// It can be used to make collision tests
///
/// \author Alessandro Ribeiro
///
class _SSE2_ALIGN_PRE OBB{

    

    public:

    vec3 center;///<Store the center of the OBB
    vec3 dimension_2;///<Store the box dimension/2 of the OBB
    quat orientation;///<Store the rotation of the OBB

    vec3 right_up_depth_vec[3];
    float center_right_up_depth_proj_min[3];
    float center_right_up_depth_proj_max[3];

    vec3 box_vertices[8];

    void computeOptimizationData();


    //--------------------------------------------------------------------------
    /// \brief Construct a ZERO OBB class
    ///
    /// The ZERO OBB class have center in the origin (0,0,0)
    /// orientation on quaternion identity (0,0,0,1)
    /// dimension_2 in zero (0,0,0)
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// OBB obb = OBB();
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    ///
    OBB();
    OBB(const vec3 &_center, const vec3 &_dimension, const quat &_orientation);
    OBB(const mat4 &_transform);

    void setOBB(const vec3 &_center, const vec3 &_dimension, const quat &_orientation);
    void setOBB(const mat4 &_transform);

    //--------------------------------------------------------------------------
    // Static methods
    //--------------------------------------------------------------------------

    /// \brief Test if a point is inside an OBB
    /// \author Alessandro Ribeiro
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// OBB obb;
    /// vec2 point;
    ///
    /// if ( OBB::pointInsideOBB( point, obb ) ){
    ///     ...
    /// }
    /// \endcode
    ///
    /// \param ptn 2D point to test against the OBB
    /// \param b The OBB
    /// \return true if the point is inside the OBB, otherwise false
    /// \warning The z is not used in the test
    ///
    static bool pointInsideOBB(const vec2& ptn,const OBB& b);
    /// \brief Test if a point is inside an OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// OBB obb;
    /// vec3 point;
    ///
    /// if ( OBB::pointInsideOBB( point, obb ) ){
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param ptn 3D point to test against the OBB
    /// \param b The OBB
    /// \return true if the point is inside the OBB, otherwise false
    ///
    static bool pointInsideOBB(const vec3& ptn,const OBB& b);

    /// \brief Test if there is some overlaped area between two AABBs
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// OBB obb_a;
    /// OBB obb_b;
    ///
    /// if ( OBB::obbOverlapsOBB( obb_a, obb_b ) ){
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param a The first OBB to test
    /// \param b The second OBB to test against
    /// \return true if there are some overlaped area between the two AABBs, otherwise false
    ///
    static bool obbOverlapsOBB(const OBB& a,const OBB& b);

    /// \brief Create an OBB that that contains the AABB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// AABB aabb;
    ///
    /// OBB obb = OBB::fromAABB( aabb );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param a The AABB
    /// \return OBB containing the AABB
    ///
    static OBB fromAABB(const AABB& a);

    /// \brief Create an OBB that that contains the triangle vertex a, b and c
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// // triangle vertex
    /// vec3 a, b, c;
    ///
    /// OBB obb = OBB::fromTriangle( a, b, c );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param a Triangle vertex a
    /// \param b Triangle vertex b
    /// \param c Triangle vertex c
    /// \return OBB containing the triangle
    ///
    static OBB fromTriangle(const vec3& a, const vec3& b, const vec3& c);

    /// \brief Create an OBB that that contains the triangle
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// Triangle triangle;
    ///
    /// OBB obb = OBB::fromTriangle( triangle );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param triangle Triangle
    /// \return OBB containing the triangle
    ///
    static OBB fromTriangle(const Triangle& triangle);

    /// \brief Create an OBB that that contains the sphere with pos as position and radius
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// // sphere parameters
    /// vec3 sphere_pos;
    /// float sphere_radius;
    ///
    /// OBB obb = OBB::fromSphere( sphere_pos, sphere_radius );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param pos sphere center position
    /// \param radius sphere radius
    /// \return OBB containing the sphere
    ///
    static OBB fromSphere(const vec3& pos, const float &radius);
    
    /// \brief Create an OBB that that contains the sphere
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// Sphere sphere;
    ///
    /// OBB obb = OBB::fromSphere( sphere );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param sphere Sphere
    /// \return OBB containing the sphere
    ///
    static OBB fromSphere(const Sphere& sphere);

    /// \brief Create an OBB that that contains the line segment
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// // segment points
    /// vec3 segment_a, segment_b;
    ///
    /// OBB obb = OBB::fromLineSegment( segment_a, segment_b );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param a Point 'a' from the line segment
    /// \param b Point 'b' from the line segment
    /// \return OBB containing the line segment
    ///
    static OBB fromLineSegment(const vec3& a, const vec3& b);

    /// \brief Create an OBB that that contains the line segment
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// LineSegment lineSegment;
    ///
    /// OBB obb = OBB::fromLineSegment( lineSegment );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param ls the line segment
    /// \return OBB containing the line segment
    ///
    static OBB fromLineSegment(const LineSegment& ls);

    /// \brief Create an OBB that that contains the frustum
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// Frustum frustum;
    ///
    /// OBB obb = OBB::fromFrustum( frustum );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param frustum the frustum
    /// \return OBB containing the line segment
    ///
    static OBB fromFrustum(const Frustum& frustum);
    
    /// \brief Raycast test against an OBB
    ///
    /// Intersect ray R(t) = p + t*d, |d| = 1, against OBB a.
    ///
    /// When intersecting it returns the intersection distance tmin and normal of intersection.
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// Ray ray;
    /// OBB obb;
    ///
    /// float tmin;
    /// vec3 normal;
    /// if ( OBB::raycastOBB(ray, obb, &tmin, &normal) ) {
    ///     vec3 collision_ptn = ray.origin + ray.dir * tmin;
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param r The Ray
    /// \param a The OBB
    /// \param outTmin **return** the distance from ray origin
    /// \param outNormal **return** the surface normal
    /// \return true, if the ray intersects the obb
    ///
    static bool raycastOBB(const Ray &r, const OBB& a, float *outTmin, vec3 *outNormal);

    /// \brief Test if a line segment intersects the OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// //line segment a and b
    /// vec3 a, b;
    /// OBB obb;
    ///
    /// if ( OBB::segmentIntersectsOBB( a, b, obb ) ) {
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param p0 line segment point a
    /// \param p1 line segment point b
    /// \param obb The OBB
    /// \return true, if the obb overlaps the line segment
    ///
    static bool segmentIntersectsOBB(const vec3& p0, const vec3& p1, const OBB &obb);

    /// \brief Test if a line segment intersects the OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// LineSegment ls;
    /// OBB obb;
    ///
    /// if ( OBB::segmentIntersectsOBB( ls, obb ) ) {
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param ls line segment
    /// \param obb The OBB
    /// \return true, if the obb overlaps the line segment
    ///
    static bool segmentIntersectsOBB(const LineSegment& ls, const OBB &obb);

    /// \brief Compute the point in the OBB surface that is closer to another specified point
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// vec3 ptn_to_input;
    /// OBB obb;
    ///
    /// vec3 closest_point = OBB::closestPointToOBB( ptn_to_input, obb );
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param p The input point to compute closest point
    /// \param obb The OBB
    /// \return The closest point in the OBB related to p
    ///
    static vec3 closestPointToOBB(const vec3 &p, const OBB &obb);


    /// \brief Test if there is some overlaped area between AABB and OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// AABB aabb_a;
    /// OBB obb_b;
    ///
    /// if ( OBB::aabbOverlapsAABB( aabb_a, obb_b ) ){
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param a The first AABB to test
    /// \param b The second OBB to test against
    /// \return true if there are some overlaped area between the AABB and OBB, otherwise false
    ///
    static bool aabbOverlapsOBB(const AABB& a,const OBB& b);

    /// \brief Test if a sphere overlaps the OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// vec3 sphere_center;
    /// float sphere_radius;
    /// OBB obb;
    ///
    /// vec3 penetration;
    /// if ( OBB::sphereOverlapsOBB( sphere_center, sphere_radius, obb, &penetration ) ) {
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param center The sphere center
    /// \param radius The sphere radius
    /// \param obb The OBB
    /// \param penetration **return** the amount one shape is inside the other
    /// \return true, if the obb overlaps the sphere
    ///
    static bool sphereOverlapsOBB(const vec3 &center, const float &radius, const OBB& obb, vec3 *penetration);

    /// \brief Test if a sphere overlaps the OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// Sphere sphere;
    /// OBB obb;
    ///
    /// vec3 penetration;
    /// if ( OBB::sphereOverlapsOBB( sphere, obb, &penetration ) ) {
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param sphere The Sphere
    /// \param obb The OBB
    /// \param penetration **return** the amount one shape is inside the other
    /// \return true, if the obb overlaps the sphere
    ///
    static bool sphereOverlapsOBB(const Sphere &sphere, const OBB& obb, vec3 *penetration);

    /// \brief Test if a plane intersects the OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// Plane plane;
    /// OBB obb;
    ///
    /// if ( OBB::planeIntersectsOBB( plane, obb ) ) {
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param plane The Plane
    /// \param obb The OBB
    /// \return true, if the obb intersects the plane
    ///
    static bool planeIntersectsOBB(const Plane &plane, const OBB &obb);

    /// \brief Test if a triangle intersects the OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// // triangle vertex a, b, c
    /// vec3 a, b, c;
    /// OBB obb;
    ///
    /// if ( OBB::triangleIntersectsOBB( a, b, c, obb ) ) {
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param v0 The triangle vertex a
    /// \param v1 The triangle vertex b
    /// \param v2 The triangle vertex c
    /// \param obb The OBB
    /// \return true, if the obb intersects the triangle
    ///
    static bool triangleIntersectsOBB(const vec3 &v0, const vec3 &v1, const vec3 &v2, const OBB &obb);

    /// \brief Test if a triangle overlaps the OBB
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// Triangle triangle;
    /// OBB obb;
    ///
    /// if ( OBB::triangleIntersectsOBB( triangle, obb ) ) {
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param t The triangle
    /// \param obb The OBB
    /// \return true, if the obb overlaps the triangle
    ///
    static bool triangleIntersectsOBB(const Triangle &t, const OBB &obb);

    //
    // Cloned methods from other collision classes
    //

    /// \brief Test if a frustum overlaps the obb
    ///
    /// Example:
    ///
    /// \code
    /// #include <aRibeiroCore/aRibeiroCore.h>
    /// using namespace aRibeiro;
    /// using namespace aRibeiro::collision;
    ///
    /// Frustum frustum;
    /// OBB obb;
    ///
    /// if ( OBB::frustumOverlapsOBB( frustum, obb ) ) {
    ///     ...
    /// }
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    /// \param f The frustum
    /// \param obb The OBB
    /// \return true, if the obb overlaps the frustum
    ///
    static bool frustumOverlapsOBB(const Frustum &f, const OBB &obb);

    SSE2_CLASS_NEW_OPERATOR
} _SSE2_ALIGN_POS;

}

}

#endif
