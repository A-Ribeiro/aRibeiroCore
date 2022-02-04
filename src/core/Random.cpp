
#include "Random.h"
#include <librandom/Random.hpp>

#include <time.h>

RandomLib::Random _aribeiro_random;

namespace aRibeiro {


    bool Random::initialized = false;
    
    void Random::setSeed(long int v) {
        _aribeiro_random.Reseed(v);
        initialized = true;
    }

    double Random::getDouble() {
        if (!initialized)
            Random::setSeed(time(NULL));
        return _aribeiro_random.FloatN();
    }
    float Random::getFloat() {
        if (!initialized)
            Random::setSeed(time(NULL));
        return (float)_aribeiro_random.FloatN();
    }
    vec2 Random::getVec2() {
        return vec2(Random::getFloat(), Random::getFloat());
    }
    vec3 Random::getVec3() {
        return vec3(Random::getFloat(), Random::getFloat(), Random::getFloat());
    }
    vec4 Random::getVec4() {
        return vec4(Random::getFloat(), Random::getFloat(), Random::getFloat(), Random::getFloat());
    }
    
    mat4 Random::getMat4() {
        return mat4(Random::getFloat(), Random::getFloat(), Random::getFloat(), Random::getFloat(),
            Random::getFloat(), Random::getFloat(), Random::getFloat(), Random::getFloat(),
            Random::getFloat(), Random::getFloat(), Random::getFloat(), Random::getFloat(),
            Random::getFloat(), Random::getFloat(), Random::getFloat(), Random::getFloat());
    }
    quat Random::getQuat(){
        if (!initialized)
            Random::setSeed(time(NULL));
        return quatFromEuler( DEG2RAD( Random::getFloat() * 360.0f ) , DEG2RAD(Random::getFloat() * 360.0f), DEG2RAD(Random::getFloat() * 360.0f));
    }
    
    vec2 Random::getVec2Direction(){
        vec2 result = Random::getVec2() * 2.0f - vec2(1.0f);
        while (sqrLength(result) < 0.001f){
            result = Random::getVec2() * 2.0f - vec2(1.0f);
        }
        result = normalize(result);
        return result;
    }
    vec3 Random::getVec3Direction(){
        vec3 result = Random::getVec3() * 2.0f - vec3(1.0f);
        while (sqrLength(result) < 0.001f){
            result = Random::getVec3() * 2.0f - vec3(1.0f);
        }
        result = normalize(result);
        return result;
    }
    vec4 Random::getVec4Ptn(){
        return vec4(Random::getVec3(), 1.0f);
    }
    vec4 Random::getVec4Vec(){
        return vec4(Random::getVec3(), 0.0f);
    }
    vec4 Random::getVec4VecDirection(){
        return vec4(Random::getVec3Direction(), 0.0f);
    }
    int Random::getRange(int min, int max) {
        int64_t delta = (int64_t)max - (int64_t)min + 1LL;
        double rand_f = Random::getDouble() * (double)delta;
        int64_t rand_i = (int64_t)rand_f;
        if (rand_i >= delta)
            rand_i = delta - 1LL;
        return (int)((int64_t)min + rand_i);
    }

    vec3 Random::getVec3_uvw() {
        float u = Random::getFloat();
        float v = Random::getFloat();
        float w = 1.0f - u - v;
        if (w < 0) {
            u = 1.0f - u;
            v = 1.0f - v;
            w = - w;
        }
        return vec3(u,v,w);
    }

    vec3 Random::getVec3PointInsideTriangle(const vec3 &a, const vec3 &b, const vec3 &c) {
        float u = Random::getFloat();
        float v = Random::getFloat();
        float w = 1.0f - u - v;
        if (w < 0) {
            u = 1.0f - u;
            v = 1.0f - v;
            w = - w;
        }
        return a * u + b * v + c * w;
    }

}

