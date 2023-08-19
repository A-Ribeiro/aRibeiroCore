
#include "Random.h"

#include <aRibeiroCore/SetNullAndDelete.h>

#include <librandom/Random.hpp>

#include <time.h>

//RandomLib::Random _aribeiro_random;

namespace aRibeiro {

	RandomInstantiable RandomInstantiable_aribeiro_random;


	//bool Random::initialized = false;

	void Random::setSeed(long int v) {
		RandomInstantiable_aribeiro_random.setSeed(v);
	}

	double Random::getDouble() {
		return RandomInstantiable_aribeiro_random.getDouble();
	}
	float Random::getFloat() {
		return RandomInstantiable_aribeiro_random.getFloat();
	}
	vec2 Random::getVec2() {
		return RandomInstantiable_aribeiro_random.getVec2();
	}
	vec3 Random::getVec3() {
		return RandomInstantiable_aribeiro_random.getVec3();
	}
	vec4 Random::getVec4() {
		return RandomInstantiable_aribeiro_random.getVec4();
	}

	mat4 Random::getMat4() {
		return RandomInstantiable_aribeiro_random.getMat4();
	}
	quat Random::getQuat() {
		return RandomInstantiable_aribeiro_random.getQuat();
	}

	vec2 Random::getVec2Direction() {
		return RandomInstantiable_aribeiro_random.getVec2Direction();
	}
	vec3 Random::getVec3Direction() {
		return RandomInstantiable_aribeiro_random.getVec3Direction();
	}
	vec4 Random::getVec4Ptn() {
		return RandomInstantiable_aribeiro_random.getVec4Ptn();
	}
	vec4 Random::getVec4Vec() {
		return RandomInstantiable_aribeiro_random.getVec4Vec();
	}
	vec4 Random::getVec4VecDirection() {
		return RandomInstantiable_aribeiro_random.getVec4VecDirection();
	}
	int Random::getRange(int min, int max) {
		return RandomInstantiable_aribeiro_random.getRange(min, max);
	}

	vec3 Random::getVec3_uvw() {
		return RandomInstantiable_aribeiro_random.getVec3_uvw();
	}

	vec3 Random::getVec3PointInsideTriangle(const vec3& a, const vec3& b, const vec3& c) {
		return RandomInstantiable_aribeiro_random.getVec3PointInsideTriangle(a, b, c);
	}


	RandomInstantiable::RandomInstantiable() {
		_aribeiro_internal_type = (void*)new RandomLib::Random();
		this->setSeed(time(NULL));
	}

	RandomInstantiable::~RandomInstantiable() {
		RandomLib::Random* lib = (RandomLib::Random*)_aribeiro_internal_type;
		_aribeiro_internal_type = NULL;
		aRibeiro::setNullAndDelete(lib);
	}

	void RandomInstantiable::setSeed(long int v) {
		RandomLib::Random* lib = (RandomLib::Random*)_aribeiro_internal_type;
		lib->Reseed(v);
	}

	double RandomInstantiable::getDouble() {
		RandomLib::Random* lib = (RandomLib::Random*)_aribeiro_internal_type;
		return lib->FloatN();
	}

	float RandomInstantiable::getFloat() {
		RandomLib::Random* lib = (RandomLib::Random*)_aribeiro_internal_type;
		return (float)lib->FloatN();
	}

	vec2 RandomInstantiable::getVec2() {
		return vec2(this->getFloat(), this->getFloat());
	}

	vec3 RandomInstantiable::getVec3() {
		return vec3(this->getFloat(), this->getFloat(), this->getFloat());
	}

	vec4 RandomInstantiable::getVec4() {
		return vec4(this->getFloat(), this->getFloat(), this->getFloat(), this->getFloat());
	}

	mat4 RandomInstantiable::getMat4() {
		return mat4(this->getFloat(), this->getFloat(), this->getFloat(), this->getFloat(),
			this->getFloat(), this->getFloat(), this->getFloat(), this->getFloat(),
			this->getFloat(), this->getFloat(), this->getFloat(), this->getFloat(),
			this->getFloat(), this->getFloat(), this->getFloat(), this->getFloat());
	}

	quat RandomInstantiable::getQuat() {
		return quatFromEuler(DEG2RAD(this->getFloat() * 360.0f), DEG2RAD(this->getFloat() * 360.0f), DEG2RAD(this->getFloat() * 360.0f));
	}

	vec2 RandomInstantiable::getVec2Direction() {
		vec2 result = this->getVec2() * 2.0f - vec2(1.0f);
		while (sqrLength(result) < 0.001f) {
			result = this->getVec2() * 2.0f - vec2(1.0f);
		}
		result = normalize(result);
		return result;
	}

	vec3 RandomInstantiable::getVec3Direction() {
		vec3 result = this->getVec3() * 2.0f - vec3(1.0f);
		while (sqrLength(result) < 0.001f) {
			result = this->getVec3() * 2.0f - vec3(1.0f);
		}
		result = normalize(result);
		return result;
	}

	vec4 RandomInstantiable::getVec4Ptn() {
		return vec4(this->getVec3(), 1.0f);
	}

	vec4 RandomInstantiable::getVec4Vec() {
		return vec4(this->getVec3(), 0.0f);
	}

	vec4 RandomInstantiable::getVec4VecDirection() {
		return vec4(this->getVec3Direction(), 0.0f);
	}

	int RandomInstantiable::getRange(int min, int max) {
		int64_t delta = (int64_t)max - (int64_t)min + 1LL;
		double rand_f = this->getDouble() * (double)delta;
		int64_t rand_i = (int64_t)rand_f;
		if (rand_i >= delta)
			rand_i = delta - 1LL;
		return (int)((int64_t)min + rand_i);
	}

	vec3 RandomInstantiable::getVec3_uvw() {
		float u = this->getFloat();
		float v = this->getFloat();
		float w = 1.0f - u - v;
		if (w < 0) {
			u = 1.0f - u;
			v = 1.0f - v;
			w = -w;
		}
		return vec3(u, v, w);
	}

	vec3 RandomInstantiable::getVec3PointInsideTriangle(const vec3& a, const vec3& b, const vec3& c) {
		float u = this->getFloat();
		float v = this->getFloat();
		float w = 1.0f - u - v;
		if (w < 0) {
			u = 1.0f - u;
			v = 1.0f - v;
			w = -w;
		}
		return a * u + b * v + c * w;
	}


}

