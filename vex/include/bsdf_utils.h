#ifndef __bsdf_utils_h__
#define __bsdf_utils_h__

vector SphericalDirection(float sintheta;
                          float costheta;
                          float phi) {
    return set(sintheta*cos(phi), sintheta*sin(phi), costheta);
}

vector SphericalDirection(float costheta;
                          float phi) {
    float sintheta = sqrt(max(0.0,1.0-costheta*costheta));
    return set(sintheta*cos(phi), sintheta*sin(phi), costheta);
}

float SphericalTheta(vector w) {
    return acos(clamp(w.z, -1.0, 1.0));
}

float SphericalPhi(vector w) {
    float p = atan2(w.y, w.x);
    return (p < 0.0) ? p + 6.28318530717958647692  : p;
}

float CosTheta    (vector w) { return w.z; }
float AbsCosTheta (vector w) { return abs(w.z); }
float CosThetaSqr (vector w) { return w.z*w.z; }
float SinTheta    (vector w) { return sqrt(SinThetaSqr(w)); }
float SinThetaSqr (vector w) { return max(0.0, 1.0 - w.z*w.z); }
float CosPhi      (vector w) { return w.x / SinTheta(w); }
float CosPhiSqr   (vector w) { float cph = CosPhi(w); return cph*cph; }
float SinPhi      (vector w) { return w.y / SinTheta(w); }
float SinPhiSqr   (vector w) { float sph = SinPhi(w); return sph*sph; }

float tanTheta(vector w) {
    float tmp = 1.0 - w.z*w.z;
    if (tmp <= 0.0)
        return 0.0;
    return sqrt(tmp)/w.z;
}
float tanThetaSqr(vector w) {
    float wSqr = w.z*w.z;
    float tmp = 1.0 - wSqr;
    if (tmp <= 0.0)
        return 0.0;
    return tmp/wSqr;
}

void ConcentricSampleDisk(float u1; 
                          float u2; 
                          export float dx; 
                          export float dy) {
    float r, theta;

    float sx = 2 * u1 - 1;
    float sy = 2 * u2 - 1;

    if (sx == 0.0 && sy == 0.0) {
        dx = 0.0;
        dy = 0.0;
        return;
    }
    if (sx >= -sy) {
        if (sx > sy) {

            r = sx;
            if (sy > 0.0) 
                theta = sy/r;
            else
                theta = 8.0f + sy/r;
        }
        else {

            r = sy;
            theta = 2.0f - sx/r;
        }
    }
    else {
        if (sx <= sy) {

            r = -sx;
            theta = 4.0f - sy/r;
        }
        else {

            r = -sy;
            theta = 6.0f + sx/r;
        }
    }
    theta *= 3.14159265358979323846  / 4.f;
    dx = r * cos(theta);
    dy = r * sin(theta);
}

struct UniformHemisphereSampler {
    vector sample(float u1; float u2;                      ) {
        float z = u1;
        float r = sqrt(max(0.f, 1.f - z*z));
        float phi = 2. * 3.14159265358979323846  * u2;
        float x = r * cos(phi);
        float y = r * sin(phi);
        return set(x,y,z);
    }
    float pdf() {
        return 0.5;
    }
    float pdf(vector w) {
        return 0.5;
    }
}

struct CosHemisphereSampler {
    float s_pdf = 0.0;
    vector sample(float u1; float u2;) {
        vector s = set(cos(u1*3.14159265358979323846 *2.),
                       sin(u1*3.14159265358979323846 *2.),
                       0.0);
        s *= sqrt(u2);
        s.z = sqrt(1. - u2);
        this.s_pdf = s.z *  2.0;
        return s;
    }
    float pdf() {
        return this.s_pdf * 2.0;
    }
    float pdf(float costheta) {

        return costheta * 2.0;
    }
}

struct CosHemisphereDiskSampler {
    float s_pdf = 0.0;
    vector sample(float u1; float u2;) {
        float s1 = 0;
        float s2 = 0;
        ConcentricSampleDisk(u1,u2,s1,s2);
        vector ret = set(s1, s2, sqrt(max(0., 1.f-s1*s1-s2*s2)));
        this.s_pdf = ret.z/3.14159265358979323846 ;
        return ret;
    }
    float pdf() {
        return this.s_pdf * 2.0;
    }
    float pdf(float costheta) {
        return costheta * 2.0;
    }
} 
#endif
