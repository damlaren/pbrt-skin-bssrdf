#ifndef PBRT_MATERIALS_BECKMANN_H
#define PBRT_MATERIALS_BECKMANN_H

#include "core/reflection.h"

class Beckmann : public MicrofacetDistribution {
public:
    Beckmann(float rough) { 
        if (rough > 1.f || isnan(rough)) rough = 0.35f;
        m = rough;
    }
    // Blinn Public Methods
    float D(const Vector &wh) const {
        float costhetah = AbsCosTheta(wh);
        float tanToM = SinTheta(wh)/(costhetah * m);
        return (expf(-tanToM * tanToM)/(M_PI * m*m*costhetah*costhetah*costhetah*costhetah));
    }
    void Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const{
        float phi= u2 * 2.f * M_PI;

        float tan_theta = m * sqrt(-log(u1));
        float cos_theta = cos(atan(tan_theta));
        float sin_theta = sqrtf(max(0.f, 1.f - cos_theta*cos_theta));
        Vector wh = SphericalDirection(sin_theta,cos_theta,phi);
        if (!SameHemisphere(wo, wh)) wh = -wh;
        *wi = -wo + 2.f * Dot(wo, wh) * wh;
        float beckPdf =  (Dot(wo, wh) <= 0.f) ? 0.f : 
            u1/(4.f * M_PI * m*m*Dot(wo, wh)*pow(cos_theta,3));
        *pdf =  beckPdf;
    }
    float Pdf(const Vector &wo, const Vector &wi) const{
        Vector wh = Normalize(wo + wi);
        float costheta = AbsCosTheta(wh);
        float sintheta = SinTheta(wh);

        float tanToM = sintheta/(costheta * m);
        float beckPdf =  (Dot(wo, wh) <= 0.f) ? 0.f : 
            exp(-tanToM * tanToM)/(4*M_PI*m*m*pow(costheta,3)*Dot(wo,wh));
        return beckPdf;
    }
private:
    float m;
};

#endif
