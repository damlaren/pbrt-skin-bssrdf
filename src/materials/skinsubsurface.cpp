#include "stdafx.h"
#include "materials/skinsubsurface.h"
#include "textures/constant.h"
#include "volume.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"
#include "iostream"

SkinSubsurfaceMaterial::SkinSubsurfaceMaterial(float C_h, float C_m, float beta, float rho,
					      Reference<Texture<float> > e, Reference<Texture<float> > bump,
					      Reference<Texture<Spectrum> > kd)
  : eta(e) , bumpMap(bump), Kd(kd)
{
    //wavelengths in nm
    float lambdaR = 650.f, lambdaG = 532.f, lambdaB = 473.f;
    float sa_em[3] = { 6.6e10 * pow(lambdaR,-3.33),  6.6e10 * pow(lambdaG,-3.33),  6.6e10 * pow(lambdaB,-3.33)};
    float sa_pm[3] = { 2.9e14 * pow(lambdaR,-4.75),  2.9e14 * pow(lambdaG,-4.75),  2.9e14 * pow(lambdaB,-4.75) };
    float sa_ba[3] = { 0.0244 + 8.53*expf((-lambdaR - 154.f)/66.2),
                        0.0244 + 8.53*expf((-lambdaG - 154.f)/66.2),
                        0.0244 + 8.53*expf((-lambdaB - 154.f)/66.2) };
    float sps[3] = {14.74 * pow(lambdaR,-0.22) + 2.2e11 * pow(lambdaR,-4),
                    14.74 * pow(lambdaG,-0.22) + 2.2e11 * pow(lambdaG,-4),
                    14.74 * pow(lambdaB,-0.22) + 2.2e11 * pow(lambdaB,-4) };

    Spectrum spec_aEm = Spectrum::FromRGB(sa_em);
    Spectrum spec_aPm = Spectrum::FromRGB(sa_pm);
    Spectrum spec_aBase = Spectrum::FromRGB(sa_ba);
    Spectrum sa_epi = C_m * (spec_aEm * beta + (1 - beta) * spec_aPm) + (1 - C_m) * spec_aBase;
    Spectrum sps_epi = Spectrum::FromRGB(sps);

    rho_s = Spectrum(rho);
    episigma_a = new ConstantTexture<Spectrum>(sa_epi);
    episigmap_s = new ConstantTexture<Spectrum>(sps_epi);
}
        

BSDF *SkinSubsurfaceMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
                const DifferentialGeometry &dgShading,
                MemoryArena &arena) const {
     DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    
    Spectrum kd = Kd->Evaluate(dgs).Clamp();
    if (!kd.IsBlack()){
        BxDF *diff = BSDF_ALLOC(arena, Lambertian)(kd);
        bsdf->Add(diff);
    }
    
    float e = eta->Evaluate(dgs);
    Fresnel *fresnel = BSDF_ALLOC(arena,FresnelDielectric)(1.,e);
    float rough = 0.35f;
    BxDF *spec = BSDF_ALLOC(arena,Microfacet)(rho_s, fresnel, BSDF_ALLOC(arena, Beckmann)(rough));
    bsdf->Add(spec);
    return bsdf;
}

BSSRDF *SkinSubsurfaceMaterial::GetBSSRDF(const DifferentialGeometry &dgGeom,
                    const DifferentialGeometry &dgShading,
                    MemoryArena &arena) const {
    float e = eta->Evaluate(dgShading);
    return BSDF_ALLOC(arena, BSSRDF)(episigma_a->Evaluate(dgShading), episigmap_s->Evaluate(dgShading),e);
}

SkinSubsurfaceMaterial *CreateSkinSubsurfaceMaterial(const Transform &xform,
        const TextureParams &mp) {
    float kdDefault[3] = {0.9725,0.764,0.3764};
    Reference<Texture<float> > ior = mp.GetFloatTexture("index", 1.4f);
    float C_h = mp.FindFloat("Ch", 0.001f);
    float C_m = mp.FindFloat("Cm", 0.001f);
    float B = mp.FindFloat("beta", 0.25);
    float rhos = mp.FindFloat("rho", 0.5);
    Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum::FromRGB(kdDefault));
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");

    return new SkinSubsurfaceMaterial(C_h, C_m, B, rhos, ior, bumpMap, Kd);
}
