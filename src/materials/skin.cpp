#include "stdafx.h"
#include "materials/skin.h"
#include "textures/constant.h"
#include "volume.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"

#include <iostream>

SkinMicrofacet::SkinMicrofacet(const Spectrum &reflectance, Fresnel *f, MicrofacetDistribution *d, Reference<Texture<Spectrum> > kdref) : Microfacet(reflectance, f, d), Kd(kdref)
{
}

Spectrum SkinMicrofacet::f(const Vector &wo,
			   const Vector &wi) const
{
  return Microfacet::f(wo, wi);
}

Spectrum SkinMicrofacet::Sample_f(const Vector &wo, Vector *wi,
				  float u1, float u2,
				  float *pdf) const
{
  return Microfacet::Sample_f(wo, wi, u1, u2, pdf);
}

float SkinMicrofacet::Pdf(const Vector &wo,
			  const Vector &wi) const
{
  // PDF is calculated only from the microfacet distribution.
  return Microfacet::Pdf(wo, wi);
}

// SkinMaterial Method Definitions
BSDF *SkinMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    Spectrum R = Kr->Evaluate(dgs).Clamp();
    float e = eta->Evaluate(dgs);

    if (!R.IsBlack()) {
      MicrofacetDistribution *md = BSDF_ALLOC(arena, Beckmann)(0.35f);
      Fresnel *fr = BSDF_ALLOC(arena, FresnelDielectric)(1., e);
      bsdf->Add(BSDF_ALLOC(arena, SkinMicrofacet)(R, fr, md, Kr));
    }
    
    return bsdf;
}

BSSRDF *SkinMaterial::GetBSSRDF(const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    float e = eta->Evaluate(dgShading);
    return BSDF_ALLOC(arena, BSSRDF)(scale * sigma_a->Evaluate(dgShading),
        scale * sigma_prime_s->Evaluate(dgShading), e);
}


SkinMaterial *CreateSkinMaterial(const Transform &xform,
        const TextureParams &mp) {
    float sa_rgb[3] = { .0011f, .0024f, .014f }, sps_rgb[3] = { 2.55f, 3.21f, 3.77f };
    Spectrum sa = Spectrum::FromRGB(sa_rgb), sps = Spectrum::FromRGB(sps_rgb);
    string name = mp.FindString("name");
    bool found = GetVolumeScatteringProperties(name, &sa, &sps);
    if (name != "" && !found)
        Warning("Named material \"%s\" not found.  Using defaults.", name.c_str());
    float scale = mp.FindFloat("scale", 1.f);

    Reference<Texture<Spectrum> > sigma_a, sigma_prime_s;
    sigma_a = mp.GetSpectrumTexture("sigma_a", sa);
    sigma_prime_s = mp.GetSpectrumTexture("sigma_prime_s", sps);
    Reference<Texture<float> > ior = mp.GetFloatTexture("index", 1.44f);
    Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");

    // Add a texture
    Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.0f));
    return new SkinMaterial(scale, Kr, sigma_a, sigma_prime_s, ior, bumpMap, Kd);
}
