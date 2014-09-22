#include "stdafx.h"
#include "materials/skindj.h"
#include "textures/constant.h"
#include "volume.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"

#include <iostream>

SkinDJMicrofacet::SkinDJMicrofacet(const Spectrum &reflectance, Fresnel *f, MicrofacetDistribution *d) : Microfacet(reflectance, f, d)
{
}

Spectrum SkinDJMicrofacet::f(const Vector &wo,
			     const Vector &wi) const
{
  return Microfacet::f(wo, wi);
}

Spectrum SkinDJMicrofacet::Sample_f(const Vector &wo,
				    Vector *wi,
				    float u1, float u2,
				    float *pdf) const
{
  return Microfacet::Sample_f(wo, wi, u1, u2, pdf);
}

float SkinDJMicrofacet::Pdf(const Vector &wo,
			    const Vector &wi) const
{
  // PDF is calculated only from the microfacet distribution.
  return Microfacet::Pdf(wo, wi);
}

// SkinMaterial Method Definitions
BSDF *SkinDJMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);

    // Integrator will add distribution with right color.
    
    return bsdf;
}

BSSRDF *SkinDJMaterial::GetBSSRDF(const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, MemoryArena &arena) const {
  //irrelevant, integrator overwrites it
    return BSDF_ALLOC(arena, BSSRDF)(Spectrum(), Spectrum(), 0);
}


SkinDJMaterial *CreateSkinDJMaterial(const Transform &xform,
        const TextureParams &mp) {
    float scale = mp.FindFloat("scale", 1.f);
    Reference<Texture<float> > ior = mp.GetFloatTexture("index", 1.4f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    Reference<Texture<Spectrum> > tex = mp.GetSpectrumTexture("skintex", Spectrum(0.5f * INV_PI));

    // Add a texture
    return new SkinDJMaterial(scale, ior, bumpMap, tex);
}
