#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_MATERIALS_SKIN_SUBSURFACE_H
#define PBRT_MATERIALS_SKIN_SUBSURFACE_H

// materials/skin.h*
#include "pbrt.h"
#include "material.h"
#include "spectrum.h"
#include "textures/constant.h"
#include "materials/beckmann.h"

// SkinSubsurfaceMaterial Declarations
class SkinSubsurfaceMaterial : public Material {
public:
    // SkinSubsurfaceMaterial Public Methods
    SkinSubsurfaceMaterial( float C_h, float C_m, float beta, float rho, 
                Reference<Texture<float> > e, 
                Reference<Texture<float> > bump,
                Reference<Texture<Spectrum> > kd);
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
    BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom,
                      const DifferentialGeometry &dgShading,
                      MemoryArena &arena) const;
private:
    Reference<Texture<Spectrum> > Kd, episigma_a, episigmap_s;
    Spectrum rho_s;
    Reference<Texture<float> > eta, bumpMap;
};

SkinSubsurfaceMaterial *CreateSkinSubsurfaceMaterial(const Transform &xform, const TextureParams &mp);

#endif // PBRT_MATERIALS_SUBSURFACE_H
