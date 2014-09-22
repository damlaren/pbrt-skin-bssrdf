/*
 * A material for the tastiest part of the human body.
 * Based on "A Spectral BSSRDF for Shading Human Skin" [Donner & Jensen 2006].
 *
 * This material is meant for use with the SkinSubsurfaceIntegrator.
 */


#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_MATERIALS_SKINDJ_H
#define PBRT_MATERIALS_SKINDJ_H

// materials/subsurface.h*
#include "pbrt.h"
#include "material.h"
#include "core/reflection.h"
#include "materials/beckmann.h"

// The BRDF defined in DF2006 is a Torrance-Sparrowe micofacet
class SkinDJMicrofacet : public Microfacet {
 public:
  SkinDJMicrofacet(const Spectrum &reflectance, Fresnel *f, MicrofacetDistribution *d);
  Spectrum f(const Vector &wo, const Vector &wi) const;
  Spectrum Sample_f(const Vector &wo, Vector *wi,
		    float u1, float u2, float *pdf) const;
  float Pdf(const Vector &wo, const Vector &wi) const;
};

// SubsurfaceMaterial Declarations
class SkinDJMaterial : public Material {
public:
    // SkinDJMaterial Public Methods
    SkinDJMaterial(float sc,
		   Reference<Texture<float> > e,
		   Reference<Texture<float> > bump,
		   Reference<Texture<Spectrum> > tex) {
        scale = sc;
        eta = e;
        bumpMap = bump;
	Kd = tex;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
    BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom,
                      const DifferentialGeometry &dgShading,
                      MemoryArena &arena) const;

    // SkinDJMaterial Private Data
    float scale;
    Reference<Texture<float> > eta, bumpMap;
    Reference<Texture<Spectrum> > Kd;
};


SkinDJMaterial *CreateSkinDJMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_SKIN_H
