/*
 * A material for the tastiest part of the human body.
 * Based on "A Spectral BSSRDF for Shading Human Skin" [Donner & Jensen 2006].
 *
 * The material is adapted from the Subsurface material, but has a different
 * BRDF and parameterization.
 */


#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_MATERIALS_SKIN_H
#define PBRT_MATERIALS_SKIN_H

// materials/subsurface.h*
#include "pbrt.h"
#include "material.h"
#include "core/reflection.h"
#include "materials/beckmann.h"
#include "materials/skindj.h"

// The BRDF defined in DF2006 is a Torrance-Sparrowe micofacet
// BRDF modulated by an additional "oiliness" parameter.
// We don't include the oiliness here, but fold it straight it into the integrator.
// The representation of skin parameters is very closely tied to the
// Skin BSSRDF model.
class SkinMicrofacet : public Microfacet {
 public:
  SkinMicrofacet(const Spectrum &reflectance, Fresnel *f, MicrofacetDistribution *d, Reference<Texture<Spectrum> > kdref);
  Spectrum f(const Vector &wo, const Vector &wi) const;
  Spectrum Sample_f(const Vector &wo, Vector *wi,
		    float u1, float u2, float *pdf) const;
  float Pdf(const Vector &wo, const Vector &wi) const;

  Reference<Texture<Spectrum> > Kd;
};

// SubsurfaceMaterial Declarations
class SkinMaterial : public Material {
public:
    // SkinMaterial Public Methods
    SkinMaterial(float sc, Reference<Texture<Spectrum> > kr,
		 Reference<Texture<Spectrum> > sa,
		 Reference<Texture<Spectrum> > sps,
		 Reference<Texture<float> > e,
		 Reference<Texture<float> > bump,
		 Reference<Texture<Spectrum> > kd) {
        scale = sc;
        Kr = kr;
        sigma_a = sa;
        sigma_prime_s = sps;
        eta = e;
        bumpMap = bump;
	Kd = kd;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
    BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom,
                      const DifferentialGeometry &dgShading,
                      MemoryArena &arena) const;
private:
    // SkinMaterial Private Data
    float scale;
    Reference<Texture<Spectrum> > Kr, sigma_a, sigma_prime_s, Kd;
    Reference<Texture<float> > eta, bumpMap;
};


SkinMaterial *CreateSkinMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_SKIN_H
