/*
 * An attempt to simulate the appearance of oil and the roughness of human skin, based on:
 * http://graphics.ucsd.edu/~henrik/papers/skin_bssrdf/skin_bssrdf.pdf [DJ2006]
 * http://graphics.ucsd.edu/papers/layered/layered.pdf [DJ2005]
 * "A Reflectance Model for Computer Graphics" [CT82]
 */

#ifndef WET_H
#define WET_H

#include "core/reflection.h"
#include "integrator.h"
#include "pbrt.h"

class BSDF;

class Wet {
 public:

  RNG rng;

  float oiliness; //*SIGH* total hack, get this from BSDF when first used...

  // Integrate a BRDF over a hemisphere as in [DJ2006], eq. 2.
  // This is different from the rho calculations 
  Spectrum integrate_BRDF(BSDF *bsdf, const Vector& wi,
			  int sqrtSamples, BxDFType bxdfType);
};

#endif
