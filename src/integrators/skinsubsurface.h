
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

/* Implements multipole model from "A Spectral BSSRDF for
   Shading Human Skin" by Donner & Jensen 2006. */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_SKINSUBSURFACE_H
#define PBRT_INTEGRATORS_SKINSUBSURFACE_H

#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
#include "renderers/surfacepoints.h"
struct SkinSubsurfaceOctreeNode;


// SkinSubsurfaceIntegrator Declarations
class SkinSubsurfaceIntegrator : public SurfaceIntegrator {
public:
    // SkinSubsurfaceIntegrator Public Methods
    SkinSubsurfaceIntegrator(int mdepth, float merror, float mindist, float vCm, float vBm, float vCh, float vrho,
			     const string &fn)
    {
        maxSpecularDepth = mdepth;
        maxError = merror;
        minSampleDist = mindist;
        filename = fn;

	Cm = vCm;
	Bm = vBm;
	Ch = vCh;
	rho = vrho;
	
        octree = NULL;
    }
    ~SkinSubsurfaceIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect, const Sample *sample,
        RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *, const Camera *, const Renderer *);
private:
    // SkinSubsurfaceIntegrator Private Data
    int maxSpecularDepth;
    float maxError, minSampleDist;
    string filename;
    vector<IrradiancePoint> irradiancePoints;
    BBox octreeBounds;
    SkinSubsurfaceOctreeNode *octree;
    MemoryArena octreeArena;

    // Skin parameters
    float Cm;
    float Bm;
    float Ch;
    float rho;

    Spectrum totalDiffuseReflectance;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
};


SkinSubsurfaceIntegrator *CreateSkinSubsurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_SKINSUBSURFACE_H
