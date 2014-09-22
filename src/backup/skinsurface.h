#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATOR_SKINSUBSURFACE_H
#define PBRT_INTEGRATOR_SKINSUBSURFACE_H

// integrators/skinsubsurface.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
#include "renderers/surfacepoints.h"

struct SkinSubsurfaceOctreeNode;
struct SkinIrradiancePoint {
    SkinIrradiancePoint() { }
    SkinIrradiancePoint(const SurfacePoint &sp, const Spectrum &ee)
        : p(sp.p), n(sp.n), E(ee), area(sp.area),
          rayEpsilon(sp.rayEpsilon) { }
    Point p;
    Normal n;
    Spectrum E;
    float area, rayEpsilon;
};


// SkinSubsurfaceIntegrator Declarations
class SkinSubsurfaceIntegrator : public SurfaceIntegrator {
public:
    // SkinSubsurfaceIntegrator Public Methods
    SkinSubsurfaceIntegrator(float et, int mdepth, float merror, float mindist, const string &fn) 
        : maxSpecularDepth(mdepth), maxError(merror), minSampleDist(mindist), filename(fn), epidermisT(et)
    {
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
    vector<SkinIrradiancePoint> irradiancePoints;
    BBox octreeBounds;
    SkinSubsurfaceOctreeNode *octree;
    MemoryArena octreeArena;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;

    float epidermisT;
};

SkinSubsurfaceIntegrator *CreateSkinSubsurfaceIntegrator(const ParamSet &params);

#endif