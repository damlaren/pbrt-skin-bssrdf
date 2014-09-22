#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PBDSUBSURFACE_H
#define PBRT_INTEGRATORS_PBDSUBSURFACE_H

// integrators/pbdsubsurface.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
#include "renderers/surfacepoints.h"

struct PBDSubsurfaceOctreeNode;

// PBDSubsurfaceIntegrator Helper Declarations
struct PBDIrradiancePoint {
    PBDIrradiancePoint() { }
    PBDIrradiancePoint(const SurfacePoint &sp, const Spectrum &ee)
    : p(sp.p), n(sp.n), E(ee), area(sp.area),
    rayEpsilon(sp.rayEpsilon) { }
    Point p;
    Normal n;
    Spectrum E;
    float area, rayEpsilon;
};



// PBDSubsurfaceIntegrator Declarations
class PBDSubsurfaceIntegrator : public SurfaceIntegrator {
public:
    // PBDSubsurfaceIntegrator Public Methods
    PBDSubsurfaceIntegrator(int mdepth, float merror, float mindist,
                            const string &fn) {
        maxSpecularDepth = mdepth;
        maxError = merror;
        minSampleDist = mindist;
        filename = fn;
        octree = NULL;
    }
    PBDSubsurfaceIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
                const RayDifferential &ray, const Intersection &isect, const Sample *sample,
                RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *, const Camera *, const Renderer *);
private:
    // PBDSubsurfaceIntegrator Private Data
    int maxSpecularDepth;
    float maxError, minSampleDist;
    string filename;
    vector<PBDIrradiancePoint> irradiancePoints;
    BBox octreeBounds;
    PBDSubsurfaceOctreeNode *octree;
    MemoryArena octreeArena;
    
    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
};


PBDSubsurfaceIntegrator *CreatePBDSubsurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_PBDSUBSURFACE_H
