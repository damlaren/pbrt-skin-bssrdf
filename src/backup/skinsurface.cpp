// integrators/skinsubsurface.cpp*

#include "stdafx.h"
#include "integrators/skinsurface.h"
#include "integrators/pbdsubsurface.cpp"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "reflection.h"
#include "octree.h"
#include "camera.h"
#include "floatfile.h"

struct SkinDiffusionReflectance;

// DipoleSubsurfaceIntegrator Local Declarations
struct SkinSubsurfaceOctreeNode {
    // SubsurfaceOctreeNode Methods
    SkinSubsurfaceOctreeNode() {
        isLeaf = true;
        sumArea = 0.f;
        for (int i = 0; i < 8; ++i)
            ips[i] = NULL;
    }
    void Insert(const BBox &nodeBound, SkinIrradiancePoint *ip,
                MemoryArena &arena) {
        Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
        if (isLeaf) {
            // Add _IrradiancePoint_ to leaf octree node
            for (int i = 0; i < 8; ++i) {
                if (!ips[i]) {
                    ips[i] = ip;
                    return;
                }
            }

            // Convert leaf node to interior node, redistribute points
            isLeaf = false;
            SkinIrradiancePoint *localIps[8];
            for (int i = 0; i < 8; ++i) {
                localIps[i] = ips[i];
                children[i] = NULL;
            }
            for (int i = 0; i < 8; ++i)  {
                SkinIrradiancePoint *ip = localIps[i];
                // Add _IrradiancePoint_ _ip_ to interior octree node
                int child = (ip->p.x > pMid.x ? 4 : 0) +
                    (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
                if (!children[child])
                    children[child] = arena.Alloc<SkinSubsurfaceOctreeNode>();
                BBox childBound = octreeChildBound(child, nodeBound, pMid);
                children[child]->Insert(childBound, ip, arena);
            }
            /* fall through to interior case to insert the new point... */
        }
        // Add _IrradiancePoint_ _ip_ to interior octree node
        int child = (ip->p.x > pMid.x ? 4 : 0) +
            (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
        if (!children[child])
            children[child] = arena.Alloc<SkinSubsurfaceOctreeNode>();
        BBox childBound = octreeChildBound(child, nodeBound, pMid);
        children[child]->Insert(childBound, ip, arena);
    }
    void InitHierarchy() {
        if (isLeaf) {
            // Init node leaf from irradiance point
            float sumWt = 0.f;
            uint32_t i;
            for (i = 0; i < 8; ++i) {
                if (!ips[i]) break;
                float wt = ips[i]->E.y();
                E += ips[i]->E;
                p += wt * ips[i]->p;
                sumWt += wt;
                sumArea += ips[i]->area;
            }
            if (sumWt > 0.f) p /= sumWt;
            E /= i;
        }
        else {
            // Init interior node
            float sumWt = 0.f;
            uint32_t nChildren = 0;
            for (uint32_t i = 0; i < 8; ++i) {
                if (!children[i]) continue;
                ++nChildren;
                children[i]->InitHierarchy();
                float wt = children[i]->E.y();
                E += children[i]->E;
                p += wt * children[i]->p;
                sumWt += wt;
                sumArea += children[i]->sumArea;
            }
            if (sumWt > 0.f) p /= sumWt;
            E /= nChildren;
        }
    }
    Spectrum Mo(const BBox &nodeBound, const Point &p, const SkinDiffusionReflectance &Rd,
                float maxError);

    // SubsurfaceOctreeNode Public Data
    Point p;
    bool isLeaf;
    Spectrum E;
    float sumArea;
    union {
        SkinSubsurfaceOctreeNode *children[8];
        SkinIrradiancePoint *ips[8];
    };
};


struct SkinDiffusionReflectance {
    // SkinDiffusionReflectance Public Methods

    //Takes epidermis scattering values
    SkinDiffusionReflectance(const Spectrum &sigma_a, const Spectrum &sigmap_s, float eta, float thickness) {
        int nSamples = 5;
        
        A0 = (1.f + threeC2(eta))/(1.f - twoC1(eta));
        Cphi = 0.25 * (1.f - twoC1(eta));
        Cphi_corrective = (1.f - twoC1(1/eta));
        C_e = 0.5 * (1.f - threeC2(eta));
        
        //Epidermis
        Ad = (1.f + threeC2(1.f))/(1.f - twoC1(1.f));
        sigmap_t_epi = sigma_a + sigmap_s;
        sigma_tr_epi = Sqrt(3.f * sigma_a * sigmap_t_epi);
        alphap_epi = sigmap_s / sigmap_t_epi;
        D_epi = (2 * sigma_a + sigmap_t_epi)/(3.f*sigmap_t_epi*sigmap_t_epi);

        //dermis - todo: incorporate into material?
        Spectrum sigmaps_derm = sigmap_s/2.f;
        float rgb_derm[3] = { 0.125f, 0.595f, 12.65f };
        Spectrum sigmaa_derm = Spectrum::FromRGB(rgb_derm);
        sigmap_t_derm = sigmaa_derm + sigmaps_derm;
        sigma_tr_derm = Sqrt(3.f * sigmaa_derm * sigmap_t_derm);
        alphap_derm = sigmaps_derm / sigmap_t_derm;
        D_derm = (2 * sigma_a + sigmap_t_epi)/(3.f*sigmap_t_epi*sigmap_t_epi);

        Spectrum l0 = Spectrum(1.f)/sigmap_t_epi;
        Spectrum ld = Spectrum(1.f)/sigmap_t_derm;
        zb0 = 2.f * A0 * eta * D_epi;
        zbd = 2.f * Ad * D_derm;

        zr1 = 2.f*(Spectrum(thickness) + 2.f*zb0) + l0;
        zv1 = 2.f*(Spectrum(thickness) + 2.f*zb0) - l0 - 2.f*zb0;
        zr2 = 4.f*(Spectrum(thickness) + 2.f*zbd) + ld;
        zv2 = 4.f*(Spectrum(thickness) + 2.f*zb0) - l0 - 2.f*zbd;
    }

    Spectrum operator()(float d2) const {
        
        Spectrum dr1 = Sqrt(Spectrum(d2) + zr1 * zr1);
        Spectrum dv1= Sqrt(Spectrum(d2) + zv1 * zv1);
        Spectrum dr2 = Sqrt(Spectrum(d2) + zr2 * zr2);
        Spectrum dv2= Sqrt(Spectrum(d2) + zv2 * zv2);

        Spectrum Rd;

        Rd /= Cphi_corrective;
        return Rd.Clamp();
    }

    // SkinDiffusionReflectance Data
    Spectrum sigmap_t_epi, sigma_tr_epi, alphap_epi;
    Spectrum sigmap_t_derm, sigma_tr_derm, alphap_derm;
    Spectrum D_epi, D_derm;
    Spectrum zb0, zbd, zr1, zr2, zv1, zv2;
    float Ad, A0, Cphi, Cphi_corrective, C_e;
};



// SkinSubsurfaceIntegrator Method Definitions
SkinSubsurfaceIntegrator::~SkinSubsurfaceIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}


void SkinSubsurfaceIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    // Allocate and request samples for sampling all lights
    uint32_t nLights = scene->lights.size();
    lightSampleOffsets = new LightSampleOffsets[nLights];
    bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
    for (uint32_t i = 0; i < nLights; ++i) {
        const Light *light = scene->lights[i];
        int nSamples = light->nSamples;
        if (sampler) nSamples = sampler->RoundSize(nSamples);
        lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
    }
}


void SkinSubsurfaceIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    vector<SurfacePoint> pts;
    // Get _SurfacePoint_s for translucent objects in scene
    if (filename != "") {
        // Initialize _SurfacePoint_s from file
        vector<float> fpts;
        if (ReadFloatFile(filename.c_str(), &fpts)) {
            if ((fpts.size() % 8) != 0)
                Error("Excess values (%d) in points file \"%s\"", int(fpts.size() % 8),
                      filename.c_str());
            for (u_int i = 0; i < fpts.size(); i += 8)
                pts.push_back(SurfacePoint(Point(fpts[i], fpts[i+1], fpts[i+2]),
                                           Normal(fpts[i+3], fpts[i+4], fpts[i+5]),
                                           fpts[i+6], fpts[i+7]));
        }
    }
    if (pts.size() == 0) {
        Point pCamera = camera->CameraToWorld(camera->shutterOpen,
                                              Point(0, 0, 0));
        FindPoissonPointDistribution(pCamera, camera->shutterOpen,
                                     minSampleDist, scene, &pts);
    }

    // Compute irradiance values at sample points
    RNG rng;
    MemoryArena arena;
    PBRT_SUBSURFACE_STARTED_COMPUTING_IRRADIANCE_VALUES();
    ProgressReporter progress(pts.size(), "Computing Irradiances");
    for (uint32_t i = 0; i < pts.size(); ++i) {
        SurfacePoint &sp = pts[i];
        Spectrum E(0.f);
        for (uint32_t j = 0; j < scene->lights.size(); ++j) {
            // Add irradiance from light at point
            const Light *light = scene->lights[j];
            Spectrum Elight = 0.f;
            int nSamples = RoundUpPow2(light->nSamples);
            uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
            uint32_t compScramble = rng.RandomUInt();
            for (int s = 0; s < nSamples; ++s) {
                float lpos[2];
                Sample02(s, scramble, lpos);
                float lcomp = VanDerCorput(s, compScramble);
                LightSample ls(lpos[0], lpos[1], lcomp);
                Vector wi;
                float lightPdf;
                VisibilityTester visibility;
                Spectrum Li = light->Sample_L(sp.p, sp.rayEpsilon,
                    ls, camera->shutterOpen, &wi, &lightPdf, &visibility);
                if (Dot(wi, sp.n) <= 0.) continue;
                if (Li.IsBlack() || lightPdf == 0.f) continue;
                Li *= visibility.Transmittance(scene, renderer, NULL, rng, arena);
                if (visibility.Unoccluded(scene))
                    Elight += Li * AbsDot(wi, sp.n) / lightPdf;
            }
            E += Elight / nSamples;
        }
        irradiancePoints.push_back(SkinIrradiancePoint(sp, E));
        PBRT_SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(&sp, &E);
        arena.FreeAll();
        progress.Update();
    }
    progress.Done();
    PBRT_SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES();

    // Create octree of clustered irradiance samples
    octree = octreeArena.Alloc<SkinSubsurfaceOctreeNode>();
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octreeBounds = Union(octreeBounds, irradiancePoints[i].p);
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octree->Insert(octreeBounds, &irradiancePoints[i], octreeArena);
    octree->InitHierarchy();
}


Spectrum SkinSubsurfaceIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    // Evaluate BSSRDF and possibly compute subsurface scattering
    BSSRDF *bssrdf = isect.GetBSSRDF(ray, arena);
    if (bssrdf && octree) {
        Spectrum sigma_a  = bssrdf->sigma_a();
        Spectrum sigmap_s = bssrdf->sigma_prime_s();
        Spectrum sigmap_t = sigmap_s + sigma_a;
        if (!sigmap_t.IsBlack()) {
            // Use hierarchical integration to evaluate reflection from dipole model
            PBRT_SUBSURFACE_STARTED_OCTREE_LOOKUP(const_cast<Point *>(&p));
            SkinDiffusionReflectance Rd(sigma_a, sigmap_s, bssrdf->eta(),epidermisT);
            Spectrum Mo = octree->Mo(octreeBounds, p, Rd, maxError);
            FresnelDielectric fresnel(1.f, bssrdf->eta());
            Spectrum Ft = Spectrum(1.f) - fresnel.Evaluate(AbsDot(wo, n));
            float Fdt = 1.f - Fdr(bssrdf->eta());
            L += (INV_PI * Ft) * (Fdt * Mo);
            PBRT_SUBSURFACE_FINISHED_OCTREE_LOOKUP();
        }
    }
    L += UniformSampleAllLights(scene, renderer, arena, p, n,
        wo, isect.rayEpsilon, ray.time, bsdf, sample, rng, lightSampleOffsets,
        bsdfSampleOffsets);
    if (ray.depth < maxSpecularDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


Spectrum SkinSubsurfaceOctreeNode::Mo(const BBox &nodeBound, const Point &pt,
        const SkinDiffusionReflectance &Rd, float maxError) {
    // Compute $M_\roman{o}$ at node if error is low enough
    float dw = sumArea / DistanceSquared(pt, p);
    if (dw < maxError && !nodeBound.Inside(pt))
    {
        PBRT_SUBSURFACE_ADDED_INTERIOR_CONTRIBUTION(const_cast<SubsurfaceOctreeNode *>(this));
        return Rd(DistanceSquared(pt, p)) * E * sumArea;
    }

    // Otherwise compute $M_\roman{o}$ from points in leaf or recursively visit children
    Spectrum Mo = 0.f;
    if (isLeaf) {
        // Accumulate $M_\roman{o}$ from leaf node
        for (int i = 0; i < 8; ++i) {
            if (!ips[i]) break;
            PBRT_SUBSURFACE_ADDED_POINT_CONTRIBUTION(const_cast<IrradiancePoint *>(ips[i]));
            Mo += Rd(DistanceSquared(pt, ips[i]->p)) * ips[i]->E * ips[i]->area;
        }
    }
    else {
        // Recursively visit children nodes to compute $M_\roman{o}$
        Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
        for (int child = 0; child < 8; ++child) {
            if (!children[child]) continue;
            BBox childBound = octreeChildBound(child, nodeBound, pMid);
            Mo += children[child]->Mo(childBound, pt, Rd, maxError);
        }
    }
    return Mo;
}


SkinSubsurfaceIntegrator *CreateSkinSubsurfaceIntegrator(const ParamSet &params) {
    float epidermis = params.FindOneFloat("epidermis", 0.25f);
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float maxError = params.FindOneFloat("maxerror", .05f);
    float minDist = params.FindOneFloat("minsampledistance", .25f);
    string pointsfile = params.FindOneFilename("pointsfile", "");
    if (PbrtOptions.quickRender) { maxError *= 4.f; minDist *= 4.f; }
    return new SkinSubsurfaceIntegrator(epidermis, maxDepth, maxError, minDist, pointsfile);
}


