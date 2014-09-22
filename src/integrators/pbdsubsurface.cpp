
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


// integrators/pbdsubsurface.cpp*
#include "stdafx.h"
#include "integrators/pbdsubsurface.h"
#include "integrators/wet.h"
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
struct PBDDiffusionReflectance;

// Deterministic random sample
float deterministic_uniform_sample(int i, int N) {
  assert(i > 0 && N > 0 && i <= N);
  float xi_i = (static_cast<float>(i) - 0.5f) / static_cast<float>(N);
  assert(xi_i >= 0 && xi_i <= 1.0f);
  return xi_i;
}

// Exponential sampling
Spectrum exponential_ti(int i, int N, Spectrum sigmap_t) {
  assert( sigmap_t != 0);
  float xi_i = deterministic_uniform_sample(i, N);
  return (Spectrum(-log(1 - xi_i)) / sigmap_t);
}

Spectrum exponential_pdf(Spectrum ti, Spectrum sigmap_t) {
  return (sigmap_t * Exp(-sigmap_t * ti));
}

inline float twoC1(float eta){
    if (eta < 1)
        return (0.919317 -
                3.4793 * eta +
                6.75335 * (eta * eta) -
                7.80989 * (eta * eta * eta) +
                4.98554 * (eta*eta*eta*eta) -
                1.36881*(eta*eta*eta*eta*eta));
    else
        return (-9.23372 +
                22.2272 * eta  -
                29.9292 * eta * eta +
                10.2291 * (eta*eta*eta) -
                2.54396 * (eta*eta*eta*eta) +
                0.254913 * (eta*eta*eta*eta*eta));
}

inline float threeC2(float eta){
    if (eta < 1)
        return (0.828421 -
                2.62051 * eta +
                3.36231 * (eta * eta) -
                1.95284 * (eta * eta * eta) +
                0.23649 * (eta*eta*eta*eta) -
                0.145787* (eta*eta*eta*eta*eta));
    else
        return (-1641.1 +
                135.926/(eta*eta*eta) -
                656.175/(eta*eta) +
                1376.53/eta +
                1213.67*eta -
                568.556 * (eta*eta) +
                164.798*(eta*eta*eta) -
                27.0181*(eta*eta*eta*eta) +
                1.91826*(eta*eta*eta*eta*eta));
}

// DipoleSubsurfaceIntegrator Local Declarations
struct PBDSubsurfaceOctreeNode {
    // SubsurfaceOctreeNode Methods
    PBDSubsurfaceOctreeNode() {
        isLeaf = true;
        sumArea = 0.f;
        for (int i = 0; i < 8; ++i)
            ips[i] = NULL;
    }
    void Insert(const BBox &nodeBound, PBDIrradiancePoint *ip,
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
            PBDIrradiancePoint *localIps[8];
            for (int i = 0; i < 8; ++i) {
                localIps[i] = ips[i];
                children[i] = NULL;
            }
            for (int i = 0; i < 8; ++i)  {
                PBDIrradiancePoint *ip = localIps[i];
                // Add _IrradiancePoint_ _ip_ to interior octree node
                int child = (ip->p.x > pMid.x ? 4 : 0) +
                (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
                if (!children[child])
                    children[child] = arena.Alloc<PBDSubsurfaceOctreeNode>();
                BBox childBound = octreeChildBound(child, nodeBound, pMid);
                children[child]->Insert(childBound, ip, arena);
            }
            /* fall through to interior case to insert the new point... */
        }
        // Add _IrradiancePoint_ _ip_ to interior octree node
        int child = (ip->p.x > pMid.x ? 4 : 0) +
        (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
        if (!children[child])
            children[child] = arena.Alloc<PBDSubsurfaceOctreeNode>();
        BBox childBound = octreeChildBound(child, nodeBound, pMid);
        children[child]->Insert(childBound, ip, arena);
    }
    void InitHierarchy() {
        if (isLeaf) {
            // Init _SubsurfaceOctreeNode_ leaf from _IrradiancePoint_s
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
            // Init interior _SubsurfaceOctreeNode_
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
    Spectrum Mo(const BBox &nodeBound, const Point &p, const PBDDiffusionReflectance &Rd,
                float maxError);
    
    // SubsurfaceOctreeNode Public Data
    Point p;
    bool isLeaf;
    Spectrum E;
    float sumArea;
    union {
        PBDSubsurfaceOctreeNode *children[8];
        PBDIrradiancePoint *ips[8];
    };
};


struct PBDDiffusionReflectance {
    // PBDDiffusionReflectance Public Methods
    PBDDiffusionReflectance(const Spectrum &sigma_a, const Spectrum &sigmap_s, float eta, Vector wo ) {
        sigmap_t = sigma_a + sigmap_s;
        alphap = sigmap_s / sigmap_t;
        
        for (int i = 1; i <= nSamples; i++){
            ti[i] = exponential_ti(i,nSamples,sigmap_t);
            pdf_ti[i] = exponential_pdf(ti[i],sigmap_t);
        }

        A = (1.f + threeC2(eta))/(1.f - twoC1(eta));
        Cphi = 0.25 * (1.f - twoC1(eta));
        Cphi_corrective = (1.f - twoC1(1/eta));
        C_e = 0.5 * (1.f - threeC2(eta));
        D_g = (sigma_a + sigmap_t)/(3.f*sigmap_t*sigmap_t);
        
        sigma_tr = Sqrt(sigma_a/D_g);
        zr = Spectrum(1.f) / sigmap_t;
        zv = zr + 4.f*A*D_g;

    }
    Spectrum operator()(float sqd) const {
        Spectrum Rd(0.f);
        for (int i = 1; i <= nSamples; i++){
            
            Spectrum Q = alphap * sigmap_t * Exp(-sigmap_t * ti[i]);
            
            Spectrum dr = Sqrt(Spectrum(sqd) + zr * zr);
            Spectrum dv = Sqrt(Spectrum(sqd) + zv * zv);
            
            Spectrum Rd_phi = (Cphi*alphap*alphap)/(4.f*M_PI*D_g) * (Exp(-sigma_tr * dr)/dr - Exp(-sigma_tr * dv)/dv);
            Spectrum Rd_E = (C_e*alphap*alphap)/(4.f*M_PI) * (zr*(dr * sigma_tr + Spectrum(1.f))* Exp(-sigma_tr * dr)/(dr*dr*dr) -
             (zv * (dv * sigma_tr + Spectrum(1.f)) * Exp(-sigma_tr * dv)/ (dv * dv * dv)));
            
            Spectrum kappa = Spectrum(1.f);// - Exp(-2.f * sigmap_t * (ti[i] + dr));
            Rd += (Rd_phi + Rd_E)* Q * kappa/pdf_ti[i];
        }
        Rd /= nSamples;
        return Rd.Clamp();
    }
    
    // DiffusionReflectance Data
    static const int nSamples = 5; //TODO: MAKE INPUT VAR
    Spectrum ti[nSamples];
    Spectrum pdf_ti[nSamples];
    Spectrum zr, zv, sigmap_t, sigma_tr, alphap, D_g;
    float A, Cphi, Cphi_corrective, C_e;

};


// PBDSubsurfaceIntegrator Method Definitions
PBDSubsurfaceIntegrator::PBDSubsurfaceIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}


void PBDSubsurfaceIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
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

void PBDSubsurfaceIntegrator::Preprocess(const Scene *scene,
                                         const Camera *camera, 
                                         const Renderer *renderer) {
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
        irradiancePoints.push_back(PBDIrradiancePoint(sp, E));
        PBRT_SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(&sp, &E);
        arena.FreeAll();
        progress.Update();
    }
    progress.Done();
    PBRT_SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES();
    
    // Create octree of clustered irradiance samples
    octree = octreeArena.Alloc<PBDSubsurfaceOctreeNode>();
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octreeBounds = Union(octreeBounds, irradiancePoints[i].p);
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octree->Insert(octreeBounds, &irradiancePoints[i], octreeArena);
    octree->InitHierarchy();
}


Spectrum PBDSubsurfaceIntegrator::Li(const Scene *scene, const Renderer *renderer,
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
            PBDDiffusionReflectance Rd(sigma_a, sigmap_s, bssrdf->eta(),wo);
            Spectrum Mo = octree->Mo(octreeBounds, p, Rd, maxError);
            FresnelDielectric fresnel(1.f, bssrdf->eta());
            Spectrum Ft = Spectrum(1.f) - fresnel.Evaluate(AbsDot(wo, n));
            float Fdt = 1.f - Fdr(bssrdf->eta());

	    // modulate SSS contribution by rho_dr
            L += (INV_PI * Ft) * (Fdt * Mo);
	    //Wet wet;
	    //Spectrum rho_dr = wet.integrate_BRDF(bsdf, ray.d, 10,
	    //BxDFType(BSDF_REFLECTION | BSDF_GLOSSY));
	//L += (INV_PI * Ft) * (Fdt * Mo) * (Spectrum(1.0f) - rho_dr);
	    //L += (INV_PI * Ft) * (Fdt * Mo) * (Spectrum(0.0f));

            //L += (INV_PI * Ft) * (Fdt * Mo);
            PBRT_SUBSURFACE_FINISHED_OCTREE_LOOKUP();
        }
    }
    L += UniformSampleAllLights(scene, renderer, arena, p, n,
            wo, isect.rayEpsilon, ray.time, bsdf, sample, rng, 
            lightSampleOffsets, bsdfSampleOffsets);
    if (ray.depth < maxSpecularDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


Spectrum PBDSubsurfaceOctreeNode::Mo(const BBox &nodeBound, const Point &pt,
                                     const PBDDiffusionReflectance &Rd, float maxError) {
    // Compute $M_\roman{o}$ at node if error is low enough
    float dw = sumArea / DistanceSquared(pt, p);
    if (dw < maxError && !nodeBound.Inside(pt))
    {
        PBRT_SUBSURFACE_ADDED_INTERIOR_CONTRIBUTION(const_cast<PBDSubsurfaceOctreeNode *>(this));
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

PBDSubsurfaceIntegrator *CreatePBDSubsurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float maxError = params.FindOneFloat("maxerror", .05f);
    float minDist = params.FindOneFloat("minsampledistance", .25f);
    string pointsfile = params.FindOneFilename("pointsfile", "");
    if (PbrtOptions.quickRender) { maxError *= 4.f; minDist *= 4.f; }
    return new PBDSubsurfaceIntegrator(maxDepth, maxError, minDist, pointsfile);
}


