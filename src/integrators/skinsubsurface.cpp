
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

#include "stdafx.h"
#include "integrators/skinsubsurface.h"
#include "scene.h"
#include "materials/beckmann.h"
#include "materials/skindj.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "reflection.h"
#include "octree.h"
#include "camera.h"
#include "floatfile.h"
#include "core/reflection.h"

#include <iostream> //sigh.. cleanup

// forward declarations
struct DiffusionReflectance;

// Sigma*

static const float OXY_WAVELENGTHS[] = {
  400,
  402,
  404,
  406,
  408,
  410,
  412,
  414,
  416,
  418,
  420,
  422,
  424,
  426,
  428,
  430,
  432,
  434,
  436,
  438,
  440,
  442,
  444,
  446,
  448,
  450,
  452,
  454,
  456,
  458,
  460,
  462,
  464,
  466,
  468,
  470,
  472,
  474,
  476,
  478,
  480,
  482,
  484,
  486,
  488,
  490,
  492,
  494,
  496,
  498,
  500,
  502,
  504,
  506,
  508,
  510,
  512,
  514,
  516,
  518,
  520,
  522,
  524,
  526,
  528,
  530,
  532,
  534,
  536,
  538,
  540,
  542,
  544,
  546,
  548,
  550,
  552,
  554,
  556,
  558,
  560,
  562,
  564,
  566,
  568,
  570,
  572,
  574,
  576,
  578,
  580,
  582,
  584,
  586,
  588,
  590,
  592,
  594,
  596,
  598,
  600,
  602,
  604,
  606,
  608,
  610,
  612,
  614,
  616,
  618,
  620,
  622,
  624,
  626,
  628,
  630,
  632,
  634,
  636,
  638,
  640,
  642,
  644,
  646,
  648,
  650,
  652,
  654,
  656,
  658,
  660,
  662,
  664,
  666,
  668,
  670,
  672,
  674,
  676,
  678,
  680,
  682,
  684,
  686,
  688,
  690,
  692,
  694,
  696,
  698,
  700
};
static const int N_OXY_WAVELENGTHS = sizeof(OXY_WAVELENGTHS) / sizeof(float);

Spectrum SigmaADeoxy(const Spectrum &l) {
  static const float VALUES[] = {
    119.5931833,
    126.4978986,
    135.699187,
    144.9004753,
    153.9025274,
    162.7931786,
    172.1058679,
    183.4880437,
    194.8702195,
    206.5630326,
    218.2815535,
    230.2357302,
    247.0101395,
    258.0645395,
    268.240586,
    283.1083256,
    295.7266233,
    295.7266233,
    292.9844465,
    268.6262047,
    221.3450791,
    194.544586,
    151.4217144,
    127.0527609,
    92.82696744,
    55.32127349,
    33.5488186,
    19.37197907,
    16.44170614,
    13.86427423,
    12.52660614,
    11.18893805,
    10.31572614,
    9.716731907,
    9.118594605,
    8.653067256,
    8.199751163,
    8.05964307,
    7.922748465,
    7.850123628,
    7.792709302,
    7.97009386,
    8.147478419,
    8.324862977,
    8.514673023,
    8.93564,
    9.356392744,
    9.777359721,
    10.19811247,
    10.65335665,
    11.17329907,
    11.69324149,
    12.21318391,
    12.73312633,
    13.25306874,
    13.80386065,
    14.42684893,
    15.04983721,
    15.67282549,
    16.29581377,
    16.91880205,
    17.59449153,
    18.4227146,
    19.25093767,
    20.07894651,
    20.90716958,
    21.73603535,
    22.54154977,
    23.34706419,
    24.15043628,
    24.95380837,
    25.78717302,
    26.62268,
    27.45818698,
    28.11588093,
    28.60647349,
    28.96424186,
    29.19989767,
    29.2106093,
    29.0092307,
    28.80785209,
    27.99805302,
    27.08542233,
    26.15136837,
    25.14447535,
    24.13972465,
    23.21209767,
    22.34231349,
    21.4725293,
    20.60253088,
    19.82722326,
    19.10761609,
    18.38800893,
    17.59470577,
    16.64329898,
    15.17002167,
    13.64125814,
    12.09064288,
    10.60451163,
    9.136161674,
    7.860835256,
    7.295904,
    6.730972744,
    6.166255721,
    5.611607628,
    5.057816465,
    4.601286884,
    4.157182791,
    3.933738233,
    3.710079442,
    3.486420651,
    3.316962698,
    3.163572186,
    3.009967442,
    2.874358233,
    2.757601488,
    2.640844744,
    2.533728465,
    2.464959814,
    2.39597693,
    2.327208279,
    2.258439628,
    2.189670977,
    2.123623079,
    2.066058791,
    2.008494502,
    1.950930214,
    1.893365926,
    1.835801637,
    1.778237349,
    1.728085507,
    1.681875544,
    1.635644158,
    1.589434195,
    1.543224233,
    1.49701427,
    1.450804307,
    1.407315098,
    1.368089116,
    1.328863135,
    1.289637153,
    1.250411172,
    1.211206614,
    1.171980633,
    1.132754651,
    1.0989916,
    1.07141987,
    1.043869563,
    1.016297833,
    0.988726102,
    0.960982986
  };
  assert((sizeof(VALUES) / sizeof(float)) == N_OXY_WAVELENGTHS);

  return InterpolateSpectrumSamples(OXY_WAVELENGTHS, VALUES, N_OXY_WAVELENGTHS, l);
}

Spectrum SigmaAOxy(const Spectrum &l) {
  static const float VALUES[] = {
    142.588906,
    152.2250865,
    165.342546,
    189.7072149,
    226.1867349,
    250.0308186,
    267.897814,
    280.794614,
    279.5092186,
    276.1029209,
    257.2718791,
    231.306893,
    201.5050019,
    174.6166735,
    151.62952,
    131.7915851,
    114.6786884,
    88.54874326,
    71.13592093,
    63.80916744,
    54.93993953,
    49.69124186,
    43.61989116,
    40.87771442,
    35.90751907,
    33.64308093,
    31.52646326,
    28.68145488,
    26.50913674,
    25.43797395,
    23.82266047,
    22.13022326,
    21.31999572,
    19.85571619,
    18.67615172,
    17.78622967,
    16.93508372,
    16.12828391,
    15.45195172,
    14.84524512,
    14.26210409,
    13.76529879,
    13.48615377,
    13.21257879,
    12.94757312,
    12.684924,
    12.36486056,
    12.02787274,
    11.70266772,
    11.38646047,
    11.21121823,
    11.03104865,
    10.93550093,
    10.68270651,
    10.70948558,
    10.73048037,
    10.79217935,
    10.94149944,
    11.24806623,
    12.05572298,
    12.96235516,
    14.16634214,
    15.67603898,
    17.40446726,
    19.27557442,
    21.4001187,
    23.4991693,
    25.1316214,
    26.64624558,
    27.69598512,
    28.51221116,
    28.54220372,
    27.90164837,
    26.70837302,
    24.99022791,
    23.0385693,
    21.24929898,
    19.71753619,
    18.46513265,
    17.91841116,
    17.46702316,
    17.47066512,
    18.16456437,
    19.54615014,
    21.51537581,
    23.83122977,
    26.33560837,
    28.55077302,
    29.7461907,
    29.3112986,
    26.83477023,
    23.19281674,
    18.5523253,
    14.24667935,
    10.58480223,
    7.712800558,
    5.606680279,
    4.112622419,
    3.044030419,
    2.412472837,
    1.713860465,
    1.426788837,
    1.139717209,
    0.958262233,
    0.882423907,
    0.806585581,
    0.730747256,
    0.65490893,
    0.594495349,
    0.549506512,
    0.504517674,
    0.459528837,
    0.41454,
    0.378977395,
    0.352841023,
    0.326704651,
    0.300568279,
    0.274431907,
    0.256436372,
    0.246581674,
    0.236726977,
    0.226872279,
    0.217017581,
    0.209090977,
    0.203092465,
    0.197093953,
    0.191095442,
    0.18509693,
    0.179526884,
    0.174385302,
    0.171171814,
    0.168172558,
    0.165173302,
    0.162174047,
    0.159603256,
    0.15746093,
    0.155318605,
    0.152962047,
    0.151033953,
    0.149534326,
    0.148677395,
    0.147820465,
    0.146963535,
    0.146106605,
    0.146963535,
    0.147820465,
    0.148677395,
    0.149534326,
    0.151033953,
    0.153176279,
    0.155318605
  };
  assert((sizeof(VALUES) / sizeof(float)) == N_OXY_WAVELENGTHS);

  return InterpolateSpectrumSamples(OXY_WAVELENGTHS, VALUES, N_OXY_WAVELENGTHS, l);
}

Spectrum SigmaABaseline(const Spectrum &l) {
  return Spectrum(0.0244f) + 8.53f * Exp(-(l - 154.0f) / 66.2f);
}

Spectrum SigmaADerm(const float Ch, const Spectrum &l) {
  static const float GAMMA = 0.75f;
  return Ch * (GAMMA * SigmaAOxy(l) + (1 - GAMMA) * SigmaADeoxy(l)) + (1 - Ch) * SigmaABaseline(l);
}

Spectrum SigmaAEm(const Spectrum &l) {
  return 6.6e10 * Pow(l, -3.33);
}

Spectrum SigmaAPm(const Spectrum &l) {
  return 2.9e14 * Pow(l, -4.75);
}

Spectrum SigmaAEpi(const float Cm, const float Bm, const Spectrum &l) {
  return Cm * (Bm * SigmaAEm(l) + (1 - Bm) * SigmaAPm(l)) + (1 - Cm) * SigmaABaseline(l);
}

Spectrum SigmaSPrime(const Spectrum &l) {
  return 14.74f * Pow(l, -.22f) + 2.2e11 * Pow(l, -4);
}

Spectrum SigmaSPrimeDerm(const Spectrum &l) {
  return 0.50f * SigmaSPrime(l);
}

Spectrum SigmaSPrimeEpi(const Spectrum &l) {
  return SigmaSPrime(l);
}

Spectrum SigmaTDerm(const float Ch, const Spectrum &l) {
  return SigmaADerm(Ch, l) + SigmaSPrimeDerm(l);
}

Spectrum SigmaTEpi(const float Cm, const float Bm, const Spectrum &l) {
  return SigmaAEpi(Cm, Bm, l) + SigmaSPrimeEpi(l);
}

Spectrum DiffusionCoefficient(const Spectrum &sigmaA, const Spectrum &sigmaSPrime) {
  //TODO: Actually don't recall where this expression for D comes from. It's usually 1 / (3 * sigma_t_prime)
  return Spectrum(1.0f) / (sigmaSPrime * 3.0f + sigmaA);
}

float CalcA(const float eta) {
  const float f = Fdr(eta);
  return (1 + f) / (1 - f);
}

// SkinSubsurfaceIntegrator Local Declarations
struct SkinSubsurfaceOctreeNode {
    // SkinSubsurfaceOctreeNode Methods
    SkinSubsurfaceOctreeNode() {
        isLeaf = true;
        sumArea = 0.f;
        for (int i = 0; i < 8; ++i)
            ips[i] = NULL;
    }
    void Insert(const BBox &nodeBound, IrradiancePoint *ip,
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
            IrradiancePoint *localIps[8];
            for (int i = 0; i < 8; ++i) {
                localIps[i] = ips[i];
                children[i] = NULL;
            }
            for (int i = 0; i < 8; ++i)  {
                IrradiancePoint *ip = localIps[i];
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
            // Init _SkinSubsurfaceOctreeNode_ leaf from _IrradiancePoint_s
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
            // Init interior _SkinSubsurfaceOctreeNode_
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
    Spectrum Mo(const BBox &nodeBound, const Point &p, float maxError);

    // SkinSubsurfaceOctreeNode Public Data
    Point p;
    bool isLeaf;
    Spectrum E;
    float sumArea;
    union {
        SkinSubsurfaceOctreeNode *children[8];
        IrradiancePoint *ips[8];
    };
};

// sigh...
#define MAX_R 25.0f
#define DELTA_R 0.1f

struct DiffusionReflectance {
  static const int N_R = MAX_R / DELTA_R + 1;
  static float distances[N_R];
  static Spectrum REpiPos[N_R];  // aka R1+
  static Spectrum REpiNeg[N_R];  // R1-
  static Spectrum TEpiPos[N_R];  // T1+
  static Spectrum TEpiNeg[N_R];  // T1-
  static Spectrum RDerm[N_R]; // R2+
  static Spectrum R[N_R]; // Final profile
  static Spectrum totalDiffuseReflectance;

  static Spectrum RValue(float d) {
    int imin = d / DELTA_R;
    int imax = imin + 1;
    assert (imin >= 0);
    if (imax >= N_R) {
      imax = N_R - 1;
    }

    float t = d / DELTA_R - imin;
    assert(t >= 0.0f && t <= 1.0f);

    return Lerp(t, R[imin], R[imax]);
  }

  static void Convolve(const Spectrum *R1, const Spectrum *R2,
		       Spectrum* Rout)
  {
    memset(Rout, 0, sizeof(Spectrum) * N_R);
    for (int ri = 0; ri < N_R; ri++) {
      assert(R1[ri] >= 0);
      assert(R2[ri] >= 0);
      Spectrum Rsum = Spectrum(0.0f);
      for (int rj = 0; rj <= ri; rj++) {
	int rk = ri - rj;
	if (rk >= 0 && rk < N_R) {
	  Rsum += R1[rj] * R2[rk];
	}
      }
      Rout[ri] = (Rsum * DELTA_R).Clamp(0, 1);
    }
  }

  // Evaluate multipole diffusion model at radius r.
  // d: slab thickness
  // etas: ratios of indices of refraction at entry and exit interfaces
  static void EvaluateMultipole(float etaTop, float etaBot, int nPoles, float slabThickness, float r, Spectrum sigmaA,
				Spectrum sigmaSPrime, Spectrum &Rsum, Spectrum &Tsum)
  {
    // sanity check: see if the dipole approximation matches... and it seems it does.
    /*
    Spectrum sigmaTPrime = sigmaA + sigmaSPrime;
    Spectrum str = Sqrt(sigmaA * sigmaTPrime * 3.0f);
    Spectrum alphaPrime = sigmaSPrime / sigmaTPrime;
    Spectrum zpos = Spectrum(1.f) / sigmaTPrime;
    Spectrum A = CalcA(etaTop);
    //Spectrum zneg = -(zpos + 4 * A * D);
    Spectrum zneg = (Spectrum(1.0f) + (4.0f / 3.0f) * A) / sigmaTPrime;
    Spectrum dpos = Sqrt(Spectrum(r * r) + zpos * zpos);
    Spectrum dneg = Sqrt(Spectrum(r * r) + zneg * zneg);
    Spectrum dppos = (str * dpos + 1.0f) * Exp(-str * dpos) / (sigmaTPrime * Pow(dpos, 3));
    Spectrum dpneg = zneg * (str * dneg + 1.0f) * Exp(-str * dneg) / (sigmaTPrime * Pow(dneg, 3));
    Rsum = alphaPrime / (4.0f * M_PI) * (dppos + dpneg);
    Tsum = 0;
    */

    Spectrum sigmaTPrime = sigmaA + sigmaSPrime;
    Spectrum str = Sqrt(3.0f * sigmaA * sigmaTPrime); // effective transport coefficient
    Spectrum alphaPrime = sigmaSPrime / sigmaTPrime;
    Spectrum D = DiffusionCoefficient(sigmaA, sigmaSPrime);
    Spectrum meanFreePath = Spectrum(1.0f) / sigmaTPrime;
    float A0 = CalcA(etaTop);
    float Ad = CalcA(etaBot);
    Spectrum zb0 = 2.0f * A0 * D;
    Spectrum zbd = 2.0f * Ad * D;

    // Evaluate multipole model for epidermis
    Rsum = 0;
    Tsum = 0;
    for (int i = -nPoles; i <= nPoles; i++) {
      Spectrum zri = (zb0 + zbd + slabThickness) * 2.0f * i - meanFreePath;
      Spectrum zvi = (zb0 + zbd + slabThickness) * 2.0f * i + meanFreePath + 2 * zb0;
      Spectrum dri = Sqrt(zri * zri + r * r);
      Spectrum dvi = Sqrt(zvi * zvi + r * r);

      Spectrum numer1 = alphaPrime * zvi * (str * dvi + 1.0f) * Exp(-str * dvi);
      Spectrum numer2 = alphaPrime * zri * (str * dri + 1.0f) * Exp(-str * dri);
      Spectrum denom1 = 4.0f * M_PI * Pow(dvi, 3);
      Spectrum denom2 = 4.0f * M_PI * Pow(dri, 3);
      Rsum += numer1 / denom1 - numer2 / denom2;

      Spectrum offsetr = zri + slabThickness;
      Spectrum offsetv = zvi + slabThickness;
      dri = Sqrt(offsetr * offsetr + r * r);
      dvi = Sqrt(offsetv * offsetv + r * r);
      numer1 = alphaPrime * offsetr * (str * dri + 1.0f) * Exp(-str * dri);
      numer2 = alphaPrime * offsetv * (str * dvi + 1.0f) * Exp(-str * dvi);
      denom1 = 4.0f * M_PI * Pow(dri, 3);
      denom2 = 4.0f * M_PI * Pow(dvi, 3);
      Tsum += numer1 / denom1 - numer2 / denom2;
    }
    Rsum = Rsum.Clamp(0, 1);
    Tsum = Tsum.Clamp(0, 1);
  }

  // Pre-build convolved profiles for all wavelengths
  static void BuildProfile(float Cm, float Bm, float Ch, float rho, float eta, const Spectrum& l, int nPoles, float epiThickness) {
      memset(REpiPos, 0, sizeof(REpiPos));
      memset(REpiNeg, 0, sizeof(REpiNeg));
      memset(TEpiPos, 0, sizeof(TEpiPos));
      memset(TEpiNeg, 0, sizeof(TEpiNeg));
      memset(RDerm, 0, sizeof(RDerm));

      // Build epidermis profile (and fill out distance values)
      //std::cout << "Building profile" << std::endl;
      for (int ri = 0; ri < N_R; ri++) {
	float r = ri * DELTA_R;
	Spectrum sa = SigmaAEpi(Cm, Bm, l);
	Spectrum ssp = SigmaSPrimeEpi(l);
	distances[ri] = r;

	EvaluateMultipole(eta, 1.0, nPoles, epiThickness, r,
			  sa, ssp, REpiPos[ri], TEpiPos[ri]);
	EvaluateMultipole(1.0, eta, nPoles, epiThickness,
			  r, sa, ssp, REpiNeg[ri], TEpiNeg[ri]);
      }

      // Evaluate dipole model for dermis
      for (int ri = 0; ri < N_R; ri++) {
	float d = distances[ri];
	Spectrum stprime = SigmaADerm(Ch, l) + SigmaSPrimeDerm(l);
	Spectrum str = Sqrt(SigmaADerm(Ch, l) * stprime * 3.0f);
	Spectrum alphaPrime = SigmaSPrimeDerm(l) / stprime;
	Spectrum D = Spectrum(1.0f) / (stprime * 3.0f);
	Spectrum meanFreePath = Spectrum(1.0f) / stprime;
	float A = (1.f + Fdr(1.0f)) / (1.f - Fdr(1.0f));
	Spectrum zb = A * D * 2.0f;
	Spectrum zr = meanFreePath;
	Spectrum zv = -meanFreePath - zb * 2.0f;
	Spectrum dr = Sqrt(zr * zr + d * d);
	Spectrum dv = Sqrt(zv * zv + d * d);

	Spectrum numer1 = alphaPrime * zr * (str * dr + 1.0f) * Exp(-str * dr);
	Spectrum denom1 = 4.0f * M_PI * Pow(dr, 3);
	Spectrum numer2 = alphaPrime * zv * (str * dv + 1.0f) * Exp(-str * dv);
	Spectrum denom2 = 4.0f * M_PI * Pow(dv, 3);
	RDerm[ri] = (numer1 / denom1 - numer2 / denom2).Clamp(0, 1);
      }

      // Convolve profiles together
      //std::cout << "Convolving" << std::endl;
      Spectrum trans1[N_R], temp1[N_R];
      Convolve(TEpiPos, RDerm, temp1);
      Convolve(temp1, TEpiNeg, trans1); // T1+ * R2+ * T1-

      Spectrum trans2[N_R], temp2[N_R];
      Convolve(RDerm, TEpiNeg, temp1);
      Convolve(REpiNeg, temp1, temp2);
      Convolve(RDerm, temp2, temp1);
      Convolve(TEpiPos, temp1, trans2); // T1+ * R2+ * R1- * R2+ * T1-

      // Add it all together
      memset(R, 0, sizeof(R));
      for (int i = 0; i < N_R; i++) {
	R[i] = REpiPos[i] + trans1[i] + trans2[i];
        //R[i] = trans1[i] + trans2[i];
        //R[i] = REpiPos[i];
	R[i] = R[i].Clamp(0, 1.0);
      }

      // Compute total diffuse reflectance
      std::cout << "Computing total diffuse reflectance=";
      totalDiffuseReflectance = Spectrum(0.f);
      for (int ri = 0; ri < N_R; ri++) {
	totalDiffuseReflectance += R[ri] * distances[ri] * DELTA_R;
      }
      totalDiffuseReflectance = totalDiffuseReflectance.Clamp();
      totalDiffuseReflectance.Print(stdout);
    }
};

float DiffusionReflectance::distances[N_R];
Spectrum DiffusionReflectance::REpiPos[N_R];
Spectrum DiffusionReflectance::REpiNeg[N_R];
Spectrum DiffusionReflectance::TEpiPos[N_R];
Spectrum DiffusionReflectance::TEpiNeg[N_R];
Spectrum DiffusionReflectance::RDerm[N_R];
Spectrum DiffusionReflectance::R[N_R];
Spectrum DiffusionReflectance::totalDiffuseReflectance;

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
        irradiancePoints.push_back(IrradiancePoint(sp, E));
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

    // Compute specular & diffuse components separately
    Spectrum LSpec(0.f);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const DifferentialGeometry &dgs = bsdf->dgShading;
    const Point &p = dgs.p;
    const Normal &n = dgs.nn;

    // Add specular highlights, using only a Torrance-Sparrow
    // microfacet model.
    //    Spectrum shine = totalDiffuseReflectance;
    Spectrum shine = Spectrum(rho);
    MicrofacetDistribution *md = BSDF_ALLOC(arena, Beckmann)(0.35f);
    Fresnel *fr = BSDF_ALLOC(arena, FresnelDielectric)(1., 1.4f);
    bsdf->Add(BSDF_ALLOC(arena, SkinDJMicrofacet)(shine, fr, md));
    LSpec = UniformSampleAllLights(scene, renderer, arena, p, n,
				   wo, isect.rayEpsilon, ray.time,
				   bsdf, sample, rng,
				   lightSampleOffsets,
				   bsdfSampleOffsets);

    // Compute rho_dr; darkens some patches of skin when more oil is present
    Spectrum rhodr = bsdf->rho(ray.d, rng, BxDFType(BSDF_REFLECTION | BSDF_GLOSSY), 10).Clamp(0, 1);

    // Scattering component. Clear out specular contribution.
    bsdf->Clear();
    const GeometricPrimitive* geoPrim = dynamic_cast<const GeometricPrimitive*>(isect.primitive);
    assert(geoPrim != NULL);
    Reference<Material> skinMaterial = geoPrim->material;
    const SkinDJMaterial *djMat = dynamic_cast<const SkinDJMaterial*>(skinMaterial.GetPtr());
    assert(djMat != NULL);    
    Spectrum texr = djMat->Kd->Evaluate(dgs).Clamp(0, 1);
    //Spectrum texr = Spectrum(0.5f / M_PI);
    bsdf->Add(BSDF_ALLOC(arena, Lambertian)(texr));

    // Evaluate BSSRDF and possibly compute subsurface scatter
    BSSRDF *bssrdf = isect.GetBSSRDF(ray, arena);
    Spectrum LScatter = Spectrum(0.0f);
    if (bssrdf && octree) {	
      // Use hierarchical integration to evaluate reflection from multipole model
      PBRT_SUBSURFACE_STARTED_OCTREE_LOOKUP(const_cast<Point *>(&p));
      Spectrum Mo = octree->Mo(octreeBounds, p, maxError);
      float eta = 1.4f;
      FresnelDielectric fresnel(1.f, eta);
      //Spectrum Ft = Spectrum(1.f) - fresnel.Evaluate(AbsDot(wo, n));
      Spectrum Ft = Spectrum(1.f) - rhodr;

      LScatter = (INV_PI * Ft) * Mo;

      PBRT_SUBSURFACE_FINISHED_OCTREE_LOOKUP();
    }

    // Diffuse component. S4 of the Skin BSSRDF paper says this is
    // used to scale the scattering component.
    Spectrum LDiffuse = UniformSampleAllLights(scene, renderer, arena, p, n,
					       wo, isect.rayEpsilon, ray.time,
					       bsdf, sample, rng,
					       lightSampleOffsets,
					       bsdfSampleOffsets);

    // Add it all together. 
    //L = LSpec + LScatter;
    L = LScatter * LDiffuse; // fraction rho_dr is reflected right away by oil film, LScatter is modulated
    L += (LSpec - L).Clamp(0, 1); // Don't let specular component overwhelm diffuse. We just want the highlights.
    //L += LSpec;

    if (ray.depth < maxSpecularDepth) {
      // Trace rays for specular reflection and refraction.
      L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene,
			     sample, arena);
      L += SpecularTransmit(ray, bsdf, rng, isect,
			      renderer, scene, sample, arena);
    }


    
    return L;
}


Spectrum SkinSubsurfaceOctreeNode::Mo(const BBox &nodeBound, const Point &pt, float maxError) {
    // Compute $M_\roman{o}$ at node if error is low enough
    float dw = sumArea / DistanceSquared(pt, p);
    if (dw < maxError && !nodeBound.Inside(pt))
    {
        PBRT_SUBSURFACE_ADDED_INTERIOR_CONTRIBUTION(const_cast<SkinSubsurfaceOctreeNode *>(this));
        //return Rd(DistanceSquared(pt, p)) * E * sumArea;
	return DiffusionReflectance::RValue(Distance(pt, p)) * E * sumArea;
    }

    // Otherwise compute $M_\roman{o}$ from points in leaf or recursively visit children
    Spectrum Mo = 0.f;
    if (isLeaf) {
        // Accumulate $M_\roman{o}$ from leaf node
        for (int i = 0; i < 8; ++i) {
            if (!ips[i]) break;
            PBRT_SUBSURFACE_ADDED_POINT_CONTRIBUTION(const_cast<IrradiancePoint *>(ips[i]));
	    Mo += DiffusionReflectance::RValue(Distance(pt, ips[i]->p)) * ips[i]->E * ips[i]->area;
            //Mo += Rd(DistanceSquared(pt, ips[i]->p)) * ips[i]->E * ips[i]->area;
        }
    }
    else {
        // Recursively visit children nodes to compute $M_\roman{o}$
        Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
        for (int child = 0; child < 8; ++child) {
            if (!children[child]) continue;
            BBox childBound = octreeChildBound(child, nodeBound, pMid);
            Mo += children[child]->Mo(childBound, pt, maxError);
        }
    }
    return Mo;
}


SkinSubsurfaceIntegrator *CreateSkinSubsurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float maxError = params.FindOneFloat("maxerror", .05f);
    float minDist = params.FindOneFloat("minsampledistance", .25f);
    string pointsfile = params.FindOneFilename("pointsfile", "");
    if (PbrtOptions.quickRender) { maxError *= 4.f; minDist *= 4.f; }

    float Cm = params.FindOneFloat("Cm", 0.0f);
    float Bm = params.FindOneFloat("Bm", 0.0f);
    float Ch = params.FindOneFloat("Ch", 0.0f);
    float rho = params.FindOneFloat("rho", 0.0f);

    float RGBWavelengths[] = {650, 510, 475};
    Spectrum c = RGBSpectrum::FromRGB(RGBWavelengths);

    

    DiffusionReflectance::BuildProfile(Cm, Bm, Ch, rho, 1.4f, c, 0, 0.25f);
    
    return new SkinSubsurfaceIntegrator(maxDepth, maxError, minDist, Cm, Bm, Ch, rho, pointsfile);
}


