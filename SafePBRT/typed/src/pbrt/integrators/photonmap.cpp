
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


// integrators/photonmap.cpp*
#include "integrators/photonmap.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"


// PhotonIntegrator Local Declarations
struct Photon {
    Photon(const Point<world_s> &pp, const Spectrum<power_d> &wt, const Direction<world_s> &w)
        : p(pp), alpha(wt), wi(w) { }
    Photon() { }
    Point<world_s> p;
    Spectrum<power_d> alpha;
    Direction<world_s> wi;
};


class PhotonShootingTask : public Task {
public:
    PhotonShootingTask(int tn, const mfloat<time_d> &ti, Mutex &m, PhotonIntegrator *in,
        ProgressReporter &prog, bool &at, int &ndp,
        std::vector<Photon> &direct, std::vector<Photon> &indir, std::vector<Photon> &caustic,
        std::vector<RadiancePhoton> &rps, std::vector<Spectrum<rho_d> > &rpR, std::vector<Spectrum<rho_d> > &rpT,
        uint32_t &ns, Distribution1D<power_d> *distrib, const Scene *sc,
        const Renderer *sr)
    : taskNum(tn), time(ti), mutex(m), integrator(in), progress(prog),
      abortTasks(at), nDirectPaths(ndp),
      directPhotons(direct), indirectPhotons(indir), causticPhotons(caustic),
      radiancePhotons(rps), rpReflectances(rpR), rpTransmittances(rpT),
      nshot(ns), lightDistribution(distrib), scene(sc), renderer (sr) { }
    void Run();

    int taskNum;
    mfloat<time_d> time;
    Mutex &mutex;
    PhotonIntegrator *integrator;
    ProgressReporter &progress;
    bool &abortTasks;
    int &nDirectPaths;
    std::vector<Photon> &directPhotons, &indirectPhotons, &causticPhotons;
    std::vector<RadiancePhoton> &radiancePhotons;
    std::vector<Spectrum<rho_d> > &rpReflectances, &rpTransmittances;
    uint32_t &nshot;
    const Distribution1D<power_d> *lightDistribution;
    const Scene *scene;
    const Renderer *renderer;
};


struct RadiancePhoton {
    RadiancePhoton(const Point<world_s> &pp, const Normal<world_s> &nn)
        : p(pp), n(nn), Lo() { }
    RadiancePhoton() { }
    Point<world_s> p;
    Normal<world_s> n;
    Spectrum<radiance_d> Lo;
};


class ComputeRadianceTask : public Task {
public:
    ComputeRadianceTask(ProgressReporter &prog, uint32_t tn, uint32_t nt,
        std::vector<RadiancePhoton> &rps, const std::vector<Spectrum<rho_d> > &rhor,
        const std::vector<Spectrum<rho_d> > &rhot,
        uint32_t nlookup, const mfloat<area_d> &md2,
        int ndirect, KdTree<Photon> *direct,
        int nindirect, KdTree<Photon> *indirect,
        int ncaus, KdTree<Photon> *caustic)
        : progress(prog), taskNum(tn), numTasks(nt), radiancePhotons(rps),
          rpReflectances(rhor), rpTransmittances(rhot), nLookup(nlookup),
          maxDistSquared(md2),
          nDirectPaths(ndirect), nIndirectPaths(nindirect), nCausticPaths(ncaus),
          directMap(direct), indirectMap(indirect), causticMap(caustic) { }
    void Run();

private:
    ProgressReporter &progress;
    uint32_t taskNum, numTasks;
    std::vector<RadiancePhoton> &radiancePhotons;
    const std::vector<Spectrum<rho_d> > &rpReflectances, &rpTransmittances;
    uint32_t nLookup;
    mfloat<area_d> maxDistSquared;
    int nDirectPaths, nIndirectPaths, nCausticPaths;
    KdTree<Photon> *directMap, *indirectMap, *causticMap;
};


struct PhotonProcess {
    // PhotonProcess Public Methods
    PhotonProcess(uint32_t mp, ClosePhoton *buf);
    void operator()(const Point<world_s> &p, const Photon &photon, const mfloat<area_d> dist2,
                    mfloat<area_d> &maxDistSquared);
    ClosePhoton *photons;
    uint32_t nLookup, nFound;
};


struct ClosePhoton {
    // ClosePhoton Public Methods
    ClosePhoton(const Photon *p = NULL, const mfloat<area_d> &md2 = mfloat<area_d>(INFINITY))
        : photon(p), distanceSquared(md2) { }
    bool operator<(const ClosePhoton &p2) const {
        return distanceSquared == p2.distanceSquared ?
            (photon < p2.photon) : (distanceSquared < p2.distanceSquared);
    }
    const Photon *photon;
    mfloat<area_d> distanceSquared;
};


PhotonProcess::PhotonProcess(uint32_t mp, ClosePhoton *buf) {
    photons = buf;
    nLookup = mp;
    nFound = 0;
}


struct RadiancePhotonProcess {
    // RadiancePhotonProcess Methods
    RadiancePhotonProcess(const Normal<world_s> &nn)
        :  n(nn) {
        photon = NULL;
    }
    void operator()(const Point<world_s> &p, const RadiancePhoton &rp,
                    const mfloat<area_d> &distSquared, mfloat<area_d> &maxDistSquared) {
        if (Dot(rp.n, n) > 0) {
            photon = &rp;
            maxDistSquared = distSquared;
        }
    }
    const Normal<world_s> &n;
    const RadiancePhoton *photon;
};


inline float kernel(const Photon *photon, const Point<world_s> &p, const mfloat<area_d> maxDist2);
static Spectrum<radiance_d> LPhoton(KdTree<Photon> *map, int nPaths, int nLookup,
    ClosePhoton *lookupBuf, BSDF *bsdf, RNG &rng, const Intersection &isect,
    const Direction<world_s> &w, mfloat<area_d> maxDistSquared);
static Spectrum<irradiance_d> EPhoton(KdTree<Photon> *map, int count, int nLookup,
    ClosePhoton *lookupBuf, const mfloat<area_d> maxDist2, const Point<world_s> &p, const Normal<world_s> &n);

// PhotonIntegrator Local Definitions
inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot) {
    return (found < needed && (found == 0 || found < shot / 1024));
}


inline void PhotonProcess::operator()(const Point<world_s> &p,
        const Photon &photon, const mfloat<area_d> distSquared, mfloat<area_d> &maxDistSquared) {
    if (nFound < nLookup) {
        // Add photon to unordered array of photons
        photons[nFound++] = ClosePhoton(&photon, distSquared);
        if (nFound == nLookup) {
            std::make_heap(&photons[0], &photons[nLookup]);
            maxDistSquared = photons[0].distanceSquared;
        }
    }
    else {
        // Remove most distant photon from heap and add new photon
        std::pop_heap(&photons[0], &photons[nLookup]);
        photons[nLookup-1] = ClosePhoton(&photon, distSquared);
        std::push_heap(&photons[0], &photons[nLookup]);
        maxDistSquared = photons[0].distanceSquared;
    }
}


inline float kernel(const Photon *photon, const Point<world_s> &p,
                    const mfloat<area_d> maxDist2) {
    float s = (1.f - DistanceSquared(photon->p, p) / maxDist2);
    return 3.f * INV_PI * s * s;
}


Spectrum<radiance_d> LPhoton(KdTree<Photon> *map, int nPaths, int nLookup,
      ClosePhoton *lookupBuf, BSDF *bsdf, RNG &rng,
      const Intersection &isect, const Direction<world_s> &wo, mfloat<area_d> maxDist2) {
    Spectrum<radiance_d> L;
    BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
        BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
    if (map && bsdf->NumComponents(nonSpecular) > 0) {
        PBRT_PHOTON_MAP_STARTED_LOOKUP(const_cast<DifferentialGeometry *>(&isect.dg));
        // Do photon map lookup at intersection point
        PhotonProcess proc(nLookup, lookupBuf);
        map->Lookup(isect.dg.p, proc, maxDist2);

        // Estimate reflected radiance due to incident photons
        ClosePhoton *photons = proc.photons;
        int nFound = proc.nFound;
        Normal<world_s> Nf = Faceforward(bsdf->dgShading.nn, wo);
        if (bsdf->NumComponents(BxDFType(BSDF_REFLECTION |
                BSDF_TRANSMISSION | BSDF_GLOSSY)) > 0) {
            // Compute exitant radiance from photons for glossy surface
            for (int i = 0; i < nFound; ++i) {
                const Photon *p = photons[i].photon;
                float k = kernel(p, isect.dg.p, maxDist2);
                L += (k / (nPaths * maxDist2)) * bsdf->f(wo, p->wi) *
                     p->alpha;
            }
        }
        else {
            // Compute exitant radiance from photons for diffuse surface
            Spectrum<irradiance_d> Lr, Lt;
            for (int i = 0; i < nFound; ++i) {
                if (Dot(Nf, photons[i].photon->wi) > 0.f) {
                    float k = kernel(photons[i].photon, isect.dg.p, maxDist2);
                    Lr += (k / (nPaths * maxDist2)) * photons[i].photon->alpha;
                }
                else {
                    float k = kernel(photons[i].photon, isect.dg.p, maxDist2);
                    Lt += (k / (nPaths * maxDist2)) * photons[i].photon->alpha;
                }
            }
            Direction<local_s> _wo = bsdf->WorldToLocal(wo);
            L += Lr * bsdf->rho(_wo, rng, BSDF_ALL_REFLECTION) / hemisphericalProjectedAngle +
                 Lt * bsdf->rho(_wo, rng, BSDF_ALL_TRANSMISSION) / hemisphericalProjectedAngle;
        }
        PBRT_PHOTON_MAP_FINISHED_LOOKUP(const_cast<DifferentialGeometry *>(&isect.dg),
            proc.nFound, proc.nLookup, &L);
    }
    return L;
}


Spectrum<irradiance_d> EPhoton(KdTree<Photon> *map, int count, int nLookup,
        ClosePhoton *lookupBuf, const mfloat<area_d> maxDist2, const Point<world_s> &p,
        const Normal<world_s> &n) {
    if (!map) return Spectrum<irradiance_d>();
    // Lookup nearby photons at irradiance computation point
    PhotonProcess proc(nLookup, lookupBuf);
    mfloat<area_d> md2 = maxDist2;
    map->Lookup(p, proc, md2);
    Assert(md2 > mfloat<area_d>(0.0f));

    // Accumulate irradiance value from nearby photons
    if (proc.nFound == 0) return Spectrum<irradiance_d>();
    ClosePhoton *photons = proc.photons;
    Spectrum<power_d> E;
    for (uint32_t i = 0; i < proc.nFound; ++i)
        if (Dot(n, photons[i].photon->wi) > 0.)
            E += photons[i].photon->alpha;
    return E / (count * md2 * pi);
}



// PhotonIntegrator Method Definitions
PhotonIntegrator::PhotonIntegrator(int ncaus, int nind,
        int nl, int mdepth, int mphodepth, const mfloat<length_d> & mdist, bool fg,
        int gs, float ga) {
    nCausticPhotonsWanted = ncaus;
    nIndirectPhotonsWanted = nind;
    nLookup = nl;
    maxSpecularDepth = mdepth;
    maxPhotonDepth = mphodepth;
    maxDistSquared = mdist * mdist;
    finalGather = fg;
    cosGatherAngle = cos(Radians(ga));
    gatherSamples = gs;
    nCausticPaths = nIndirectPaths = 0;
    causticMap = indirectMap = NULL;
    radianceMap = NULL;
    lightSampleOffsets = NULL;
    bsdfSampleOffsets = NULL;
}


PhotonIntegrator::~PhotonIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
    delete causticMap;
    delete indirectMap;
    delete radianceMap;
}


void PhotonIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
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

    // Request samples for final gathering
    if (finalGather) {
        gatherSamples = max(1, gatherSamples/2);
        if (sampler) gatherSamples = sampler->RoundSize(gatherSamples);
        bsdfGatherSampleOffsets = BSDFSampleOffsets(gatherSamples, sample);
        indirGatherSampleOffsets = BSDFSampleOffsets(gatherSamples, sample);
    }
}


void PhotonIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    // Declare shared variables for photon shooting
    Mutex *mutex = Mutex::Create();
    int nDirectPaths = 0;
    std::vector<Photon> causticPhotons, directPhotons, indirectPhotons;
    std::vector<RadiancePhoton> radiancePhotons;
    bool abortTasks = false;
    causticPhotons.reserve(nCausticPhotonsWanted);
    indirectPhotons.reserve(nIndirectPhotonsWanted);
    uint32_t nshot = 0;
    std::vector<Spectrum<rho_d> > rpReflectances, rpTransmittances;

    // Compute light power CDF for photon shooting
    Distribution1D<power_d> *lightDistribution = ComputeLightSamplingCDF(scene);

    // Run parallel tasks for photon shooting
    ProgressReporter progress(nCausticPhotonsWanted+nIndirectPhotonsWanted, "Shooting photons");
    std::vector<Task *> photonShootingTasks;
    int nTasks = NumSystemCores();
    for (int i = 0; i < nTasks; ++i)
        photonShootingTasks.push_back(new PhotonShootingTask(
            i, camera ? camera->shutterOpen : mfloat<time_d>(0.0f), *mutex, this, progress, abortTasks, nDirectPaths,
            directPhotons, indirectPhotons, causticPhotons, radiancePhotons,
            rpReflectances, rpTransmittances,
            nshot, lightDistribution, scene, renderer));
    EnqueueTasks(photonShootingTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < photonShootingTasks.size(); ++i)
        delete photonShootingTasks[i];
    Mutex::Destroy(mutex);
    progress.Done();

    // Build kd-trees for indirect and caustic photons
    KdTree<Photon> *directMap = NULL;
    if (directPhotons.size() > 0)
        directMap = new KdTree<Photon>(directPhotons);
    if (causticPhotons.size() > 0)
        causticMap = new KdTree<Photon>(causticPhotons);
    if (indirectPhotons.size() > 0)
        indirectMap = new KdTree<Photon>(indirectPhotons);

    // Precompute radiance at a subset of the photons
    if (finalGather && radiancePhotons.size()) {
        // Launch tasks to compute photon radiances
        std::vector<Task *> radianceTasks;
        uint32_t numTasks = 64;
        ProgressReporter progRadiance(numTasks, "Computing photon radiances");
        for (uint32_t i = 0; i < numTasks; ++i)
            radianceTasks.push_back(new ComputeRadianceTask(progRadiance,
                i, numTasks, radiancePhotons, rpReflectances, rpTransmittances,
                nLookup, maxDistSquared, nDirectPaths, directMap,
                nIndirectPaths, indirectMap,
                nCausticPaths, causticMap));
        EnqueueTasks(radianceTasks);
        WaitForAllTasks();
        for (uint32_t i = 0; i < radianceTasks.size(); ++i)
            delete radianceTasks[i];
        progRadiance.Done();
        radianceMap = new KdTree<RadiancePhoton>(radiancePhotons);
    }
    delete directMap;
}


void PhotonShootingTask::Run() {
    // Declare local variables for _PhotonShootingTask_
    MemoryArena arena;
    RNG rng(31 * taskNum);
    std::vector<Photon> localDirectPhotons, localIndirectPhotons, localCausticPhotons;
    std::vector<RadiancePhoton> localRadiancePhotons;
    uint32_t totalPaths = 0;
    bool causticDone = (integrator->nCausticPhotonsWanted == 0);
    bool indirectDone = (integrator->nIndirectPhotonsWanted == 0);
    PermutedHalton halton(6, rng);
    std::vector<Spectrum<rho_d> > localRpReflectances, localRpTransmittances;
    while (true) {
        // Follow photon paths for a block of samples
        const uint32_t blockSize = 4096;
        for (uint32_t i = 0; i < blockSize; ++i) {
            float u[6];
            halton.Sample(++totalPaths, u);
            // Choose light to shoot photon from
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum];

            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential<world_s> photonRay;
            mfloat<invsolidanglearea_d> pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal<world_s> Nl;
            Spectrum<radiance_d> Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);
            if (pdf == mfloat<invsolidanglearea_d>(0.0f) || Le.IsBlack()) continue;
            Spectrum<power_d> alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
            if (!alpha.IsBlack()) {
                // Follow photon path through scene and record intersections
                PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);
                bool specularPath = true;
                Intersection photonIsect;
                int nIntersections = 0;
                while (scene->Intersect(photonRay, &photonIsect)) {
                    ++nIntersections;
                    // Handle photon/surface intersection
                    alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
                    BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, arena);
                    BxDFType specularType = BxDFType(BSDF_REFLECTION |
                        BSDF_TRANSMISSION | BSDF_SPECULAR);
                    bool hasNonSpecular = (photonBSDF->NumComponents() >
                        photonBSDF->NumComponents(specularType));
                    Direction<world_s> wo = -photonRay.d;
                    if (hasNonSpecular) {
                        // Deposit photon at surface
                        Photon photon(photonIsect.dg.p, alpha, wo);
                        bool depositedPhoton = false;
                        if (specularPath && nIntersections > 1) {
                            if (!causticDone) {
                                PBRT_PHOTON_MAP_DEPOSITED_CAUSTIC_PHOTON(&photonIsect.dg, &alpha, &wo);
                                depositedPhoton = true;
                                localCausticPhotons.push_back(photon);
                            }
                        }
                        else {
                            // Deposit either direct or indirect photon
                            // stop depositing direct photons once indirectDone is true; don't
                            // want to waste memory storing too many if we're going a long time
                            // trying to get enough caustic photons desposited.
                            if (nIntersections == 1 && !indirectDone) {
                                PBRT_PHOTON_MAP_DEPOSITED_DIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                                depositedPhoton = true;
                                localDirectPhotons.push_back(photon);
                            }
                            else if (nIntersections > 1 && !indirectDone) {
                                PBRT_PHOTON_MAP_DEPOSITED_INDIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                                depositedPhoton = true;
                                localIndirectPhotons.push_back(photon);
                            }
                        }

                        // Possibly create radiance photon at photon intersection point
                        if (depositedPhoton && integrator->finalGather &&
                                rng.RandomFloat() < .125f) {
                            Normal<world_s> n = photonIsect.dg.nn;
                            n = Faceforward(n, -photonRay.d);
                            localRadiancePhotons.push_back(RadiancePhoton(photonIsect.dg.p, n));
                            Spectrum<rho_d> rho_r = photonBSDF->rho(rng, BSDF_ALL_REFLECTION);
                            localRpReflectances.push_back(rho_r);
                            Spectrum<rho_d> rho_t = photonBSDF->rho(rng, BSDF_ALL_TRANSMISSION);
                            localRpTransmittances.push_back(rho_t);
                        }
                    }
                    if (nIntersections >= integrator->maxPhotonDepth) break;

                    // Sample new photon ray direction
                    Direction<world_s> wi;
                    mfloat<invsolidangle_d> pdf;
                    BxDFType flags;
                    Spectrum<brdf_d> fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng),
                                                       &pdf, BSDF_ALL, &flags);
                    if (fr.IsBlack() || pdf == mfloat<invsolidangle_d>(0.0f))
                        break;
                    Spectrum<power_d> anew = alpha * fr *
                        AbsDot(wi, photonBSDF->dgShading.nn) / pdf;

                    // Possibly terminate photon path with Russian roulette
                    float continueProb = min(1.f, anew.y() / alpha.y());
                    if (rng.RandomFloat() > continueProb)
                        break;
                    alpha = anew / continueProb;
                    specularPath &= ((flags & BSDF_SPECULAR) != 0);
                    if (indirectDone && !specularPath) break;
                    photonRay = RayDifferential<world_s>(photonIsect.dg.p, wi, photonRay,
                                                photonIsect.rayEpsilon);
                }
                PBRT_PHOTON_MAP_FINISHED_RAY_PATH(&photonRay, &alpha);
            }
            arena.FreeAll();
        }

        // Merge local photon data with data in _PhotonIntegrator_
        { MutexLock lock(mutex);

        // Give up if we're not storing enough photons
        if (abortTasks)
            return;
        if (nshot > 500000 &&
            (unsuccessful(integrator->nCausticPhotonsWanted,
                                      causticPhotons.size(), blockSize) ||
             unsuccessful(integrator->nIndirectPhotonsWanted,
                                      indirectPhotons.size(), blockSize))) {
            Error("Unable to store enough photons.  Giving up.\n");
            abortTasks = true;
            return;
        }
        progress.Update(localIndirectPhotons.size() + localCausticPhotons.size());
        nshot += blockSize;

        // Merge indirect photons into shared array
        if (!indirectDone) {
            integrator->nIndirectPaths += blockSize;
            for (uint32_t i = 0; i < localIndirectPhotons.size(); ++i)
                indirectPhotons.push_back(localIndirectPhotons[i]);
            localIndirectPhotons.erase(localIndirectPhotons.begin(),
                                       localIndirectPhotons.end());
            if (indirectPhotons.size() >= integrator->nIndirectPhotonsWanted)
                indirectDone = true;
        }

        // Merge direct, caustic, and radiance photons into shared array
        nDirectPaths += blockSize;
        for (uint32_t i = 0; i < localDirectPhotons.size(); ++i)
            directPhotons.push_back(localDirectPhotons[i]);
        localDirectPhotons.erase(localDirectPhotons.begin(),
                                 localDirectPhotons.end());
        
        if (!causticDone) {
            integrator->nCausticPaths += blockSize;
            for (uint32_t i = 0; i < localCausticPhotons.size(); ++i)
                causticPhotons.push_back(localCausticPhotons[i]);
            localCausticPhotons.erase(localCausticPhotons.begin(), localCausticPhotons.end());
            if (causticPhotons.size() >= integrator->nCausticPhotonsWanted)
                causticDone = true;
        }
        
        for (uint32_t i = 0; i < localRadiancePhotons.size(); ++i)
            radiancePhotons.push_back(localRadiancePhotons[i]);
        localRadiancePhotons.erase(localRadiancePhotons.begin(), localRadiancePhotons.end());
        for (uint32_t i = 0; i < localRpReflectances.size(); ++i)
            rpReflectances.push_back(localRpReflectances[i]);
        localRpReflectances.erase(localRpReflectances.begin(), localRpReflectances.end());
        for (uint32_t i = 0; i < localRpTransmittances.size(); ++i)
            rpTransmittances.push_back(localRpTransmittances[i]);
        localRpTransmittances.erase(localRpTransmittances.begin(), localRpTransmittances.end());
        }

        // Exit task if enough photons have been found
        if (indirectDone && causticDone)
            break;
    }
}


void ComputeRadianceTask::Run() {
    // Compute range of radiance photons to process in task
    uint32_t taskSize = radiancePhotons.size() / numTasks;
    uint32_t excess = radiancePhotons.size() % numTasks;
    uint32_t rpStart = min(taskNum, excess) * (taskSize+1) +
                       max(0, (int)taskNum-(int)excess) * taskSize;
    uint32_t rpEnd = rpStart + taskSize + (taskNum < excess ? 1 : 0);
    if (taskNum == numTasks-1) Assert(rpEnd == radiancePhotons.size());
    ClosePhoton *lookupBuf = new ClosePhoton[nLookup];
    for (uint32_t i = rpStart; i < rpEnd; ++i) {
        // Compute radiance for radiance photon _i_
        RadiancePhoton &rp = radiancePhotons[i];
        const Spectrum<rho_d> &rho_r = rpReflectances[i], &rho_t = rpTransmittances[i];
        if (!rho_r.IsBlack()) {
            // Accumulate outgoing radiance due to reflected irradiance
            Spectrum<irradiance_d> E = EPhoton(directMap, nDirectPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, rp.n) +
                         EPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, rp.n) +
                         EPhoton(causticMap, nCausticPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, rp.n);
            rp.Lo += rho_r * E / hemisphericalProjectedAngle;
        }
        if (!rho_t.IsBlack()) {
            // Accumulate outgoing radiance due to transmitted irradiance
            Spectrum<irradiance_d> E = EPhoton(directMap, nDirectPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, -rp.n) +
                         EPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, -rp.n) +
                         EPhoton(causticMap, nCausticPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, -rp.n);
            rp.Lo += rho_t * E / hemisphericalProjectedAngle;
        }
    }
    delete[] lookupBuf;
    progress.Update();
}


Spectrum<radiance_d> PhotonIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential<world_s> &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum<radiance_d> L;
    Direction<world_s> wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point<world_s> &p = bsdf->dgShading.p;
    const Normal<world_s> &n = bsdf->dgShading.nn;
    L += UniformSampleAllLights(scene, renderer, arena, p, n,
        wo, isect.rayEpsilon, ray.time, bsdf, sample, rng,
        lightSampleOffsets, bsdfSampleOffsets);
    // Compute caustic lighting for photon map integrator
    ClosePhoton *lookupBuf = arena.Alloc<ClosePhoton>(nLookup);
    L += LPhoton(causticMap, nCausticPaths, nLookup, lookupBuf, bsdf,
                 rng, isect, wo, maxDistSquared);

    // Compute indirect lighting for photon map integrator
    if (finalGather) {
    #if 1
        // Do one-bounce final gather for photon map
        BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
            BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
        if (bsdf->NumComponents(nonSpecular) > 0) {
            // Find indirect photons around point for importance sampling
            const uint32_t nIndirSamplePhotons = 50;
            PhotonProcess proc(nIndirSamplePhotons,
                               arena.Alloc<ClosePhoton>(nIndirSamplePhotons));
            mfloat<area_d> searchDist2 = maxDistSquared;
            while (proc.nFound < nIndirSamplePhotons) {
                mfloat<area_d> md2 = searchDist2;
                proc.nFound = 0;
                indirectMap->Lookup(p, proc, md2);
                searchDist2 *= 2.f;
            }

            // Copy photon directions to local array
            Direction<world_s> *photonDirs = arena.Alloc<Direction<world_s> >(nIndirSamplePhotons);
            for (uint32_t i = 0; i < nIndirSamplePhotons; ++i)
                photonDirs[i] = proc.photons[i].photon->wi;

            // Use BSDF to do final gathering
            Spectrum<radiance_d> Li;
            for (int i = 0; i < gatherSamples; ++i) {
                // Sample random direction from BSDF for final gather ray
                Direction<world_s> wi;
                mfloat<invsolidangle_d> pdf;
                BSDFSample bsdfSample(sample, bsdfGatherSampleOffsets, i);
                Spectrum<brdf_d> fr = bsdf->Sample_f(wo, &wi, bsdfSample,
                                             &pdf, BxDFType(BSDF_ALL & ~BSDF_SPECULAR));
                if (fr.IsBlack() || pdf == mfloat<invsolidangle_d>(0.0f)) continue;
                Assert(pdf >= mfloat<invsolidangle_d>(0.0f));

                // Trace BSDF final gather ray and accumulate radiance
                RayDifferential<world_s> bounceRay(p, wi, ray, isect.rayEpsilon);
                Intersection gatherIsect;
                if (scene->Intersect(bounceRay, &gatherIsect)) {
                    // Compute exitant radiance _Lindir_ using radiance photons
                    Spectrum<radiance_d> Lindir;
                    Normal<world_s> nGather = gatherIsect.dg.nn;
                    nGather = Faceforward(nGather, -bounceRay.d);
                    RadiancePhotonProcess proc(nGather);
                    mfloat<area_d> md2 = mfloat<area_d>(INFINITY);
                    radianceMap->Lookup(gatherIsect.dg.p, proc, md2);
                    if (proc.photon != NULL)
                        Lindir = proc.photon->Lo;
                    Lindir *= renderer->Transmittance(scene, bounceRay, NULL, rng, arena);

                    // Compute MIS weight for BSDF-sampled gather ray

                    // Compute PDF for photon-sampling of direction _wi_
                    mfloat<invsolidangle_d> photonPdf(0.0f);
                    mfloat<invsolidangle_d> conePdf = UniformConePdf(cosGatherAngle);
                    for (uint32_t j = 0; j < nIndirSamplePhotons; ++j)
                        if (Dot(photonDirs[j], wi) > .999f * cosGatherAngle)
                            photonPdf += conePdf;
                    photonPdf /= nIndirSamplePhotons;
                    float wt = PowerHeuristic(gatherSamples, pdf, gatherSamples, photonPdf);
                    Li += fr * Lindir * AbsDot(wi, n) * wt / pdf;
                }
            }
            L += Li / gatherSamples;

            // Use nearby photons to do final gathering
            Li = Spectrum<radiance_d>();
            for (int i = 0; i < gatherSamples; ++i) {
                // Sample random direction using photons for final gather ray
                BSDFSample gatherSample(sample, indirGatherSampleOffsets, i);
                int photonNum = min((int)nIndirSamplePhotons - 1,
                    Floor2Int(gatherSample.uComponent * nIndirSamplePhotons));

                // Sample gather ray direction from _photonNum_
                Direction<world_s> vx, vy;
                CoordinateSystem(photonDirs[photonNum], &vx, &vy);
                Direction<world_s> wi = UniformSampleCone(gatherSample.uDir[0], gatherSample.uDir[1],
                                              cosGatherAngle, vx, vy, photonDirs[photonNum]);

                // Trace photon-sampled final gather ray and accumulate radiance
                Spectrum<brdf_d> fr = bsdf->f(wo, wi);
                if (fr.IsBlack()) continue;
                RayDifferential<world_s> bounceRay(p, wi, ray, isect.rayEpsilon);
                Intersection gatherIsect;
                PBRT_PHOTON_MAP_STARTED_GATHER_RAY(&bounceRay);
                if (scene->Intersect(bounceRay, &gatherIsect)) {
                    // Compute exitant radiance _Lindir_ using radiance photons
                    Spectrum<radiance_d> Lindir;
                    Normal<world_s> nGather = gatherIsect.dg.nn;
                    nGather = Faceforward(nGather, -bounceRay.d);
                    RadiancePhotonProcess proc(nGather);
                    mfloat<area_d> md2(INFINITY);
                    radianceMap->Lookup(gatherIsect.dg.p, proc, md2);
                    if (proc.photon != NULL)
                        Lindir = proc.photon->Lo;
                    Lindir *= renderer->Transmittance(scene, bounceRay, NULL, rng, arena);

                    // Compute PDF for photon-sampling of direction _wi_
                    mfloat<invsolidangle_d> photonPdf = mfloat<invsolidangle_d>(0.0f);
                    mfloat<invsolidangle_d> conePdf = UniformConePdf(cosGatherAngle);
                    for (uint32_t j = 0; j < nIndirSamplePhotons; ++j)
                        if (Dot(photonDirs[j], wi) > .999f * cosGatherAngle)
                            photonPdf += conePdf;
                    photonPdf /= nIndirSamplePhotons;

                    // Compute MIS weight for photon-sampled gather ray
                    mfloat<invsolidangle_d> bsdfPdf = bsdf->Pdf(wo, wi);
                    float wt = PowerHeuristic(gatherSamples, photonPdf, gatherSamples, bsdfPdf);
                    Li += fr * Lindir * AbsDot(wi, n) * wt / photonPdf;
                }
                PBRT_PHOTON_MAP_FINISHED_GATHER_RAY(&bounceRay);
            }
            L += Li / gatherSamples;
        }
    #else
        // for debugging / examples: use the photon map directly
        Normal<world_s> nn = Faceforward(n, -ray.d);
        RadiancePhotonProcess proc(nn);
        mfloat<area_d> md2 = mfloat<area_d>(INFINITY);
        radianceMap->Lookup(p, proc, md2);
        if (proc.photon)
            L += proc.photon->Lo;
    #endif
    }
    else
        L += LPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                     bsdf, rng, isect, wo, maxDistSquared);
    if (ray.depth+1 < maxSpecularDepth) {
        Direction<world_s> wi;
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample, arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample, arena);
    }
    return L;
}


PhotonIntegrator *CreatePhotonMapSurfaceIntegrator(const ParamSet &params) {
    int nCaustic = params.FindOneInt("causticphotons", 20000);
    int nIndirect = params.FindOneInt("indirectphotons", 100000);
    int nUsed = params.FindOneInt("nused", 50);
    if (PbrtOptions.quickRender) nCaustic = nCaustic / 10;
    if (PbrtOptions.quickRender) nIndirect = nIndirect / 10;
    if (PbrtOptions.quickRender) nUsed = max(1, nUsed / 10);
    int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
    int maxPhotonDepth = params.FindOneInt("maxphotondepth", 5);
    bool finalGather = params.FindOneBool("finalgather", true);
    int gatherSamples = params.FindOneInt("finalgathersamples", 32);
    if (PbrtOptions.quickRender) gatherSamples = max(1, gatherSamples / 4);
    mfloat<length_d> maxDist = params.FindOneFloat("maxdist", mfloat<length_d>(0.1f));
    float gatherAngle = params.FindOneFloat("gatherangle", 10.f);
    return new PhotonIntegrator(nCaustic, nIndirect,
        nUsed, maxSpecularDepth, maxPhotonDepth, maxDist, finalGather, gatherSamples,
        gatherAngle);
}


