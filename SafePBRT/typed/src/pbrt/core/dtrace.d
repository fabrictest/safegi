
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

provider PBRT {
probe started_ray_intersection(const struct Ray *);
probe finished_ray_intersection(const struct Ray *,
                                const struct Intersection *, int hit);
probe started_ray_intersectionp(const struct Ray *);
probe finished_ray_intersectionp(const struct Ray *, int hit);
probe atomic_memory_op();

probe started_parsing();
probe finished_parsing();
probe allocated_cached_transform();

probe started_preprocessing();
probe finished_preprocessing();

probe started_task(const struct Task *);
probe finished_task(const struct Task *);

probe started_generating_camera_ray(const struct Sample *);
probe finished_generating_camera_ray(const struct Sample *, const struct RayDifferential *, float weight);

probe sample_outside_image_extent(const struct CameraSample *);

probe started_rendering();
probe finished_rendering();
probe started_rendertask(int num);
probe finished_rendertask(int num);
probe started_camera_ray_integration(const struct RayDifferential *, const struct Sample *);
probe finished_camera_ray_integration(const struct RayDifferential *, const struct Sample *, const void *L);

probe started_specular_reflection_ray(const struct RayDifferential *);
probe finished_specular_reflection_ray(const struct RayDifferential *);
probe started_specular_refraction_ray(const struct RayDifferential *);
probe finished_specular_refraction_ray(const struct RayDifferential *);

probe started_adding_image_sample(const struct Sample *, const struct RayDifferential *, const void *L, const void *T);
probe finished_adding_image_sample();
probe created_shape(struct Shape *);
probe created_triangle(struct Triangle *);

probe ray_triangle_intersection_test(const struct Ray *, const struct Triangle *);
probe ray_triangle_intersection_hit(const struct Ray *, float t);
probe ray_triangle_intersectionp_test(const struct Ray *, const struct Triangle *);
probe ray_triangle_intersectionp_hit(const struct Ray *, float t);

probe started_trilinear_texture_lookup(float s, float t);
probe finished_trilinear_texture_lookup();
probe started_ewa_texture_lookup(float s, float t);
probe finished_ewa_texture_lookup();

probe loaded_image_map(const char *filename, int width, int height, int elementSize, void *mipPtr);
probe mipmap_trilinear_filter(void *mipPtr, float s, float t, float width, float level, int nLevels);
probe mipmap_ewa_filter(void *mipPtr, float s, float t, float ds0, float ds1, float dt0, float dt1, float minorLength, float majorLength, float level, int nLevels);
probe accessed_texel(void *mipPtr, int level, int s, int t);

probe started_bsdf_shading(const struct Ray *);
probe finished_bsdf_shading(const struct Ray *, const struct BSDF *);
probe started_bssrdf_shading(const struct Ray *);
probe finished_bssrdf_shading(const struct Ray *, const struct BSSRDF *);

probe grid_started_construction(struct GridAccel *, unsigned int numPrimitives);
probe grid_finished_construction(struct GridAccel *);
probe grid_voxelized_primitive(const int *vmin, const int *vmax);
probe grid_bounds_and_resolution(const struct BBox *bbox, const int *nvoxels);
probe grid_intersection_test(const struct GridAccel *, const struct Ray *);
probe grid_intersectionp_test(const struct GridAccel *, const struct Ray *);
probe grid_ray_missed_bounds();
probe grid_ray_traversed_voxel(const int v[3], int nprims);
probe grid_ray_primitive_intersection_test(const struct Primitive *);
probe grid_ray_primitive_intersectionp_test(const struct Primitive *);
probe grid_ray_primitive_hit(const struct Primitive *);

probe kdtree_started_construction(struct KdTreeAccel *, int nprims);
probe kdtree_finished_construction(struct KdTreeAccel *);
probe kdtree_created_leaf(int nprims, int depth);
probe kdtree_created_interior_node(int axis, float split);
probe kdtree_intersection_test(const struct KdTreeAccel *, const struct Ray *);
probe kdtree_intersectionp_test(const struct KdTreeAccel *, const struct Ray *);
probe kdtree_ray_missed_bounds();
probe kdtree_intersection_traversed_interior_node(const struct KdAccelNode *);
probe kdtree_intersection_traversed_leaf_node(const struct KdAccelNode *, int nprims);
probe kdtree_intersectionp_traversed_interior_node(const struct KdAccelNode *);
probe kdtree_intersectionp_traversed_leaf_node(const struct KdAccelNode *, int nprims);
probe kdtree_intersection_hit(const struct Primitive *);
probe kdtree_intersectionp_missed();
probe kdtree_intersectionp_hit(const struct Primitive *);
probe kdtree_intersection_finished();
probe kdtree_intersectionp_primitive_test(const struct Primitive *);
probe kdtree_intersection_primitive_test(const struct Primitive *);

probe bvh_started_construction(struct BVHAccel *, int nprims);
probe bvh_finished_construction(struct BVHAccel *);
probe bvh_intersection_started(const struct BVHAccel *, const struct Ray *);
probe bvh_intersection_traversed_interior_node(const struct LinearBVHNode *);
probe bvh_intersection_traversed_leaf_node(const struct LinearBVHNode *);
probe bvh_intersection_primitive_test(const struct Primitive *);
probe bvh_intersection_primitive_hit(const struct Primitive *);
probe bvh_intersection_primitive_missed(const struct Primitive *);
probe bvh_intersection_finished();
probe bvh_intersectionp_started(const struct BVHAccel *, const struct Ray *);
probe bvh_intersectionp_traversed_interior_node(const struct LinearBVHNode *);
probe bvh_intersectionp_traversed_leaf_node(const struct LinearBVHNode *);
probe bvh_intersectionp_primitive_test(const struct Primitive *);
probe bvh_intersectionp_primitive_hit(const struct Primitive *);
probe bvh_intersectionp_primitive_missed(const struct Primitive *);
probe bvh_intersectionp_finished();

probe supersample_pixel_yes(int xpos, int ypos);
probe supersample_pixel_no(int xpos, int ypos);

probe irradiance_cache_started_ray(const struct RayDifferential *);
probe irradiance_cache_finished_ray(const struct RayDifferential *, float dist, const void *L);
probe irradiance_cache_added_new_sample(const struct Point *, const struct Normal *, float maxDist, const void *E, const struct Vector *primaryDir, float pixelSpacing);
probe irradiance_cache_started_interpolation(const struct Point *p, const struct Normal *n);
probe irradiance_cache_finished_interpolation(const struct Point *p, const struct Normal *n, int successful, int nfound);
probe irradiance_cache_checked_sample(const struct IrradianceSample *, float perr, float nerr);
probe irradiance_cache_started_computing_irradiance(const struct Point *P, const struct Normal *N);
probe irradiance_cache_finished_computing_irradiance(const struct Point *P, const struct Normal *N);

probe photon_map_started_ray_path(const struct RayDifferential *, const void *alpha);
probe photon_map_finished_ray_path(const struct RayDifferential *, const void *alpha);
probe photon_map_deposited_direct_photon(const struct DifferentialGeometry *, const void *alpha, const struct Vector *wo);
probe photon_map_deposited_indirect_photon(const struct DifferentialGeometry *, const void *alpha, const struct Vector *wo);
probe photon_map_deposited_caustic_photon(const struct DifferentialGeometry *, const void *alpha, const struct Vector *wo);
probe photon_map_started_gather_ray(const struct RayDifferential *);
probe photon_map_finished_gather_ray(const struct RayDifferential *);
probe photon_map_started_lookup(const struct DifferentialGeometry *);
probe photon_map_finished_lookup(const struct DifferentialGeometry *, int nFound, int nWanted, const void *L);

probe subsurface_started_rays_for_points();
probe subsurface_finished_rays_for_points(int totalRaysTraced, int numPointsAdded);
probe subsurface_added_point_to_octree(const struct SurfacePoint *, float minSampleDist);
probe subsurface_computed_irradiance_at_point(const struct SurfacePoint *, const void *E);
probe subsurface_added_interior_contribution(const struct SubsurfaceOctreeNode *node);
probe subsurface_added_point_contribution(const struct IrradiancePoint *node);
probe subsurface_started_computing_irradiance_values();
probe subsurface_finished_computing_irradiance_values();
probe subsurface_started_octree_lookup(const struct Point *);
probe subsurface_finished_octree_lookup();

probe mlt_accepted_mutation(float a, const struct MLTSample *current, const struct MLTSample *proposed);
probe mlt_rejected_mutation(float a, const struct MLTSample *current, const struct MLTSample *proposed);
probe mlt_started_mlt_task(struct MLTTask *);
probe mlt_finished_mlt_task(struct MLTTask *);
};


