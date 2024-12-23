#ifndef _SPACES_H_
#define _SPACES_H_

struct unknown_s { }; //unknown space;
struct mirror_s { }; // normal along mirror direction, x,y whatever
struct local_s { };  // normal along z, x,y over tangents, o on the surface location
struct world_s { };  // world
struct object_s { };  // shape
struct camera_s { };   // camera_s
struct light_s { }; // source
struct raster_s { }; // raster space
struct screen_s { }; // screen space
struct texture_s { };
struct material_s{ }; // material space
struct volume_s { }; // volume space

#endif
