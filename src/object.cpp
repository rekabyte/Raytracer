#include "object.h"

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection d'une sphère.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Sphere::local_intersect(Ray ray, 
							  double t_min, double t_max, 
							  Intersection *hit) 
{
    //on verifie si le ray commence dans la sphere
    if (length(ray.origin) <= radius + EPSILON) return false;

    //calculs necessaires pour resoudre l'equation de la sphere:
    auto a = dot(ray.direction, ray.direction);
    auto b = 2.0 * dot(ray.origin, ray.direction);
    auto c = dot(ray.origin, ray.origin) - this->radius * this->radius;
    auto discriminant = b * b - 4 * a * c;

    //si discriminant negatif, alors logiquement il n'intersecte pas
    if (discriminant < 0) return false;

    //si discriminant positif
    //on resout l'equation de la sphere
    auto sqrt_discriminant = sqrt(discriminant);
    auto root = (-b - sqrt_discriminant) / (2.0 * a);
    if (root <= t_min || root >= t_max) {
        root = (-b + sqrt_discriminant) / (2.0 * a);
        if (root < t_min || root > t_max) {
            return false;
        }
    }

    //mis a jour des infos de l'intersection:
    hit->depth = root;
    hit->position = ray.origin + root * ray.direction;
    hit->normal = normalize(hit->position);

    //Pour les UVs:
    //on convertit le point d'intersection en coord spheriques:
    double theta = atan2(hit->position.z, hit->position.x);
    double phi = acos(hit->position.y / radius);

    //on calcule u et v
    double u = (theta + PI) / (2 * PI);
    double v = phi / PI;

    //Pour une raison que je ne veux pas savoir (ni debug, plus de temps :( )), les textures
    //sont decalees de 90 degre sur l'axe des U
    double offset = -0.25;
    u = u - offset;

    //il faut que le u reste dans l'interval [0,1]
    if (u < 0) u += 1;
    if (u > 1) u -= 1;

    //les textures sont inversees par rapport a l'axe x, donc je les corrige:
    u = 1.0 - u;

    hit->uv = double2{u, v};

    return true;
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
	//les coins min et max LOCALS de l'AABB
    double3 min_local = double3{-radius, -radius, -radius};
    double3 max_local = double3{radius, radius, radius};

    //calcul de la position local des 8 points:
    std::vector<double3> corners_local = retrieve_corners(AABB{min_local, max_local});

    //on transforme tout ces points en espace global avec la formule ci-dessous:
    std::vector<double3> corners_global;
    for (const auto& corner_local : corners_local) {
        corners_global.push_back(mul(transform, {corner_local, 1}).xyz());
    }

    //finalement, on construit l'aabb a partir des ces coins:
    return construct_aabb(corners_global);
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un quad (rectangle).
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Quad::local_intersect(Ray const ray, double t_min, double t_max, Intersection* hit) {
    //equation du plan
    double t = -(ray.origin.z) / (ray.direction.z);

    //signifie que le rayon ne frappe pas le plan du quad:
    if (t < t_min || t > t_max) {
        return false;
    }

    //coordonnees du point d'intersection sur le plan du quad
    double3 intersection_point = ray.origin + ray.direction * t;

    //verification si le point d'intersection est a l'interieur du quad
    if (abs(intersection_point.x) <= half_size && abs(intersection_point.y) <= half_size) {
        // Le point d'intersection est à l'interieur du quad
        hit->depth = t;
        hit->position = intersection_point; //position de l'intersection

        //normale du plan du quad (z+ dans l'espace local)
        double3 local_normal = double3(0.0, 0.0, 1.0);

        //inverser la normale si dans la meme direction du rayon
        if (dot(local_normal, ray.direction) > 0) {
            local_normal = -local_normal;
        }

        hit->normal = local_normal;

        //normalizer x et y pour l'utiliser dans u et v:
        double min_x = -half_size;
        double max_x = half_size;
        double min_y = -half_size;
        double max_y = half_size;

        double u = (hit->position.x + half_size) / (half_size * 2);
        double v = (hit->position.y + half_size) / (half_size * 2);

        //(1-v) parce que j'inverse j'inverse la texture par rapport a l'axe v
        hit->uv = double2{u, 1 - v};

        return true;
    }

    return false; //le point d'intersection n'est pas à l'interieur du quad
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {

    //return Object::compute_aabb();

	//determinez les coins du quad dans le repere local
    double half_width = half_size;
    double3 min_local = double3{-half_width, -half_width, 0}; // le quad est sur le plan XY dans l'espace local
    double3 max_local = double3{half_width, half_width, 0};

    //ajouter un petit epsilon pour garantir que le AABB ait un volume
    min_local -= double3{EPSILON, EPSILON, EPSILON};
    max_local += double3{EPSILON, EPSILON, EPSILON};

    //transformer les coins du AABB local dans le repere global
    std::vector<double3> corners_local = retrieve_corners(AABB{min_local, max_local});
    std::vector<double3> corners_global;
    for (const auto& corner_local : corners_local) {
        corners_global.push_back(mul(transform, {corner_local, 1}).xyz());
    }

    //construir le AABB global a partir des coins transformees
    return construct_aabb(corners_global);
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un cylindre.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Cylinder::local_intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    double a = ray.direction.x * ray.direction.x + ray.direction.z * ray.direction.z;
    double b = 2 * (ray.origin.x * ray.direction.x + ray.origin.z * ray.direction.z);
    double c = ray.origin.x * ray.origin.x + ray.origin.z * ray.origin.z - radius * radius;

    //resolution de l'equation quadratique pour trouver les points d'intersection
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return false;
    }

    auto sqrt_discriminant = sqrt(discriminant);
    auto root = (-b - sqrt_discriminant) / (2.0 * a);
    auto root2 = (-b + sqrt_discriminant) / (2.0 * a);

    if (root <= t_min || root >= t_max) {
        if (root2 < t_min || root2 > t_max) {
            return false;
        }
    }

    double3 point = ray.origin + root * ray.direction;
    //verification de l'appartenance du point au cylindre (si le premier point est trop "loin", check le deuxième)
    if (point.y < -half_height || point.y > half_height) {
        root = root2;
        point = ray.origin + root * ray.direction;
        if (point.y < -half_height || point.y > half_height) {
            return false;
        }
    }

    //calcul de la normale du cylindre au point d'intersection
    double3 normal = normalize(double3(point.x, 0.0, point.z));

    //inverser la normale si elle est dans la meme direction que le rayon incident
    if (dot(normal, ray.direction) > 0) {
        normal = -normal;
    }

    hit->depth = root;
    hit->position = point;
    hit->normal = normal;

    //calculer les coordonnees cylindriques
    double theta = atan2(point.z, point.x);
    double u = (theta + PI) / (2 * PI);
    double v = (point.y + half_height) / (2 * half_height);

    //offset pour les textures:
    u = 0.5 - u;
    v = 1 - v;

    if (u < 0) u += 1;

    hit->uv = double2{u, v};

    return true;
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
    return Object::compute_aabb();

	//le centre du cylindre est à l'origine dans le repere local
    double3 center = double3{0, 0, 0};

    //les points minimaux et maximaux du cylindre dans le repere local
    double3 min_local = double3{center.x - radius, center.y - half_height, center.z - radius};
    double3 max_local = double3{center.x + radius, center.y + half_height, center.z + radius};

    //les coins du AABB dans le repere local
    std::vector<double3> corners_local = retrieve_corners(AABB{min_local, max_local});

    //transformer les coins dans le repere global
    std::vector<double3> corners_global;
    for (const auto& corner_local : corners_local) {
        corners_global.push_back(mul(transform, {corner_local, 1}).xyz());
    }

    //construire le AABB global à partir des coins transformes
    return construct_aabb(corners_global);
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un mesh.
//
// Référez-vous au PDF pour la paramétrisation pour les coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
//
bool Mesh::local_intersect(Ray ray,  
						   double t_min, double t_max, 
						   Intersection* hit)
{
	bool hit_anything = false;
    double closest_so_far = t_max;
    Intersection temp_hit;

    for (const auto& tri : triangles) {
        if (intersect_triangle(ray, t_min, closest_so_far, tri, &temp_hit)) {
            //si une intersection est trouvee, met a jour les informations de l'intersection
            hit_anything = true;
            closest_so_far = temp_hit.depth;
            *hit = temp_hit;
        }
    }

    return hit_anything;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un triangle.
// S'il y a intersection, remplissez hit avec l'information sur la normale et les coordonnées texture.
bool Mesh::intersect_triangle(Ray ray, 
                              double t_min, double t_max,
                              Triangle const tri,
                              Intersection *hit)
{
    //recuperation des positions des sommets du triangle
    double3 const &p0 = positions[tri[0].pi];
    double3 const &p1 = positions[tri[1].pi];
    double3 const &p2 = positions[tri[2].pi];

    //calcul des aretes du triangle
    double3 edge1 = p1 - p0;
    double3 edge2 = p2 - p0;

    //calcul du produit vectoriel entre la direction du rayon et edge2
    double3 h = cross(ray.direction, edge2);
    double a = dot(edge1, h);

    //si le produit vectoriel est proche de zero, le rayon est parallele au triangle
    if (a > -0.00001 && a < 0.00001)
        return false;

    double f = 1.0 / a;
    double3 s = ray.origin - p0;
    double u = f * dot(s, h);

    //verifie si le point d'intersection est à l'interieur du triangle
    if (u < 0.0 || u > 1.0)
        return false;

    double3 q = cross(s, edge1);
    double v = f * dot(ray.direction, q);

    //verifie si le point d'intersection est à l'interieur du triangle
    if (v < 0.0 || u + v > 1.0)
        return false;

    //a cette etape, on calcule t pour trouver ou se situe le point d'intersection sur la ligne
    double t = f * dot(edge2, q);

    //verifie si le point d'intersection se trouve à l'interieur du segment de rayon spécifiee
    if (t > t_min && t < t_max) {
        hit->depth = t;
        hit->position = ray.origin + t * ray.direction;
        hit->normal = cross(edge1, edge2);

        //calcul des coordonnees barycentriques
        double w0 = 1.0 - u - v;
        double w1 = u;
        double w2 = v;

        //interpolation des coordonnees UV
        double2 uv0 = tex_coords[tri[0].ti];
        double2 uv1 = tex_coords[tri[1].ti];
        double2 uv2 = tex_coords[tri[2].ti];
        double2 interpolated_uv = w0 * uv0 + w1 * uv1 + w2 * uv2;

        hit->uv = interpolated_uv;

        return true;
    }

    return false;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
//les points minimaux et maximaux de la mesh dans le repere local
    double3 min_local;
    double3 max_local;

    //parcours de tous les triangles de la mesh
    for (const auto& triangle : triangles) {
        //parcours des sommets de chaque triangle
        for (int i = 0; i < 3; ++i) {
            int pos_index = triangle[i].pi;

            //verifie si l'indice de position est valide
            if (pos_index != -1) {
                //acces a la position du sommet
                double3 pos = positions[pos_index];

                //met a jour des valeurs minimales et maximales pour chaque coordonnees
                min_local = min(min_local, pos);
                max_local = max(max_local, pos);
            }
        }
    }

    //les coins du AABB dans le repere local
    std::vector<double3> corners_local = retrieve_corners(AABB{min_local, max_local});

    //les coins du AABB dans le repere global
    std::vector<double3> corners_global;
    for (const auto& corner_local : corners_local) {
        corners_global.push_back(mul(transform, {corner_local, 1}).xyz());
    }

    //construire le AABB global a partir des coins transformes
    return construct_aabb(corners_global);
    
}