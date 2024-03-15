#include "object.h"

// Fonction retournant soit la valeur v0 ou v1 selon le signe.
int rsign(double value, double v0, double v1) {
	return (int(std::signbit(value)) * (v1-v0)) + v0;
}

//multiplication matricielle:
double3 multiply(const double4x4& m, const double3& v)
{
    double3 result;

    result.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3];
    result.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3];
    result.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3];

    double w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3];

    if (w != 1.0 && w != 0.0) {
        result.x /= w;
        result.y /= w;
        result.z /= w;
    }

    return result;
}

//fonction auxiliaire:
Ray transform_ray(const Ray& ray, const double4x4& inverse_transform)
{
    double3 transformed_origin = multiply(inverse_transform, ray.origin);
    double3 transformed_direction = multiply(inverse_transform, ray.direction);

    return Ray(transformed_origin, transformed_direction);
}

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
    auto a = dot(ray.direction, ray.direction);
    auto b = 2.0 * dot(ray.origin, ray.direction);
    auto c = dot(ray.origin, ray.origin) - this->radius * this->radius;
    auto discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return false;

    auto sqrt_discriminant = sqrt(discriminant);
    auto root = (-b - sqrt_discriminant) / (2.0 * a);
    if (root <= t_min || root >= t_max) {
        root = (-b + sqrt_discriminant) / (2.0 * a);
        if (root < t_min || root > t_max) {
            return false;
        }
    }

    hit->depth = root;
    hit->position = ray.origin + root*ray.direction;
    hit->normal = normalize(hit->position);
    //std::cout << "sphere hit at : (" << hit->position.x << "," << hit->position.y << "," << hit->position.z << ")" << std::endl;
    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
	// Les coins minimaux et maximaux du AABB dans l'espace local
    double3 min_local = double3{-radius, -radius, -radius};
    double3 max_local = double3{radius, radius, radius};

    // Les coins du AABB dans l'espace local
    std::vector<double3> corners_local = retrieve_corners(AABB{min_local, max_local});

    // Transformer les coins dans le repère global
    std::vector<double3> corners_global;
    for (const auto& corner_local : corners_local) {
        corners_global.push_back(mul(transform, {corner_local, 1}).xyz());
    }

    // Construire le AABB final à partir des coins transformés
    return construct_aabb(corners_global);
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un quad (rectangle).
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Quad::local_intersect(Ray const ray, double t_min, double t_max, Intersection* hit) {
    // Équation du plan du quad dans l'espace local
    double t = -(ray.origin.z) / (ray.direction.z);
    if (t < t_min || t > t_max) {
        return false; // Le rayon ne frappe pas le plan du quad
    }

    // Coordonnées du point d'intersection sur le plan du quad
    double3 intersection_point = ray.origin + ray.direction * t;

    // Vérification si le point d'intersection est à l'intérieur du quad
    if (abs(intersection_point.x) <= half_size && abs(intersection_point.y) <= half_size) {
        // Le point d'intersection est à l'intérieur du quad
        hit->depth = t; // Profondeur de l'intersection
        hit->position = intersection_point; // Position de l'intersection

        // Normale du plan du quad (Z+ dans l'espace local)
        double3 local_normal = double3(0.0, 0.0, 1.0);

        // Inverse the normal if it is in the same direction as the ray
        if (dot(local_normal, ray.direction) > 0) {
            local_normal = -local_normal;
        }

        hit->normal = local_normal;

        return true;
    }

    return false; // Le point d'intersection n'est pas à l'intérieur du quad
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {

    //return Object::compute_aabb();

	// Déterminez les coins du quad dans le repère local
    double half_width = half_size;
    double3 min_local = double3{-half_width, -half_width, 0}; // le quad est sur le plan XY dans l'espace local
    double3 max_local = double3{half_width, half_width, 0};

    // Ajoutez un petit epsilon pour garantir que le AABB ait un volume
    min_local -= double3{EPSILON, EPSILON, EPSILON};
    max_local += double3{EPSILON, EPSILON, EPSILON};

    // Transformez les coins du AABB local dans le repère global
    std::vector<double3> corners_local = retrieve_corners(AABB{min_local, max_local});
    std::vector<double3> corners_global;
    for (const auto& corner_local : corners_local) {
        corners_global.push_back(mul(transform, {corner_local, 1}).xyz());
    }

    // Construisez le AABB global à partir des coins transformés
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

    // Résolution de l'équation quadratique pour trouver les points d'intersection
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
    // Vérification de l'appartenance du point au cylindre (si le premier point est trop "loin", check le deuxième)
    // évite des problèmes qui pourrait s'apparenter à du backface culling (ce n'en est pas en réalité, l'effet est trompeur)
    if (point.y < -half_height || point.y > half_height) {
        root = root2;
        point = ray.origin + root * ray.direction;
        if (point.y < -half_height || point.y > half_height) {
            return false;
        }
    }

    // Calcul de la normale du cylindre au point d'intersection
    double3 normal = normalize(double3(point.x, 0.0, point.z));

    // Inverser la normale si elle est dans la même direction que le rayon incident
    if (dot(normal, ray.direction) > 0) {
        normal = -normal;
    }

    hit->depth = root;
    hit->position = point;
    hit->normal = normal;

    return true;
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
    return Object::compute_aabb();

	// Le centre du cylindre est à l'origine dans le repère local
    double3 center = double3{0, 0, 0};

    // Les points minimaux et maximaux du cylindre dans le repère local
    double3 min_local = double3{center.x - radius, center.y - half_height, center.z - radius};
    double3 max_local = double3{center.x + radius, center.y + half_height, center.z + radius};

    // Les coins du AABB dans le repère local
    std::vector<double3> corners_local = retrieve_corners(AABB{min_local, max_local});

    // Transformer les coins dans le repère global
    std::vector<double3> corners_global;
    for (const auto& corner_local : corners_local) {
        corners_global.push_back(mul(transform, {corner_local, 1}).xyz());
    }

    // Construire le AABB global à partir des coins transformés
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
bool Mesh::intersect_triangle(Ray  ray, 
							  double t_min, double t_max,
							  Triangle const tri,
							  Intersection *hit)
{
	// Extrait chaque position de sommet des données du maillage.
	double3 const &p0 = positions[tri[0].pi]; // ou Sommet A (Pour faciliter les explications)
	double3 const &p1 = positions[tri[1].pi]; // ou Sommet B
	double3 const &p2 = positions[tri[2].pi]; // ou Sommet C

	// Triangle en question. Respectez la convention suivante pour vos variables.
	//
	//     A
	//    / \
	//   /   \
	//  B --> C
	//
	// Respectez la règle de la main droite pour la normale.

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Pour plus de d'informations sur la géométrie, référez-vous à la classe dans object.hpp.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.

	double3 edge1 = p1 - p0;
    double3 edge2 = p2 - p0;

    double3 h = cross(ray.direction, edge2);
    double a = dot(edge1, h);

    if (a > -0.00001 && a < 0.00001)
        return false;    // This ray is parallel to this triangle.

    double f = 1.0 / a;
    double3 s = ray.origin - p0;
    double u = f * dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    double3 q = cross(s, edge1);
    double v = f * dot(ray.direction, q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * dot(edge2, q);

    if (t > t_min && t < t_max) {
        hit->depth = t;
        hit->position = ray.origin + t * ray.direction;
        hit->normal = cross(edge1, edge2); // assuming counter clockwise winding
        return true;
    }

    return false;

}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
// Les points minimaux et maximaux de la mesh dans le repère local
    double3 min_local;// = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    double3 max_local;// = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};

    // Parcours de tous les triangles de la mesh
    for (const auto& triangle : triangles) {
        // Parcours des sommets de chaque triangle
        for (int i = 0; i < 3; ++i) {
            int pos_index = triangle[i].pi;

            // Vérifie si l'indice de position est valide
            if (pos_index != -1) {
                // Accès à la position du sommet
                double3 pos = positions[pos_index];

                // Mise à jour des valeurs minimales et maximales pour chaque coordonnée
                min_local = min(min_local, pos);
                max_local = max(max_local, pos);
            }
        }
    }

    // Les coins du AABB dans le repère local
    std::vector<double3> corners_local = retrieve_corners(AABB{min_local, max_local});

    // Les coins du AABB dans le repère global
    std::vector<double3> corners_global;
    for (const auto& corner_local : corners_local) {
        corners_global.push_back(mul(transform, {corner_local, 1}).xyz());
    }

    // Construire le AABB global à partir des coins transformés
    return construct_aabb(corners_global);
    
}