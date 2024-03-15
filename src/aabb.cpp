#include "aabb.h" 

// @@@@@@ VOTRE CODE ICI
// Implémenter l'intersection d'un rayon avec un AABB dans l'intervalle décrit.
bool AABB::intersect(Ray ray, double t_min, double t_max)  {
    //calcul des valeurs de t pour les intersections avec les plans x
    double tmin = (min.x - ray.origin.x) / ray.direction.x;
    double tmax = (max.x - ray.origin.x) / ray.direction.x;

    if (tmin > tmax) std::swap(tmin, tmax);

    //calcul des valeurs de t pour les intersections avec les plans y
    double tymin = (min.y - ray.origin.y) / ray.direction.y;
    double tymax = (max.y - ray.origin.y) / ray.direction.y;

    if (tymin > tymax) std::swap(tymin, tymax);

    //verification de l'intersection entre le segment de rayon et le volume englobant sur l'axe y
    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    //calcul des valeurs de t pour les intersections avec les plans z
    double tzmin = (min.z - ray.origin.z) / ray.direction.z;
    double tzmax = (max.z - ray.origin.z) / ray.direction.z;

    if (tzmin > tzmax) std::swap(tzmin, tzmax);

    //verification de l'intersection entre le segment de rayon et le volume englobant sur l'axe z
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    return true;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction qui permet de trouver les 8 coins de notre AABB.
std::vector<double3> retrieve_corners(AABB aabb) {
    std::vector<double3> corners(8);

    //calcul des 8 coins de l'AABB
    corners[0] = aabb.min; // min x, min y, min z
    corners[1] = {aabb.min.x, aabb.min.y, aabb.max.z}; // min x, min y, max z
    corners[2] = {aabb.min.x, aabb.max.y, aabb.min.z}; // min x, max y, min z
    corners[3] = {aabb.min.x, aabb.max.y, aabb.max.z}; // min x, max y, max z
    corners[4] = {aabb.max.x, aabb.min.y, aabb.min.z}; // max x, min y, min z
    corners[5] = {aabb.max.x, aabb.min.y, aabb.max.z}; // max x, min y, max z
    corners[6] = {aabb.max.x, aabb.max.y, aabb.min.z}; // max x, max y, min z
    corners[7] = aabb.max; // max x, max y, max z

    return corners;
}

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.
AABB construct_aabb(std::vector<double3> points) {
    //si aucun point n'est provide:
    if (points.empty()) {
        return AABB{};
    }

     //init une AABB avec les coordonnwes du premier point
    AABB aabb;
    aabb.min = points[0];
    aabb.max = points[0];

    //parcourt tous les points pour trouver les coordonnees min et max
    for (const auto& point : points) {
        aabb.min.x = std::min(aabb.min.x, point.x);
        aabb.min.y = std::min(aabb.min.y, point.y);
        aabb.min.z = std::min(aabb.min.z, point.z);
        aabb.max.x = std::max(aabb.max.x, point.x);
        aabb.max.y = std::max(aabb.max.y, point.y);
        aabb.max.z = std::max(aabb.max.z, point.z);
    }

    return aabb;
};

AABB combine(AABB a, AABB b) {
	return AABB{min(a.min,b.min),max(a.max,b.max)};
};

bool compare(AABB a, AABB b, int axis){
	return a.min[axis] < b.min[axis];
};