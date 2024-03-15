#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
	bool hasIntersection = false;
    double closestIntersection = t_max;

    // Pile de nœuds à parcourir
    std::vector<BVHNode*> nodesToVisit;
    nodesToVisit.push_back(root);

    while (!nodesToVisit.empty()) {
        // Prend le dernier nœud de la pile
        BVHNode* node = nodesToVisit.back();
        nodesToVisit.pop_back();

        // Vérifie si le rayon intersecte le AABB du nœud
        if (node->aabb.intersect(ray, t_min, t_max)) {
            // Si le nœud est une feuille, effectue l'intersection avec la géométrie
            if (!node->left && !node->right) {
                Intersection tempHit;
                if (objects[node->idx]->intersect(ray, t_min, t_max, &tempHit)) {
                    if (tempHit.depth < closestIntersection) {
                        hasIntersection = true;
                        closestIntersection = tempHit.depth;
                        *hit = tempHit; // Update hit with the intersection data
                    }
                }
            } else {
                // Si le nœud n'est pas une feuille, ajoute ses enfants à la pile
                if (node->left) {
                    nodesToVisit.push_back(node->left);
                }
                if (node->right) {
                    nodesToVisit.push_back(node->right);
                }
            }
        }
    }

    return hasIntersection;
}

// @@@@@@ VOTRE CODE ICI
// - Parcourir tous les objets
// 		- Détecter l'intersection avec l'AABB
//			- Si intersection, détecter l'intersection avec la géométrie.
//				- Si intersection, mettre à jour les paramètres.
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.

bool Naive::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    bool hasIntersection = false;
    double closestIntersection = t_max;

    for(size_t i = 0; i < objects.size(); i++) {
        Intersection tempHit;

        // First, intersect with the AABB
        if(aabbs[i].intersect(ray, t_min, t_max)) {
            // If the ray intersects the AABB, then intersect with the geometry
            if(objects[i]->intersect(ray, t_min, t_max, &tempHit)) {
                if(tempHit.depth < closestIntersection) {
                    hasIntersection = true;
                    closestIntersection = tempHit.depth;
                    *hit = tempHit; // Update hit with the intersection data
                }
            }
        }
    }

    return hasIntersection;
}

//Sans AABB:
//bool Naive::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
//    bool hasIntersection = false;
//    double closestIntersection = t_max;
//
//    for(auto obj : objects) {
//        Intersection tempHit;
//
//        if(obj->intersect(ray, t_min, t_max, &tempHit)) {
//            if(tempHit.depth < closestIntersection) {
//                hasIntersection = true;
//                closestIntersection = tempHit.depth;
//                *hit = tempHit; // Update hit with the intersection data
//            }
//        }
//    }
//
//    return hasIntersection;
//}