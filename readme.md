# RayTracer en C++

## Description du Projet

Ce projet est un RayTracer développé en C++, conçu pour créer des images tridimensionnelles réalistes en simulant le comportement des rayons lumineux. Un RayTracer fonctionne en traçant le chemin des rayons de lumière à travers les pixels d'une image et en simulant leurs interactions avec les objets dans la scène virtuelle. Ce processus permet de produire des effets de lumière réalistes tels que les ombres, les reflets et les réfractions.

Le RayTracer utilise des algorithmes de calcul d'intersection et des modèles de matériaux pour déterminer la couleur de chaque pixel de l'image finale. Les capacités de ce projet incluent la gestion de diverses formes géométriques, matériaux et sources lumineuses, ainsi que des techniques avancées telles que l'anticrénelage et les ombres douces.

## Fonctionnalités Clés

- **Géométrie 3D** : Support pour les sphères, les plans, les triangles et autres formes géométriques.
- **Matériaux** : Implémentation de matériaux diffus, réfléchissants et réfractifs.
- **Sources Lumineuses** : Gestion de différentes sources de lumière, y compris ambiantes, directionnelles et ponctuelles.
- **Effets Visuels** : Anticrénelage pour réduire les effets d'aliasing, et ombres douces pour un rendu plus réaliste.
- **Performance** : Optimisation des calculs pour des rendus efficaces même avec des scènes complexes.

## Technologies Utilisées

- **Langage** : C++11 pour des performances élevées et une gestion avancée des ressources.
- **Construction** : Utilisation de CMake pour une génération multiplateforme des fichiers de construction.
- **Imagerie** : Intégration de bibliothèques d'images comme stb_image pour le chargement et la sauvegarde d'images.

## Structure du Projet

Le projet est structuré en plusieurs composants principaux :

- **Core** : Contient les classes de base pour les vecteurs, les rayons, les intersections, etc.
- **Geometry** : Définit les formes géométriques prises en charge (sphères, plans, triangles).
- **Materials** : Gère les propriétés des matériaux et leurs interactions avec la lumière.
- **Renderer** : Implémente l'algorithme de traçage des rayons et gère le rendu final de l'image.
- **Scene** : Gère la composition des objets, des lumières et de la caméra dans la scène.
