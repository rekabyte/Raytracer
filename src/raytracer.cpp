#include "raytracer.h"

// Inspiration de mon code de ce merveilleux article ( qui a ete d'une grande aide <3 ):
// https://raytracing.github.io/books/RayTracingInOneWeekend.html

void Raytracer::render(const Scene& scene, Frame* output)
{
    // Crée le z_buffer.
    double *z_buffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        z_buffer[i] = scene.camera.z_far; //Anciennement DBL_MAX. À remplacer avec la valeur de scene.camera.z_far
    }

    // @@@@@@ VOTRE CODE ICI
	// Calculez les paramètres de la caméra pour les rayons.

    //calcul des vecteurs essentiels pour la camera:
    //tout les vecteurs de direction sont sont normalizees
    double3 cameraPos = scene.camera.position;
    double3 cameraDir = normalize(scene.camera.center - cameraPos);
    double3 right = normalize(cross(cameraDir, scene.camera.up));   // right = (direction cam . vector up de la cam)
    double3 up = normalize(cross(right, cameraDir));                // up = (camera right .  camera dir)
    double jitter_radius = scene.jitter_radius;

    //calcul des parametres du viewport:
    double fovy = scene.camera.fovy;
    double ratio = scene.camera.aspect;

    //** la hauteur du viewport est calcule selon la formule suivante: h = 2 * znear * tan(fov / 2)
    double vw_height = 2 * scene.camera.z_near * tan(deg2rad(fovy / 2));
    double vp_width = vw_height * ratio;

    // vw_bottomLeft = point de depart des rayons
    //** calculee a partir de la formule: cameraDirection * znear + cameraPosition - cameraRight * (vw half width) - up * (vw half height)
    double3 vw_bottomLeft = cameraDir * scene.camera.z_near + cameraPos - right * (vp_width / 2.0) - up * (vw_height / 2.0);


    // Itère sur tous les pixels de l'image.
    for (int y = 0; y < scene.resolution[1]; y++) {
        if (y % 40){
            std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
        }

        for(int x = 0; x < scene.resolution[0]; x++) {

            int avg_z_depth = 0;
            double3 avg_ray_color{0,0,0};
            
            for(int iray = 0; iray < scene.samples_per_pixel; iray++) {
                // Génère le rayon approprié pour ce pixel.
                Ray ray;
                // Initialise la profondeur de récursivité du rayon.
                int ray_depth = 0;
                // Initialize la couleur du rayon
                double3 ray_color{0,0,0};

                // @@@@@@ VOTRE CODE ICI
				// Mettez en place le rayon primaire en utilisant les paramètres de la caméra.
				// Lancez le rayon de manière uniformément aléatoire à l'intérieur du pixel dans la zone délimité par jitter_radius. 
				//Faites la moyenne des différentes couleurs obtenues suite à la récursion.

                //genere vecteur 2d du jitter avec x et y qui sont randomisees pour avoir un bel echantillonage:
                double2 jitter_offset = random_in_unit_disk() * jitter_radius;

                //calcul de la largeur/hauteur d'un pixel
                double pixelWidth = vp_width / scene.resolution[0];
                double pixelHeight = vw_height / scene.resolution[1];

                //coordonnes du pixel = vw_bottomLeft + (cameraRight * x * pixelWidth) + (cameraUp * y * pixelHeight)
                // pour integrer la randomisation du jitter,
                // on remplace x et y respectivement par (x + jitter_offset.x) et (y + jitter_offset.y)
                double3 viewportPixelCoord = vw_bottomLeft 
                                            + right * (x + jitter_offset.x) * pixelWidth 
                                            + up * (y + jitter_offset.y) * pixelHeight;


                //Parametres des rayons
                //(origine = cameraPosition ET direction = vecteur qui passe par la pos de la camera et traverse le pixel actuel):
                ray.origin = cameraPos;
                ray.direction = normalize(viewportPixelCoord - cameraPos);

                //recursion
                double z_depth = scene.camera.z_far; //depth max de la scene
                trace(scene, ray, ray_depth, &ray_color, &z_depth);

                //moyenne des différentes couleurs obtenues suite à la récursion
                avg_ray_color += ray_color;
                avg_z_depth += (int)z_depth;
            }

            avg_z_depth = avg_z_depth / scene.samples_per_pixel;
            avg_ray_color = avg_ray_color / scene.samples_per_pixel;

            // Test de profondeur
            if(avg_z_depth >= scene.camera.z_near && avg_z_depth <= scene.camera.z_far && 
                avg_z_depth < z_buffer[x + y*scene.resolution[0]]) {
                z_buffer[x + y*scene.resolution[0]] = avg_z_depth;

                // Met à jour la couleur de l'image (et sa profondeur)
                output->set_color_pixel(x, y, avg_ray_color);
                output->set_depth_pixel(x, y, (avg_z_depth - scene.camera.z_near) / (scene.camera.z_far-scene.camera.z_near));
            }
        }
    }

    delete[] z_buffer;
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		- Détermine si le rayon intersecte la géométrie.
//      	- Calculer la contribution associée à la réflexion.
//			- Calculer la contribution associée à la réfraction.
//			- Mettre à jour la couleur avec le shading + 
//			  Ajouter réflexion selon material.reflection +
//			  Ajouter réfraction selon material.refraction 
//            pour la couleur de sortie.
//          - Mettre à jour la nouvelle profondeur.
void Raytracer::trace(const Scene& scene, Ray ray, int ray_depth, double3* out_color, double* out_z_depth)
{
    Intersection hit;

    // Fait appel à l'un des containers spécifiées.
    if(scene.container->intersect(ray,EPSILON,*out_z_depth,&hit)) {       
        Material& material = ResourceManager::Instance()->materials[hit.key_material];

        // @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réflection d'un rayon de manière récursive
        bool isRayDepth = (ray_depth <= 1) ? true : false;

        double3 reflection_color = {0,0,0};
        if(isRayDepth) {
            Ray reflection;

            //calcul des parametres du ray de reflection:
            //reflection calculee comme suit: reflection = incident - 2 * (normal . incident) * normal
            //avec incident = ray.direction
            reflection.origin = hit.position;
            reflection.direction = normalize(ray.direction - 2 * (dot(hit.normal, ray.direction) * hit.normal));

            //recursion:
            double reflection_ray_z = scene.camera.z_far;
            trace(scene, reflection, ray_depth+1, &reflection_color, &reflection_ray_z);
        }

        // @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réfraction d'un rayon de manière récursive.
		// 
		// Assumez que l'extérieur/l'air a un indice de réfraction de 1.
		//
		// Toutes les géométries sont des surfaces et non pas de volumes.
        
        
        *out_color = shade(scene, hit) + reflection_color * material.k_reflection;
        *out_z_depth = hit.depth;
    }
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		* Calculer la contribution des lumières dans la scène.
//			- Itérer sur toutes les lumières.
//				- Inclure la contribution spéculaire selon le modèle de Blinn en incluant la composante métallique.
//	          	- Inclure la contribution diffuse. (Faites attention au produit scalare. >= 0)
//   	  	- Inclure la contribution ambiante
//      * Calculer si le point est dans l'ombre
//			- Itérer sur tous les objets et détecter si le rayon entre l'intersection et la lumière est occludé.
//				- Ne pas considérer les points plus loins que la lumière.
//			- Par la suite, intégrer la pénombre dans votre calcul
//		* Déterminer la couleur du point d'intersection.
//        	- Si texture est présente, prende la couleur à la coordonnées uv
//			- Si aucune texture, prendre la couleur associé au matériel.

double3 Raytracer::shade(const Scene& scene, Intersection hit)
{
    Material& material = ResourceManager::Instance()->materials[hit.key_material]; //lorsque vous serez rendu à la partie texture.

    //toutes les variables pour le calcul de la couleur finale du pixel:
    double3 lightDir;
    double3 diffuseComponent = {0.0,0.0,0.0};
    double3 specularComponent = {0.0,0.0,0.0};
    double3 lightIntensity = {0.0,0.0,0.0};
    double3 color = material.color_albedo;

    //y-a-t-il une texture?
    bool isTextured = (material.texture_albedo.width() > 0 && material.texture_albedo.height() > 0) ? true : false;

    // si une texture est présente, on récupère la couleur aux coordonées UV
    if (isTextured) {
        //calcul des coord UV:
        int tex_width = material.texture_albedo.width()-1;
        int tex_height = material.texture_albedo.height()-1;
        int U = (hit.uv.x * tex_width);
        int V = (hit.uv.y * tex_height);

        //echantillonnage de la texture:
        unsigned char R, G, B;
        material.texture_albedo.get_pixel(U, V, R, G, B);

        //normalisation de couleur, pour passer de [0..255] a [0..1]
        double3 textureColor = { R / 255.0, G / 255.0, B / 255.0 };
        color = textureColor;
    }

    //Calculer la contribution des lumières dans la scène.
    //			- Itérer sur toutes les lumières.

    double3 finalLight = {0.0,0.0,0.0}; // lumiere finale
    

    for(SphericalLight light : scene.lights) {

        //Calcule des parametres de la lumiere:
        double radius = light.radius;

        //direction de la lumiere = vecteur qui commencer par la position de l'intersection 
        // et qui passe par la position de la lumiere
        lightDir = light.position - hit.position;
        double light_depth = length(lightDir);

        //on peut normalizer la direction mtn que light_depth a ete calculee
        lightDir = normalize(lightDir);

        //contribution initiale:
        double light_contribution = 1.0;
        

        //*** Calculer si le point est dans l'ombre
        //shawowRay = commencer a la position du hit et se dirige vers la lumiere
        Ray shadowRay;
        shadowRay.origin = hit.position;
        shadowRay.direction = lightDir;

        //intersection entre le vecteur de shadow et l'objet stockee:
        Intersection shadow_hit;

        //est-ce que la lumiere a un rayon?
        bool isLightDefined = (light.radius > 0.0) ? true : false;

        //si oui:
        if(isLightDefined) {
            //nbre de shadow samples, nbre plus grand => shadows plus realistes mais bcp plus de calculs
            //todo a changer plus tard
            int shadow_ray_samples = 5;
            double hitCount = 0;

            for(int i = 0; i < shadow_ray_samples; i++) {
                
                //calcul d'un offset en x,y comme pour le jitter_radius
                double2 jitter_offset = random_in_unit_disk() * light.radius;

                //generer des points d'echantillonnage à l'interieur de la source lumineuse afin de calculer les soft shadows:
                // avec: cross(lightDir, {0, 1, 0}) et cross(lightDir, {1, 0, 0})
                // etant les vecteurs perpendiculaires a l'axe x et y respectivement
                // ils sont ensuites multipliee par le offset du jitter pour avoir une randomization uniforme
                double3 sampledPoint = light.position 
                                    + jitter_offset.x * normalize(cross(lightDir, {0, 1, 0}))
                                    + jitter_offset.y * normalize(cross(lightDir, {1, 0, 0}));

                //distance entre l'intersection et le point echantillonnee:
                double samplePointDistance = length(sampledPoint - hit.position);

                //on verifie si l'ombre intersecte avec un objet/geometrie, si oui, on augmente le nombre de hitCount,
                // qui va servir plus tard pour calculer la penombre
                shadowRay.origin = hit.position;
                shadowRay.direction = normalize(sampledPoint - hit.position);
                Intersection shadow_hit;

                if(scene.container->intersect(shadowRay, EPSILON, samplePointDistance, &shadow_hit)) {
                    hitCount++;
                }
            }

            //Pour la penombre:
            //(Au lieu de simplement mettre la contribution à 0, remplacez ce dernier avec le facteur d'occlusion.)
            //(inclure un facteur d'occlusion en fonction du nombre de rayons qui atteignent la lumière)
            light_contribution = 1.0 - hitCount / shadow_ray_samples;
        }
        //si non:
        // Si le point s'avère être dans l'ombre, le terme d'éclairage devrait être de 0 pour cette source de lumière.
        else{
            if(scene.container->intersect(shadowRay,EPSILON,light_depth,&shadow_hit)) light_contribution = 0.0;
        }

        //On commence a calculer la lumiere finale ici:
        double3 viewDir = normalize(scene.camera.position - hit.position);
        double3 halfVec = normalize(viewDir + lightDir);

        //le coeff de Lambert est calcule selon la formule: k_lambert = max(0, N . L)
        double k_lambert = std::max(0.0, dot(hit.normal, lightDir));

        //le coeff de specularite est calcule selon la formule: k_specular = max(0, N . H)^(brillance)
        //avec H = halfway vector
        double specularCoef = pow(std::max(0.0, dot(hit.normal, halfVec)), material.shininess);

        //lumiere diffuse = coeff_diffuse * color * coeff_lambert
        diffuseComponent = material.k_diffuse * color * k_lambert;

        //specularLight = coeff_specular × (materialReflectance × color + (1 − materialReflectance ) ) × coeff_specular
        specularComponent = material.k_specular * (material.metallic * color + (1 - material.metallic))* specularCoef;

        //lightIntensity = (diffuse + specular) * (emissionCoeff / (lightLength)^2)
        lightIntensity = (diffuseComponent + specularComponent)* light.emission/pow(light_depth,2);

        //la lumiere finale correspond a l'intesition multiplie par le facteur de contribution
        finalLight += lightIntensity * light_contribution;
    }
    
    //Apres le calcul de toutes les lumieres, on peut calculer la lumiere ambiente avec la formule:
    // ambientIntensity = ambientLight * materialReflectivity * surfaceColor
    double3 ambientIntensity = scene.ambient_light * material.k_ambient * color;

    finalLight += ambientIntensity;

    return finalLight;
}