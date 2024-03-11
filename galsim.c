#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define EPSILON_O 1e-3

// Global variables
double* pos_x, *pos_y, *mass, *vx, *vy, *fx, *fy;
double G, theta_0;

static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

// Structure for the tree
typedef struct node_t {
    int particle_index;
    int has_particle;
    int has_children;
    double min_x, max_x, min_y, max_y, node_mass, cm_x, cm_y;
    struct node_t* children;
} node_t;

void update_mass_and_cm(node_t* node, int particle);
node_t* get_child(node_t* node, int particle);
void subdivide(node_t* node);

/* ------------------------------
            Tree assembly
   ------------------------------ */

void insert(node_t* node, int particle) {
    // If the node is empty, insert the particle and set the mass and center of mass
    if (!node->has_particle) {
        node->particle_index = particle;
        node->cm_x = pos_x[particle];
        node->cm_y = pos_y[particle];
        node->node_mass = mass[particle];
        node->has_particle = 1;

    // If the node has children, insert the particle in the appropriate child
    } else if (node->has_children) {
        update_mass_and_cm(node, particle);
        insert(get_child(node, particle), particle);

    // If the node is a leaf, subdivide it and insert the old and new particles
    } else {
        node->has_children = 1;
        subdivide(node);
        insert(get_child(node, node->particle_index), node->particle_index); // Old
        insert(get_child(node, particle), particle); // New
        update_mass_and_cm(node, particle);
    }
}

void subdivide(node_t* node) {
    // A few variable names for readability
    double min_x = node->min_x;
    double max_x = node->max_x;
    double min_y = node->min_y;
    double max_y = node->max_y;
    double mid_x = (min_x + max_x) / 2;
    double mid_y = (min_y + max_y) / 2;

    // Allocate memory for the children
    node->children = malloc(4 * sizeof(node_t));

    // Initialize the children
    for (int i = 0; i < 4; i++) { 
        node->children[i].has_particle = 0;
        node->children[i].has_children = 0;
    }

    node->children[0].min_x = min_x;
    node->children[0].max_x = mid_x;
    node->children[0].min_y = min_y;
    node->children[0].max_y = mid_y;

    node->children[1].min_x = min_x;
    node->children[1].max_x = mid_x;
    node->children[1].min_y = mid_y;
    node->children[1].max_y = max_y;

    node->children[2].min_x = mid_x;
    node->children[2].max_x = max_x;
    node->children[2].min_y = min_y;
    node->children[2].max_y = mid_y;

    node->children[3].min_x = mid_x;
    node->children[3].max_x = max_x;
    node->children[3].min_y = mid_y;
    node->children[3].max_y = max_y;
}

// Function to get the child of a node that a particle belongs to
node_t* get_child(node_t* node, int particle) {
    double mid_x = (node->min_x + node->max_x) / 2;
    double mid_y = (node->min_y + node->max_y) / 2;
    // Clever way to get the index of the child: it compares pairs of bits of left and right outcomes
    // so that 0 | 0 = 0,   0 | 1 = 1,   2 | 0 = 2,    2 | 1 = 3
    int index = (pos_x[particle] >= mid_x) *2 | (pos_y[particle] >= mid_y);
    return &node->children[index];
}

void update_mass_and_cm(node_t* node, int particle) {
    double total_mass = node->node_mass + mass[particle];
    node->cm_x = (node->node_mass * node->cm_x + mass[particle] * pos_x[particle]) / total_mass;
    node->cm_y = (node->node_mass * node->cm_y + mass[particle] * pos_y[particle]) / total_mass;
    node->node_mass = total_mass;
}

/* ------------------------------
         Force calculation
   ------------------------------ */

void calculate_force(int particle_index, node_t* node) {
    double r_x = pos_x[particle_index] - node->cm_x;
    double r_y = pos_y[particle_index] - node->cm_y;
    double r_squared = r_x * r_x + r_y * r_y;
    double r_plummer = sqrt(r_squared) + EPSILON_O;
    double r_cubed = r_plummer * r_plummer * r_plummer;
    double force_factor = -G * mass[particle_index] * node->node_mass / r_cubed;

    fx[particle_index] += force_factor * r_x;
    fy[particle_index] += force_factor * r_y;
}

void update_forces(int particle, node_t* node) {
    // If the node has no children and has a particle that is not the one we are calculating the force on, calculate the force
    if (!node->has_children && node->has_particle && node->particle_index != particle) {
        calculate_force(particle, node);
    } 
    // If the node has children, check theta and decide whether to use it or its children
    else if (node->has_children) {
        double r_x = pos_x[particle] - node->cm_x;
        double r_y = pos_y[particle] - node->cm_y;
        double r_squared = r_x * r_x + r_y * r_y;
        double theta = (node->max_x - node->min_x) / sqrt(r_squared);

        if (theta < theta_0) {
            calculate_force(particle, node);
        } 
        else {
            for (int i = 0; i < 4; i++) {
                update_forces(particle, &node->children[i]);
            }
        }
    }
}

// Function to free the tree's children
void free_tree(node_t* node) {
    if (node->has_children) {
        for (int i = 0; i < 4; i++) {
            free_tree(&node->children[i]);
        }
        free(node->children);
    }
}

/* ------------------------------
            Main function
   ------------------------------ */

int main(int argc, char* argv[]) {
    double time_tol = get_wall_seconds();

    // Set number of threads
    omp_set_num_threads(8);

    if (argc != 6) {
        printf("You should enter the following parameters in order:\n");
        printf("N filename nsteps delta_t theta_0\n");
        return 1;
    }

    int N = atoi(argv[1]);
    char* filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    theta_0 = atof(argv[5]);

    FILE* data_file = fopen(filename, "rb");
    if (data_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    // Allocate memory
    pos_x = malloc(N * sizeof(double));
    pos_y = malloc(N * sizeof(double));
    mass = malloc(N * sizeof(double));
    vx = malloc(N * sizeof(double));
    vy = malloc(N * sizeof(double));
    fx = malloc(N * sizeof(double));
    fy = malloc(N * sizeof(double));

    double mass_inver[N];
    double ax[N];
    double ay[N];
    double brightness[N];
    G = 100.0 / N;

    // Read data from file
    for (int i = 0; i < N; i++) {
        fread(&pos_x[i], sizeof(double), 1 , data_file);
        fread(&pos_y[i], sizeof(double), 1 , data_file);
        fread(&mass[i], sizeof(double), 1 , data_file);
        mass_inver[i] = 1.0 / mass[i];
        fread(&vx[i], sizeof(double), 1 , data_file);
        fread(&vy[i], sizeof(double), 1 , data_file);
        fread(&brightness[i], sizeof(double), 1, data_file);
    }

    fclose(data_file);
    
    /* Calculate first step acceleration for Verlet Velocity */

    // Create root node
    node_t* root = malloc(sizeof(node_t));
    root->has_particle = 0;
    root->has_children = 0;
    root->min_x = 0.0;
    root->max_x = 1.0;
    root->min_y = 0.0;
    root->max_y = 1.0;

    // Initialize forces and insert particles
    for (int i = 0; i < N; i++) {
        fx[i] = 0;
        fy[i] = 0;
        insert(root, i);
    }

    // Calculate forces
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        update_forces(i, root);
    }

    // Free tree
    free_tree(root);
    free(root);
    
    // Time iterations
    for (int step = 0; step < nsteps; step++) {
        
        // Update position
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            ax[i] = fx[i] * mass_inver[i];
            ay[i] = fy[i] * mass_inver[i];
            pos_x[i] += delta_t * vx[i] + 0.5 * delta_t * delta_t * ax[i];
            pos_y[i] += delta_t * vy[i] + 0.5 * delta_t * delta_t * ay[i];
        }

        // Create root, initialize, insert particles and update forces
        root = malloc(sizeof(node_t));
        root->has_particle = 0;
        root->has_children = 0;
        root->min_x = 0.0;
        root->max_x = 1.0;
        root->min_y = 0.0;
        root->max_y = 1.0;

        for (int i = 0; i < N; i++) {
            fx[i] = 0;
            fy[i] = 0;
            insert(root, i);
        }

        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            update_forces(i, root);
        }

        // Update velocities
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            vx[i] += 0.5 * delta_t * (fx[i] * mass_inver[i] + ax[i]);
            vy[i] += 0.5 * delta_t * (fy[i] * mass_inver[i] + ay[i]);
        }

        free_tree(root);
        free(root);

    }

    FILE* rfile = fopen("result.gal", "w");
    if (rfile == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    
    // Write data to file
    for (int i = 0; i < N; i++) {
        fwrite(&pos_x[i], sizeof(double), 1 , rfile);
        fwrite(&pos_y[i], sizeof(double), 1 , rfile);
        fwrite(&mass[i], sizeof(double), 1 , rfile);
        fwrite(&vx[i], sizeof(double), 1 , rfile);
        fwrite(&vy[i], sizeof(double), 1 , rfile);
        fwrite(&brightness[i], sizeof(double), 1, rfile);
    }

    fclose(rfile);

    // Free memory
    free(pos_x);
    free(pos_y);
    free(mass);
    free(vx);
    free(vy);
    free(fx);
    free(fy);

    printf("The execution took %7.8f wall seconds.\n", get_wall_seconds() - time_tol);

    return 0;
}
