#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <stdbool.h>

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
    bool is_leaf;
    double x, y, width, node_mass, cm_x, cm_y;
    struct node_t* children[4];
} node_t;

void update_mass_and_cm(node_t* node, int particle);
void get_size(node_t* node, int particle, double* dimensions, int* index);

/* ------------------------------
            Tree assembly
   ------------------------------ */

node_t* create_node(double x, double y, double width) {
    node_t* node = malloc(sizeof(node_t));
    node->x = x;
    node->y = y;
    node->width = width;
    for (int i = 0; i < 4; i++) {
        node->children[i] = NULL;
    }
    return node;
}

void insert(node_t** node, int particle, double x, double y, double width) {
    // If there's no node, create it and add particle
    if (*node == NULL) {
        *node = create_node(x, y, width);
        (*node)->particle_index = particle;
        (*node)->cm_x = pos_x[particle];
        (*node)->cm_y = pos_y[particle];
        (*node)->node_mass = mass[particle];
        (*node)->is_leaf = true;

    // If there is, set leaf to false and add old and new particles
    } else {
        if ((*node)->is_leaf) {
            (*node)->is_leaf = false;
            int index;
            double dimensions[3]; // x, y, width
            get_size(*node, (*node)->particle_index, dimensions, &index);
            insert(&(*node)->children[index], (*node)->particle_index, dimensions[0], dimensions[1], dimensions[2]); // Old
        } 
        int index;
        double dimensions[3]; // x, y, width
        get_size(*node, particle, dimensions, &index);
        insert(&(*node)->children[index], particle, dimensions[0], dimensions[1], dimensions[2]); // New
        update_mass_and_cm(*node, particle);
    }
}

void get_size(node_t* node, int particle, double* dimensions, int* index) {
    double mid_x = node->x + node->width / 2;
    double mid_y = node->y + node->width / 2;
    // Clever way to get the index of the child
    *index = (pos_x[particle] >= mid_x) *2 | (pos_y[particle] >= mid_y);
    
    dimensions[2] = node->width / 2;
    dimensions[0] = (*index < 2) ? node->x : mid_x;
    dimensions[1] = (*index % 2 == 0) ? node->y : mid_y;
}


void update_mass_and_cm(node_t* node, int particle) {
    double total_mass = node->node_mass + mass[particle];
    node->cm_x = (node->node_mass * node->cm_x + mass[particle] * pos_x[particle]) / total_mass;
    node->cm_y = (node->node_mass * node->cm_y + mass[particle] * pos_y[particle]) / total_mass;
    node->node_mass = total_mass;
}

// Clean memory
void free_tree(node_t** node) {
    if (*node == NULL) {
        return;
    }
    if ((*node)->is_leaf == false) {
        for (int i = 0; i < 4; i++) {
            free_tree(&(*node)->children[i]);
        }
    }
    free(*node);
    *node = NULL;
}

/* ------------------------------
         Force calculation
   ------------------------------ */

void get_r(node_t* node, int particle, double* distances) {
    double r_x = pos_x[particle] - node->cm_x;
    double r_y = pos_y[particle] - node->cm_y;
    double r_squared = r_x * r_x + r_y * r_y;
    distances[0] = r_x;
    distances[1] = r_y;
    distances[2] = sqrt(r_squared);
}

void calculate_force(int particle_index, node_t* node, double distances[3]) {
    
    double r_plummer = distances[2] + EPSILON_O;
    double r_cubed = r_plummer * r_plummer * r_plummer;
    double force_factor = -G * mass[particle_index] * node->node_mass / r_cubed;

    fx[particle_index] += force_factor * distances[0];
    fy[particle_index] += force_factor * distances[1];
}

void update_forces(int particle, node_t* node) {
    // If the node has a particle that is not the one we are calculating the force on, calculate the force
    double distances[3]; // x, y, r
    get_r(node, particle, distances);
    
    if (node->is_leaf) {
        calculate_force(particle, node, distances);
    } else {
        double theta = node->width / distances[2];

        if (theta < theta_0) {
            calculate_force(particle, node, distances);
        } else {
            for (int i = 0; i < 4; i++) {
                if (node->children[i] != NULL) {
                    update_forces(particle, node->children[i]);
                }
            }
        }
    
    }
}

/* ------------------------------
            Main function
   ------------------------------ */

int main(int argc, char* argv[]) {
    double time_tol = get_wall_seconds();

    if (argc != 7) {
        printf("You should enter the following parameters in order:\n");
        printf("N filename nsteps delta_t theta_0 processes\n");
        return 1;
    }

    int N = atoi(argv[1]);
    char* filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    theta_0 = atof(argv[5]);
    int processes = atoi(argv[6]);

    omp_set_num_threads(processes);

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
    double half_dt_squared = 0.5 * delta_t * delta_t;

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
    node_t* root = NULL;

    // Initialize forces and insert particles
    for (int i = 0; i < N; i++) {
        fx[i] = 0;
        fy[i] = 0;
        insert(&root, i, 0, 0, 1);
    }

    // Calculate forces
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        update_forces(i, root);
    }

    // Free tree
    free_tree(&root);
    
    // Time iterations
    for (int step = 0; step < nsteps; step++) {
        
        // Update position
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            ax[i] = fx[i] * mass_inver[i];
            ay[i] = fy[i] * mass_inver[i];
            pos_x[i] += delta_t * vx[i] + half_dt_squared * ax[i];
            pos_y[i] += delta_t * vy[i] + half_dt_squared * ay[i];
        }

        // Create root node
        node_t* root = NULL;

        // Initialize forces and insert particles
        for (int i = 0; i < N; i++) {
            fx[i] = 0;
            fy[i] = 0;
            insert(&root, i, 0, 0, 1);
        }

        // Calculate forces
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

        free_tree(&root);

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
