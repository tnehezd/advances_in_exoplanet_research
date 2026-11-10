#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <yaml.h>

/*
 * A deliberately small 1D, two-fluid teaching model.
 *
 * The gas surface density evolves by viscous diffusion. The dust surface
 * density is advected by the radial drift speed calculated from the current
 * gas state. Both fields use the same time loop and the same final time.
 *
 * The dust feels drag, but its drag force is not fed back into the gas
 * velocity. The Stokes number is prescribed and constant. The YAML reader
 * below is configuration plumbing; the main physics fields are the arrays
 * sigma_g, pressure, dpdr, sigma_d, v_g, and v_d.
 *
 * Radius is stored as r/r0, time is in years, and radial speeds are in
 * r0/year. The reference radius in AU and stellar mass in solar masses set
 * Omega0 through Kepler's third law. Sigma and pressure use relative teaching
 * units; the pressure normalization cancels in (H/P)*(dP/dr).
 */

typedef struct {
    double r_inner;
    double r_outer;
    int n_grid;
    double reference_radius_au;
    double stellar_mass_msun;
    double omega0_per_year;
    double h;
    double alpha_active;
    double dead_inner;
    double dead_outer;
    double dead_reduction;
    double edge_width;
    double t_end;
    double stokes_number;
    double initial_dust_to_gas;
} Config;

typedef struct {
    double *x;
    double *r;
    double *r_edge;
    double *alpha;
    double *nu;              /* Viscosity in (r0)^2/year. */
    double *sigma_g;
    double *sigma_g_next;
    double *sigma_g_initial;
    double *sigma_g_middle;
    double *pressure;        /* Relative midplane-pressure profile. */
    double *dpdr;            /* Pressure proxy gradient per unit r/r0. */
    double *gas_velocity;    /* Radial gas speed in r0/year. */
    double *sigma_d;
    double *sigma_d_next;
    double *sigma_d_initial;
    double *sigma_d_middle;
    double *dust_velocity;   /* Radial dust speed in r0/year. */
    double *dust_face_velocity;
    double *q;
    double *dust_flux;
} Model;

/* Find a named child in a YAML mapping, for example grid.r_inner. */
static yaml_node_t *mapping_value(yaml_document_t *document,
                                  yaml_node_t *mapping, const char *key) {
    if (mapping == NULL || mapping->type != YAML_MAPPING_NODE) {
        return NULL;
    }
    /* Check each key in this YAML mapping until the requested name is found. */
    for (yaml_node_pair_t *pair = mapping->data.mapping.pairs.start;
         pair < mapping->data.mapping.pairs.top; ++pair) {
        yaml_node_t *key_node = yaml_document_get_node(document, pair->key);
        if (key_node != NULL && key_node->type == YAML_SCALAR_NODE &&
            strcmp((const char *)key_node->data.scalar.value, key) == 0) {
            return yaml_document_get_node(document, pair->value);
        }
    }
    return NULL;
}

/* Read one numeric YAML scalar and report its key if it is absent or invalid. */
static int read_number(yaml_document_t *document, yaml_node_t *section,
                       const char *key, double *value) {
    yaml_node_t *node = mapping_value(document, section, key);
    if (node == NULL || node->type != YAML_SCALAR_NODE) {
        fprintf(stderr, "Missing numeric YAML setting: %s\n", key);
        return 0;
    }

    /* strtod converts the text and also lets us detect invalid trailing text. */
    errno = 0;
    char *end = NULL;
    const double parsed = strtod((const char *)node->data.scalar.value, &end);
    if (errno != 0 || end == (char *)node->data.scalar.value ||
        *end != '\0' || !isfinite(parsed)) {
        fprintf(stderr, "Invalid numeric YAML setting: %s\n", key);
        return 0;
    }
    *value = parsed;
    return 1;
}

/* Grid size is an integer; reject fractional or unreasonable array lengths. */
static int read_integer(yaml_document_t *document, yaml_node_t *section,
                        const char *key, int *value) {
    double number;
    if (!read_number(document, section, key, &number) ||
        number < 1.0 || number > 10000000.0 || floor(number) != number) {
        fprintf(stderr, "YAML setting %s must be a positive integer.\n", key);
        return 0;
    }
    *value = (int)number;
    return 1;
}

/* Read all nested settings, then check physical and numerical constraints. */
static int read_config(const char *path, Config *config) {
    FILE *input = fopen(path, "r");
    if (input == NULL) {
        perror(path);
        return 0;
    }

    /* Ask libyaml to read the whole YAML file into a document tree. */
    yaml_parser_t parser;
    yaml_document_t document;
    if (!yaml_parser_initialize(&parser)) {
        fprintf(stderr, "Could not initialize the YAML parser.\n");
        fclose(input);
        return 0;
    }
    yaml_parser_set_input_file(&parser, input);
    if (!yaml_parser_load(&parser, &document)) {
        fprintf(stderr, "YAML parse error at line %zu, column %zu: %s\n",
                parser.problem_mark.line + 1, parser.problem_mark.column + 1,
                parser.problem != NULL ? parser.problem : "unknown error");
        yaml_parser_delete(&parser);
        fclose(input);
        return 0;
    }
    yaml_parser_delete(&parser);
    fclose(input);

    /*
     * The settings are grouped into mappings such as "grid", "disk", and
     * "dust". A missing section or value makes its read_number/read_integer
     * call fail, so the final result is false if any required setting is bad.
     */
    yaml_node_t *root = yaml_document_get_root_node(&document);
    yaml_node_t *grid = mapping_value(&document, root, "grid");
    yaml_node_t *disk = mapping_value(&document, root, "disk");
    yaml_node_t *dead_zone = mapping_value(&document, root, "dead_zone");
    yaml_node_t *evolution = mapping_value(&document, root, "evolution");
    yaml_node_t *dust = mapping_value(&document, root, "dust");

    const int ok =
        read_number(&document, grid, "r_inner", &config->r_inner) &&
        read_number(&document, grid, "r_outer", &config->r_outer) &&
        read_integer(&document, grid, "n_grid", &config->n_grid) &&
        read_number(&document, disk, "reference_radius_au",
                    &config->reference_radius_au) &&
        read_number(&document, disk, "stellar_mass_msun",
                    &config->stellar_mass_msun) &&
        read_number(&document, disk, "h", &config->h) &&
        read_number(&document, disk, "alpha_active", &config->alpha_active) &&
        read_number(&document, dead_zone, "r_inner", &config->dead_inner) &&
        read_number(&document, dead_zone, "r_outer", &config->dead_outer) &&
        read_number(&document, dead_zone, "reduction",
                    &config->dead_reduction) &&
        read_number(&document, dead_zone, "edge_width",
                    &config->edge_width) &&
        read_number(&document, evolution, "t_end", &config->t_end) &&
        read_number(&document, dust, "stokes_number",
                    &config->stokes_number) &&
        read_number(&document, dust, "initial_dust_to_gas",
                    &config->initial_dust_to_gas);
    yaml_document_delete(&document);
    if (!ok) {
        return 0;
    }

    /* Reject physically or numerically invalid values before allocating arrays. */
    if (config->n_grid < 3 || config->r_inner <= 0.0 ||
        config->r_outer <= config->r_inner || config->h <= 0.0 ||
        config->reference_radius_au <= 0.0 ||
        config->stellar_mass_msun <= 0.0 || config->alpha_active <= 0.0 ||
        config->dead_inner <= 0.0 ||
        config->dead_outer <= config->dead_inner ||
        config->dead_reduction <= 0.0 || config->dead_reduction > 1.0 ||
        config->edge_width <= 0.0 || config->t_end <= 0.0 ||
        config->stokes_number < 0.0 ||
        config->initial_dust_to_gas < 0.0) {
        fprintf(stderr, "YAML settings are outside their valid ranges.\n");
        return 0;
    }
    /*
     * Kepler's third law in AU, solar-mass, and year units:
     * P0 = sqrt(r0^3/Mstar) years, so Omega0 = 2*pi/P0 per year.
     */
    config->omega0_per_year =
        6.2831853071795864769 *
        sqrt(config->stellar_mass_msun /
             (config->reference_radius_au * config->reference_radius_au *
              config->reference_radius_au));
    return 1;
}

/*
 * Return the viscosity parameter alpha at radius r.
 * The two tanh factors turn the dead zone on at its inner edge and off at
 * its outer edge. Their product is close to 1 inside the dead zone and close
 * to 0 outside it.
 */
static double dead_zone_alpha(double r, const Config *config) {
    const double inner_edge = 0.5 *
        (1.0 + tanh((r - config->dead_inner) / config->edge_width));
    const double outer_edge = 0.5 *
        (1.0 + tanh((config->dead_outer - r) / config->edge_width));
    const double inside_dead_zone = inner_edge * outer_edge;
    const double alpha_drop =
        (1.0 - config->dead_reduction) * inside_dead_zone;
    return config->alpha_active * (1.0 - alpha_drop);
}

/*
 * Allocate the arrays used by the model.
 * Each array has one value per radial grid cell, except r_edge and the
 * dust face arrays, which have one extra value because they live between cells.
 */
static int allocate_model(Model *model, int n) {
    model->x = calloc((size_t)n, sizeof(double));
    model->r = calloc((size_t)n, sizeof(double));
    model->r_edge = calloc((size_t)n + 1, sizeof(double));
    model->alpha = calloc((size_t)n, sizeof(double));
    model->nu = calloc((size_t)n, sizeof(double));
    model->sigma_g = calloc((size_t)n, sizeof(double));
    model->sigma_g_next = calloc((size_t)n, sizeof(double));
    model->sigma_g_initial = calloc((size_t)n, sizeof(double));
    model->sigma_g_middle = calloc((size_t)n, sizeof(double));
    model->pressure = calloc((size_t)n, sizeof(double));
    model->dpdr = calloc((size_t)n, sizeof(double));
    model->gas_velocity = calloc((size_t)n, sizeof(double));
    model->sigma_d = calloc((size_t)n, sizeof(double));
    model->sigma_d_next = calloc((size_t)n, sizeof(double));
    model->sigma_d_initial = calloc((size_t)n, sizeof(double));
    model->sigma_d_middle = calloc((size_t)n, sizeof(double));
    model->dust_velocity = calloc((size_t)n, sizeof(double));
    model->dust_face_velocity = calloc((size_t)n + 1, sizeof(double));
    model->q = calloc((size_t)n, sizeof(double));
    model->dust_flux = calloc((size_t)n + 1, sizeof(double));

    return model->x != NULL && model->r != NULL && model->r_edge != NULL &&
           model->alpha != NULL && model->nu != NULL &&
           model->sigma_g != NULL && model->sigma_g_next != NULL &&
           model->sigma_g_initial != NULL && model->sigma_g_middle != NULL &&
           model->pressure != NULL && model->dpdr != NULL &&
           model->gas_velocity != NULL && model->sigma_d != NULL &&
           model->sigma_d_next != NULL && model->sigma_d_initial != NULL &&
           model->sigma_d_middle != NULL && model->dust_velocity != NULL &&
           model->dust_face_velocity != NULL && model->q != NULL &&
           model->dust_flux != NULL;
}

/* Free all arrays; free(NULL) is safe if allocation stopped part-way through. */
static void free_model(Model *model) {
    free(model->x);
    free(model->r);
    free(model->r_edge);
    free(model->alpha);
    free(model->nu);
    free(model->sigma_g);
    free(model->sigma_g_next);
    free(model->sigma_g_initial);
    free(model->sigma_g_middle);
    free(model->pressure);
    free(model->dpdr);
    free(model->gas_velocity);
    free(model->sigma_d);
    free(model->sigma_d_next);
    free(model->sigma_d_initial);
    free(model->sigma_d_middle);
    free(model->dust_velocity);
    free(model->dust_face_velocity);
    free(model->q);
    free(model->dust_flux);
}

/*
 * Set up the grid and initial disk.
 * The grid is evenly spaced in x = sqrt(r), not in r itself. The initial gas
 * surface density is Sigma_g = 1/r; the initial dust is a fixed fraction of it.
 */
static void initialize_model(Model *model, const Config *config) {
    const int n = config->n_grid;
    const double x_inner = sqrt(config->r_inner);
    const double x_outer = sqrt(config->r_outer);
    const double dx = (x_outer - x_inner) / (n - 1);

    for (int i = 0; i < n; ++i) {
        const double x = x_inner + i * dx;
        const double r = x * x;
        const double alpha = dead_zone_alpha(r, config);
        const double omega_k =
            config->omega0_per_year * pow(r, -1.5);
        const double viscosity =
            alpha * config->h * config->h * r * r * omega_k;
        const double sigma_g = 1.0 / r;

        model->x[i] = x;
        model->r[i] = r;
        model->alpha[i] = alpha;
        model->nu[i] = viscosity;
        model->sigma_g[i] = sigma_g;
        model->sigma_d[i] = config->initial_dust_to_gas * sigma_g;
    }

    /*
     * Copy the nearest interior value to each boundary.
     * This approximates zero radial gradient at the edge of the modeled disk.
     */
    model->sigma_g[0] = model->sigma_g[1];
    model->sigma_g[n - 1] = model->sigma_g[n - 2];
    model->sigma_d[0] = model->sigma_d[1];
    model->sigma_d[n - 1] = model->sigma_d[n - 2];

    /*
     * Dust is stored at cell centers, but its flux crosses cell edges.
     * These edges are also used to calculate each cell's annular area.
     */
    model->r_edge[0] = model->r[0];
    model->r_edge[n] = model->r[n - 1];
    for (int i = 1; i < n; ++i) {
        model->r_edge[i] = 0.5 * (model->r[i - 1] + model->r[i]);
    }

    /* Save the starting profiles; these arrays are not changed by the solver. */
    for (int i = 0; i < n; ++i) {
        model->sigma_g_initial[i] = model->sigma_g[i];
        model->sigma_d_initial[i] = model->sigma_d[i];
        model->sigma_g_middle[i] = model->sigma_g[i];
        model->sigma_d_middle[i] = model->sigma_d[i];
    }
}

/*
 * Calculate pressure, its radial gradient, and the gas radial velocity.
 * For this flat, locally isothermal disk, P is proportional to Sigma_g*h/r^2.
 * The helper q = nu*Sigma_g*sqrt(r) is used in the viscous gas velocity.
 */
static void calculate_gas_fields(Model *model, const Config *config) {
    const int n = config->n_grid;
    for (int i = 0; i < n; ++i) {
        const double r = model->r[i];
        const double sigma_g = model->sigma_g[i];
        model->pressure[i] = sigma_g * config->h / (r * r);
        model->q[i] = model->nu[i] * sigma_g * model->x[i];
    }

    /*
     * Use centered differences for dP/dr and dq/dr:
     * (f[i+1] - f[i-1]) / (r[i+1] - r[i-1]).
     */
    for (int i = 1; i < n - 1; ++i) {
        const double dr = model->r[i + 1] - model->r[i - 1];
        const double pressure_difference =
            model->pressure[i + 1] - model->pressure[i - 1];
        const double dqdr = (model->q[i + 1] - model->q[i - 1]) / dr;
        model->dpdr[i] = pressure_difference / dr;
        model->gas_velocity[i] = -3.0 * dqdr /
                                 (model->sigma_g[i] * model->x[i]);
    }

    /* At the two edges, reuse the nearest interior gradient and velocity. */
    model->dpdr[0] = model->dpdr[1];
    model->dpdr[n - 1] = model->dpdr[n - 2];
    model->gas_velocity[0] = model->gas_velocity[1];
    model->gas_velocity[n - 1] = model->gas_velocity[n - 2];
}

/*
 * Calculate the dust radial velocity from the current gas fields.
 * This RAPID-style expression uses directly visible P and dP/dr:
 * v_d = v_g/(1+St^2) + St/(1+St^2) * (H/P)*(dP/dr)*c_s.
 * Since dP/dr is normally negative, the pressure term points inward.
 */
static void calculate_dust_velocity(Model *model, const Config *config) {
    const double st = config->stokes_number;
    const double denominator = 1.0 + st * st;
    for (int i = 0; i < config->n_grid; ++i) {
        const double r = model->r[i];
        const double scale_height = config->h * r;
        const double sound_speed =
            config->h * config->omega0_per_year / sqrt(r);
        const double gas_part = model->gas_velocity[i] / denominator;
        const double pressure_part =
            st / denominator * (scale_height / model->pressure[i]) *
            model->dpdr[i] * sound_speed;
        model->dust_velocity[i] = gas_part + pressure_part;
    }
}

/*
 * Estimate the largest stable timestep for the explicit gas-diffusion update.
 * A finer grid (smaller dx) or a larger viscosity requires a smaller timestep.
 * The factor 0.35 is a conservative safety factor for this FTCS scheme.
 */
static double gas_stable_timestep(const Model *model, const Config *config) {
    double max_nu = 0.0;
    for (int i = 0; i < config->n_grid; ++i) {
        if (model->nu[i] > max_nu) {
            max_nu = model->nu[i];
        }
    }

    const double dx = model->x[1] - model->x[0];
    const double x_inner = model->x[0];
    const double safety_factor = 0.35;
    const double numerator = safety_factor * dx * dx * x_inner * x_inner;
    const double denominator = 3.0 * max_nu;
    return numerator / denominator;
}

/*
 * Estimate the largest stable timestep for dust advection.
 * First estimate the dust velocity at each cell edge. Then calculate how
 * quickly dust can leave each cell and keep the most restrictive timestep.
 */
static double dust_stable_timestep(Model *model, const Config *config) {
    const int n = config->n_grid;

    /* At the two outer edges, use the velocity of the nearest cell. */
    model->dust_face_velocity[0] = model->dust_velocity[0];
    model->dust_face_velocity[n] = model->dust_velocity[n - 1];

    /* Between cells, take the average of the two neighboring velocities. */
    for (int i = 1; i < n; ++i) {
        model->dust_face_velocity[i] =
            0.5 * (model->dust_velocity[i - 1] + model->dust_velocity[i]);
    }

    double dt_limit = INFINITY;
    for (int i = 0; i < n; ++i) {
        /* Area of this annular cell, omitting the common factor 2*pi. */
        const double r_left = model->r_edge[i];
        const double r_right = model->r_edge[i + 1];
        const double annular_area = 0.5 * (r_right * r_right - r_left * r_left);

        /*
         * Add the outward flow through the cell's two edges.
         * A positive velocity at the outer edge carries dust outward.
         * A negative velocity at the inner edge also carries dust outward.
         */
        double outward_flux_rate = 0.0;
        if (model->dust_face_velocity[i + 1] > 0.0) {
            outward_flux_rate +=
                r_right * model->dust_face_velocity[i + 1];
        }
        if (model->dust_face_velocity[i] < 0.0) {
            outward_flux_rate +=
                r_left * -model->dust_face_velocity[i];
        }

        if (outward_flux_rate > 0.0) {
            /* Keep this step to at most 40% of the cell's outflow time. */
            const double cell_dt = 0.4 * annular_area / outward_flux_rate;
            if (cell_dt < dt_limit) {
                dt_limit = cell_dt;
            }
        }
    }
    return dt_limit;
}

/*
 * Advance the gas surface density by one explicit FTCS step.
 * Calculate every new value into sigma_g_next first; only after the whole
 * profile is valid do we copy it back to sigma_g.
 */
static int advance_gas(Model *model, const Config *config, double dt) {
    const int n = config->n_grid;
    const double dx = model->x[1] - model->x[0];

    /* The boundary cells are set separately, so update interior cells first. */
    for (int i = 1; i < n - 1; ++i) {
        const double second_derivative =
            (model->q[i + 1] - 2.0 * model->q[i] + model->q[i - 1]) /
            (dx * dx);

        /* This is the discretized gas diffusion equation for this model. */
        const double x = model->x[i];
        const double x_cubed = x * x * x;
        const double change =
            dt * 3.0 * second_derivative / (4.0 * x_cubed);
        model->sigma_g_next[i] = model->sigma_g[i] + change;

        if (!isfinite(model->sigma_g_next[i]) ||
            model->sigma_g_next[i] <= 0.0) {
            fprintf(stderr, "Gas Sigma became non-positive at cell %d.\n", i);
            return 0;
        }
    }

    /* Copy the nearest interior values to the two boundaries. */
    model->sigma_g_next[0] = model->sigma_g_next[1];
    model->sigma_g_next[n - 1] = model->sigma_g_next[n - 2];

    /* The update is complete, so replace the old gas profile with the new one. */
    for (int i = 0; i < n; ++i) {
        model->sigma_g[i] = model->sigma_g_next[i];
    }
    return 1;
}

/*
 * Advance Sigma_d using conservative upwind fluxes:
 * dSigma_d/dt = -(1/r) d(r*Sigma_d*v_d)/dr.
 * Inflow from outside the modeled domain is set to zero; outward drift may
 * leave through either boundary. This update conserves dust mass up to those
 * explicitly allowed boundary losses.
 */
static int advance_dust(Model *model, const Config *config, double dt) {
    const int n = config->n_grid;

    /*
     * Calculate fluxes through internal cell edges.
     * Use the dust density from the side the dust is coming from (upwind).
     */
    for (int face = 1; face < n; ++face) {
        const double velocity = model->dust_face_velocity[face];
        double sigma_upwind;
        if (velocity >= 0.0) {
            sigma_upwind = model->sigma_d[face - 1];
        } else {
            sigma_upwind = model->sigma_d[face];
        }
        model->dust_flux[face] =
            model->r_edge[face] * sigma_upwind * velocity;
    }

    /*
     * At the inner edge, allow dust to leave inward, but do not add dust
     * flowing in from outside the grid.
     */
    model->dust_flux[0] = 0.0;
    if (model->dust_face_velocity[0] < 0.0) {
        model->dust_flux[0] =
            model->r_edge[0] * model->sigma_d[0] *
            model->dust_face_velocity[0];
    }

    /*
     * At the outer edge, allow dust to leave outward, but do not add dust
     * flowing in from beyond the grid.
     */
    model->dust_flux[n] = 0.0;
    if (model->dust_face_velocity[n] > 0.0) {
        model->dust_flux[n] =
            model->r_edge[n] * model->sigma_d[n - 1] *
            model->dust_face_velocity[n];
    }

    /* Update each cell from the difference between its outgoing and incoming flux. */
    for (int i = 0; i < n; ++i) {
        const double r_left = model->r_edge[i];
        const double r_right = model->r_edge[i + 1];
        const double annular_area = 0.5 * (r_right * r_right - r_left * r_left);
        const double flux_difference =
            model->dust_flux[i + 1] - model->dust_flux[i];
        const double change = dt * flux_difference / annular_area;
        model->sigma_d_next[i] = model->sigma_d[i] - change;

        if (!isfinite(model->sigma_d_next[i]) ||
            model->sigma_d_next[i] < -1e-12) {
            fprintf(stderr, "Dust Sigma became negative at cell %d.\n", i);
            return 0;
        }
        if (model->sigma_d_next[i] < 0.0) {
            model->sigma_d_next[i] = 0.0;
        }
    }

    /* All cells have been updated; now replace the old dust profile. */
    for (int i = 0; i < n; ++i) {
        model->sigma_d[i] = model->sigma_d_next[i];
    }
    return 1;
}

/*
 * Save the initial, halfway, and final profile.
 * The dust flag selects Sigma_d; otherwise the function writes Sigma_g.
 */
static int write_snapshots(const char *filename, const Model *model,
                           const Config *config, int dust) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror(filename);
        return 0;
    }
    fprintf(file, "r,sigma_initial,sigma_middle,sigma_final\n");
    const double *initial =
        dust ? model->sigma_d_initial : model->sigma_g_initial;
    const double *middle =
        dust ? model->sigma_d_middle : model->sigma_g_middle;
    const double *final = dust ? model->sigma_d : model->sigma_g;
    for (int i = 0; i < config->n_grid; ++i) {
        fprintf(file, "%.12g,%.12g,%.12g,%.12g\n",
                model->r[i], initial[i], middle[i], final[i]);
    }
    if (fclose(file) != 0) {
        perror(filename);
        return 0;
    }
    return 1;
}

/* Write final gas values and the derived quantities used in the plots. */
static int write_gas_output(const Model *model, const Config *config) {
    FILE *file = fopen("gas.csv", "w");
    if (file == NULL) {
        perror("gas.csv");
        return 0;
    }
    fprintf(file, "r,sigma_g,pressure,dpdr,v_g,alpha\n");
    for (int i = 0; i < config->n_grid; ++i) {
        fprintf(file, "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                model->r[i], model->sigma_g[i], model->pressure[i],
                model->dpdr[i], model->gas_velocity[i], model->alpha[i]);
    }
    if (fclose(file) != 0) {
        perror("gas.csv");
        return 0;
    }
    return 1;
}

/* Write the initial/final dust density and its final radial drift speed. */
static int write_dust_output(const Model *model, const Config *config) {
    FILE *file = fopen("dust.csv", "w");
    if (file == NULL) {
        perror("dust.csv");
        return 0;
    }
    fprintf(file, "r,sigma_d_initial,sigma_d_final,v_d\n");
    for (int i = 0; i < config->n_grid; ++i) {
        fprintf(file, "%.12g,%.12g,%.12g,%.12g\n",
                model->r[i], model->sigma_d_initial[i], model->sigma_d[i],
                model->dust_velocity[i]);
    }
    if (fclose(file) != 0) {
        perror("dust.csv");
        return 0;
    }
    return 1;
}

/* Record inputs and the stability-limited steps for interpretation in Python. */
static int write_parameters(const Config *config, double gas_dt,
                            double minimum_dt) {
    FILE *file = fopen("parameters.csv", "w");
    if (file == NULL) {
        perror("parameters.csv");
        return 0;
    }
    fprintf(file, "name,value\n");
    fprintf(file, "r_inner,%.12g\nr_outer,%.12g\nn_grid,%d\n",
            config->r_inner, config->r_outer, config->n_grid);
    fprintf(file, "reference_radius_au,%.12g\nstellar_mass_msun,%.12g\n",
            config->reference_radius_au, config->stellar_mass_msun);
    fprintf(file, "omega0_per_year,%.12g\n", config->omega0_per_year);
    fprintf(file, "h,%.12g\nalpha_active,%.12g\n",
            config->h, config->alpha_active);
    fprintf(file, "dead_inner,%.12g\ndead_outer,%.12g\n",
            config->dead_inner, config->dead_outer);
    fprintf(file, "dead_reduction,%.12g\nedge_width,%.12g\n",
            config->dead_reduction, config->edge_width);
    fprintf(file, "t_end,%.12g\nstokes_number,%.12g\n",
            config->t_end, config->stokes_number);
    fprintf(file, "initial_dust_to_gas,%.12g\n",
            config->initial_dust_to_gas);
    fprintf(file, "gas_stable_dt,%.12g\nminimum_used_dt,%.12g\n",
            gas_dt, minimum_dt);
    if (fclose(file) != 0) {
        perror("parameters.csv");
        return 0;
    }
    return 1;
}

/*
 * Evolve gas and dust together to the requested final time.
 * At each step: calculate the current velocities, choose a stable shared dt,
 * update both densities, and save the halfway profiles when we reach them.
 */
static int evolve_two_fluids(Model *model, const Config *config,
                             double *gas_dt, double *minimum_dt) {
    const int n = config->n_grid;
    const double middle_time = 0.5 * config->t_end;
    double time = 0.0;
    int middle_saved = 0;

    /* The gas viscosity is fixed, so its diffusion limit is calculated once. */
    *gas_dt = gas_stable_timestep(model, config);
    *minimum_dt = INFINITY;

    while (time < config->t_end) {
        /* Velocities depend on the gas density at the start of this step. */
        calculate_gas_fields(model, config);
        calculate_dust_velocity(model, config);

        /* The dust limit can change as the gas and dust velocities change. */
        const double dust_dt = dust_stable_timestep(model, config);
        double dt = *gas_dt;
        if (dust_dt < dt) {
            dt = dust_dt;
        }

        /* Do not step past the requested end time. */
        const double time_left = config->t_end - time;
        if (time_left < dt) {
            dt = time_left;
        }

        /* Shorten this step if needed, so the middle snapshot is exact. */
        if (!middle_saved && time < middle_time && time + dt > middle_time) {
            dt = middle_time - time;
        }
        if (!isfinite(dt) || dt <= 0.0) {
            fprintf(stderr, "Could not choose a positive common timestep.\n");
            return 0;
        }
        if (dt < *minimum_dt) {
            *minimum_dt = dt;
        }

        /*
         * Both updates use the same old-time state and the same dt.
         * This keeps gas and dust synchronized at every step.
         */
        if (!advance_gas(model, config, dt) ||
            !advance_dust(model, config, dt)) {
            return 0;
        }
        time += dt;

        if (!middle_saved && time >= middle_time) {
            /* Copy the synchronized halfway profiles for the snapshot files. */
            for (int i = 0; i < n; ++i) {
                model->sigma_g_middle[i] = model->sigma_g[i];
                model->sigma_d_middle[i] = model->sigma_d[i];
            }
            middle_saved = 1;
        }
    }

    /* Refresh diagnostics from the final, synchronized gas and dust state. */
    calculate_gas_fields(model, config);
    calculate_dust_velocity(model, config);
    return middle_saved;
}

int main(int argc, char **argv) {
    /* The program expects exactly one command-line argument: the YAML file. */
    if (argc != 2) {
        fprintf(stderr, "Usage: %s CONFIG.yaml\n", argv[0]);
        return EXIT_FAILURE;
    }

    Config config;
    if (!read_config(argv[1], &config)) {
        return EXIT_FAILURE;
    }

    /* Allocate and initialize all gas/dust arrays before starting the solver. */
    Model model = {0};
    if (!allocate_model(&model, config.n_grid)) {
        fprintf(stderr, "Could not allocate model arrays.\n");
        free_model(&model);
        return EXIT_FAILURE;
    }

    initialize_model(&model, &config);
    double gas_dt;
    double minimum_dt;

    /* Stop at the first failed step so we do not write incomplete output files. */
    int ok = evolve_two_fluids(&model, &config, &gas_dt, &minimum_dt);
    if (ok) {
        ok = write_gas_output(&model, &config);
    }
    if (ok) {
        ok = write_dust_output(&model, &config);
    }
    if (ok) {
        ok = write_snapshots("gas_snapshots.csv", &model, &config, 0);
    }
    if (ok) {
        ok = write_snapshots("dust_snapshots.csv", &model, &config, 1);
    }
    if (ok) {
        ok = write_parameters(&config, gas_dt, minimum_dt);
    }

    if (ok) {
        printf("Gas and dust evolved together to t = %.6g years.\n",
               config.t_end);
        printf("St = %.6g; initial dust-to-gas ratio = %.6g.\n",
               config.stokes_number, config.initial_dust_to_gas);
        printf("Gas stability limit = %.4g; smallest common timestep = %.4g.\n",
               gas_dt, minimum_dt);
    }

    free_model(&model);
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
