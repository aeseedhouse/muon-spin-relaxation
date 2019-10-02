#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

// keeps variables in one place so no need to redefine throughout the code
// config - holds variables put in from an external file
typedef struct {

	size_t num_dim;
	double mu_const;

	bool prnt_posit, use_lods;

	size_t num_lods; double lod_reduce; double lod_rad;

	size_t vec_size_db, vec_size_zu;

	size_t prg_chck_int, mag_num; bool prg_chck;

	size_t bar_len; clock_t str_t, end_t;

	double* mag_mom; double* mag_dim; size_t* diple_pnts; size_t num_diples;
	double* m_grid_max; double* m_grid_min; size_t* mag_pnts; size_t num_mags;
	size_t total_diples;

	double* c_grid_max; double* c_grid_min; size_t* c_grid_pnts; size_t c_grid_num_pnts;

	double *mean, *stdv;

}config_s;


// grid - holds variables for any grid structure e.g muon grid or magnet grid, also "mag_field" stores the mag field
typedef struct {

	double* min; double* max; double* int_arr;
	size_t* pnts;

	size_t axes_len; double* axes;
	size_t num_pnts; double* posit;

	double* mag_field;
	size_t* dim_com_1; size_t* dim_com_2; size_t com_idx;

}grid_s;

// mag - holds variables for the position of the magnets and how many dipoles in each
typedef struct {

	size_t* pnts; size_t num_diples;

	double* centre; double* dim; double* mag_mom;

	double* posit; double* diple_moms;

}mag_s;

// mag_grid - holds variables for the number of magnets in the sample film
typedef struct {

	size_t num_mags; mag_s* mags;
	grid_s m_grid;

	double* centre; double* dim; double* grid_mom;
	double* posit; double* mag_moms;

	size_t num_lods; double lod_reduce;

}mag_grid_s;


// reading file info
char* split() { return strtok(NULL, ","); } // reads string until comma
void get_argf(double *val) { *val = atof(split()); }
void get_args(char **val) { *val = split(); }
void get_argzu(size_t *val) { *val = atoi(split()); }
void get_argb(bool* val) { *val = atoi(split()); }

// skips a line in a file
void lineSkip(FILE* f, char* buffer, size_t size) {
	getline(&buffer, &size, f);
	char* t = strtok(buffer, ",");
}

// Reads data from config file into config
void readConfig(char* filename, config_s* config)
{
	printf("Loading config file... \n");
	FILE* f = fopen(filename, "r");

	if (f == NULL) { printf("Could not find config file: \"%s\"! Using defaults.\n", filename); return; }
	else {
		char* buffer = NULL;
		size_t lines_read = 0, i = 0, size;

		lineSkip(f, buffer, size); get_argzu(&config->num_dim);
		lineSkip(f, buffer, size); get_argf(&config->mu_const);
		getline(&buffer, &size, f);
		lineSkip(f, buffer, size); get_argb(&config->prnt_posit);
		lineSkip(f, buffer, size); get_argzu(&config->bar_len);
		lineSkip(f, buffer, size); get_argzu(&config->prg_chck_int);
		getline(&buffer, &size, f);
		lineSkip(f, buffer, size); get_argf(&config->lod_reduce);
		lineSkip(f, buffer, size); get_argzu(&config->num_lods);
		lineSkip(f, buffer, size); get_argf(&config->lod_rad);
		getline(&buffer, &size, f);
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argf(&config->mag_mom[dim_idx]); }
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argf(&config->mag_dim[dim_idx]); }
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argzu(&config->diple_pnts[dim_idx]); }
		getline(&buffer, &size, f);
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argf(&config->m_grid_max[dim_idx]); }
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argf(&config->m_grid_min[dim_idx]); }
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argzu(&config->mag_pnts[dim_idx]); }
		getline(&buffer, &size, f);
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argf(&config->c_grid_max[dim_idx]); }
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argf(&config->c_grid_min[dim_idx]); }
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argzu(&config->c_grid_pnts[dim_idx]); }
		getline(&buffer, &size, f);
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argf(&config->mean[dim_idx]); }
		lineSkip(f, buffer, size); for (size_t dim_idx = 0; dim_idx < config->num_dim; dim_idx++) { get_argf(&config->stdv[dim_idx]); }

		fflush(stdout);

		printf("Config loaded.\n");

	}
	fclose(f);
}


// progress of code - michael did this
void printProgress(size_t count, double time_to_cmplt, size_t max, config_s config)
{
	double percent = ((double)count / (double)max);
	count = floor(percent*config.bar_len);
	max = config.bar_len;

	char prefix[] = "Progress: [";
	char suffix[] = "]";
	size_t prefix_length = sizeof(prefix) - 1;
	size_t suffix_length = sizeof(suffix) - 1;
	char *buffer = calloc(max + prefix_length + suffix_length + 1, 1); // +1 for \0
	size_t i = 0;

	strcpy(buffer, prefix);
	for (; i < max; ++i) {
		buffer[prefix_length + i] = i < count ? '#' : ' ';
	}

	strcpy(&buffer[prefix_length + i], suffix);
	printf("\b%c[2K\r%s", 5, buffer);
	printf(" %.2f%%. ", percent * 100);

	printf("Mag: %zu/%zu. ", config.mag_num, config.num_mags);

	int n, day, hr, min, sec;
	n = (size_t)time_to_cmplt;

	if (n>86400) {
		;min = n / 60; sec = n % 60; hr = min / 60;  min = min % 60; day = hr / 24; hr = day % 24;
		printf("ETC: %d:%d:%d:%d days \n", day, hr, min, sec);
	}
	else if (n>3600) {
		min = n / 60; sec = n % 60; hr = min / 60; min = min % 60;
		printf("ETC: %d:%d:%d hours \n", hr, min, sec);
	}
	else if (n>60) {
		min = n / 60; sec = n % 60;
		printf("ETC: %d:%d mins \n", min, sec);
	}
	else {
		printf("ETC: %d secs \n", n);
	}

	fflush(stdout);
	free(buffer);

}

size_t calcArrProduct(size_t* arr, size_t arr_len) {
	size_t arr_prod = 1;
	for (size_t pnt_idx = 0; pnt_idx < arr_len; pnt_idx++) {
		arr_prod *= arr[pnt_idx];
	}
	return arr_prod;
}

size_t calcMax(size_t* arr, size_t arr_len) {
	size_t max = 0;
	for (size_t pnt_idx = 0; pnt_idx < arr_len; pnt_idx++) {
		if (arr[pnt_idx] > max) { max = arr[pnt_idx]; }
	}
	return max;
}

double* calcVectSum(double* vect_1, double* vect_2, size_t num_dim) {
	double* sum = malloc(sizeof(double)*num_dim);
	for (size_t dim_idx = 0; dim_idx < num_dim; dim_idx++) {
		sum[dim_idx] = vect_1[dim_idx] + vect_2[dim_idx];
	}
	return sum;

}

double* calcVectDiff(double* vect_1, double* vect_2, size_t num_dim) {
	double* diff = malloc(sizeof(double)*num_dim);
	for (size_t dim_idx = 0; dim_idx < num_dim; dim_idx++) {
		diff[dim_idx] = vect_1[dim_idx] - vect_2[dim_idx];
	}
	return diff;
}

double calcVectDot(double* vect_1, double* vect_2, size_t num_dim) {
	double dot = 0;
	for (size_t dim_idx = 0; dim_idx < num_dim; dim_idx++) {
		dot += vect_1[dim_idx] * vect_2[dim_idx];
	}
	return dot;
}

double calcVectNorm(double* vect, size_t num_dim) {
	double norm = 0;
	for (size_t dim_idx = 0; dim_idx < num_dim; dim_idx++) {
		norm += vect[dim_idx] * vect[dim_idx];
	}
	norm = sqrt(norm);
	return norm;
}

double randn()
{
	double U1, U2, W, mult;
	static double X1, X2;
	static int call = 0;

	if (call == 1)
	{
		call = !call;
		return ((double)X2);
	}

	do
	{
		U1 = -1 + ((double)rand() / RAND_MAX) * 2;
		U2 = -1 + ((double)rand() / RAND_MAX) * 2;
		W = pow(U1, 2) + pow(U2, 2);
	} while (W >= 1 || W == 0);

	mult = sqrt((-2 * log(W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return  ((double)X1);
}

double* randomGauss(config_s config) {
	double u, v;
	double z0, z1, z2;

	u = (double)rand() / (double)RAND_MAX;
	v = (double)rand() / (double)RAND_MAX;
	z2 = randn();
	z0 = (-2.*log(u)) * cos(2.*3.14159*v);
	z1 = (-2.*log(u)) * sin(2.*3.14159*v);

	double* z_comp = malloc(sizeof(double) * 3);
	z_comp[0] = z0 * config.stdv[0] + config.mean[0];
	z_comp[1] = z1 * config.stdv[1] + config.mean[1];
	z_comp[2] = z2 * config.stdv[2] + config.mean[2];

	return z_comp;
}

double* randomArr(config_s config) {
	double a, b, c, norm;
	double* arr = malloc(sizeof(double) * 3);

	a = (double)rand() / (double)RAND_MAX;
	b = (double)rand() / (double)RAND_MAX;
	c = (double)rand() / (double)RAND_MAX;

	arr[0] = a; arr[1] = b; arr[2] = c;

	norm = calcVectNorm(arr, config.num_dim);

	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
		arr[dim_idx] = norm*arr[dim_idx];
	}

	return arr;
}


// creates the axes array for each dimension
void makeAxes(grid_s grid, config_s config, grid_s* ret_grid) {
	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
		if (grid.pnts[dim_idx] > 1) {
			if (grid.max[dim_idx] < grid.min[dim_idx]) { printf("Grid maximum must be greater than grid minimum! Exiting.\n"); exit(1); }
			else if (grid.max[dim_idx] == grid.min[dim_idx]) { printf("Grid maximum cannot be equal to the grid minimum for non singular grids! Exiting.\n"); exit(1); }
			else { grid.int_arr[dim_idx] = (grid.max[dim_idx] - grid.min[dim_idx]) / (double)(grid.pnts[dim_idx] - 1); }
		}
		else if (grid.pnts[dim_idx] == 1) {
			if (grid.max[dim_idx] < grid.min[dim_idx]) { printf("Grid maximum must be greater than grid minimum! Exiting.\n"); exit(1); }
			else if (grid.max[dim_idx] == grid.min[dim_idx]) { grid.int_arr[dim_idx] = 0; }
			else { printf("Differing minumum and maximum detected for singular array, assuming midpoint.\n"); grid.min[dim_idx] += (grid.max[dim_idx] - grid.min[dim_idx]) / 2.0; }
		}
		else { printf("Points grid cannot have 0 value! Exiting.\n"); exit(1); }
	}
	for (size_t pnt_idx = 0; pnt_idx < calcMax(grid.pnts, config.num_dim); pnt_idx++) {
		for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
			grid.axes[config.num_dim*pnt_idx + dim_idx] = grid.min[dim_idx] + pnt_idx*grid.int_arr[dim_idx];
		}
	}

	*ret_grid = grid;
}

// 
void makeGridLayer(grid_s* grid, size_t layer, config_s config) {

	for (size_t idx_arr = 0; idx_arr < grid->pnts[layer]; idx_arr++) {
		grid->dim_com_1[layer] = idx_arr*config.num_dim;

		grid->dim_com_2[layer] = 1;
		for (size_t idx = layer + 1; idx < config.num_dim; idx++) { grid->dim_com_2[layer] *= grid->pnts[idx]; }
		grid->dim_com_2[layer] *= grid->dim_com_1[layer];
		grid->com_idx += grid->dim_com_2[layer];

		if (layer >= config.num_dim - 1) {
			for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) { grid->posit[grid->com_idx + dim_idx] = grid->axes[grid->dim_com_1[dim_idx] + dim_idx]; }
			grid->com_idx -= grid->dim_com_2[layer];
			continue;
		}

		makeGridLayer(grid, layer + 1, config);
		grid->com_idx -= grid->dim_com_2[layer];
	}
}

// same as grid layer but for gaussian distribution
void makeGaussGridLayer(grid_s* grid, size_t layer, config_s config) {

	double* rand_num = malloc(sizeof(double) * 3);
	for (size_t idx_arr = 0; idx_arr < grid->pnts[layer]; idx_arr++) {
		grid->dim_com_1[layer] = idx_arr*config.num_dim;

		grid->dim_com_2[layer] = 1;
		for (size_t idx = layer + 1; idx < config.num_dim; idx++) { grid->dim_com_2[layer] *= grid->pnts[idx]; }
		grid->dim_com_2[layer] *= grid->dim_com_1[layer];
		grid->com_idx += grid->dim_com_2[layer];

		if (layer >= config.num_dim - 1) {
			rand_num = randomGauss(config);
			for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) { grid->posit[grid->com_idx + dim_idx] = rand_num[dim_idx]; }
			grid->com_idx -= grid->dim_com_2[layer];
			continue;
		}

		makeGaussGridLayer(grid, layer + 1, config);
		grid->com_idx -= grid->dim_com_2[layer];
	}
}

// creates a grid using the layers from makeGridLayer
void makeGrid(grid_s grid, config_s config, grid_s* ret_grid) {

	grid.dim_com_1 = malloc(config.vec_size_zu); grid.dim_com_2 = malloc(config.vec_size_zu);
	grid.com_idx = 0;

	makeGridLayer(&grid, 0, config);

	free(grid.dim_com_1);
	free(grid.dim_com_2);

	*ret_grid = grid;
}

// allocates space for variables
config_s initConfig() {

	config_s config;
	config.str_t = clock();

	//Calculates vector size:
	config.vec_size_db = config.num_dim * sizeof(double);
	config.vec_size_zu = config.num_dim * sizeof(size_t);

	config.mag_mom = malloc(config.vec_size_db); config.mag_dim = malloc(config.vec_size_db); config.diple_pnts = malloc(config.vec_size_db);
	config.m_grid_max = malloc(config.vec_size_db); config.m_grid_min = malloc(config.vec_size_db); config.mag_pnts = malloc(config.vec_size_db);
	config.c_grid_max = malloc(config.vec_size_db); config.c_grid_min = malloc(config.vec_size_db); config.c_grid_pnts = malloc(config.vec_size_db);
	config.mean = malloc(config.vec_size_db); config.stdv = malloc(config.vec_size_db);

	return config;
}

// same as make grid but uses makeGaussGridLayer instead
void makeGaussGrid(grid_s grid, config_s config, grid_s* ret_grid) {

	grid.dim_com_1 = malloc(config.vec_size_zu); grid.dim_com_2 = malloc(config.vec_size_zu);
	grid.com_idx = 0;

	makeGaussGridLayer(&grid, 0, config);

	free(grid.dim_com_1);
	free(grid.dim_com_2);

	*ret_grid = grid;
}

// initialise the grid, i.e. allocate all space
grid_s initGrid(double* max, double* min, size_t* pnts, config_s config) {

	grid_s grid;
	grid.max = malloc(config.vec_size_zu), grid.min = malloc(config.vec_size_zu), grid.pnts = malloc(config.vec_size_zu);

	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
		grid.max[dim_idx] = max[dim_idx];
		grid.min[dim_idx] = min[dim_idx];
		grid.pnts[dim_idx] = pnts[dim_idx];
	}

	grid.num_pnts = calcArrProduct(grid.pnts, config.num_dim);
	grid.mag_field = calloc(grid.num_pnts*config.num_dim, sizeof(double));
	grid.posit = malloc(sizeof(double)*grid.num_pnts*config.num_dim);
	grid.axes_len = calcMax(grid.pnts, config.num_dim)*config.num_dim;
	grid.axes = calloc(grid.axes_len, sizeof(double));
	grid.int_arr = malloc(config.vec_size_db);

	makeAxes(grid, config, &grid);
	makeGrid(grid, config, &grid);

	return grid;
}

// same as initialise grid but for gauss, again.
grid_s initGaussGrid(double* max, double* min, size_t* pnts, config_s config) {

	grid_s grid;
	grid.max = malloc(config.vec_size_zu), grid.min = malloc(config.vec_size_zu), grid.pnts = malloc(config.vec_size_zu);

	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
		grid.max[dim_idx] = max[dim_idx];
		grid.min[dim_idx] = min[dim_idx];
		grid.pnts[dim_idx] = pnts[dim_idx];
	}

	grid.num_pnts = calcArrProduct(grid.pnts, config.num_dim);
	grid.mag_field = calloc(grid.num_pnts*config.num_dim, sizeof(double));
	grid.posit = malloc(sizeof(double)*grid.num_pnts*config.num_dim);
	grid.axes_len = calcMax(grid.pnts, config.num_dim)*config.num_dim;
	grid.axes = calloc(grid.axes_len, sizeof(double));
	grid.int_arr = malloc(config.vec_size_db);

	makeAxes(grid, config, &grid);
	makeGaussGrid(grid, config, &grid);

	return grid;
}

// deallocates memory
void freeGrid(grid_s* grid) {
	free(grid->max); free(grid->min); free(grid->pnts);
	free(grid->mag_field); free(grid->posit); free(grid->axes); free(grid->int_arr);
}

// initialise the magnetic grid, saving the variables to mag structure
mag_s initMag(double* centre, double* dim, size_t* pnts, double* mag_mom, config_s config) {
	mag_s mag;

	mag.centre = malloc(config.vec_size_db), mag.dim = malloc(config.vec_size_zu); mag.pnts = malloc(config.vec_size_db); mag.mag_mom = malloc(config.vec_size_db);
	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
		mag.centre[dim_idx] = centre[dim_idx];
		mag.dim[dim_idx] = dim[dim_idx];
		mag.pnts[dim_idx] = pnts[dim_idx];
		mag.mag_mom[dim_idx] = mag_mom[dim_idx];
	}

	mag.num_diples = 1;
	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
		mag.num_diples *= mag.pnts[dim_idx];
	}
	mag.posit = malloc(config.vec_size_db*mag.num_diples), mag.diple_moms = malloc(config.vec_size_db*mag.num_diples);

	double* d_grid_max = malloc(config.vec_size_zu); double* d_grid_min = malloc(config.vec_size_zu); size_t* d_grid_pnts = malloc(config.vec_size_zu);
	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
		d_grid_max[dim_idx] = mag.dim[dim_idx] / 2. + mag.centre[dim_idx];
		d_grid_min[dim_idx] = -mag.dim[dim_idx] / 2. + mag.centre[dim_idx];
		d_grid_pnts[dim_idx] = mag.pnts[dim_idx];
	}

	grid_s d_grid = initGrid(d_grid_max, d_grid_min, d_grid_pnts, config);

	free(d_grid_max); free(d_grid_min); free(d_grid_pnts);

	double* rand = malloc(sizeof(double) * 3);

	for (size_t diple_idx = 0; diple_idx < mag.num_diples; diple_idx++) {
		rand = randomArr(config);
		for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
			mag.posit[config.num_dim*diple_idx + dim_idx] = d_grid.posit[config.num_dim*diple_idx + dim_idx];
			mag.diple_moms[config.num_dim*diple_idx + dim_idx] = rand[dim_idx] * mag.mag_mom[0]; //randomise the magnetic moment!!
		}
	}

	freeGrid(&d_grid);

	return mag;
}


// equation used to calculate the magnetic field from  a dipole
void eqBField(double* int_b_field, size_t num_dim, double mu_const, double r_dot_m, double norm_r, double* diff, double* mag_mo, double* b_field) {
	double const_1 = mu_const*(3.0 * r_dot_m / pow(norm_r, 5));
	double const_2 = mu_const*(1.0 / pow(norm_r, 3));
	for (size_t dim_idx = 0; dim_idx < num_dim; dim_idx++) {
		b_field[dim_idx] = int_b_field[dim_idx] + const_1*diff[dim_idx] - const_2*mag_mo[dim_idx];
	}
}

// calculating the magnetic field
void calcBField(double* int_b_field, config_s config, double* mag_posit, double* field_posit, double* mag_mo, double* ret_b_field) {
	double* diff = malloc(config.vec_size_db);
	diff = calcVectDiff(field_posit, mag_posit, config.num_dim);
	double r_dot_m = calcVectDot(diff, mag_mo, config.num_dim);
	double norm_r = calcVectNorm(diff, config.num_dim);
	double* b_field = malloc(sizeof(double)*config.num_dim);
	eqBField(int_b_field, config.num_dim, config.mu_const, r_dot_m, norm_r, diff, mag_mo, b_field);
	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) { ret_b_field[dim_idx] = b_field[dim_idx]; }
	free(diff);
	free(b_field);
}


// calculating the magnetic field at each point in the grid defined, meatyyyyyy
void calcGridField(grid_s c_grid, mag_s* mags, config_s config, grid_s* ret_c_grid)
{
	config.prg_chck = 0;
	clock_t srt_t, end_t;
	double loop_t = 0, total_t = 0, avg_t = 0;
	size_t chck_idx = 0;
	double time_to_cmplt = 0;
	size_t lod_idx;

	//printf("Centre: %f,%f,%f \n", mags[0].centre[0], mags[0].centre[1], mags[0].centre[2]);


	for (size_t pnt_idx = 0; pnt_idx < c_grid.num_pnts; pnt_idx++) {

		if (pnt_idx % config.prg_chck_int == 0) { config.prg_chck = 1; }

		size_t p_idx = config.num_dim*pnt_idx;

		lod_idx = 0;
		if (config.use_lods == 1) {

			double* pnt_diff = calcVectDiff(&c_grid.posit[pnt_idx], mags[0].centre, config.num_dim);
			double diff_mag = calcVectNorm(pnt_diff, config.num_dim);
			free(pnt_diff);

			lod_idx = floor(diff_mag / config.lod_rad);
			if (lod_idx >= config.num_lods) { lod_idx = (config.num_lods - 1); }
		}

		for (size_t mag_idx = 0; mag_idx < mags[lod_idx].num_diples; mag_idx++) {
			size_t m_idx = config.num_dim*mag_idx;
			calcBField(&c_grid.mag_field[p_idx], config, &mags[lod_idx].posit[m_idx], &c_grid.posit[p_idx], &mags[lod_idx].diple_moms[m_idx], &c_grid.mag_field[p_idx]);
		}

		if (config.prg_chck) {
			if (chck_idx != 0) {
				end_t = clock();
				loop_t = (double)(end_t - srt_t) / CLOCKS_PER_SEC;
				total_t += loop_t;
				avg_t = total_t / (config.prg_chck_int*chck_idx);
				time_to_cmplt = avg_t*((double)c_grid.num_pnts*(config.num_mags - config.mag_num + 1) - (double)pnt_idx);
				//printProgress(pnt_idx, time_to_cmplt, c_grid.num_pnts, config);
			}
			config.prg_chck = 0;
			srt_t = clock();

			if (config.prnt_posit) {
				printf("Poisition: %f, %f, %f. \n", c_grid.posit[p_idx + 0], c_grid.posit[p_idx + 1], c_grid.posit[p_idx + 2]);
				printf("B Field: %.10e, %.10e, %.10e.\n", c_grid.mag_field[p_idx + 0], c_grid.mag_field[p_idx + 1], c_grid.mag_field[p_idx + 2]);
			}

			chck_idx++;
		}

	}

	*ret_c_grid = c_grid;
}

mag_s randMoms(mag_s mag, config_s config, mag_s* ret_mag) {
	double* rand = malloc(config.vec_size_db);
	for (size_t diple_idx = 0; diple_idx < mag.num_diples; diple_idx++) {
		rand = randomArr(config);
		for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
			mag.diple_moms[config.num_dim*diple_idx + dim_idx] = rand[dim_idx] * mag.mag_mom[0]; //randomise the magnetic moment!!
		}
	}

	*ret_mag = mag;
}


// initialise the magnetic grid - when using more than one dipole for the magnet then chooses how to approximate
// the magnets dependiong on how far they are from the position of the muon
mag_grid_s initMagGrid(double* max, double* min, size_t* pnts, double* mag_dim, size_t* diple_pnts, double* mag_mom, size_t num_lods, double lod_reduce, config_s config) {
	mag_grid_s mag_grid;
	mag_grid.m_grid = initGrid(max, min, pnts, config);

	mag_grid.lod_reduce = lod_reduce;
	mag_grid.num_lods = num_lods;

	mag_grid.num_mags = mag_grid.m_grid.num_pnts;
	mag_grid.mags = malloc(sizeof(mag_s)*config.num_lods);

	for (size_t lod_idx = 0; lod_idx < mag_grid.num_lods; lod_idx++) {
		for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
			size_t new_pnts;
			if (lod_idx != 0) { new_pnts = floor(diple_pnts[dim_idx] * mag_grid.lod_reduce); }
			else { new_pnts = diple_pnts[dim_idx]; }
			if (new_pnts > 0) { diple_pnts[dim_idx] = new_pnts; }
			else { printf("Lod devision too fine, in dimension %zu Try increasing magnet pnts, or reducing reduce or num_lods.", dim_idx); }
		}
		mag_grid.mags[lod_idx] = initMag(&mag_grid.m_grid.posit[0], mag_dim, diple_pnts, mag_mom, config);

	}

	return mag_grid;
}

// saving data in files
void printDoubles(char* output_filename, size_t srt_line, size_t end_line, size_t prec, double* data, double* data_2, size_t num_cols, size_t print_number)
{
	size_t num_lines = (end_line - srt_line);

	if (output_filename != NULL) {
		FILE* f = fopen(output_filename, "w");
		for (size_t line_idx = srt_line; line_idx < (srt_line + num_lines); ++line_idx) {
			for (size_t col_idx = 0; col_idx < num_cols; ++col_idx) {
				fprintf(f, "%.*e,", prec, data[line_idx*num_cols + col_idx]);

			}
			for (size_t col_idx = 0; col_idx < num_cols; ++col_idx) {
				fprintf(f, "%.*e", prec, data_2[line_idx*num_cols + col_idx]);
				if (col_idx != num_cols - 1) { fprintf(f, ","); }
			}

			fprintf(f, "%s\n");
		}
		fclose(f);

	}
	else { return; }
}

void printAr(char* output_filename, double* data, size_t ar_len)
{
	FILE* f = fopen(output_filename, "w");
	for (size_t ar_idx = 0; ar_idx < ar_len-1; ++ar_idx) {
		fprintf(f, "%.*e", 7, data[ar_idx]);
		fprintf(f, ",");
	}

	fclose(f);

}

// saves grid points and magnetic field with it
void printGridField(char* output_filename, grid_s grid, config_s config, size_t prec, size_t print_number)
{
	if (output_filename != NULL) {

		char* filename_mag_field = malloc(sizeof(output_filename) + sizeof(char) * 12 + sizeof(size_t));
		sprintf(filename_mag_field, "%s_mag_field_%zu.csv", output_filename, print_number);

		printf("Printing magnetic field to file: %s\n", filename_mag_field);
		printDoubles(filename_mag_field, 0, grid.num_pnts, prec, grid.posit, grid.mag_field, config.num_dim, print_number);
		printf("Printing completed. \n");
	}
	else printf("No output selected. Canceling printing. \n");
}

// checks if reading the config file, then you choose to run the code or not
void chckContinue() {
	printf("Do you wish to continue? Yes (y). No (n).\n");
	char enter = 0;
	while (enter != 'y') {
		enter = getchar();
		if (enter == 'n') { printf("Exiting.\n"); exit(1); }
	}
	printf("Continuing.\n");
}

//Checks Dimension against default to ward off errors.
void chckDim(config_s config, const size_t num_dim_default)
{
	if (config.num_dim != num_dim_default) {
		printf("Warning inputted number of dimesions does not match default! \n");
		printf("This will lead to unexpected errors if other inputs are not ajusted acordingly. \n");
		printf("Number of dimensions is set to: %zu, Default is %zu.", config.num_dim, num_dim_default);
		chckContinue();
	}
}

// estimates the amount of memory needed to run/save 
void chckMem(config_s config) {
	printf("Number of co-ordinate grid points: ||%zu points.\n", config.c_grid_num_pnts);
	printf("Number of dipoles per magnet:      ||%zu dipoles.\n", config.num_diples);
	printf("Number of magnets:                 ||%zu magnets.\n", config.num_mags);
	printf("Total number of dipoles:           ||%zu dipoles.\n", config.total_diples);
	printf("Number of LODs:                    ||%zu lods.\n", config.num_lods);

	size_t min_mem = config.num_diples*config.vec_size_db * 2 + config.c_grid_num_pnts*config.vec_size_db * 2;

	if (min_mem < 1e3) { printf("Minumum Memory required:           ||%zu bytes.\n", min_mem); }
	else if (min_mem < 1e6) { printf("Minumum Memory required:           ||%.2f KB.\n", min_mem / 1.e3); }
	else if (min_mem < 1e9) { printf("Minumum Memory required:           ||%.2f MB.\n", min_mem / 1.e6); }
	else if (min_mem < 1e12) { printf("Minumum Memory required:           ||%.2f GB.\n", min_mem / 1.e9); }
	chckContinue();
}

// changing the centre of the magnet as the code moves onto the next one in the grid
void moveMag(mag_s* mag, double* old_centre, double* new_centre, config_s config) {
	double* centre_diff = malloc(config.vec_size_db);
	centre_diff = calcVectDiff(new_centre, old_centre, config.num_dim);

	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) { mag->centre[dim_idx] += centre_diff[dim_idx]; }

	for (size_t pnt_idx = 0; pnt_idx < mag->num_diples; pnt_idx++) {
		for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) { mag->posit[pnt_idx*config.num_dim + dim_idx] += centre_diff[dim_idx]; }
	}

	free(centre_diff);
}

int main()
{
	// So that random numbers can be generated
	srand(time(NULL));

	//Config Defaults (Program will revert to these if no config file is found, or if readConfig function is commented out):

	const size_t num_dim = 3;

	config_s config = initConfig();

	config.num_dim = 3;
	config.mu_const = 1.0E-7;

	config.prnt_posit = 0;
	config.bar_len = 45;
	config.prg_chck_int = 500;

	config.use_lods = 0;
	config.lod_reduce = 0.5; // how much reduced by each step
	config.num_lods = 2; // how many layers
	config.lod_rad = 0.01; // how far away needs to be before approximation comes in

	double mag_mom[3] = { 1.,0.,0. };
	double mag_dim[3] = { 0.01,0.05,0.002 };
	size_t diple_pnts[3] = { 4,4,1 };

	double m_grid_max[3] = { 0.5,0.5,0. };
	double m_grid_min[3] = { -0.5,-0.5,0. };
	size_t mag_pnts[3] = { 10, 10, 1 };

	double c_grid_max[3] = { 1,1,0 };
	double c_grid_min[3] = { -1, -1, 0 };
	size_t c_grid_pnts[3] = { 100, 100, 1 };

	config.vec_size_db = config.num_dim * sizeof(double);
	config.vec_size_zu = config.num_dim * sizeof(size_t);

	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
		config.mag_mom[dim_idx] = mag_mom[dim_idx]; config.mag_dim[dim_idx] = mag_dim[dim_idx]; config.diple_pnts[dim_idx] = diple_pnts[dim_idx];
		config.m_grid_max[dim_idx] = m_grid_max[dim_idx]; config.m_grid_min[dim_idx] = m_grid_min[dim_idx]; config.mag_pnts[dim_idx] = mag_pnts[dim_idx];
		config.c_grid_max[dim_idx] = c_grid_max[dim_idx]; config.c_grid_min[dim_idx] = c_grid_min[dim_idx]; config.c_grid_pnts[dim_idx] = c_grid_pnts[dim_idx];
	}


	printf("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* Program Start *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* \n");
	printf("Running Magnificient Magnetic Field Program, now with Gaussian muons! \n");

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initilising ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	printf("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~* Initilize Config *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* \n");

	readConfig("config.csv", &config); //<- Read config from file, overwrites defaults.

	config.num_diples = calcArrProduct(config.diple_pnts, config.num_dim);
	config.num_mags = calcArrProduct(config.mag_pnts, config.num_dim);
	config.c_grid_num_pnts = calcArrProduct(config.c_grid_pnts, config.num_dim);
	config.total_diples = config.num_mags*config.num_diples;
	config.mag_num = 1;

	printf("*~~~~~~~~~~~~~~~~~~~~~~~~~* Check Memory Required *~~~~~~~~~~~~~~~~~~~~~~~~~~~* \n");
	chckMem(config);
	//chckDim(config, num_dim);

	printf("*~~~~~~~~~~~~~~~~~~~~~~~~* Initillise Coordinate Grid *~~~~~~~~~~~~~~~~~~~~~~~* \n");

	//Initillise Coordinate Grid:
	grid_s c_grid = initGaussGrid(config.c_grid_max, config.c_grid_min, config.c_grid_pnts, config);

	printf("Complete.\n");

	printf("*~~~~~~~~~~~~~~~~~~~~~~~~~* Initillise Magnet Grid *~~~~~~~~~~~~~~~~~~~~~~~~~~* \n");

	grid_s m_grid = initGrid(config.m_grid_max, config.m_grid_min, config.mag_pnts, config);
	mag_s* mags = malloc(sizeof(mag_s)*config.num_lods);

	size_t* new_pnts = malloc(config.vec_size_zu);
	double* new_mom = malloc(config.vec_size_db);

	for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) { new_mom[dim_idx] = config.mag_mom[dim_idx]; }

	for (size_t lod_idx = 0; lod_idx < config.num_lods; lod_idx++) {
		for (size_t dim_idx = 0; dim_idx < config.num_dim; dim_idx++) {
			if (lod_idx != 0) {
				size_t prev_pnts = new_pnts[dim_idx];
				new_pnts[dim_idx] = floor(new_pnts[dim_idx] * config.lod_reduce);
				if (new_pnts[dim_idx] > 0 && prev_pnts > 0) {
					for (size_t dim_idx_2 = 0; dim_idx_2 < config.num_dim; dim_idx_2++) { new_mom[dim_idx_2] = new_mom[dim_idx_2] * ((double)prev_pnts / (double)new_pnts[dim_idx]); }

				}

			}
			else {
				new_pnts[dim_idx] = config.diple_pnts[dim_idx];
				for (size_t dim_idx_2 = 0; dim_idx_2 < config.num_dim; dim_idx_2++) { new_mom[dim_idx_2] = config.mag_mom[dim_idx_2]; }

			}
			if (new_pnts[dim_idx] > 0) {
				config.diple_pnts[dim_idx] = new_pnts[dim_idx];
				for (size_t dim_idx_2 = 0; dim_idx_2 < config.num_dim; dim_idx_2++) { config.mag_mom[dim_idx_2] = new_mom[dim_idx_2]; }

			}
			else { printf("Lod devision too fine, in dimension %zu Try increasing magnet pnts, or reducing reduce or num_lods. \n", dim_idx); }
		}
		mags[lod_idx] = initMag(&m_grid.posit[0], config.mag_dim, config.diple_pnts, config.mag_mom, config);
	}

	printf("Complete.\n");

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Map B Field ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	printf("*~~~~~~~~~~~~~~~* Mapping Magnetic Field to Co-ordinate Grid *~~~~~~~~~~~~~~~~* \n");

	double* angle = malloc(config.num_mags);

	for (size_t lod_idx = 0; lod_idx < config.num_lods; lod_idx++) {
		moveMag(&mags[lod_idx], mags[lod_idx].centre, &m_grid.posit[0], config);
	}
	calcGridField(c_grid, mags, config, &c_grid);
	//angle[0] = atan(mags[0].diple_moms[0] / mags[0].diple_moms[2]);

	config.mag_num += 1;
	
	srand(time(NULL));
	for (size_t mag_idx = 1; mag_idx < m_grid.num_pnts; mag_idx++) {
		//printf("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* \n");
		for (size_t lod_idx = 0; lod_idx < config.num_lods; lod_idx++) {
			moveMag(&mags[lod_idx], &m_grid.posit[(mag_idx - 1)*config.num_dim], &m_grid.posit[mag_idx*config.num_dim], config);
			randMoms(mags[lod_idx], config, &mags[lod_idx]);
		}
		//angle[mag_idx] = atan(mags[0].diple_moms[0] / mags[0].diple_moms[2]);
		calcGridField(c_grid, mags, config, &c_grid);

		config.mag_num += 1;
	}

	//printAr("angle.dat", angle, config.num_mags);
	printf("Complete.\n");

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Print B Field ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	printf("*~~~~~~~~~~~~~~~~~~~~~~~~* Printing Co-ordinate Grid *~~~~~~~~~~~~~~~~~~~~~~~~* \n");

	printGridField("c_grid", c_grid, config, 20, 1);

	freeGrid(&c_grid);

	printf("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* End Program *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* \n");

	return 0;

}