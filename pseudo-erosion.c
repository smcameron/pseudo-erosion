/*
	Copyright (C) 2017 Stephen M. Cameron
	Author: Stephen M. Cameron

	This file is part of pseudo-erosion.

	pseudo-erosion is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	pseudo-erosion is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with pseudo-erosion; if not, write to the Free Software
	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>

#include "open-simplex-noise.h"
#include "png_utils.h"

#define DEFAULT_IMAGE_SIZE 1024
#define DEFAULT_FEATURE_SIZE 512
#define DEFAULT_GRID_SIZE 4

static char *output_file = "output.png";
static int image_size = DEFAULT_IMAGE_SIZE;
static int feature_size = DEFAULT_FEATURE_SIZE;
static int grid_size = DEFAULT_GRID_SIZE;
static int seed = 123456;
static char *input_image = NULL;

static struct option long_options[] = {
	{ "featuresize", required_argument, NULL, 'f' },
	{ "gridsize", required_argument, NULL, 'g' },
	{ "size", required_argument, NULL, 's' },
	{ "seed", required_argument, NULL, 'S' },
	{ "outputfile", required_argument, NULL, 'o' },
	{ "input", required_argument, NULL, 'i' },
	{ 0, 0, 0, 0 },
};

static void usage(void)
{
	fprintf(stderr, "pseudo_erosion: Usage:\n\n");
	fprintf(stderr, "	pseudo_erosion [-g gridsize] [-o outputfile] [-s imagesize] \\\n");
	fprintf(stderr, "		[-i inputfile] [-f featuresize]\n");
	fprintf(stderr, "\n");
	exit(1);
}

struct grid_point {
	double x, y;
	int cx, cy; /* connected to gridpoint(grid, cx, cy) */
};

struct grid {
	struct grid_point *g;
	int dim;
};

static struct grid *allocate_grid(int dim)
{
	struct grid_point *gp;
	struct grid *g;

	gp = malloc(sizeof(*gp) * (dim + 1) * (dim + 1));
	memset(gp, 0,  sizeof(*gp) * (dim + 1) * (dim + 1));
	g = malloc(sizeof(*g));
	g->g = gp;
	g->dim = dim;
	return g;
}

static void free_grid(struct grid *grid)
{
	free(grid->g);
	grid->g = NULL;
	free(grid);
}

static inline struct grid_point *gridpoint(struct grid *grid, int x, int y)
{
	return &grid->g[(grid->dim + 1) * y + x];
}

static double *allocate_heightmap(int dim)
{
	double *hm;

	hm = calloc((size_t) dim * (size_t) dim, sizeof(*hm));
	if (!hm) {
		fprintf(stderr, "pseudo-erosion: Out of memory allocating %dx%d heightmap\n", dim, dim);
		exit(1);
	}
	return hm;
}

static int write_heightmap_png(const char *filename, const double *hm, int dim)
{
	size_t i, n = (size_t) dim * (size_t) dim;
	double minv, maxv, range;
	uint32_t *pixels;
	int rc;

	if (n == 0)
		return -1;

	minv = hm[0];
	maxv = hm[0];
	for (i = 1; i < n; i++) {
		if (hm[i] < minv)
			minv = hm[i];
		if (hm[i] > maxv)
			maxv = hm[i];
	}
	range = maxv - minv;

	pixels = malloc(n * sizeof(*pixels));
	if (!pixels) {
		fprintf(stderr, "pseudo-erosion: Out of memory allocating PNG buffer\n");
		exit(1);
	}

	for (i = 0; i < n; i++) {
		double t = (range > 0.0) ? (hm[i] - minv) / range : 0.0;
		uint32_t v;
		if (t < 0.0)
			t = 0.0;
		else if (t > 1.0)
			t = 1.0;
		v = (uint32_t) (t * 255.0 + 0.5);
		pixels[i] = (0xffu << 24) | (v * 0x010101u);
	}

	rc = png_utils_write_png_image(filename, (unsigned char *) pixels, dim, dim, 1, 0);
	free(pixels);
	return rc;
}

static double *read_heightmap_png(const char *filename, int *out_dim)
{
	int w = 0, h = 0, a = 0, dim, bpp, row_bytes, x, y;
	char whynot[256];
	unsigned char *raw;
	double *hm;

	raw = (unsigned char *) png_utils_read_png_image(filename, 0, 0, 0, &w, &h, &a, whynot, sizeof(whynot));
	if (!raw) {
		fprintf(stderr, "pseudo-erosion: Failed to read '%s': %s\n", filename, whynot);
		exit(1);
	}

	dim = (w < h) ? w : h;
	bpp = a ? 4 : 3;
	row_bytes = w * bpp;
	if (row_bytes & 0x03)
		row_bytes += 4 - (row_bytes & 0x03);

	hm = allocate_heightmap(dim);
	for (y = 0; y < dim; y++) {
		const unsigned char *row = raw + (size_t) y * (size_t) row_bytes;
		for (x = 0; x < dim; x++) {
			const unsigned char *p = row + (size_t) x * (size_t) bpp;
			hm[(size_t) y * (size_t) dim + (size_t) x] =
				((double) p[0] + (double) p[1] + (double) p[2]) / (3.0 * 255.0);
		}
	}

	free(raw);
	*out_dim = dim;
	return hm;
}

/* Offsets for Moore neighborhood, including self */
static const int xo[] = { -1, 0, 1, 1, 1, 0, -1, -1, 0 };
static const int yo[] = { -1, -1, -1, 0, 1, 1, 1, 0, 0 };

static void setup_grid_points(struct osn_context *ctx, struct grid *grid, const double dim, const double feature_size)
{
	int i, x, y;
	double xoffset, yoffset;

	for (y = 0; y < grid->dim + 1; y++) {
		for (x = 0; x < grid->dim + 1; x++) {
			double ox, oy;
			ox = ((double) x * dim / (double) grid->dim / feature_size);
			oy = ((double) y * dim / (double) grid->dim / feature_size);
			xoffset = 0.5 * open_simplex_noise3(ctx, ox, oy, 25.7) * dim / grid->dim / feature_size;
			yoffset = 0.5 * open_simplex_noise3(ctx, ox, oy, 95.9) * dim / grid->dim / feature_size;
			gridpoint(grid, x, y)->x = ox + xoffset;
			gridpoint(grid, x, y)->y = oy + yoffset;
		}
	}
	/* Set up connections. Each grid point is "connected to" it's lowest neighbor,
	 * (possibly itself)
	 */
	for (y = 0; y < grid->dim + 1; y++) {
		for (x = 0; x < grid->dim + 1; x++) {
			int lown = -1;
			double lowest_value = 100000.0;
			/* Find the lowest neighbor, lown (index into xo[], yo[]) */
			for (i = 0; i < 9; i++) { /* Check Moore neighborhood */
				int nx, ny;
				double value;
				double px, py;
				nx = x + xo[i];
				ny = y + yo[i];
				if (nx < 0 || nx > grid->dim || ny < 0 || ny > grid->dim)
					continue;
				px = gridpoint(grid, nx, ny)->x;
				py = gridpoint(grid, nx, ny)->y;
				value = open_simplex_noise4(ctx, px, py, 0.0, 0.0);
				if (value < lowest_value) {
					lown = i;
					lowest_value = value;
				}
			}
			/* Set the connection to lowest neighbor */
			gridpoint(grid, x, y)->cx = x + xo[lown];
			gridpoint(grid, x, y)->cy = y + yo[lown];
		}
	}
}

static void setup_grid_points_from_heightmap(struct osn_context *ctx, struct grid *grid,
		const double dim, const double feature_size, const double *hm)
{
	int i, x, y;
	double xoffset, yoffset;
	int max_idx = (int) dim * (int) dim - 1;

	for (y = 0; y < grid->dim + 1; y++) {
		for (x = 0; x < grid->dim + 1; x++) {
			double ox, oy;
			ox = ((double) x * dim / (double) grid->dim / feature_size);
			oy = ((double) y * dim / (double) grid->dim / feature_size);
			xoffset = 0.5 * open_simplex_noise3(ctx, ox, oy, 25.7) * dim / grid->dim / feature_size;
			yoffset = 0.5 * open_simplex_noise3(ctx, ox, oy, 95.9) * dim / grid->dim / feature_size;
			gridpoint(grid, x, y)->x = ox + xoffset;
			gridpoint(grid, x, y)->y = oy + yoffset;
		}
	}
	/* Set up connections. Each grid point is "connected to" it's lowest neighbor,
	 * (possibly itself)
	 */
	for (y = 0; y < grid->dim + 1; y++) {
		for (x = 0; x < grid->dim + 1; x++) {
			int lown = -1;
			double lowest_value = 100000.0;
			/* Find the lowest neighbor, lown (index into xo[], yo[]) */
			for (i = 0; i < 9; i++) { /* Check Moore neighborhood */
				int nx, ny, idx;
				double value;
				double px, py;
				nx = x + xo[i];
				ny = y + yo[i];
				if (nx < 0 || nx > grid->dim || ny < 0 || ny > grid->dim)
					continue;
				px = gridpoint(grid, nx, ny)->x;
				py = gridpoint(grid, nx, ny)->y;
				idx = (int) (py * dim + px);
				if (idx < 0)
					idx = 0;
				else if (idx > max_idx)
					idx = max_idx;
				value = hm[idx];
				if (value < lowest_value) {
					lown = i;
					lowest_value = value;
				}
			}
			/* Set the connection to lowest neighbor */
			gridpoint(grid, x, y)->cx = x + xo[lown];
			gridpoint(grid, x, y)->cy = y + yo[lown];
		}
	}
}

static inline double sqr(double x)
{
	return x * x;
}

static void pseudo_erosion(double *hm, struct osn_context *ctx, struct grid *grid, int dim, float feature_size)
{
	int i, x, y, gx, gy, cx, cy, ngx, ngy;
	double f1, f2, x1, y1, x2, y2, px, py, h;

	for (y = 0; y < dim; y++) {
		ngy = grid->dim * y / dim;
		for (x = 0; x < dim; x++) { /* For each pixel... */
			double minh = 10000.0;
			ngx = grid->dim * x / dim;
			for (i = 0; i < 9; i++) {
				gx = ngx + xo[i];
				gy = ngy + yo[i];
				if (gx < 0 || gy < 0 || gx > grid->dim || gy > grid->dim)
					continue;
				px = (double) x / feature_size;
				py = (double) y / feature_size;
				x1 = gridpoint(grid, gx, gy)->x;
				y1 = gridpoint(grid, gx, gy)->y;
				cx = gridpoint(grid, gx, gy)->cx;
				cy = gridpoint(grid, gx, gy)->cy;
				x2 = gridpoint(grid, cx, cy)->x;
				y2 = gridpoint(grid, cx, cy)->y;
				f1 = ((y1 - y2) * (py - y1) + (x1 - x2) * (px - x1)) / (sqr(y1 - y2) + sqr(x1 - x2));
				if (f1 > 0.0) {
					h = sqrt(sqr(px - x1) + sqr(py - y1));
				} else if (f1 < -1.0) {
					h = sqrt(sqr(px - x2) + sqr(py - y2));
				} else {
					f2 = fabs(((y1 - y2) * (px - x1) - (x1 - x2) * (py - y1)) /
							sqrt(sqr(x1 - x2) + sqr(y1 - y2)));
					h = f2;
				}
				if (h < minh)
					minh = h;
			}
			hm[y * dim + x] = minh;
		}
		printf(".");
		fflush(stdout);
	}
	printf("\n");
	fflush(stdout);
}

static void process_int_option(char *option_name, char *option_value, int *value)
{
	int tmp;

	if (sscanf(option_value, "%d", &tmp) == 1) {
		*value = tmp;
	} else {
		fprintf(stderr, "Bad %s option '%s'\n", option_name, option_value);
		usage();
	}
}

static void process_options(int argc, char *argv[])
{
	int c;

	while (1) {
		int option_index = 0;
		c = getopt_long(argc, argv, "f:g:i:o:s:S:", long_options, &option_index);
		if (c == -1)
			break;
		switch (c) {
		case 'f':
			process_int_option("size", optarg, &feature_size);
			break;
		case 'g':
			process_int_option("size", optarg, &grid_size);
			break;
		case 'i':
			input_image = optarg;
			break;
		case 'o':
			output_file = optarg;
			break;
		case 's':
			process_int_option("size", optarg, &image_size);
			break;
		case 'S':
			process_int_option("seed", optarg, &seed);
			break;
		default:
			fprintf(stderr, "pseudo_erosion: Unknown option '%s'\n",
				option_index > 0 && option_index < argc &&
				argv[option_index] ? argv[option_index] : "(none)");
			usage();
			break;
		}
	}
	return;
}

/* Combine heightmaps a,b as 0.5*a + 0.25*b */
static void combine_heightmaps_f1(double *hm1, const double *hm2, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			int idx = y * imsize + x;
			hm1[idx] = 0.25 * hm2[idx] + 0.5 * hm1[idx];
		}
	}
}

/* Combine heightmaps a,b as a + sqr(b) */
static void combine_heightmaps_f2(double *hm1, const double *hm2, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			int idx = y * imsize + x;
			hm1[idx] = hm2[idx] * hm2[idx] + hm1[idx];
		}
	}
}

/* Combine heightmaps a,b,c as a + b * 0.5 * c */
static void combine_heightmaps_f3(double *hm1, const double *hm2, const double *hm3, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			int idx = y * imsize + x;
			hm1[idx] = hm1[idx] + hm2[idx] * 0.5 * hm3[idx];
		}
	}
}

/* Combine heightmaps a,b,c,d as a + sqrt(b * c) * 0.3333 * d */
static void combine_heightmaps_f4(double *hm1, const double *hm2, const double *hm3, const double *hm4, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			int idx = y * imsize + x;
			double prod = hm2[idx] * hm3[idx];
			hm1[idx] = hm1[idx] + sqrt(prod > 0.0 ? prod : 0.0) * 0.3333 * hm4[idx];
		}
	}
}

int main(int argc, char *argv[])
{
	double *img, *img2, *img3, *img4, *img5 = NULL;
	struct osn_context *ctx;
	struct grid *g, *g2, *g3, *g4, *g5;

	process_options(argc, argv);

	open_simplex_noise(seed, &ctx);
	printf("pseudo-erosion: Generating %d x %d heightmap image '%s'\n",
		image_size, image_size, output_file);
	g = allocate_grid(grid_size);
	/* First iteration, or input image */
	if (input_image) {
		img = read_heightmap_png(input_image, &image_size);
	} else {
		img = allocate_heightmap(image_size);
		setup_grid_points(ctx, g, image_size, feature_size);
		pseudo_erosion(img, ctx, g, image_size, feature_size);
	}

	write_heightmap_png("img-a.png", img, image_size);

	/* 2nd iteration */
	img2 = allocate_heightmap(image_size);
	g2 = allocate_grid(grid_size * 2);
	setup_grid_points(ctx, g2, image_size, feature_size / 2);
	pseudo_erosion(img2, ctx, g2, image_size, feature_size / 2);
	combine_heightmaps_f1(img, img2, image_size);

	write_heightmap_png("img2.png", img2, image_size);
	write_heightmap_png("img-b.png", img, image_size);

	/* 3rd iteration */
	img3 = allocate_heightmap(image_size);
	g3 = allocate_grid(grid_size * 4);
	setup_grid_points_from_heightmap(ctx, g3, image_size, feature_size / 4, img);
	pseudo_erosion(img3, ctx, g3, image_size, feature_size / 4);
	combine_heightmaps_f2(img, img3, image_size);

	write_heightmap_png("img3.png", img3, image_size);
	write_heightmap_png("img-c.png", img, image_size);

	/* 4th iteration */
	img4 = allocate_heightmap(image_size);
	g4 = allocate_grid(grid_size * 8);
	setup_grid_points_from_heightmap(ctx, g4, image_size, feature_size / 8, img);
	pseudo_erosion(img4, ctx, g4, image_size, feature_size / 8);
	combine_heightmaps_f3(img, img3, img4, image_size);

	write_heightmap_png("img4.png", img4, image_size);
	write_heightmap_png("img-d.png", img, image_size);

	/* 5th iteration */
	img5 = allocate_heightmap(image_size);
	g5 = allocate_grid(grid_size * 16);
	setup_grid_points_from_heightmap(ctx, g5, image_size, feature_size / 16, img);
	pseudo_erosion(img5, ctx, g5, image_size, feature_size / 16);
	combine_heightmaps_f4(img, img3, img4, img5, image_size);

	write_heightmap_png("img5.png", img5, image_size);
	write_heightmap_png("img-e.png", img, image_size);

	write_heightmap_png(output_file, img, image_size);
	open_simplex_noise_free(ctx);
	free_grid(g);
	free_grid(g2);
	free_grid(g3);
	free_grid(g4);
	free_grid(g5);
	free(img);
	free(img2);
	free(img3);
	free(img4);
	free(img5);
	return 0;
}
