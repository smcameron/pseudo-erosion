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
#define DEFAULT_FEATURE_SIZE 256
#define DEFAULT_GRID_SIZE 8
#define DEFAULT_WARP 0.25

static char *output_file = "output.png";
static int image_size = DEFAULT_IMAGE_SIZE;
static int feature_size = DEFAULT_FEATURE_SIZE;
static int grid_size = DEFAULT_GRID_SIZE;
static int seed = 123456;
static double warp_amount = DEFAULT_WARP;
static char *input_image = NULL;

static struct option long_options[] = {
	{ "featuresize", required_argument, NULL, 'f' },
	{ "gridsize", required_argument, NULL, 'g' },
	{ "size", required_argument, NULL, 's' },
	{ "seed", required_argument, NULL, 'S' },
	{ "warp", required_argument, NULL, 'w' },
	{ "outputfile", required_argument, NULL, 'o' },
	{ "input", required_argument, NULL, 'i' },
	{ 0, 0, 0, 0 },
};

static void usage(void)
{
	fprintf(stderr, "pseudo_erosion: Usage:\n\n");
	fprintf(stderr, "	pseudo_erosion [-g gridsize] [-o outputfile] [-s imagesize] \\\n");
	fprintf(stderr, "		[-i inputfile] [-f featuresize] [-S seed] [-w warp]\n");
	fprintf(stderr, "\n");
	exit(1);
}

#define GRID_PAD 2

struct grid_point {
	double x, y, h;
	int cx, cy; /* connected to gridpoint(grid, cx, cy) */
};

struct grid {
	struct grid_point *g;
	int dim;
	int stride;
};

static struct grid *allocate_grid(int dim)
{
	struct grid_point *gp;
	struct grid *g;
	int stride = dim + 1 + 2 * GRID_PAD;

	gp = calloc((size_t) stride * (size_t) stride, sizeof(*gp));
	if (!gp) {
		fprintf(stderr, "pseudo-erosion: Out of memory allocating %dx%d grid\n", dim, dim);
		exit(1);
	}
	g = malloc(sizeof(*g));
	if (!g) {
		fprintf(stderr, "pseudo-erosion: Out of memory allocating grid struct\n");
		exit(1);
	}
	g->g = gp;
	g->dim = dim;
	g->stride = stride;
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
	return &grid->g[(y + GRID_PAD) * grid->stride + (x + GRID_PAD)];
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

static double hash_offset(int seed_val, int layer, int gx, int gy, int axis)
{
	uint64_t h = (uint64_t) (uint32_t) seed_val;
	h ^= (uint64_t) (uint32_t) layer * 0x9e3779b97f4a7c15ULL;
	h ^= (uint64_t) (uint32_t) gx * 0xbf58476d1ce4e5b9ULL;
	h ^= (uint64_t) (uint32_t) gy * 0x94d049bb133111ebULL;
	h ^= (uint64_t) (uint32_t) axis * 0x517cc1b727220a95ULL;
	h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
	h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
	h = h ^ (h >> 31);
	return ((double) (h >> 11) * (1.0 / 9007199254740992.0) - 0.5) * 0.9;
}

static double sample_heightmap(const double *hm, int img_dim, double fx, double fy)
{
	int x0, y0, x1, y1;
	double tx, ty, v00, v10, v01, v11;

	if (img_dim <= 0)
		return 0.0;
	if (fx < 0.0)
		fx = 0.0;
	else if (fx > (double) (img_dim - 1))
		fx = (double) (img_dim - 1);
	if (fy < 0.0)
		fy = 0.0;
	else if (fy > (double) (img_dim - 1))
		fy = (double) (img_dim - 1);

	x0 = (int) floor(fx);
	y0 = (int) floor(fy);
	x1 = (x0 + 1 < img_dim) ? x0 + 1 : x0;
	y1 = (y0 + 1 < img_dim) ? y0 + 1 : y0;
	tx = fx - (double) x0;
	ty = fy - (double) y0;

	v00 = hm[(size_t) y0 * (size_t) img_dim + (size_t) x0];
	v10 = hm[(size_t) y0 * (size_t) img_dim + (size_t) x1];
	v01 = hm[(size_t) y1 * (size_t) img_dim + (size_t) x0];
	v11 = hm[(size_t) y1 * (size_t) img_dim + (size_t) x1];

	return (v00 * (1.0 - tx) + v10 * tx) * (1.0 - ty) +
	       (v01 * (1.0 - tx) + v11 * tx) * ty;
}

/* Set up connections. Each grid point is "connected to" its lowest neighbor,
 * (possibly itself).
 */
static void connect_lowest_neighbors(struct grid *grid)
{
	int i, x, y;

	for (y = -GRID_PAD; y <= grid->dim + GRID_PAD; y++) {
		for (x = -GRID_PAD; x <= grid->dim + GRID_PAD; x++) {
			int lown = 8;
			double lowest_value = HUGE_VAL;
			struct grid_point *gp = gridpoint(grid, x, y);

			for (i = 0; i < 9; i++) {
				int nx = x + xo[i];
				int ny = y + yo[i];
				double value;

				if (nx < -GRID_PAD || nx > grid->dim + GRID_PAD ||
				    ny < -GRID_PAD || ny > grid->dim + GRID_PAD)
					continue;
				value = gridpoint(grid, nx, ny)->h;
				if (value < lowest_value) {
					lown = i;
					lowest_value = value;
				}
			}
			gp->cx = x + xo[lown];
			gp->cy = y + yo[lown];
		}
	}
}

static void generate_noise_heightmap(struct osn_context *ctx, double *hm,
		int dim, double feat_size, double z_offset)
{
	int x, y;
	double minv = HUGE_VAL, maxv = -HUGE_VAL, range;
	size_t i, n = (size_t) dim * (size_t) dim;

	for (y = 0; y < dim; y++) {
		for (x = 0; x < dim; x++) {
			double v = open_simplex_noise3(ctx, (double) x / feat_size,
					(double) y / feat_size, z_offset);
			hm[(size_t) y * (size_t) dim + (size_t) x] = v;
			if (v < minv)
				minv = v;
			if (v > maxv)
				maxv = v;
		}
	}

	range = maxv - minv;
	if (range > 0.0) {
		for (i = 0; i < n; i++)
			hm[i] = (hm[i] - minv) / range;
	} else {
		for (i = 0; i < n; i++)
			hm[i] = 0.0;
	}
}

static void setup_grid_points_from_heightmap(struct grid *grid, int img_dim,
		const double *hm, int layer)
{
	int x, y;
	double cell_to_px = (double) img_dim / (double) grid->dim;

	for (y = -GRID_PAD; y <= grid->dim + GRID_PAD; y++) {
		for (x = -GRID_PAD; x <= grid->dim + GRID_PAD; x++) {
			struct grid_point *gp = gridpoint(grid, x, y);
			double fx, fy;
			gp->x = (double) x + hash_offset(seed, layer, x, y, 0);
			gp->y = (double) y + hash_offset(seed, layer, x, y, 1);
			fx = gp->x * cell_to_px - 0.5;
			fy = gp->y * cell_to_px - 0.5;
			gp->h = sample_heightmap(hm, img_dim, fx, fy);
		}
	}
	connect_lowest_neighbors(grid);
}

static void pseudo_erosion(double *hm, struct osn_context *ctx, struct grid *grid,
		int dim, int layer, double warp)
{
	int x, y, dx, dy, gx, gy, ngx, ngy;
	double maxh = 0.0;
	double z_base = (double) layer * 31.7;
	size_t i, n = (size_t) dim * (size_t) dim;

	for (y = 0; y < dim; y++) {
		double base_py = ((double) y + 0.5) * (double) grid->dim / (double) dim;
		for (x = 0; x < dim; x++) { /* For each pixel... */
			double px = ((double) x + 0.5) * (double) grid->dim / (double) dim;
			double py = base_py;
			double minh2 = HUGE_VAL;
			double h;

			if (warp > 0.0) {
				double wx = open_simplex_noise3(ctx, px * 0.5, py * 0.5, z_base + 11.3);
				double wy = open_simplex_noise3(ctx, px * 0.5, py * 0.5, z_base + 57.9);
				px += warp * wx;
				py += warp * wy;
			}

			ngx = (int) floor(px + 0.5);
			ngy = (int) floor(py + 0.5);
			for (dy = -2; dy <= 2; dy++) {
				gy = ngy + dy;
				if (gy < -GRID_PAD || gy > grid->dim + GRID_PAD)
					continue;
				for (dx = -2; dx <= 2; dx++) {
					struct grid_point *gp, *cp;
					double x1, y1, x2, y2, vx, vy, len2, ex, ey, d2;

					gx = ngx + dx;
					if (gx < -GRID_PAD || gx > grid->dim + GRID_PAD)
						continue;
					gp = gridpoint(grid, gx, gy);
					cp = gridpoint(grid, gp->cx, gp->cy);
					x1 = gp->x;
					y1 = gp->y;
					x2 = cp->x;
					y2 = cp->y;
					vx = x2 - x1;
					vy = y2 - y1;
					len2 = vx * vx + vy * vy;
					if (len2 <= 1e-12) {
						ex = px - x1;
						ey = py - y1;
					} else {
						double t = ((px - x1) * vx + (py - y1) * vy) / len2;
						if (t < 0.0)
							t = 0.0;
						else if (t > 1.0)
							t = 1.0;
						ex = px - (x1 + t * vx);
						ey = py - (y1 + t * vy);
					}
					d2 = ex * ex + ey * ey;
					if (d2 < minh2)
						minh2 = d2;
				}
			}
			h = sqrt(minh2);
			hm[(size_t) y * (size_t) dim + (size_t) x] = h;
			if (h > maxh)
				maxh = h;
		}
		printf(".");
		fflush(stdout);
	}
	printf("\n");
	fflush(stdout);

	if (maxh > 0.0) {
		for (i = 0; i < n; i++)
			hm[i] /= maxh;
	}
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

static void process_double_option(char *option_name, char *option_value, double *value)
{
	double tmp;

	if (sscanf(option_value, "%lf", &tmp) == 1) {
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
		c = getopt_long(argc, argv, "f:g:i:o:s:S:w:", long_options, &option_index);
		if (c == -1)
			break;
		switch (c) {
		case 'f':
			process_int_option("featuresize", optarg, &feature_size);
			break;
		case 'g':
			process_int_option("gridsize", optarg, &grid_size);
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
		case 'w':
			process_double_option("warp", optarg, &warp_amount);
			break;
		default:
			fprintf(stderr, "pseudo_erosion: Unknown option '%s'\n",
				option_index > 0 && option_index < argc &&
				argv[option_index] ? argv[option_index] : "(none)");
			usage();
			break;
		}
	}
	if (image_size <= 0 || feature_size <= 0 || grid_size <= 0 || warp_amount < 0.0) {
		fprintf(stderr, "pseudo_erosion: size, featuresize, and gridsize must be positive and warp >= 0\n");
		usage();
	}
	return;
}

/* Combine heightmaps a,b as a + 0.5*b */
static void combine_heightmaps_f1(double *hm1, const double *hm2, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			int idx = y * imsize + x;
			hm1[idx] = hm1[idx] + 0.5 * hm2[idx];
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

/* Combine heightmaps a,b,c,d as a + sqrt(b * c) * d / 3 */
static void combine_heightmaps_f4(double *hm1, const double *hm2, const double *hm3, const double *hm4, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			int idx = y * imsize + x;
			double prod = hm2[idx] * hm3[idx];
			hm1[idx] = hm1[idx] + sqrt(prod > 0.0 ? prod : 0.0) * hm4[idx] / 3.0;
		}
	}
}

int main(int argc, char *argv[])
{
	double *img = NULL, *img2 = NULL, *img3 = NULL, *img4 = NULL, *img5 = NULL;
	struct osn_context *ctx;
	struct grid *g1, *g2, *g3;

	process_options(argc, argv);

	open_simplex_noise(seed, &ctx);
	printf("pseudo-erosion: Generating %d x %d heightmap image '%s'\n",
		image_size, image_size, output_file);

	/* Base heightmap: either input image, or p1 + p2 / 2 */
	if (input_image) {
		img = read_heightmap_png(input_image, &image_size);
		write_heightmap_png("img-a.png", img, image_size);
		write_heightmap_png("img-b.png", img, image_size);
	} else {
		img = allocate_heightmap(image_size);
		generate_noise_heightmap(ctx, img, image_size, feature_size, 0.0);
		write_heightmap_png("img-a.png", img, image_size);

		img2 = allocate_heightmap(image_size);
		generate_noise_heightmap(ctx, img2, image_size, (double) feature_size / 2.0, 17.3);
		write_heightmap_png("img2.png", img2, image_size);

		combine_heightmaps_f1(img, img2, image_size);
		write_heightmap_png("img-b.png", img, image_size);
	}

	/* 1st erosion pass (r1, default 8x8): map = map + sqr(r1) */
	img3 = allocate_heightmap(image_size);
	g1 = allocate_grid(grid_size);
	setup_grid_points_from_heightmap(g1, image_size, img, 1);
	pseudo_erosion(img3, ctx, g1, image_size, 1, warp_amount);
	combine_heightmaps_f2(img, img3, image_size);

	write_heightmap_png("img3.png", img3, image_size);
	write_heightmap_png("img-c.png", img, image_size);

	/* 2nd erosion pass (r2, default 32x32): map = map + r1 * r2 / 2 */
	img4 = allocate_heightmap(image_size);
	g2 = allocate_grid(grid_size * 4);
	setup_grid_points_from_heightmap(g2, image_size, img, 2);
	pseudo_erosion(img4, ctx, g2, image_size, 2, warp_amount);
	combine_heightmaps_f3(img, img3, img4, image_size);

	write_heightmap_png("img4.png", img4, image_size);
	write_heightmap_png("img-d.png", img, image_size);

	/* 3rd erosion pass (r3, default 128x128): map = map + sqrt(r1 * r2) * r3 / 3 */
	img5 = allocate_heightmap(image_size);
	g3 = allocate_grid(grid_size * 16);
	setup_grid_points_from_heightmap(g3, image_size, img, 3);
	pseudo_erosion(img5, ctx, g3, image_size, 3, warp_amount);
	combine_heightmaps_f4(img, img3, img4, img5, image_size);

	write_heightmap_png("img5.png", img5, image_size);
	write_heightmap_png("img-e.png", img, image_size);

	write_heightmap_png(output_file, img, image_size);
	open_simplex_noise_free(ctx);
	free_grid(g1);
	free_grid(g2);
	free_grid(g3);
	free(img);
	free(img2);
	free(img3);
	free(img4);
	free(img5);
	return 0;
}
