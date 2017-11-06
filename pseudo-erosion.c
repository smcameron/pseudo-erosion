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
#define DEFAULT_FEATURE_SIZE 64
#define DEFAULT_GRID_SIZE 30

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

static uint32_t *allocate_image(int dim)
{
	unsigned char *image;

	image = malloc(4 * dim * dim);
	memset(image, 0, 4 * dim * dim);
	return (uint32_t *) image;
}

static inline uint32_t noise_to_color(double noise)
{
	uint32_t rgb = (0x0ff << 24) | (0x010101 * (uint32_t) ((noise + 1) * 127.5));
	return rgb;
}

static inline double color_to_noise(uint32_t color)
{
	int value = color & 0x0ff;
	return ((double) value / 127.5) - 1.0;
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

static void setup_grid_points_from_image(struct osn_context *ctx, struct grid *grid,
		const double dim, const double feature_size, uint32_t *image)
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
				value = color_to_noise(image[(int) (py * dim + px)]);
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

static void pseudo_erosion(uint32_t *image, struct osn_context *ctx, struct grid *grid, int dim, float feature_size)
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
			image[y * dim + x] = noise_to_color(minh);
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
		int option_index;
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

/* Combine images a,b as a + 0.5*b */
static void combine_images_f1(uint32_t *im1, uint32_t *im2, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			double n1, n2;
			uint32_t c1 = im1[y * imsize + x];
			uint32_t c2 = im2[y * imsize + x];
			n1 = color_to_noise(c1);
			n2 = color_to_noise(c2);
			n1 = 0.25 * n2 + 0.5 * n1;
			im1[y * imsize + x] = noise_to_color(n1);
		}
	}
}

/* Combine images a,b as a + sqr(b) */
static void combine_images_f2(uint32_t *im1, uint32_t *im2, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			double n1, n2;
			uint32_t c1 = im1[y * imsize + x];
			uint32_t c2 = im2[y * imsize + x];
			n1 = color_to_noise(c1);
			n2 = color_to_noise(c2);
			n1 = n2 * n2 + n1;
			im1[y * imsize + x] = noise_to_color(n1);
		}
	}
}

/* Combine images a,b,c as a + b * 0.5 * c */
static void combine_images_f3(uint32_t *im1, uint32_t *im2, uint32_t *im3, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			double n1, n2, n3;
			uint32_t c1 = im1[y * imsize + x];
			uint32_t c2 = im2[y * imsize + x];
			uint32_t c3 = im3[y * imsize + x];
			n1 = color_to_noise(c1);
			n2 = color_to_noise(c2);
			n3 = color_to_noise(c3);
			n1 = n1 + n2 * 0.5 * n3;
			im1[y * imsize + x] = noise_to_color(n1);
		}
	}
}

/* Combine images a,b,c,d as a + sqrt(b * c) * 0.3333 * d */
static void combine_images_f4(uint32_t *im1, uint32_t *im2, uint32_t *im3, uint32_t *im4, int imsize)
{
	int x, y;

	for (y = 0; y < imsize; y++) {
		for (x = 0; x < imsize; x++) {
			double n1, n2, n3, n4;
			uint32_t c1 = im1[y * imsize + x];
			uint32_t c2 = im2[y * imsize + x];
			uint32_t c3 = im3[y * imsize + x];
			uint32_t c4 = im4[y * imsize + x];
			n1 = color_to_noise(c1);
			n2 = color_to_noise(c2);
			n3 = color_to_noise(c3);
			n4 = color_to_noise(c4);
			n1 = n1 + sqrt(n2 * n3) * 0.3333 * n4;
			im1[y * imsize + x] = noise_to_color(n1);
		}
	}
}

int main(int argc, char *argv[])
{
	unsigned char *img, *img2, *img3, *img4, *img5 = NULL;
	struct osn_context *ctx;
	struct grid *g, *g2, *g3, *g4, *g5;

	process_options(argc, argv);

	open_simplex_noise(seed, &ctx);
	printf("pseudo-erosion: Generating %d x %d heightmap image '%s'\n",
		image_size, image_size, output_file);
	g = allocate_grid(grid_size);
	/* First iteration, or input image */
	if (input_image) {
		int w, h, a;
		char whynot[100];
		img = (unsigned char *) png_utils_read_png_image(input_image, 0, 0, 0, &w, &h, &a, whynot, 100);
		if (w < h)
			image_size = w;
		else
			image_size = h;
	} else {
		img = (unsigned char *) allocate_image(image_size);
		setup_grid_points(ctx, g, image_size, feature_size);
		pseudo_erosion((uint32_t *) img, ctx, g, image_size, feature_size);
	}

	/* 2nd iteration */
	img2 = (unsigned char *) allocate_image(image_size);
	g2 = allocate_grid(grid_size * 2);
	setup_grid_points(ctx, g2, image_size, feature_size);
	pseudo_erosion((uint32_t *) img2, ctx, g2, image_size, feature_size);
	combine_images_f1((uint32_t *) img, (uint32_t *) img2, image_size);

	/* 3rd iteration */
	img3 = (unsigned char *) allocate_image(image_size);
	g3 = allocate_grid(grid_size * 4);
	setup_grid_points_from_image(ctx, g3, image_size, feature_size, (uint32_t *) img2);
	pseudo_erosion((uint32_t *) img3, ctx, g3, image_size, feature_size);
	combine_images_f2((uint32_t *) img, (uint32_t *) img3, image_size);

	/* 4th iteration */
	img4 = (unsigned char *) allocate_image(image_size);
	g4 = allocate_grid(grid_size * 8);
	setup_grid_points_from_image(ctx, g4, image_size, feature_size, (uint32_t *) img3);
	pseudo_erosion((uint32_t *) img4, ctx, g4, image_size, feature_size);
	combine_images_f3((uint32_t *) img, (uint32_t *) img3, (uint32_t *) img4, image_size);

	/* 5th iteration */
	img5 = (unsigned char *) allocate_image(image_size);
	g5 = allocate_grid(grid_size * 16);
	setup_grid_points_from_image(ctx, g5, image_size, feature_size, (uint32_t *) img3);
	pseudo_erosion((uint32_t *) img5, ctx, g5, image_size, feature_size);
	combine_images_f4((uint32_t *) img, (uint32_t *) img3, (uint32_t *) img4, (uint32_t *) img5, image_size);

	png_utils_write_png_image(output_file, (unsigned char *) img, image_size, image_size, 1, 0);
	open_simplex_noise_free(ctx);
	free_grid(g);
	return 0;
}
