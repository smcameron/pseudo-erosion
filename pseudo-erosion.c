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

#include "open-simplex-noise.h"
#include "png_utils.h"

#define DIM 1024
#define FEATURE_SIZE 64
#define GRIDDIM 30

struct grid_point {
	double x, y;
	int cx, cy; /* connected to grid[cx][cy] */
} grid[GRIDDIM + 1][GRIDDIM + 1];

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

/* Offsets for Moore neighborhood, including self */
static const int xo[] = { -1, 0, 1, 1, 1, 0, -1, -1, 0 };
static const int yo[] = { -1, -1, -1, 0, 1, 1, 1, 0, 0 };

static void setup_grid_points(struct osn_context *ctx, const double dim, const double feature_size)
{
	int i, x, y;
	double xoffset, yoffset;

	for (y = 0; y < GRIDDIM + 1; y++) {
		for (x = 0; x < GRIDDIM + 1; x++) {
			xoffset = 0.7 * (double) rand() / (double) RAND_MAX * dim / GRIDDIM / feature_size;
			yoffset = 0.7 * (double) rand() / (double) RAND_MAX * dim / GRIDDIM / feature_size;
			grid[x][y].x = ((double) x * dim / (double) GRIDDIM / feature_size) + xoffset;
			grid[x][y].y = ((double) y * dim / (double) GRIDDIM / feature_size) + yoffset;
		}
	}
	/* Set up connections. Each grid point is "connected to" it's lowest neighbor,
	 * (possibly itself)
	 */
	for (y = 0; y < GRIDDIM + 1; y++) {
		for (x = 0; x < GRIDDIM + 1; x++) {
			int lown = -1;
			double lowest_value = 100000.0;
			/* Find the lowest neighbor, lown (index into xo[], yo[]) */
			for (i = 0; i < 9; i++) { /* Check Moore neighborhood */
				int nx, ny;
				double value;
				double px, py;
				nx = x + xo[i];
				ny = y + yo[i];
				if (nx < 0 || nx > GRIDDIM || ny < 0 || ny > GRIDDIM)
					continue;
				px = grid[nx][ny].x;
				py = grid[nx][ny].y;
				value = open_simplex_noise4(ctx, px, py, 0.0, 0.0);
				if (value < lowest_value) {
					lown = i;
					lowest_value = value;
				}
			}
			/* Set the connection to lowest neighbor */
			grid[x][y].cx = x + xo[lown];
			grid[x][y].cy = y + yo[lown];
		}
	}
}

static inline double sqr(double x)
{
	return x * x;
}

static void pseudo_erosion(uint32_t *image, struct osn_context *ctx, int dim, float feature_size)
{
	int i, x, y, gx, gy, cx, cy, ngx, ngy;
	double f1, f2, x1, y1, x2, y2, px, py, h;

	for (y = 0; y < dim; y++) {
		ngy = GRIDDIM * y / dim;
		for (x = 0; x < dim; x++) { /* For each pixel... */
			double minh = 10000.0;
			ngx = GRIDDIM * x / dim;
			for (i = 0; i < 9; i++) {
				gx = ngx + xo[i];
				gy = ngy + yo[i];
				if (gx < 0 || gy < 0 || gx > GRIDDIM || gy > GRIDDIM)
					continue;
				px = (double) x / feature_size;
				py = (double) y / feature_size;
				x1 = grid[gx][gy].x;
				y1 = grid[gx][gy].y;
				cx = grid[gx][gy].cx;
				cy = grid[gx][gy].cy;
				x2 = grid[cx][cy].x;
				y2 = grid[cx][cy].y;
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
}

int main(int argc, char *argv[])
{
	char *filename = "output.png";
	int seed = 123456;
	unsigned char *image = NULL;
	struct osn_context *ctx;
	open_simplex_noise(seed, &ctx);

	printf("pseudo-erosion: Generating %d x %d heightmap image '%s'\n", DIM, DIM, filename);
	image = (unsigned char *) allocate_image(DIM);
	setup_grid_points(ctx, DIM, FEATURE_SIZE);
	pseudo_erosion((uint32_t *) image, ctx, DIM, FEATURE_SIZE);
	png_utils_write_png_image(filename, (unsigned char *) image, DIM, DIM, 1, 0);
	open_simplex_noise_free(ctx);
	return 0;
}
