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

#include "open-simplex-noise.h"
#include "png_utils.h"

#define DIM 1024
#define FEATURE_SIZE 64

static uint32_t *allocate_image(int dim)
{
	unsigned char *image;

	image = malloc(4 * dim * dim);
	memset(image, 0, 4 * dim * dim);
	return (uint32_t *) image;
}

static unsigned char *generate_initial_noise_image(struct osn_context *ctx, int dim, double feature_size)
{
	uint32_t *image = allocate_image(dim);
	int x, y;
	double value;
	uint32_t rgb;

	for (y = 0; y < dim; y++) {
		for (x = 0; x < dim; x++) {
			value = open_simplex_noise4(ctx,
				(double) x / feature_size, (double) y / feature_size, 0.0, 0.0);
			rgb = 0x010101 * (uint32_t) ((value + 1) * 127.5);
			image[y * dim + x] = (0x0ff << 24) | (rgb);
		}
	}
	return (unsigned char *) image;
}

int main(int argc, char *argv[])
{
	char *filename = "output.png";
	int seed = 123456;
	unsigned char *image = NULL;
	struct osn_context *ctx;
	open_simplex_noise(seed, &ctx);

	printf("pseudo-erosion: Generating %d x %d heightmap image '%s'\n", DIM, DIM, filename);
	image = generate_initial_noise_image(ctx, DIM, FEATURE_SIZE);
	png_utils_write_png_image(filename, (unsigned char *) image, DIM, DIM, 1, 0);
	open_simplex_noise_free(ctx);
	return 0;
}
