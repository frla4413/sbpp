//============================================================
// Name        : mesh.hpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : hpp file for a few meshes.
//============================================================

/*
 * A Block contains two Arrays: x and y.
 * x and y represent the coordinates of the block.
 * A single-block mesh consists of a single Block.
 *
 * A multi-block mesh is a std::vector with several
 * connected blocks, i.e. std::vector<Block>.
 *
 * Single-block meshes
 *  Call:
 *  where the right-hand side is one of the following:
 *  o Block CartesianGrid(Nx, Ny)
 *  o CircleSector(Nr, Nth, r_in, r_out, angle_0, angl1)
 *
 *  Ex:  auto block = CartesianGrid(Nx, Ny)
 *
 * Multi-block meshes
 *  o Annulus(N, r_in, r_out)
 *
 *  Ex: auto blocks = MultiBlockMesh(param),
 */

#pragma once
#include "array.hpp"

struct Block {
  Array x, y;
};

Block CartesianGrid(int Nx, int Ny);
Block CircleSector(int Nr, int Nth, double r_in, double r_out,
                   double angle_0, double andgle_1);

//Multiblock meshes
std::vector<Block> Annulus(int N, double r_in, double r_out);
