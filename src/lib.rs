#![allow(clippy::needless_range_loop)]

//! Implements the Marching Cubes algorithm for generating Isosurfaces from volume data.
//! See the readme for example uses.

mod tables;

use std::{io, io::ErrorKind, time::Instant};

// use std::time::Instant;
use lin_alg::f32::Vec3;

// #[cfg(target_arch = "x86_64")]
// use lin_alg::f32::{f32x8, Vec3x8};
use crate::tables::{CUBE_CORNER_OFFSETS, EDGE_TABLE, EDGE_VERTEX_PAIRS, TRI_TABLE};

pub trait GridPoint {
    fn value(&self) -> f64;
}

/// Represents a Vertex output by the algorithm.
#[derive(Debug)]
pub struct Vertex {
    pub posit: Vec3,
    pub normal: Vec3,
}

/// Represents a mesh output by the algorithm. Convert this to your downstream application's Mesh.
#[derive(Debug)]
pub struct Mesh {
    pub vertices: Vec<Vertex>,
    /// Grouped into triangles, each 3 indices.
    pub indices: Vec<usize>,
}

#[derive(Debug)]
pub struct MarchingCubes {
    pub dims: (usize, usize, usize),
    pub values: Vec<f32>,
    pub iso_level: f32,
    pub offset: Vec3,
    scale: [f32; 3],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MeshSide {
    Both,
    OutsideOnly,
    InsideOnly,
}

impl MarchingCubes {
    /// Initialize with grid data.
    ///  `hdr` supplies the grid dimensions and some basic statistics.
    /// `values` **must** be in the same `x-fast, y-medium, z-slow` order that
    /// `read_map_data()` produces (it is, because the nested `for k { for j { for i { … }}}`
    /// loop matches our own indexing scheme).
    pub fn new(
        dims: (usize, usize, usize),
        size: (f32, f32, f32),
        sampling_interval: (f32, f32, f32),
        offset: Vec3,
        values: Vec<f32>,
        iso_level: f32,
    ) -> io::Result<Self> {
        let expected_len = dims.0 * dims.1 * dims.2;

        if values.len() != expected_len {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!(
                    "Density array has {} points, but header implies {expected_len}.",
                    values.len()
                ),
            ));
        }

        let scale = [
            size.0 / sampling_interval.0,
            size.1 / sampling_interval.1,
            size.2 / sampling_interval.2,
        ];

        Ok(Self {
            dims,
            values,
            iso_level,
            offset,
            scale,
        })
    }

    /// Constructor that uses an `impl GridPoint` vice `Vec<f32>`. This may be desired depending
    /// on your grid representation. Same as `new` otherwise.
    pub fn from_gridpoints<T: GridPoint>(
        dims: (usize, usize, usize),
        size: (f32, f32, f32),
        sampling_interval: (f32, f32, f32),
        offset: Vec3,
        values: &[T],
        iso_level: f32,
    ) -> io::Result<Self> {
        let values_vec = values.iter().map(|d| d.value() as f32).collect();
        Self::new(dims, size, sampling_interval, offset, values_vec, iso_level)
    }

    /// Central-difference gradient of the scalar field at integer voxel coords.
    /// Assumes the volume is at least 2×2×2; falls back to fwd/bwd diff at borders.
    /// We use this to compute normal vectors.
    fn gradient(&self, x: usize, y: usize, z: usize) -> Vec3 {
        let (nx, ny, nz) = self.dims;

        let gx = match (x > 0, x + 1 < nx) {
            (true, true) => (self.get_value(x + 1, y, z) - self.get_value(x - 1, y, z)) * 0.5,
            (false, true) => self.get_value(x + 1, y, z) - self.get_value(x, y, z),
            (true, false) => self.get_value(x, y, z) - self.get_value(x - 1, y, z),
            (false, false) => 0.0,
        };
        let gy = match (y > 0, y + 1 < ny) {
            (true, true) => (self.get_value(x, y + 1, z) - self.get_value(x, y - 1, z)) * 0.5,
            (false, true) => self.get_value(x, y + 1, z) - self.get_value(x, y, z),
            (true, false) => self.get_value(x, y, z) - self.get_value(x, y - 1, z),
            (false, false) => 0.0,
        };
        let gz = match (z > 0, z + 1 < nz) {
            (true, true) => (self.get_value(x, y, z + 1) - self.get_value(x, y, z - 1)) * 0.5,
            (false, true) => self.get_value(x, y, z + 1) - self.get_value(x, y, z),
            (true, false) => self.get_value(x, y, z) - self.get_value(x, y, z - 1),
            (false, false) => 0.0,
        };

        Vec3::new(gx, gy, gz)
    }

    /// Run this to generate the mesh; this function contains the primary algorithm.
    pub fn generate(&self, mesh_side: MeshSide) -> Mesh {
        // Note: We observed slowdowns vice speedups with Rayon, for electron density.
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        let (nx, ny, nz) = self.dims;

        let start = Instant::now();

        for x in 0..(nx - 1) {
            for y in 0..(ny - 1) {
                for z in 0..(nz - 1) {
                    // Corner densities
                    let mut cube = [0.; 8];

                    // #[cfg(target_arch = "x86_64")]
                    // let mut cube = f32x8::splat(0.);

                    for i in 0..8 {
                        let (dx, dy, dz) = CUBE_CORNER_OFFSETS[i];
                        cube[i] = self.get_value(x + dx, y + dy, z + dz);
                    }

                    let cube_index = compute_cube_index(&cube, self.iso_level);
                    if EDGE_TABLE[cube_index] == 0 {
                        continue;
                    }

                    // Edge vertices and normals
                    let mut pos_list = [Vec3::new_zero(); 12];
                    let mut norm_list = [Vec3::new_zero(); 12];

                    for i in 0..12 {
                        if (EDGE_TABLE[cube_index] & (1 << i)) == 0 {
                            continue;
                        }

                        let (a, b) = EDGE_VERTEX_PAIRS[i];

                        // Corner positions (voxel coords, not Å)
                        let pa = self.corner_pos(x, y, z, a);
                        let pb = self.corner_pos(x, y, z, b);

                        let va = cube[a];
                        let vb = cube[b];

                        // Interp factor μ ∈ [0,1]
                        let mu = (self.iso_level - va) / (vb - va);

                        pos_list[i] = pa + (pb - pa) * mu;

                        // Gradients (central difference) at corners a & b
                        let (ax, ay, az) = {
                            let (dx, dy, dz) = CUBE_CORNER_OFFSETS[a];
                            (x + dx, y + dy, z + dz)
                        };
                        let (bx, by, bz) = {
                            let (dx, dy, dz) = CUBE_CORNER_OFFSETS[b];
                            (x + dx, y + dy, z + dz)
                        };
                        let ga = self.gradient(ax, ay, az);
                        let gb = self.gradient(bx, by, bz);

                        // Interpolate gradient, then normalise & flip (-∇ρ points out)
                        let g = ga + (gb - ga) * mu;
                        norm_list[i] = (-g).to_normalized();
                    }

                    for tri_inside in TRI_TABLE[cube_index].chunks(3) {
                        if tri_inside[0] == -1 {
                            break;
                        }

                        let tri_inside = tri_inside.to_vec();

                        // Flip the index order from the table, to get correct-oriented faces.
                        let tri_outside = if mesh_side == MeshSide::InsideOnly {
                            Vec::new()
                        } else {
                            let mut t = tri_inside.clone();
                            let orig_0 = tri_inside[0];
                            t[0] = tri_inside[1];
                            t[1] = orig_0;
                            t
                        };

                        let tri_sets = match mesh_side {
                            MeshSide::Both => vec![&tri_outside, &tri_inside],
                            MeshSide::InsideOnly => vec![&tri_inside],
                            MeshSide::OutsideOnly => vec![&tri_outside],
                        };

                        for set in tri_sets {
                            for &edge_id in set {
                                let posit = pos_list[edge_id as usize] + self.offset;
                                let normal = norm_list[edge_id as usize];

                                vertices.push(Vertex { posit, normal });
                                indices.push(vertices.len() - 1);
                            }
                        }
                    }
                }
            }
        }

        // let elapsed = start.elapsed();
        // println!("Time taken for cubes: {:?}μs", elapsed.as_micros());

        Mesh { vertices, indices }
    }

    fn get_value(&self, x: usize, y: usize, z: usize) -> f32 {
        let (nx, ny, _) = self.dims;
        self.values[x + y * nx + z * nx * ny]
    }

    fn corner_pos(&self, x: usize, y: usize, z: usize, corner: usize) -> Vec3 {
        let (dx, dy, dz) = CUBE_CORNER_OFFSETS[corner];
        Vec3::new(
            (x + dx) as f32 * self.scale[0],
            (y + dy) as f32 * self.scale[1],
            (z + dz) as f32 * self.scale[2],
        )
    }
}

fn compute_cube_index(cube: &[f32; 8], iso: f32) -> usize {
    let mut idx = 0;
    for i in 0..8 {
        if cube[i] < iso {
            idx |= 1 << i;
        }
    }
    idx
}
