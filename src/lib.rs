//! Implements the Marching Cubes algorithm for generating Isosurfaces from volume data.

mod tables;

use graphics::{Mesh, Vertex};
use lin_alg::f32::Vec3;

use crate::tables::{CUBE_CORNER_OFFSETS, EDGE_TABLE, EDGE_VERTEX_PAIRS, TRI_TABLE};

pub trait GridPoint {
    fn value(&self) -> f64;
}

// todo: Use this custom type, and make a custom Vertex type, to make this library more accessible, i.e.
// todo: to not depend on `graphics`.

#[derive(Debug)]
pub struct Vertex_ {
    pub posit: Vec3,
    pub normal: Vec3,
}

/// Represents a mesh output by the algorithm.
#[derive(Debug)]
pub struct Mesh_ {
    pub vertices: Vec<Vertex_>,
    /// Grouped into triangles.
    pub indices: Vec<usize>,
}

/// Quantise a Vec3 to integer bins so we can use it as a HashMap key.
/// 1 Å equals one “voxel unit” in the current grid, so 1-to-1 quantisation
/// is fine.  If you change the coordinate scaling later, update this!
fn q(v: Vec3) -> (i32, i32, i32) {
    (v.x.round() as i32, v.y.round() as i32, v.z.round() as i32)
}

#[derive(Debug)]
pub struct MarchingCubes {
    pub dims: (usize, usize, usize),
    pub data: Vec<f32>,
    pub iso_level: f32,
    scale: [f32; 3],
}

impl MarchingCubes {
    /// Build a `MarchingCubes` grid from the voxel list returned by `read_map_data`.
    ///
    /// * `hdr` supplies the grid dimensions and some basic statistics.
    /// * `density` **must** be in the same `x-fast, y-medium, z-slow` order that
    ///   `read_map_data()` produces (it is, because the nested `for k { for j { for i { … }}}`
    ///   loop matches our own indexing scheme).
    ///
    /// The default iso-level is the mean density stored in the header.  You can
    /// pick another threshold later by doing `mc.iso_level = …;`.
    pub fn new<T: GridPoint>(
        dims: (usize, usize, usize),
        size: (f32, f32, f32),
        sampling_interval: (f32, f32, f32),
        density: &[T],
        iso_level: f32,
    ) -> Self {
        let expected_len = dims.0 * dims.1 * dims.2;
        assert_eq!(
            density.len(),
            expected_len,
            "Density array has {} points, but header implies {}",
            density.len(),
            expected_len
        );

        // Strip the coordinate part – Marching Cubes needs only the scalar field.
        let data: Vec<f32> = density.iter().map(|d| d.value() as f32).collect();

        let scale = [
            size.0 / sampling_interval.0,
            size.1 / sampling_interval.1,
            size.2 / sampling_interval.2,
        ];

        Self {
            dims,
            data,
            iso_level,
            scale,
        }
    }

    /// Central-difference gradient of the scalar field at integer voxel coords.
    /// Assumes the volume is at least 2×2×2; falls back to fwd/bwd diff at borders.
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

    pub fn generate(&self) -> Mesh {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        let (nx, ny, nz) = self.dims;

        for x in 0..(nx - 1) {
            for y in 0..(ny - 1) {
                for z in 0..(nz - 1) {
                    // Corner densities
                    let mut cube = [0f32; 8];
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

                        // 2.b gradients (central diff) at corners a & b   ↓↓↓
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

                        // Interpolate gradient, then normalise & flip (-∇ρ points “out”)
                        let g = ga + (gb - ga) * mu;
                        norm_list[i] = (-g).to_normalized();
                    }

                    for tri in TRI_TABLE[cube_index].chunks(3) {
                        if tri[0] == -1 {
                            break;
                        }

                        for &edge_id in tri {
                            let p = pos_list[edge_id as usize];
                            let n = norm_list[edge_id as usize];

                            vertices.push(Vertex::new([p.x, p.y, p.z], n));
                            indices.push(vertices.len() - 1);
                        }
                    }
                }
            }
        }

        Mesh {
            vertices,
            indices,
            material: 0,
        }
    }

    fn get_value(&self, x: usize, y: usize, z: usize) -> f32 {
        let (nx, ny, _) = self.dims;
        self.data[x + y * nx + z * nx * ny]
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

fn interpolate(p1: Vec3, p2: Vec3, valp1: f32, valp2: f32, iso: f32) -> Vec3 {
    if (iso - valp1).abs() < 1e-6 {
        return p1;
    }
    if (iso - valp2).abs() < 1e-6 {
        return p2;
    }
    let mu = (iso - valp1) / (valp2 - valp1);
    p1 + (p2 - p1) * mu
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
