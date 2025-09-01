# Mcubes: Marching Cubes Isosurface library

[![Crate](https://img.shields.io/crates/v/mcubes.svg)](https://crates.io/crates/mcubes)
[![Docs](https://docs.rs/mcubes/badge.svg)](https://docs.rs/mcubes)
[![PyPI](https://img.shields.io/pypi/v/mcubes.svg)](https://pypi.org/project/mcubes)


Uses the Marching Cubes algorithm to create isosurfaces in Rust or Python, from volume data. Designed to be easy to integrate
in applications.

![Electron density demo](screenshots/daedalus_iso_a.png)

Based loosely on [PyMarchingCubes](https://github.com/JustusThies/PyMarchingCubes).

Outputs a native `Mesh`, which contains native `Vertex`s. In practice, you will convert these to whatever
mesh struct your application uses. For example, [graphics::Mesh](https://docs.rs/graphics/latest/graphics/struct.Mesh.html).

Uses [lin-alg](https://github.com/david-oconnor/lin-alg) for its `Vec3` type; this is the library's only dependency.

![Surface demo](screenshots/surface_mesh_transparent.png)

Used by the [Daedalus molecule viewer](https://github.com/David-OConnor/daedalus) to view experimentally-derived electron density from protein crystals.

The grid must be  regularly spaced, along 3 orthogonal axes. Values are either a `Vec<f32>`, or points which impl
`mcubes::GridPoint`. This trait contains a single method: To get the value at that point.

Example creating a solvent-accessible-surface mesh by setting the ISO level to 0, and 
rendering only the vertices.
![Surface mesh](screenshots/surface_a.png)

Example use:
```rust
use mcubes::{GridPoint, MarchingCubes, MeshSide};

/// An example data struct from your application. If you use something like this 
/// to represent a point, you may with to use the `from_gridpoints` constructor.
pub struct ElectronDensity {
    pub coords: Vec3,
    pub density: f64,
}

impl GridPoint for ElectronDensity {
    fn value(&self) -> f64 { self.density }
}

fn create_mesh(hdr: &MapHeader, density: &[ElectronDensity], iso_level: f32) {
    let mc = MarchingCubes::from_gridpoints(
        // Number of grid points along each axis
        (hdr.nx as usize, hdr.ny as usize, hdr.nz as usize),
        // Grid dimensions per unit cell along each axis
        (hdr.cell[0], hdr.cell[1], hdr.cell[2]),
        // Sampling interval along each axis. (Usually the same as grid point number of grid points.)
        (hdr.mx as f32, hdr.my as f32, hdr.mz as f32),
        density,
        // The value to draw the isosurface at.
        iso_level,
    );
    
    // If you are using a `Vec<f32` for data instead if `GridPoint`, use this constructor:
    let mc = MarchingCubes::new(
        (hdr.nx as usize, hdr.ny as usize, hdr.nz as usize),
        (hdr.cell[0], hdr.cell[1], hdr.cell[2]),
        (hdr.mx as f32, hdr.my as f32, hdr.mz as f32),
        values,  // Vec<f32>
        iso_level,
    );


    // Use MeshSide::Inside or MeshSide::Outside as required.
    let mesh = mc.generate(MeshSide::Both);
    
    // Example of converting the generic output mesh to a graphic engine's:
    let vertices = mesh
        .vertices
        .iter()
        .map(|v| graphics::Vertex::new(v.posit.to_arr(), v.normal))
        .collect();

    scene.meshes[MESH_DENSITY_SURFACE] = graphics::Mesh {
        vertices,
        indices: mesh.indices,
        material: 0,
    };
}
```

Why another Marching Cubes library? I couldn't figure out how to use the existing ones.


[//]: # (todo: Screenshot[s])
