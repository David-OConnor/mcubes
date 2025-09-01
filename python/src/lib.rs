use lin_alg::f32::Vec3;
use mcubes_rs::{
    MarchingCubes as MarchingCubesRs, Mesh as MeshRs, MeshSide as MeshSideRs, Vertex as VertexRs,
};
use pyo3::{exceptions::PyValueError, prelude::*};

#[pyclass]
pub struct MeshSide {
    value: u8,
}

// todo: Replace this workflow with how you handle it in Na-seq? Or vice-versa.
#[pymethods]
impl MeshSide {
    #[classattr]
    pub const Both: Self = Self { value: 0 };
    #[classattr]
    pub const OutsideOnly: Self = Self { value: 1 };
    #[classattr]
    pub const InsideOnly: Self = Self { value: 2 };

    fn __repr__(&self) -> String {
        format!("Mesh Side: {:?}", self.value)
    }

    // fn __repr__(&self) -> String {
    //     match self.value {
    //         0 => "MeshSide.Both",
    //         1 => "MeshSide.OutsideOnly",
    //         2 => "MeshSide.InsideOnly",
    //         _ => "MeshSide<?>",
    //     }
    //         .to_string()
    // }
}

impl MeshSide {
    fn to_rust(&self) -> MeshSideRs {
        match self.value {
            0 => MeshSideRs::Both,
            1 => MeshSideRs::OutsideOnly,
            2 => MeshSideRs::InsideOnly,
            _ => MeshSideRs::Both,
        }
    }
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct Vertex {
    #[pyo3(get)]
    posit: (f32, f32, f32),
    #[pyo3(get)]
    normal: (f32, f32, f32),
}

#[pymethods]
impl Vertex {
    fn __repr__(&self) -> String {
        format!("{self:?}")
    }
}

impl From<&VertexRs> for Vertex {
    fn from(v: &VertexRs) -> Self {
        Vertex {
            posit: (v.posit.x, v.posit.y, v.posit.z),
            normal: (v.normal.x, v.normal.y, v.normal.z),
        }
    }
}

#[pyclass]
#[derive(Debug)]
pub struct Mesh {
    vertices: Vec<Vertex>,
    #[pyo3(get)]
    indices: Vec<usize>,
}

#[pymethods]
impl Mesh {
    #[getter]
    fn vertices(&self) -> Vec<Vertex> {
        self.vertices.clone()
    }

    fn __repr__(&self) -> String {
        format!("{self:?}",)
    }
}

impl From<MeshRs> for Mesh {
    fn from(m: MeshRs) -> Self {
        let vertices = m.vertices.iter().map(Vertex::from).collect();
        Mesh {
            vertices,
            indices: m.indices,
        }
    }
}

#[pyclass]
pub struct MarchingCubes {
    inner: MarchingCubesRs,
}

#[pymethods]
impl MarchingCubes {
    #[new]
    fn new(
        dims: (usize, usize, usize),
        size: (f32, f32, f32),
        sampling_interval: (f32, f32, f32),
        offset: (f32, f32, f32),
        values: Vec<f32>,
        iso_level: f32,
    ) -> PyResult<Self> {
        let offset_v = Vec3::new(offset.0, offset.1, offset.2);
        let inner =
            MarchingCubesRs::new(dims, size, sampling_interval, offset_v, values, iso_level)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn generate(&self, side: &MeshSide) -> Mesh {
        self.inner.generate(side.to_rust()).into()
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pymodule]
fn mcubes(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<MarchingCubes>()?;
    m.add_class::<MeshSide>()?;
    m.add_class::<Mesh>()?;
    m.add_class::<Vertex>()?;
    Ok(())
}
