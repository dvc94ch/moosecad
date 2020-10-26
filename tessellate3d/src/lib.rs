pub use bbox::BoundingBox;
use mesh::{Block, Element, Mesh, Tet4};
use nalgebra::{Point3, RealField, Vector3};
use num_traits::Float;

/// Trait to be implemented by functions that should be tesselated.
pub trait ImplicitFunction<S: RealField> {
    /// Return a `BoundingBox`, which is essential, so the algorithm knows where
    /// to search for surfaces.
    fn bbox(&self) -> &BoundingBox<S>;

    /// Evaluate the function on p and return the value. A value of zero signifies
    /// that p is on the surface to be tessellated. A negative value means p is
    /// inside the object. A positive value means p is outside the object.
    /// The magnitude of value must be continuous. Furthermore value has to be
    /// equal or greater than the euclidean distance between p and the surface.
    fn value(&self, p: &Point3<S>) -> S;

    /// Compute the normal of the function at p.
    fn normal(&self, p: &Point3<S>) -> Vector3<S>;
}

pub struct IsosurfaceStuffing<'a, S: RealField> {
    function: &'a dyn ImplicitFunction<S>,
    origin: Point3<S>,
    resolution: S,
    relative_error: S,
    dim: [usize; 3],
}

impl<'a, S: Float + RealField + alga::general::RealField + From<f64>> IsosurfaceStuffing<'a, S> {
    pub fn new(function: &'a dyn ImplicitFunction<S>, resolution: S, relative_error: S) -> Self {
        let bbox = function.bbox();
        Self {
            function,
            resolution,
            relative_error,
            origin: bbox.min,
            dim: [
                Float::ceil(bbox.dim()[0] / resolution).to_usize().unwrap(),
                Float::ceil(bbox.dim()[1] / resolution).to_usize().unwrap(),
                Float::ceil(bbox.dim()[2] / resolution).to_usize().unwrap(),
            ],
        }
    }

    pub fn tessellate(self) -> Mesh<Tet4, S> {
        let mut mesh = Mesh::new();
        let mut block = Block::new();
        //for x in 0..self.dim[0] {
        //    for y in 0..self.dim[1] {
        //        for z in 0..self.dim[2] {
        /*let offset = [
            self.origin[0] + (x as f64).into(),
            self.origin[1] + (y as f64).into(),
            self.origin[2] + (z as f64).into(),
        ];*/
        let offset: [S; 3] = [0.0.into(), 0.0.into(), 0.0.into()];
        let i_offset = mesh.vertices().len() as u32;
        //let i_offset = [0; 4]; TODO compact mesh
        for coord in Tile::COORDS.iter() {
            let pos = [
                offset[0] + coord[0].into(),
                offset[1] + coord[1].into(),
                offset[2] + coord[2].into(),
            ];
            mesh.add_vertex(Point3::new(pos[0], pos[1], pos[2]));
        }
        for tet in Tile::TETS.iter() {
            block.add_elem(Tet4::new([
                tet.node(0) + i_offset,
                tet.node(1) + i_offset,
                tet.node(2) + i_offset,
                tet.node(3) + i_offset,
            ]));
        }
        //        }
        //    }
        //}
        mesh.add_block(block);
        mesh
    }
}

pub struct Tile;

impl Tile {
    const COORDS: [[f64; 3]; 8] = [
        [-0.5, -0.5, -0.5],
        [0.5, -0.5, -0.5],
        [0.5, 0.5, -0.5],
        [-0.5, 0.5, -0.5],
        [-0.5, -0.5, 0.5],
        [0.5, -0.5, 0.5],
        [0.5, 0.5, 0.5],
        [-0.5, 0.5, 0.5],
    ];
    const TETS: [Tet4; 5] = [
        Tet4::new([0, 1, 2, 5]),
        Tet4::new([0, 2, 7, 5]),
        Tet4::new([0, 7, 5, 4]),
        Tet4::new([2, 6, 7, 5]),
        Tet4::new([0, 2, 3, 7]),
    ];
}

#[cfg(test)]
mod tests {
    use super::*;

    struct UnitSphere {
        bbox: BoundingBox<f64>,
    }

    impl UnitSphere {
        fn new() -> Self {
            Self {
                bbox: BoundingBox::new(&Point3::new(-1., -1., -1.), &Point3::new(1., 1., 1.)),
            }
        }
    }

    impl ImplicitFunction<f64> for UnitSphere {
        fn bbox(&self) -> &BoundingBox<f64> {
            &self.bbox
        }

        fn value(&self, p: &Point3<f64>) -> f64 {
            Vector3::new(p.x, p.y, p.z).norm() - 1.0
        }

        fn normal(&self, p: &Point3<f64>) -> Vector3<f64> {
            Vector3::new(p.x, p.y, p.z).normalize()
        }
    }

    #[test]
    fn test_convert_mesh() {
        let f = UnitSphere::new();
        let mesh = IsosurfaceStuffing::new(&f, 0.1, 0.2).tessellate();
    }
}
