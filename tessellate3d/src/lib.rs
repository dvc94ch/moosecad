pub use bbox::BoundingBox;
use mesh::{Element, Mesh, Tet4};
use nalgebra::{Point3, RealField, Vector3};
use num_traits::Float;
use std::collections::HashMap;

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
    resolution: S,
    warp_threshold: S,
    dim: [usize; 3],
    mesh: Mesh<Tet4, S>,
    phi: Vec<S>,
    cut_map: HashMap<(usize, usize), usize>,
}

impl<'a, S: Float + RealField + alga::general::RealField + From<f64>> IsosurfaceStuffing<'a, S> {
    pub fn new(function: &'a dyn ImplicitFunction<S>, resolution: S) -> Self {
        let bbox = function.bbox();
        Self {
            function,
            resolution,
            warp_threshold: S::from_f64(0.3).unwrap(),
            dim: [
                Float::ceil(bbox.dim()[0] / resolution).to_usize().unwrap(),
                Float::ceil(bbox.dim()[1] / resolution).to_usize().unwrap(),
                Float::ceil(bbox.dim()[2] / resolution).to_usize().unwrap(),
            ],
            mesh: Mesh::new(),
            phi: Default::default(),
            cut_map: Default::default(),
        }
    }

    /// Cut a chunk of the lattice out to roughly match the shape.
    fn cut_lattice(&mut self) {
        let bbox = self.function.bbox();
        let origin = Vector3::new(bbox.min[0], bbox.min[1], bbox.min[2]);
        let mut lattice_map = HashMap::<Vector3<usize>, usize>::new();
        let (mut fx, mut fy, mut fz) = (true, true, true);
        fn transform(coord: [u32; 3], flipped: [bool; 3]) -> Vector3<usize> {
            let x = if flipped[0] { 1 - coord[0] } else { coord[0] };
            let y = if flipped[1] { 1 - coord[1] } else { coord[1] };
            let z = if flipped[2] { 1 - coord[2] } else { coord[2] };
            Vector3::new(x as usize, y as usize, z as usize)
        }
        for x in 0..self.dim[0] {
            fx = !fx;
            for y in 0..self.dim[1] {
                fy = !fy;
                for z in 0..self.dim[2] {
                    fz = !fz;
                    let flipped = [fx, fy, fz];
                    let offset = Vector3::new(x, y, z);
                    for t in Tile::TETS.iter() {
                        let coordsi = [
                            transform(Tile::LATTICE[t.node(0)], flipped) + offset,
                            transform(Tile::LATTICE[t.node(1)], flipped) + offset,
                            transform(Tile::LATTICE[t.node(2)], flipped) + offset,
                            transform(Tile::LATTICE[t.node(3)], flipped) + offset,
                        ];
                        let coordss = [
                            Point3::from(
                                coordsi[0].map(|i| S::from_f64(i as f64).unwrap())
                                    * self.resolution
                                    + origin,
                            ),
                            Point3::from(
                                coordsi[1].map(|i| S::from_f64(i as f64).unwrap())
                                    * self.resolution
                                    + origin,
                            ),
                            Point3::from(
                                coordsi[2].map(|i| S::from_f64(i as f64).unwrap())
                                    * self.resolution
                                    + origin,
                            ),
                            Point3::from(
                                coordsi[3].map(|i| S::from_f64(i as f64).unwrap())
                                    * self.resolution
                                    + origin,
                            ),
                        ];
                        let mut intersects_f = false;
                        for c in coordss.iter() {
                            if self.function.value(c) <= S::from_f64(0.0).unwrap() {
                                intersects_f = true;
                                break;
                            }
                        }

                        if intersects_f {
                            let mut tet = [0; 4];
                            for (i, c) in coordsi.iter().enumerate() {
                                if let Some(n) = lattice_map.get(c) {
                                    tet[i] = *n;
                                } else {
                                    tet[i] = self.mesh.vertices().len();
                                    lattice_map.insert(*c, tet[i]);
                                    self.mesh.add_vertex(coordss[i]);
                                    self.phi.push(self.function.value(&coordss[i]));
                                }
                            }
                            let mut tet = Tet4::new(tet);
                            if fx ^ fy ^ fz {
                                tet.flip();
                            }
                            self.mesh.add_elem(tet, &[]);
                        }
                    }
                }
            }
        }
    }

    /// Fix some edge crossings by warping vertices if it's admissible.
    fn warp_vertices(&mut self) {
        //assert!(self.warp_threshold >= 0.0 && self.warp_threshold <= 0.5);
        let zero = S::from_f64(0.0).unwrap();
        let one = S::from_f64(1.0).unwrap();
        let size = self.mesh.vertices().len();
        let mut warp = vec![None; size];
        //let mut warp_nbr = vec![None; size];
        let mut delta = vec![None; size];
        for tet in self.mesh.elems() {
            for edge in tet.edges() {
                let n0 = edge.node(0);
                let n1 = edge.node(1);
                let phi0 = self.phi[n0];
                let phi1 = self.phi[n1];
                let x0 = self.mesh.vertex(n0);
                let x1 = self.mesh.vertex(n1);
                if phi0 > zero && phi1 > zero || phi0 < zero && phi1 < zero {
                    continue;
                }
                let alpha = phi0 / (phi0 - phi1);
                if alpha < self.warp_threshold {
                    let d = alpha * nalgebra::distance(x0, x1);
                    if warp[n0].is_none() || d < warp[n0].unwrap() {
                        warp[n0] = Some(d);
                        //warp_nbr[n0] = Some(n1);
                        delta[n0] = Some((x1 - x0) * alpha);
                    }
                } else if alpha > one - self.warp_threshold {
                    let d = (one - alpha) * nalgebra::distance(x0, x1);
                    if warp[n1].is_none() || d < warp[n1].unwrap() {
                        warp[n1] = Some(d);
                        //warp_nbr[n1] = Some(n0);
                        delta[n1] = Some((x0 - x1) * (one - alpha));
                    }
                }
            }
        }
        for (i, v) in self.mesh.vertices_mut().iter_mut().enumerate() {
            if let Some(d) = delta[i] {
                *v += d;
                self.phi[i] = zero;
            }
        }
    }

    /// Interpolates the point where two coordinates intersect the boundary of the function.
    fn cut_edge(&mut self, a: usize, b: usize) -> usize {
        let zero = S::from_f64(0.0).unwrap();
        let one = S::from_f64(1.0).unwrap();
        debug_assert!(
            self.phi[a] > zero && self.phi[b] < zero || self.phi[a] < zero && self.phi[b] > zero
        );
        let edge = (usize::min(a, b), usize::max(a, b));
        if let Some(i) = self.cut_map.get(&edge) {
            return *i;
        }
        let alpha = self.phi[a] / self.phi[a] - self.phi[b];
        let vertex =
            self.mesh.vertex(a).coords * (one - alpha) + self.mesh.vertex(b).coords * alpha;
        let index = self.mesh.vertices().len();
        self.mesh.add_vertex(Point3::from(vertex));
        self.phi.push(zero);
        self.cut_map.insert(edge, index);
        index
    }

    /// Find all tets with vertices where a corner has vphi > 0 and trim them back.
    fn trim_spikes(&mut self) {
        let zero = S::from_f64(0.0).unwrap();
        let mut new_tets = Vec::new();
        let mut ti = 0;
        let mut stats = [0; 7];
        println!("num tets {}", self.mesh.elems().len());
        while ti < self.mesh.elems().len() {
            let mut a = self.mesh.elem(ti).node(0);
            let mut b = self.mesh.elem(ti).node(1);
            let mut c = self.mesh.elem(ti).node(2);
            let mut d = self.mesh.elem(ti).node(3);
            if self.phi[a] <= zero
                && self.phi[b] <= zero
                && self.phi[c] <= zero
                && self.phi[d] <= zero
            {
                ti += 1;
                continue;
            }

            let mut flipped = false;
            if self.phi[a] < self.phi[b] || (self.phi[a] == self.phi[b] && a < b) {
                std::mem::swap(&mut a, &mut b);
                flipped = !flipped;
            }
            if self.phi[c] < self.phi[d] || (self.phi[c] == self.phi[d] && c < d) {
                std::mem::swap(&mut c, &mut d);
                flipped = !flipped;
            }
            if self.phi[a] < self.phi[c] || (self.phi[a] == self.phi[c] && a < c) {
                std::mem::swap(&mut a, &mut c);
                flipped = !flipped;
            }
            if self.phi[b] < self.phi[d] || (self.phi[b] == self.phi[d] && b < d) {
                std::mem::swap(&mut b, &mut d);
                flipped = !flipped;
            }
            if self.phi[b] < self.phi[c] || (self.phi[b] == self.phi[c] && b < c) {
                std::mem::swap(&mut b, &mut c);
                flipped = !flipped;
            }

            debug_assert!(
                self.phi[a] >= self.phi[b]
                    && self.phi[b] >= self.phi[c]
                    && self.phi[c] >= self.phi[d]
            );
            debug_assert!(self.phi[a] > zero);
            debug_assert!(self.phi[d] <= zero);

            if self.phi[d] == zero {
                // tet outside of the volume
                self.mesh.elem_swap_remove(ti);
                stats[0] += 1;
            } else if self.phi[c] > zero {
                // only d inside the volume
                self.mesh.elem_swap_remove(ti);
                stats[1] += 1;

                let ad = self.cut_edge(a, d);
                let bd = self.cut_edge(b, d);
                let cd = self.cut_edge(c, d);
                new_tets.push((Tet4::new([ad, bd, cd, d]), !flipped, 6));
            } else if self.phi[b] < zero {
                // only a outside of the volume
                self.mesh.elem_swap_remove(ti);
                stats[2] += 1;

                let ab = self.cut_edge(a, b);
                let ac = self.cut_edge(a, c);
                let ad = self.cut_edge(a, d);
                new_tets.push((Tet4::new([ab, b, c, d]), !flipped, 6));
                new_tets.push((Tet4::new([c, ab, ac, d]), flipped, 6));
                new_tets.push((Tet4::new([ab, d, ac, ad]), !flipped, 6));
            } else if self.phi[b] > zero && self.phi[c] < zero {
                // 1 out, 2 in
                self.mesh.elem_swap_remove(ti);
                stats[3] += 1;

                let ac = self.cut_edge(a, c);
                let ad = self.cut_edge(a, d);
                let bc = self.cut_edge(b, c);
                let bd = self.cut_edge(b, d);
                new_tets.push((Tet4::new([ac, bd, c, d]), flipped, 6));
                new_tets.push((Tet4::new([ac, bd, bc, d]), flipped, 6));
                new_tets.push((Tet4::new([ac, ad, bd, d]), flipped, 6));
            } else if self.phi[b] == zero && self.phi[c] == zero {
                // 1 out, 1 in, 2 surface
                self.mesh.elem_swap_remove(ti);
                stats[4] += 1;

                let ad = self.cut_edge(a, d);
                new_tets.push((Tet4::new([ad, b, c, d]), !flipped, 6));
            } else if self.phi[b] == zero {
                // 1 out, 2 in, 1 surface
                self.mesh.elem_swap_remove(ti);
                stats[5] += 1;

                let ac = self.cut_edge(a, c);
                let ad = self.cut_edge(a, d);
                new_tets.push((Tet4::new([b, ac, c, d]), flipped, 6));
                new_tets.push((Tet4::new([ad, b, ac, d]), !flipped, 6));
            } else {
                // two out, one in, one surface
                self.mesh.elem_swap_remove(ti);
                stats[6] += 1;

                let ad = self.cut_edge(a, d);
                let bd = self.cut_edge(b, d);
                new_tets.push((Tet4::new([ad, bd, c, d]), flipped, 6));
            }
        }

        println!("stats {:?}", stats);
        self.mesh.add_elem_var("tet_type");
        for (mut tet, flipped, ty) in new_tets.drain(..) {
            if flipped {
                let a = tet.node(0);
                *tet.node_mut(0) = tet.node(1);
                *tet.node_mut(1) = a;
            }
            self.mesh
                .add_elem(tet, &[S::from_f64(ty as f64 / 6.0).unwrap()]);
        }
        println!("num tets {}", self.mesh.elems().len());
    }

    /// Remove tets with corners where phi = 0 but are outside the volume.
    fn remove_exterior_tets(&mut self) {
        let zero = S::from_f64(0.0).unwrap();
        let mut i = 0;
        while i < self.mesh.elems().len() {
            let a = self.mesh.elem(i).node(0);
            let b = self.mesh.elem(i).node(1);
            let c = self.mesh.elem(i).node(2);
            let d = self.mesh.elem(i).node(3);
            let phia = self.phi[a];
            let phib = self.phi[b];
            let phic = self.phi[c];
            let phid = self.phi[d];
            debug_assert!(phia <= zero && phib <= zero && phic <= zero && phid <= zero);
            if phia == zero && phib == zero && phic == zero && phid == zero {
                let p = self.mesh.vertex(a).coords
                    + self.mesh.vertex(b).coords
                    + self.mesh.vertex(c).coords
                    + self.mesh.vertex(d).coords;
                if self
                    .function
                    .value(&Point3::from(p / S::from_f64(4.0).unwrap()))
                    > zero
                {
                    self.mesh.elem_swap_remove(i);
                }
            }
            i += 1;
        }
    }

    fn check_for_inverted_tets(&self) {
        println!("checking for inverted tets");
        let zero = S::from_f64(0.0).unwrap();
        let one = S::from_f64(1.0).unwrap();
        let mut aspect_ratio = (one, zero);
        for (i, tet) in self.mesh.elems().iter().enumerate() {
            let v = self.mesh.volume(tet);
            let ar = self.mesh.aspect_ratio(tet);
            if v <= zero {
                println!(
                    "Tet #{} is inverted! volume = {} aspect = {}",
                    i, v, ar,
                );
            }
            aspect_ratio.0 = Float::min(aspect_ratio.0, ar);
            aspect_ratio.1 = Float::max(aspect_ratio.1, ar);
        }
        println!("aspect_ratio: min = {} max = {}", aspect_ratio.0, aspect_ratio.1);
    }

    pub fn tessellate(mut self) -> Mesh<Tet4, S> {
        self.cut_lattice();
        println!("num tets {}", self.mesh.elems().len());
        self.warp_vertices();
        println!("num tets {}", self.mesh.elems().len());
        self.trim_spikes();
        self.remove_exterior_tets();
        self.mesh.compact();
        self.check_for_inverted_tets();
        self.mesh
    }
}

pub struct Tile;

impl Tile {
    pub const LATTICE: [[u32; 3]; 8] = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ];
    pub const TETS: [Tet4; 5] = [
        Tet4::new([0, 1, 5, 2]),
        Tet4::new([0, 2, 5, 7]),
        Tet4::new([0, 7, 5, 4]),
        Tet4::new([2, 6, 5, 7]),
        Tet4::new([0, 2, 7, 3]),
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
            p.coords.norm() - 1.0
        }

        fn normal(&self, p: &Point3<f64>) -> Vector3<f64> {
            p.coords.normalize()
        }
    }

    #[test]
    fn test_convert_mesh() {
        let f = UnitSphere::new();
        IsosurfaceStuffing::new(&f, 0.1).tessellate();
    }
}
