use nalgebra::{Matrix3, Point3, RealField, Vector3};
use num_traits::Float;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::marker::PhantomData;

pub use nalgebra;
pub use num_traits;

pub trait Element: Copy + Debug + Eq {
    const N_NODES: usize;
    fn node(&self, i: usize) -> usize;
    fn node_mut(&mut self, i: usize) -> &mut usize;
    fn nodes(&self) -> NodeIter<'_, Self> {
        NodeIter::new(self)
    }
    fn nodes_mut(&mut self) -> NodeIterMut<'_, Self> {
        NodeIterMut::new(self)
    }
}

pub trait BoundedElement: Element {
    const N_SIDES: usize;
    type Side: Element;
    fn side(&self, i: usize) -> Self::Side;
    fn sides(&self) -> SideIter<'_, Self> {
        SideIter::new(self)
    }
}

/// Mesh representing an exodus data model.
///
/// TODO: multiple timesteps, multiple blocks, node sets, node vars.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Mesh<E: BoundedElement, S: RealField> {
    name: Option<Name>,
    coords: Vec<Point3<S>>,
    elems: Vec<E>,
    elem_var_names: Vec<Name>,
    elem_var_values: Vec<Vec<S>>,
    side_sets: Vec<SideSet<E>>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SideSet<E: BoundedElement> {
    _marker: PhantomData<E::Side>,
    name: Option<Name>,
    sides: HashSet<(usize, usize)>,
}

impl<E: BoundedElement> SideSet<E> {
    pub fn new() -> Self {
        Self {
            _marker: PhantomData,
            name: Default::default(),
            sides: Default::default(),
        }
    }

    pub fn name(&self) -> Option<&str> {
        self.name.as_ref().map(|n| n.as_ref())
    }

    pub fn set_name(&mut self, name: &str) {
        self.name = Some(Name::from(name));
    }

    pub fn add_side(&mut self, elem: usize, side: usize) -> bool {
        self.sides.insert((elem, side))
    }

    pub fn remove_side(&mut self, elem: usize, side: usize) -> bool {
        self.sides.remove(&(elem, side))
    }

    pub fn sides(&self) -> &HashSet<(usize, usize)> {
        &self.sides
    }
}

impl<E: BoundedElement, S: Float + RealField> Mesh<E, S> {
    pub fn new() -> Self {
        Self {
            name: Default::default(),
            coords: Default::default(),
            elems: Default::default(),
            elem_var_names: Default::default(),
            elem_var_values: Default::default(),
            side_sets: Default::default(),
        }
    }

    pub fn name(&self) -> Option<&str> {
        self.name.as_ref().map(|n| n.as_ref())
    }

    pub fn set_name(&mut self, name: &str) {
        self.name = Some(Name::from(name));
    }

    pub fn add_vertex(&mut self, point: Point3<S>) -> usize {
        let i = self.coords.len();
        self.coords.push(point);
        i
    }

    pub fn vertex(&self, i: usize) -> &Point3<S> {
        &self.coords[i]
    }

    pub fn vertex_mut(&mut self, i: usize) -> &mut Point3<S> {
        &mut self.coords[i]
    }

    pub fn vertices(&self) -> &[Point3<S>] {
        &self.coords
    }

    pub fn vertices_mut(&mut self) -> &mut [Point3<S>] {
        &mut self.coords
    }

    pub fn add_elem_var(&mut self, name: &str) -> usize {
        let i = self.elem_var_names.len();
        self.elem_var_names.push(Name::from(name));
        let zero = S::from_f64(0.0).unwrap();
        let values = vec![zero; self.elems.len()];
        self.elem_var_values.push(values);
        i
    }

    pub fn elem_vars(&self) -> &[Name] {
        &self.elem_var_names
    }

    pub fn add_elem(&mut self, elem: E, vars: &[S]) {
        debug_assert_eq!(vars.len(), self.elem_var_names.len());
        self.elems.push(elem);
        for (var, val) in self.elem_var_values.iter_mut().zip(vars) {
            var.push(*val);
        }
    }

    pub fn elem(&self, i: usize) -> &E {
        &self.elems[i]
    }

    pub fn elem_var(&self, i: usize) -> &[S] {
        &self.elem_var_values[i]
    }

    pub fn elem_mut(&mut self, i: usize) -> &mut E {
        &mut self.elems[i]
    }

    pub fn elem_var_mut(&mut self, i: usize) -> &mut [S] {
        &mut self.elem_var_values[i]
    }

    pub fn elems(&self) -> &[E] {
        &self.elems
    }

    pub fn elems_mut(&mut self) -> &mut [E] {
        &mut self.elems
    }

    pub fn elem_swap_remove(&mut self, i: usize) {
        self.elems.swap_remove(i);
        for var in &mut self.elem_var_values {
            var.swap_remove(i);
        }
    }

    pub fn add_side_set(&mut self, side: SideSet<E>) -> usize {
        let i = self.side_sets.len();
        self.side_sets.push(side);
        i
    }

    pub fn side_set(&self, i: usize) -> &SideSet<E> {
        &self.side_sets[i]
    }

    pub fn side_set_mut(&mut self, i: usize) -> &mut SideSet<E> {
        &mut self.side_sets[i]
    }

    pub fn side_sets(&self) -> &[SideSet<E>] {
        &self.side_sets
    }

    pub fn side_sets_mut(&mut self) -> &mut [SideSet<E>] {
        &mut self.side_sets
    }

    pub fn diff(&self, bar2: &Bar2) -> Vector3<S> {
        let a = self.coords[bar2.node(0)];
        let b = self.coords[bar2.node(1)];
        b - a
    }

    pub fn length(&self, bar2: &Bar2) -> S {
        self.diff(bar2).magnitude()
    }

    pub fn normal(&self, tri3: &Tri3) -> Vector3<S> {
        let a = self.coords[tri3.node(0)];
        let b = self.coords[tri3.node(1)];
        let c = self.coords[tri3.node(2)];
        (b - a).cross(&(c - a))
    }

    pub fn volume(&self, tet4: &Tet4) -> S {
        let a = self.coords[tet4.node(0)];
        let b = self.coords[tet4.node(1)];
        let c = self.coords[tet4.node(2)];
        let d = self.coords[tet4.node(3)];
        Matrix3::from_columns(&[a - d, b - d, c - d]).determinant()
            / S::from_f64(6.0).unwrap()
    }

    pub fn aspect_ratio(&self, tet4: &Tet4) -> S {
        let mut max_length2 = S::from_f64(0.0).unwrap();
        let mut sum_normals = S::from_f64(0.0).unwrap();
        let ab = self.diff(&tet4.edge(0));
        let ac = self.diff(&tet4.edge(1));
        let ad = self.diff(&tet4.edge(2));
        let alpha = Float::abs(ab.dot(&ac.cross(&ad)));
        for tri in tet4.sides() {
            sum_normals += self.normal(&tri).magnitude();
        }
        for edge in tet4.edges() {
            max_length2 = Float::max(max_length2, self.diff(&edge).magnitude_squared());
        }
        let max_length = Float::sqrt(max_length2);
        let radius = alpha / sum_normals;
        S::from_f64(1.0).unwrap() / (max_length / radius / S::from_f64(2.0 * 6.0.sqrt()).unwrap())
    }

    pub fn compact(&mut self) {
        let mut remap = vec![None; self.coords.len()];
        let mut nv = 0;
        for elem in self.elems_mut() {
            for node in elem.nodes_mut() {
                let i = *node as usize;
                if remap[i].is_none() {
                    remap[i] = Some(nv);
                    nv += 1;
                }
                *node = remap[i].unwrap();
            }
        }
        let zero = S::from_f64(0.0).unwrap();
        let mut coords = vec![Point3::new(zero, zero, zero); nv as usize];
        for (i, v) in self.coords.iter().enumerate() {
            if let Some(Some(j)) = remap.get(i) {
                coords[*j as usize] = *v;
            }
        }
        self.coords = coords;
    }

    pub fn boundary(&self) -> SideSet<E> {
        let mut side_set = HashMap::new();
        for (e, elem) in self.elems().iter().enumerate() {
            for (s, side) in elem.sides().enumerate() {
                let mut nodes = side.nodes().collect::<Vec<_>>();
                nodes.sort();
                if side_set.remove(&nodes).is_none() {
                    side_set.insert(nodes, (e, s));
                }
            }
        }
        let mut sides = SideSet::new();
        for (_, (e, s)) in side_set {
            sides.add_side(e, s);
        }
        sides
    }

    pub fn to_boundary(&self) -> Mesh<E::Side, S>
    where
        E::Side: BoundedElement,
    {
        let zero = S::from_f64(0.0).unwrap();
        let boundary = self.boundary();
        let mut mesh = Mesh {
            name: boundary.name,
            coords: self.coords.clone(),
            elems: Vec::with_capacity(boundary.sides().len()),
            elem_var_names: self.elem_var_names.clone(),
            elem_var_values: vec![Vec::with_capacity(boundary.sides().len()); self.elem_var_values.len()],
            side_sets: Default::default(),
        };
        let mut values = vec![zero; self.elem_var_values.len()];
        for (e, s) in boundary.sides() {
            for (i, var) in self.elem_var_values.iter().enumerate() {
                values[i] = var[*e];
            }
            mesh.add_elem(self.elem(*e).side(*s), &values);
        }
        mesh.compact();
        mesh
    }
}

pub struct NodeIter<'a, E> {
    elem: &'a E,
    i: usize,
}

impl<'a, E: Element> NodeIter<'a, E> {
    pub fn new(elem: &'a E) -> Self {
        Self { elem, i: 0 }
    }
}

impl<'a, E: Element> Iterator for NodeIter<'a, E> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < E::N_NODES {
            let res = Some(self.elem.node(self.i));
            self.i += 1;
            res
        } else {
            None
        }
    }
}

pub struct NodeIterMut<'a, E> {
    elem: &'a mut E,
    i: usize,
}

impl<'a, E: Element> NodeIterMut<'a, E> {
    pub fn new(elem: &'a mut E) -> Self {
        Self { elem, i: 0 }
    }
}

impl<'a, E: Element> Iterator for NodeIterMut<'a, E> {
    type Item = &'a mut usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < E::N_NODES {
            let i = self.i;
            self.i += 1;
            Some(unsafe { &mut *(self.elem.node_mut(i) as *mut _) })
        } else {
            None
        }
    }
}

pub struct SideIter<'a, E> {
    elem: &'a E,
    i: usize,
}

impl<'a, E: BoundedElement> SideIter<'a, E> {
    pub fn new(elem: &'a E) -> Self {
        Self { elem, i: 0 }
    }
}

impl<'a, E: BoundedElement> Iterator for SideIter<'a, E> {
    type Item = E::Side;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < E::N_SIDES {
            let res = Some(self.elem.side(self.i));
            self.i += 1;
            res
        } else {
            None
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Node1([usize; 1]);

impl Node1 {
    pub const fn new(i: [usize; 1]) -> Self {
        Self(i)
    }
}

impl Element for Node1 {
    const N_NODES: usize = 1;

    fn node(&self, i: usize) -> usize {
        self.0[i]
    }

    fn node_mut(&mut self, i: usize) -> &mut usize {
        &mut self.0[i]
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Bar2([usize; 2]);

impl Bar2 {
    pub const fn new(i: [usize; 2]) -> Self {
        Self(i)
    }
}

impl Element for Bar2 {
    const N_NODES: usize = 2;

    fn node(&self, i: usize) -> usize {
        self.0[i]
    }

    fn node_mut(&mut self, i: usize) -> &mut usize {
        &mut self.0[i]
    }
}

impl BoundedElement for Bar2 {
    const N_SIDES: usize = 2;

    type Side = Node1;

    fn side(&self, i: usize) -> Self::Side {
        Node1([self.0[i]])
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Tri3([usize; 3]);

impl Tri3 {
    pub const fn new(i: [usize; 3]) -> Self {
        Self(i)
    }
}

impl Element for Tri3 {
    const N_NODES: usize = 3;

    fn node(&self, i: usize) -> usize {
        self.0[i]
    }

    fn node_mut(&mut self, i: usize) -> &mut usize {
        &mut self.0[i]
    }
}

impl BoundedElement for Tri3 {
    const N_SIDES: usize = 3;

    type Side = Bar2;

    fn side(&self, i: usize) -> Self::Side {
        let indices = [[0, 1], [1, 2], [2, 0]][i];
        Bar2([self.node(indices[0]), self.node(indices[1])])
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Tet4([usize; 4]);

impl Tet4 {
    pub const fn new(i: [usize; 4]) -> Self {
        Self(i)
    }

    pub fn edge(&self, i: usize) -> Bar2 {
        let indices = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]][i];
        Bar2([self.node(indices[0]), self.node(indices[1])])
    }

    pub fn edges(&self) -> impl Iterator<Item = Bar2> + '_ {
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
            .iter()
            .map(move |indices| Bar2([self.node(indices[0]), self.node(indices[1])]))
    }

    pub fn flip(&mut self) {
        let a = self.0[0];
        self.0[0] = self.0[1];
        self.0[1] = a;
    }
}

impl Element for Tet4 {
    const N_NODES: usize = 4;

    fn node(&self, i: usize) -> usize {
        self.0[i]
    }

    fn node_mut(&mut self, i: usize) -> &mut usize {
        &mut self.0[i]
    }
}

impl BoundedElement for Tet4 {
    const N_SIDES: usize = 4;

    type Side = Tri3;

    fn side(&self, i: usize) -> Self::Side {
        let indices = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]][i];
        Tri3([
            self.node(indices[0]),
            self.node(indices[1]),
            self.node(indices[2]),
        ])
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Name([u8; 32]);

impl<'a> From<&'a str> for Name {
    fn from(s: &'a str) -> Self {
        let utf8: &[u8] = s.as_ref();
        let slice_len = if utf8.len() > 31 { 31 } else { utf8.len() };
        let mut bytes = [0u8; 32];
        bytes[..slice_len].copy_from_slice(&utf8[..slice_len]);
        Self(bytes)
    }
}

impl AsRef<str> for Name {
    fn as_ref(&self) -> &str {
        std::str::from_utf8(&self.0).unwrap()
    }
}

impl AsRef<[u8]> for Name {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quad_mesh_boundary() {
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(0.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 0.0));
        mesh.add_vertex(Point3::new(0.0, 1.0, 0.0));
        mesh.add_elem(Tri3([0, 1, 2]), &[]);
        mesh.add_elem(Tri3([0, 2, 3]), &[]);
        let side = mesh.boundary();
        println!("{:?}", side);
        assert_eq!(side.sides().len(), 4);
        assert!(mesh.normal(mesh.elem(0)).magnitude() > 0.0);
        assert!(mesh.normal(mesh.elem(1)).magnitude() > 0.0);
    }

    #[test]
    fn test_bcc_tet4_boundary() {
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(0.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(1.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(1.0, 0.0, 1.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 1.0));
        let tet4 = Tet4([0, 1, 2, 3]);
        mesh.add_elem(tet4, &[]);
        let side = mesh.boundary();
        println!("{:?}", side);
        assert_eq!(side.sides().len(), 4);

        let sum_normals = mesh.normal(&tet4.side(0))
            + mesh.normal(&tet4.side(1))
            + mesh.normal(&tet4.side(2))
            + mesh.normal(&tet4.side(3));

        assert_eq!(sum_normals, Vector3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_tile_volume() {
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(-0.5, -0.5, -0.5));
        mesh.add_vertex(Point3::new(0.5, -0.5, -0.5));
        mesh.add_vertex(Point3::new(0.5, 0.5, -0.5));
        mesh.add_vertex(Point3::new(-0.5, 0.5, -0.5));
        mesh.add_vertex(Point3::new(-0.5, -0.5, 0.5));
        mesh.add_vertex(Point3::new(0.5, -0.5, 0.5));
        mesh.add_vertex(Point3::new(0.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(-0.5, 0.5, 0.5));

        mesh.add_elem(Tet4::new([0, 1, 5, 2]), &[]);
        mesh.add_elem(Tet4::new([0, 2, 5, 7]), &[]);
        mesh.add_elem(Tet4::new([0, 7, 5, 4]), &[]);
        mesh.add_elem(Tet4::new([2, 6, 5, 7]), &[]);
        mesh.add_elem(Tet4::new([0, 2, 7, 3]), &[]);

        let side = mesh.boundary();
        println!("{:?}", side);
        assert_eq!(side.sides().len(), 12);

        let sum_volumes = mesh.volume(&mesh.elem(0))
            + mesh.volume(&mesh.elem(1))
            + mesh.volume(&mesh.elem(2))
            + mesh.volume(&mesh.elem(3))
            + mesh.volume(&mesh.elem(4));

        assert!(sum_volumes < 1.0);
        assert!(sum_volumes > 0.999999);
    }

    #[test]
    fn test_compact_mesh() {
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(0.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(0.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(0.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(0.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.0, 0.0, 1.0));
        mesh.add_vertex(Point3::new(0.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 1.0));
        mesh.add_vertex(Point3::new(0.0, 0.0, 0.0));
        let tet4 = Tet4([1, 3, 5, 7]);
        mesh.add_elem(tet4, &[]);
        mesh.compact();
        assert_eq!(mesh.vertices().len(), 4);
    }

    #[test]
    fn test_vars() {
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(0.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 0.0));
        mesh.add_vertex(Point3::new(0.0, 1.0, 0.0));
        mesh.add_elem_var("color");
        mesh.add_elem(Tri3([0, 1, 2]), &[1.0]);
        mesh.add_elem_var("uv_x");
        mesh.add_elem(Tri3([0, 2, 3]), &[2.0, 2.0]);
        assert_eq!(mesh.elem_var(0)[0], 1.0);
        assert_eq!(mesh.elem_var(1)[0], 0.0);
        assert_eq!(mesh.elem_var(0)[1], 2.0);
        assert_eq!(mesh.elem_var(1)[1], 2.0);
        mesh.elem_swap_remove(0);
        assert_eq!(mesh.elem_var(0)[0], 2.0);
        assert_eq!(mesh.elem_var(1)[0], 2.0);
    }

    #[test]
    fn test_aspect_ratio() {
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(1.0, 1.0, 1.0));
        mesh.add_vertex(Point3::new(-1.0, -1.0, 1.0));
        mesh.add_vertex(Point3::new(-1.0, 1.0, -1.0));
        mesh.add_vertex(Point3::new(1.0, -1.0, -1.0));
        mesh.add_elem(Tet4::new([0, 1, 3, 2]), &[]);
        assert_eq!(mesh.aspect_ratio(mesh.elem(0)), 1.0);
    }
}
