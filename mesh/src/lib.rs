use nalgebra::{Matrix3, Point3, RealField, Vector3};
use num_traits::Float;
use std::collections::HashMap;
use std::fmt::Debug;

pub trait Element: Copy + Debug + Eq {
    const N_NODES: usize;
    fn node(&self, i: usize) -> u32;
    fn node_mut(&mut self, i: usize) -> &mut u32;
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

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Mesh<E: BoundedElement, S: RealField> {
    coords: Vec<Point3<S>>,
    blocks: Vec<Block<E>>,
    sides: Vec<Side<E>>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Block<E: BoundedElement> {
    elems: Vec<E>,
}

impl<E: BoundedElement> Block<E> {
    pub fn new() -> Self {
        Self {
            elems: Default::default(),
        }
    }

    pub fn add_elem(&mut self, elem: E) {
        self.elems.push(elem);
    }

    pub fn elem(&self, i: usize) -> &E {
        &self.elems[i]
    }

    pub fn elem_mut(&mut self, i: usize) -> &mut E {
        &mut self.elems[i]
    }

    pub fn elems(&self) -> &[E] {
        &self.elems
    }

    pub fn elems_mut(&mut self) -> &mut [E] {
        &mut self.elems
    }

    pub fn boundary(&self) -> Side<E> {
        let mut side_set = HashMap::new();
        for elem in &self.elems {
            for side in elem.sides() {
                let mut nodes = side.nodes().collect::<Vec<_>>();
                nodes.sort();
                if side_set.remove(&nodes).is_none() {
                    side_set.insert(nodes, side);
                }
            }
        }
        let mut sides = Side::new();
        for (_, side) in side_set {
            sides.add_side(side);
        }
        sides
    }

    pub fn swap_remove(&mut self, i: usize) {
        self.elems.swap_remove(i);
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Side<E: BoundedElement> {
    sides: Vec<E::Side>,
}

impl<E: BoundedElement> Side<E> {
    pub fn new() -> Self {
        Self {
            sides: Default::default(),
        }
    }

    pub fn add_side(&mut self, side: E::Side) {
        self.sides.push(side);
    }

    pub fn side(&self, i: usize) -> &E::Side {
        &self.sides[i]
    }

    pub fn side_mut(&mut self, i: usize) -> &mut E::Side {
        &mut self.sides[i]
    }

    pub fn sides(&self) -> &[E::Side] {
        &self.sides
    }

    pub fn sides_mut(&mut self) -> &mut [E::Side] {
        &mut self.sides
    }
}

impl<E: BoundedElement, S: Float + RealField> Mesh<E, S> {
    pub fn new() -> Self {
        Self {
            coords: Default::default(),
            blocks: Default::default(),
            sides: Default::default(),
        }
    }

    pub fn add_vertex(&mut self, point: Point3<S>) {
        self.coords.push(point);
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

    pub fn add_block(&mut self, block: Block<E>) {
        self.blocks.push(block);
    }

    pub fn block(&self, i: usize) -> &Block<E> {
        &self.blocks[i]
    }

    pub fn block_mut(&mut self, i: usize) -> &mut Block<E> {
        &mut self.blocks[i]
    }

    pub fn blocks(&self) -> &[Block<E>] {
        &self.blocks
    }

    pub fn blocks_mut(&mut self) -> &mut [Block<E>] {
        &mut self.blocks
    }

    pub fn add_side(&mut self, side: Side<E>) {
        self.sides.push(side);
    }

    pub fn side(&self, i: usize) -> &Side<E> {
        &self.sides[i]
    }

    pub fn side_mut(&mut self, i: usize) -> &mut Side<E> {
        &mut self.sides[i]
    }

    pub fn sides(&self) -> &[Side<E>] {
        &self.sides
    }

    pub fn sides_mut(&mut self) -> &mut [Side<E>] {
        &mut self.sides
    }

    pub fn length(&self, bar2: &Bar2) -> S {
        let a = self.coords[bar2.node(0) as usize];
        let b = self.coords[bar2.node(1) as usize];
        (b - a).magnitude()
    }

    pub fn normal(&self, tri3: &Tri3) -> Vector3<S> {
        let a = self.coords[tri3.node(0) as usize];
        let b = self.coords[tri3.node(1) as usize];
        let c = self.coords[tri3.node(2) as usize];
        (b - a).cross(&(c - a))
    }

    pub fn volume(&self, tet4: &Tet4) -> S {
        let a = self.coords[tet4.node(0) as usize];
        let b = self.coords[tet4.node(1) as usize];
        let c = self.coords[tet4.node(2) as usize];
        let d = self.coords[tet4.node(3) as usize];
        Float::abs(Matrix3::from_columns(&[a - d, b - d, c - d]).determinant())
            / S::from_f64(6.0).unwrap()
    }

    /*pub fn aspect_ratio(&self, tet4: &Tet4) -> S {
        let mut max_length = 0.0;
        let mut sum_normals = 0.0;
        for tri3 in tet4.sides() {
            sum_normals += self.normal(tri3).length();
            for bar2 in tri3.sides() {
                max_length = f32::max(max_length, self.length(bar2));
            }
        }
    }*/

    pub fn compact(&mut self) {
        let mut remap = vec![None; self.coords.len()];
        let mut nv = 0;
        for block in self.blocks.iter_mut() {
            for elem in block.elems_mut() {
                for node in elem.nodes_mut() {
                    let i = *node as usize;
                    if remap[i].is_none() {
                        remap[i] = Some(nv);
                        nv += 1;
                    }
                    *node = remap[i].unwrap();
                }
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
    type Item = u32;

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
    type Item = &'a mut u32;

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
pub struct Node1([u32; 1]);

impl Node1 {
    pub const fn new(i: [u32; 1]) -> Self {
        Self(i)
    }
}

impl Element for Node1 {
    const N_NODES: usize = 1;

    fn node(&self, i: usize) -> u32 {
        self.0[i]
    }

    fn node_mut(&mut self, i: usize) -> &mut u32 {
        &mut self.0[i]
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Bar2([u32; 2]);

impl Bar2 {
    pub const fn new(i: [u32; 2]) -> Self {
        Self(i)
    }
}

impl Element for Bar2 {
    const N_NODES: usize = 2;

    fn node(&self, i: usize) -> u32 {
        self.0[i]
    }

    fn node_mut(&mut self, i: usize) -> &mut u32 {
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
pub struct Tri3([u32; 3]);

impl Tri3 {
    pub const fn new(i: [u32; 3]) -> Self {
        Self(i)
    }
}

impl Element for Tri3 {
    const N_NODES: usize = 3;

    fn node(&self, i: usize) -> u32 {
        self.0[i]
    }

    fn node_mut(&mut self, i: usize) -> &mut u32 {
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
pub struct Tet4([u32; 4]);

impl Tet4 {
    pub const fn new(i: [u32; 4]) -> Self {
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
}

impl Element for Tet4 {
    const N_NODES: usize = 4;

    fn node(&self, i: usize) -> u32 {
        self.0[i]
    }

    fn node_mut(&mut self, i: usize) -> &mut u32 {
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
        let slice_len = usize::max(32, utf8.len());
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
        let mut block = Block::new();
        block.add_elem(Tri3([0, 1, 2]));
        block.add_elem(Tri3([0, 2, 3]));
        let side = block.boundary();
        println!("{:?}", side);
        assert_eq!(side.sides().len(), 4);
        mesh.add_block(block);
        mesh.add_side(side);
    }

    #[test]
    fn test_bcc_tet4_boundary() {
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(0.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(1.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(1.0, 0.0, 1.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 1.0));
        let tet4 = Tet4([0, 1, 2, 3]);
        let mut block = Block::new();
        block.add_elem(tet4);
        let side = block.boundary();
        println!("{:?}", side);
        assert_eq!(side.sides().len(), 4);
        mesh.add_block(block);
        mesh.add_side(side);

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

        let mut block = Block::new();
        block.add_elem(Tet4::new([0, 1, 2, 5]));
        block.add_elem(Tet4::new([0, 2, 7, 5]));
        block.add_elem(Tet4::new([0, 7, 5, 4]));
        block.add_elem(Tet4::new([2, 6, 7, 5]));
        block.add_elem(Tet4::new([0, 2, 3, 7]));

        let side = block.boundary();
        println!("{:?}", side);
        assert_eq!(side.sides().len(), 12);
        mesh.add_side(side);

        let sum_volumes = mesh.volume(&block.elem(0))
            + mesh.volume(&block.elem(1))
            + mesh.volume(&block.elem(2))
            + mesh.volume(&block.elem(3))
            + mesh.volume(&block.elem(4));

        assert!(sum_volumes < 1.0);
        assert!(sum_volumes > 0.999999);
        mesh.add_block(block);
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
        let mut block = Block::new();
        block.add_elem(tet4);
        mesh.add_block(block);
        mesh.compact();
        assert_eq!(mesh.vertices().len(), 4);
    }
}
