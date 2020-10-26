use nalgebra::{Point3, RealField, Vector3};
use std::collections::HashMap;
use std::fmt::Debug;

pub trait Element: Copy + Debug + Eq {
    const N_NODES: usize;
    fn node(&self, i: usize) -> u32;
    fn nodes(&self) -> NodeIter<'_, Self> {
        NodeIter::new(self)
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

impl<E: BoundedElement> AsRef<[E]> for Block<E> {
    fn as_ref(&self) -> &[E] {
        self.elems.as_ref()
    }
}

impl<E: BoundedElement> AsMut<[E]> for Block<E> {
    fn as_mut(&mut self) -> &mut [E] {
        self.elems.as_mut()
    }
}

impl<E: BoundedElement> Block<E> {
    pub fn new() -> Self {
        Self { elems: Default::default() }
    }

    pub fn add_elem(&mut self, elem: E) {
        self.elems.push(elem);
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
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Side<E: BoundedElement> {
    sides: Vec<E::Side>,
}

impl<E: BoundedElement> AsRef<[E::Side]> for Side<E> {
    fn as_ref(&self) -> &[E::Side] {
        self.sides.as_ref()
    }
}

impl<E: BoundedElement> AsMut<[E::Side]> for Side<E> {
    fn as_mut(&mut self) -> &mut [E::Side] {
        self.sides.as_mut()
    }
}

impl<E: BoundedElement> Side<E> {
    pub fn new() -> Self {
        Self { sides: Default::default() }
    }

    pub fn add_side(&mut self, side: E::Side) {
        self.sides.push(side);
    }
}

impl<E: BoundedElement, S: RealField> Mesh<E, S> {
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

    pub fn vertex(&self, i: usize) -> Point3<S> {
        self.coords[i]
    }

    pub fn vertices(&self) -> &[Point3<S>] {
        &self.coords
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

    pub fn add_side(&mut self, side: Side<E>) {
        self.sides.push(side);
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
        todo!()
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
}

impl BoundedElement for Tri3 {
    const N_SIDES: usize = 3;

    type Side = Bar2;

    fn side(&self, i: usize) -> Self::Side {
        let indices = [[0, 1], [1, 2], [2, 0]][i];
        Bar2([
            self.node(indices[0]),
            self.node(indices[1]),
        ])
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Tet4([u32; 4]);

impl Tet4 {
    pub const fn new(i: [u32; 4]) -> Self {
        Self(i)
    }
}

impl Element for Tet4 {
    const N_NODES: usize = 4;

    fn node(&self, i: usize) -> u32 {
        self.0[i]
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
        assert_eq!(side.as_ref().len(), 4);
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
        let mut side = block.boundary();
        println!("{:?}", side);
        assert_eq!(side.as_ref().len(), 4);
        mesh.add_block(block);
        mesh.add_side(side);

        //assert!(mesh.normal(&tet4.side(0))[2] < 0.0);
        //assert!(mesh.normal(&tet4.side(1))[2] > 0.0);
        //assert!(mesh.normal(&tet4.side(2))[2] > 0.0);
        //assert!(mesh.normal(&tet4.side(3))[2] < 0.0);
    }
    /*
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(0.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.0, 0.0, 0.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 0.0));
        mesh.add_vertex(Point3::new(0.0, 1.0, 0.0));
        mesh.add_vertex(Point3::new(0.0, 0.0, 1.0));
        mesh.add_vertex(Point3::new(1.0, 0.0, 1.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 1.0));
        mesh.add_vertex(Point3::new(0.0, 1.0, 1.0));
        mesh.add_vertex(Point3::new(0.5, 0.5, 0.5));
        let mut block = Block::new();
        block.add_elem(Tet4([1, 2, 3, 0]));
        block.add_elem(Tet4([1, 3, 4, 0]));
        block.add_elem(Tet4([2, 3, 6, 0]));
        block.add_elem(Tet4([3, 6, 7, 0]));
        block.add_elem(Tet4([5, 6, 7, 0]));
        block.add_elem(Tet4([5, 8, 7, 0]));
        block.add_elem(Tet4([1, 6, 2, 0]));
        block.add_elem(Tet4([1, 5, 6, 0]));
        block.add_elem(Tet4([3, 7, 4, 0]));
        block.add_elem(Tet4([4, 7, 8, 0]));
        block.add_elem(Tet4([1, 4, 5, 0]));
        block.add_elem(Tet4([5, 8, 6, 0]));
        let mut side = block.boundary();
        println!("{:?}", side);
        assert_eq!(side.as_ref().len(), 12);
        mesh.add_block(block);
        mesh.add_side(side);


    }*/
}
