use anyhow::Error;
use mesh::{Block, BoundedElement, Element, Mesh, Tet4};
use nalgebra::{Point3, RealField};
use netcdf::attribute::AttrValue;
use netcdf::types::{BasicType, VariableType};
use num_traits::Float;
use std::path::Path;
use thiserror::Error;

// Defined in exodus spec.
pub const NC_MAX_DIMS: usize = 65536;
pub const NC_MAX_VARS: usize = 524288;
pub const NC_MAX_VAR_DIMS: usize = 8;

pub fn write<E: BoundedElement, S: Float + RealField>(
    mesh: &Mesh<E, S>,
    path: &Path,
) -> Result<(), Error> {
    let mut file = netcdf::create(path)?;
    file.add_attribute("title", "created with moosecad")?;
    file.add_attribute("version", 5.1f32)?;
    file.add_attribute("api_version", 5.1f32)?;
    file.add_attribute("floating_point_word_size", 8i64)?;

    // dimensions
    let num_dim = 3; // TODO don't hardcode
    let num_nodes = mesh.vertices().len();
    let num_blocks = mesh.blocks().len();
    file.add_dimension("num_dim", num_dim)?;
    file.add_dimension("num_nodes", num_nodes)?;
    file.add_dimension("num_el_blk", num_blocks)?;
    file.add_dimension("len_string", 2)?;
    file.add_unlimited_dimension("time_step")?;

    // coordinate names
    let mut coor_names = file.add_variable_with_type(
        "coor_names",
        &["num_dim", "len_string"],
        &VariableType::Basic(BasicType::Char),
    )?;
    coor_names.put_chars(
        unsafe { std::mem::transmute(&b"x\0y\0z\0"[..]) },
        None,
        None,
    )?;

    // time step
    let mut time_whole = file.add_variable::<f32>("time_whole", &["time_step"])?;
    time_whole.put_values(&[0.0], None, None)?;

    // vertices
    let mut vertices = Vec::with_capacity(num_nodes * 3);
    vertices.resize(vertices.capacity(), 0.0);
    for (i, v) in mesh.vertices().iter().enumerate() {
        for j in 0..num_dim {
            vertices[i + num_nodes * j] = v[j].to_f64().unwrap();
        }
    }
    let mut coord = file.add_variable::<f64>("coord", &["num_dim", "num_nodes"])?;
    coord.put_values(&vertices, None, None)?;

    for (i, block) in mesh.blocks().iter().enumerate() {
        let num_elems = block.elems().len();
        let mut elems = Vec::with_capacity(num_elems * E::N_NODES);
        elems.resize(elems.capacity(), 0);
        for (i, e) in block.elems().iter().enumerate() {
            for j in 0..E::N_NODES {
                elems[i + num_elems * j] = e.node(j) as i64;
            }
        }
        let num_el_in_blk = format!("num_el_in_blk{}", i + 1);
        let num_nod_per_el = format!("num_nod_per_el{}", i + 1);
        let connect = format!("connect{}", i + 1);
        file.add_dimension(&num_el_in_blk, num_elems)?;
        file.add_dimension(&num_nod_per_el, E::N_NODES)?;
        let mut connect = file.add_variable::<i64>(&connect, &[&num_el_in_blk, &num_nod_per_el])?;
        connect.add_attribute("elem_type", "TETRA")?; // TODO don't hardcode "TETRA"
        connect.put_values(&elems, None, None)?;
    }

    Ok(())
}

pub fn read(path: &Path) -> Result<Mesh<Tet4, f64>, Error> {
    let file = netcdf::open(path)?;
    let num_dim = file.dimension("num_dim").ok_or(ExodusError)?.len();
    let num_nodes = file.dimension("num_nodes").ok_or(ExodusError)?.len();
    let num_blocks = file.dimension("num_el_blk").ok_or(ExodusError)?.len();
    if num_dim != 3 {
        return Err(ExodusError.into());
    }

    let mut mesh = Mesh::new();
    let coord = file.variable("coord").ok_or(ExodusError)?;
    if coord.dimensions().len() != 2 {
        return Err(ExodusError.into());
    }
    if coord.dimensions()[0].len() != num_dim {
        return Err(ExodusError.into());
    }
    if coord.dimensions()[1].len() != num_nodes {
        return Err(ExodusError.into());
    }
    let vertices = coord.values(None, None)?.into_raw_vec();
    for i in 0..num_nodes {
        let x = vertices[i];
        let y = vertices[i + num_nodes];
        let z = vertices[i + num_nodes * 2];
        mesh.add_vertex(Point3::new(x, y, z));
    }

    for i in 0..num_blocks {
        let num_el_in_blk = format!("num_el_in_blk{}", i + 1);
        let num_el_in_blk = file.dimension(&num_el_in_blk).ok_or(ExodusError)?.len();
        let num_nod_per_el = format!("num_nod_per_el{}", i + 1);
        let num_nod_per_el = file.dimension(&num_nod_per_el).ok_or(ExodusError)?.len();
        if num_nod_per_el != Tet4::N_NODES {
            return Err(ExodusError.into());
        }
        let connect = format!("connect{}", i + 1);
        let connect = file.variable(&connect).ok_or(ExodusError)?;
        let elem_type = connect.attribute("elem_type").ok_or(ExodusError)?.value()?;
        if let AttrValue::Str(elem_type) = elem_type {
            if elem_type != "TETRA" {
                return Err(ExodusError.into());
            }
        } else {
            return Err(ExodusError.into());
        }
        let elems = connect.values::<i64>(None, None)?.into_raw_vec();
        let mut block = Block::new();
        for i in 0..num_el_in_blk {
            let a = elems[i] as u32;
            let b = elems[i + num_el_in_blk] as u32;
            let c = elems[i + num_el_in_blk * 2] as u32;
            let d = elems[i + num_el_in_blk * 3] as u32;
            block.add_elem(Tet4::new([a, b, c, d]));
        }
        mesh.add_block(block);
    }
    Ok(mesh)
}

#[derive(Debug, Error)]
#[error("error parsing exodus file")]
pub struct ExodusError;

#[cfg(test)]
mod tests {
    use super::*;
    use tempdir::TempDir;

    #[test]
    fn test_roundtrip() {
        let mut mesh = Mesh::new();
        mesh.add_vertex(Point3::new(0.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(1.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(1.0, 0.0, 1.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 1.0));
        let tet4 = Tet4::new([0, 1, 2, 3]);
        let mut block = Block::new();
        block.add_elem(tet4);
        let side = block.boundary();
        println!("{:?}", side);
        assert_eq!(side.sides().len(), 4);
        mesh.add_block(block);
        //mesh.add_side(side);

        let tmp = TempDir::new("roundtrip").unwrap();
        let p = tmp.path().join("out.e");
        write(&mesh, &p).unwrap();
        let mesh2 = read(&p).unwrap();
        assert_eq!(mesh, mesh2);
    }
}
