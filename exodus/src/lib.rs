use anyhow::Error;
use mesh::{Element, Mesh, Name, Tet4};
use nalgebra::Point3;
use netcdf::attribute::AttrValue;
use netcdf::types::{BasicType, VariableType};
use std::path::Path;
use thiserror::Error;

// Defined in exodus spec.
pub const NC_MAX_DIMS: usize = 65536;
pub const NC_MAX_VARS: usize = 524288;
pub const NC_MAX_VAR_DIMS: usize = 8;

pub fn write(mesh: &Mesh<Tet4, f64>, path: &Path) -> Result<(), Error> {
    let mut file = netcdf::create(path)?;
    file.add_attribute("title", "created with moosecad")?;
    file.add_attribute("version", 5.1f32)?;
    file.add_attribute("api_version", 5.1f32)?;
    file.add_attribute("floating_point_word_size", 8i64)?;

    // dimensions
    const NUM_DIM: usize = 3;
    const LEN_STRING: usize = 32;
    let num_nodes = mesh.vertices().len();
    let num_blocks = 1;
    let num_elem = mesh.elems().len();
    let num_elem_var = mesh.elem_vars().len();
    file.add_dimension("num_dim", NUM_DIM)?;
    file.add_dimension("num_nodes", num_nodes)?;
    file.add_dimension("num_elem", num_elem)?;
    file.add_dimension("num_el_blk", num_blocks)?;
    file.add_dimension("num_elem_var", num_elem_var)?;
    file.add_dimension("len_string", LEN_STRING)?;
    file.add_unlimited_dimension("time_step")?;

    // coordinate names
    let mut coor_names = file.add_variable_with_type(
        "coor_names",
        &["num_dim", "len_string"],
        &VariableType::Basic(BasicType::Char),
    )?;
    let mut coor_names_data = [0; NUM_DIM * LEN_STRING];
    coor_names_data[LEN_STRING * 0] = 'x' as u8;
    coor_names_data[LEN_STRING * 1] = 'y' as u8;
    coor_names_data[LEN_STRING * 2] = 'z' as u8;
    coor_names.put_chars(&coor_names_data, None, None)?;

    // time step
    let mut time_whole = file.add_variable::<f32>("time_whole", &["time_step"])?;
    time_whole.put_values(&[0.0], None, None)?;

    // vertices
    let mut vertices = Vec::with_capacity(num_nodes * 3);
    vertices.resize(vertices.capacity(), 0.0);
    for (i, v) in mesh.vertices().iter().enumerate() {
        for j in 0..NUM_DIM {
            vertices[i + num_nodes * j] = v[j];
        }
    }
    let mut coord = file.add_variable::<f64>("coord", &["num_dim", "num_nodes"])?;
    coord.put_values(&vertices, None, None)?;

    // block
    let mut eb_names = file.add_variable_with_type(
        "eb_names",
        &["num_el_blk", "len_string"],
        &VariableType::Basic(BasicType::Char),
    )?;
    eb_names.put_chars(
        Name::from(mesh.name().unwrap_or_default()).as_ref(),
        None,
        None,
    )?;

    let mut eb_prop = file.add_variable::<i32>("eb_prop1", &["num_el_blk"])?;
    eb_prop.put_values(&[0], None, None)?;

    let num_elems = mesh.elems().len();
    let mut elems = vec![0; num_elems * Tet4::N_NODES];
    for (i, e) in mesh.elems().iter().enumerate() {
        for j in 0..Tet4::N_NODES {
            elems[i * Tet4::N_NODES + j] = e.node(j) as i64 + 1;
        }
    }
    file.add_dimension("num_el_in_blk1", num_elems)?;
    file.add_dimension("num_nod_per_el1", Tet4::N_NODES)?;
    let mut connect =
        file.add_variable::<i64>("connect1", &["num_el_in_blk1", "num_nod_per_el1"])?;
    connect.add_attribute("elem_type", "TETRA4")?;
    connect.put_values(&elems, None, None)?;

    // TODO side sets

    // elem vars
    let mut name_elem_var = file.add_variable_with_type(
        "name_elem_var",
        &["num_elem_var", "len_string"],
        &VariableType::Basic(BasicType::Char),
    )?;
    let mut name_elem_var_data = Vec::with_capacity(mesh.elem_vars().len() * LEN_STRING);
    for name in mesh.elem_vars() {
        name_elem_var_data.extend_from_slice(name.as_ref());
    }
    name_elem_var.put_chars(&name_elem_var_data, None, None)?;

    for i in 0..num_elem_var {
        let mut vals_elem_var = file.add_variable::<f64>(
            &format!("vals_elem_var{}eb1", i),
            &["time_step", "num_el_in_blk1"],
        )?;
        vals_elem_var.put_values(mesh.elem_var(i), None, None)?;
    }

    Ok(())
}

pub fn read(path: &Path) -> Result<Mesh<Tet4, f64>, Error> {
    let file = netcdf::open(path)?;
    let num_dim = file.dimension("num_dim").ok_or(ExodusError)?.len();
    let num_nodes = file.dimension("num_nodes").ok_or(ExodusError)?.len();
    let num_elem_var = file
        .dimension("num_elem_var")
        .map(|d| d.len())
        .unwrap_or_default();
    let len_string = file.dimension("len_string").ok_or(ExodusError)?.len();
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

    let num_el_in_blk = file.dimension("num_el_in_blk1").ok_or(ExodusError)?.len();
    let num_nod_per_el = file.dimension("num_nod_per_el1").ok_or(ExodusError)?.len();
    if num_nod_per_el != Tet4::N_NODES {
        return Err(ExodusError.into());
    }
    let connect = file.variable("connect1").ok_or(ExodusError)?;
    let elem_type = connect.attribute("elem_type").ok_or(ExodusError)?.value()?;
    if let AttrValue::Str(elem_type) = elem_type {
        if elem_type != "TETRA" && elem_type != "TETRA4" {
            return Err(ExodusError.into());
        }
    } else {
        return Err(ExodusError.into());
    }
    let elems = connect.values::<i64>(None, None)?.into_raw_vec();
    for i in 0..num_el_in_blk {
        let a = elems[i * 4] as usize - 1;
        let b = elems[i * 4 + 1] as usize - 1;
        let c = elems[i * 4 + 2] as usize - 1;
        let d = elems[i * 4 + 3] as usize - 1;
        mesh.add_elem(Tet4::new([a, b, c, d]), &[]);
    }
    let eb_names = file.variable("eb_names").ok_or(ExodusError)?;
    let mut name = vec![0; len_string];
    eb_names.raw_values(&mut name, &[0, 0], &[1, len_string])?;
    mesh.set_name(std::str::from_utf8(&name)?);

    if num_elem_var > 0 {
        let name_elem_var = file.variable("name_elem_var").ok_or(ExodusError)?;
        let mut names = vec![0; num_elem_var * len_string];
        name_elem_var.raw_values(&mut names, &[0, 0], &[num_elem_var, len_string])?;
        for i in 0..num_elem_var {
            let offset = i * len_string;
            mesh.add_elem_var(std::str::from_utf8(&names[offset..(offset + len_string)])?);
            let vals_elem_var = file
                .variable(&format!("vals_elem_var{}eb1", i))
                .ok_or(ExodusError)?;
            let values = vals_elem_var.values(None, None)?.into_raw_vec();
            mesh.elem_var_mut(i).copy_from_slice(&values);
        }
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
        mesh.set_name("tet");
        mesh.add_vertex(Point3::new(0.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(1.5, 0.5, 0.5));
        mesh.add_vertex(Point3::new(1.0, 0.0, 1.0));
        mesh.add_vertex(Point3::new(1.0, 1.0, 1.0));
        mesh.add_elem_var("phi");
        mesh.add_elem(Tet4::new([0, 1, 2, 3]), &[1.0]);
        let tmp = TempDir::new("roundtrip").unwrap();
        let p = tmp.path().join("out.e");
        write(&mesh, &p).unwrap();
        let mesh2 = read(&p).unwrap();
        assert_eq!(mesh, mesh2);
    }
}
