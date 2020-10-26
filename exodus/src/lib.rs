use mesh::{BoundedElement, Mesh};
use nalgebra::RealField;
use std::io;
use std::path::Path;

// Defined in exodus spec.
pub const NC_MAX_DIMS: usize = 65536;
pub const NC_MAX_VARS: usize = 524288;
pub const NC_MAX_VAR_DIMS: usize = 8;

pub fn write<E: BoundedElement, S: RealField>(mesh: &Mesh<E, S>, path: &Path) -> io::Result<()> {
    Ok(())
}
