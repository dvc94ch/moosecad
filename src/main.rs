use anyhow::Error;
use clap::Clap;
use mesh::{BoundedElement, Element1, SideSet};
use nalgebra as na;
use std::path::PathBuf;

#[derive(Clap)]
struct Opts {
    #[clap(short = 'i', long = "input")]
    input: PathBuf,
    #[clap(short = 'o', long = "output")]
    output: PathBuf,
    #[clap(long = "fade-range")]
    fade_range: Option<f64>,
    #[clap(long = "r-multiplier")]
    r_multiplier: Option<f64>,
    #[clap(long = "tessellation-resolution")]
    tessellation_resolution: Option<f64>,
    #[clap(long = "tessellation-error")]
    tessellation_error: Option<f64>,
}

struct ObjectAdaptor<'a> {
    implicit: &'a dyn luascad::implicit3d::Object<f64>,
    tessellation_resolution: f64,
}

impl<'a> tessellate3d::ImplicitFunction<f64> for ObjectAdaptor<'a> {
    fn bbox(&self) -> &tessellate3d::BoundingBox<f64> {
        self.implicit.bbox()
    }

    fn value(&self, p: &na::Point3<f64>) -> f64 {
        self.implicit.approx_value(p, self.tessellation_resolution)
    }

    fn normal(&self, p: &na::Point3<f64>) -> na::Vector3<f64> {
        self.implicit.normal(p)
    }
}

fn main() -> Result<(), Error> {
    let opts = Opts::parse();

    let fade_range = opts.fade_range.unwrap_or(0.1);
    let r_multiplier = opts.r_multiplier.unwrap_or(1.0);
    let tessellation_resolution = opts.tessellation_resolution.unwrap_or(0.12);
    let tessellation_error = opts.tessellation_error.unwrap_or(0.2);

    let bytes = std::fs::read(opts.input)?;
    let utf8 = String::from_utf8(bytes)?;
    let mut output = luascad::eval(&utf8)?;
    for (_, out) in &mut output {
        out.object_mut().set_parameters(&luascad::implicit3d::PrimitiveParameters {
            fade_range,
            r_multiplier,
        });
    }
    let volumes = output
        .iter()
        .filter_map(|(n, o)| o.volume().map(|o| (n, o)))
        .collect::<Vec<_>>();
    let boundaries = output
        .iter()
        .filter_map(|(n, o)| o.boundary().map(|o| (n, o)))
        .collect::<Vec<_>>();

    if volumes.len() < 1 {
        println!("no volumes");
        return Ok(());
    }

    println!("tessellating {}", volumes[0].0);
    let mut mesh = tessellate3d::IsosurfaceStuffing::new(
        &ObjectAdaptor {
            implicit: volumes[0].1,
            tessellation_resolution,
        },
        tessellation_resolution,
        tessellation_error,
    )
    .tessellate();
    mesh.set_name(volumes[0].0);
    println!("tessellation complete");

    let surface = mesh.boundary();
    for (name, boundary) in boundaries {
        println!("adding boundary {}", name);
        let mut side_set = SideSet::new();
        side_set.set_name(name);
        for (elem, side) in surface.sides() {
            let mut add = true;
            for node in mesh.elem(*elem).side(*side).nodes() {
                let point = mesh.vertex(node);
                let phi = boundary.approx_value(point, f64::EPSILON);
                if phi > f64::EPSILON || phi < -f64::EPSILON {
                    add = false;
                    break;
                }
            }
            if add {
                side_set.add_side(*elem, *side);
            }
        }
        mesh.add_side_set(side_set);
    }

    exodus::write(&mesh, &opts.output)?;

    Ok(())
}
