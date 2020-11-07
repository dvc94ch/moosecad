use anyhow::Error;
use clap::Clap;
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

struct ObjectAdaptor {
    implicit: Box<dyn luascad::implicit3d::Object<f64>>,
    tessellation_resolution: f64,
}

impl tessellate3d::ImplicitFunction<f64> for ObjectAdaptor {
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
    let mut obj = luascad::eval(&utf8)?;
    obj.set_parameters(&luascad::implicit3d::PrimitiveParameters {
        fade_range,
        r_multiplier,
    });
    println!("starting tessellation");
    let mesh = tessellate3d::IsosurfaceStuffing::new(
        &ObjectAdaptor {
            implicit: obj,
            tessellation_resolution,
        },
        tessellation_resolution,
        tessellation_error,
    )
    .tessellate();
    println!("tessellation complete");
    exodus::write(&mesh, &opts.output)?;

    Ok(())
}
