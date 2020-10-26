use anyhow::Error;
use clap::Clap;
use inotify::{Inotify, WatchMask};
use kiss3d::light::Light;
use kiss3d::resource::Mesh;
use kiss3d::window::Window;
use mesh::{BoundedElement, Element};
use nalgebra as na;
use std::cell::RefCell;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::time::Duration;

#[derive(Clap)]
struct Opts {
    #[clap(short = 'i', long = "input")]
    input: PathBuf,
}

fn main() -> Result<(), Error> {
    let opts = Opts::parse();

    let mut inotify = Inotify::init()?;
    inotify.add_watch(&opts.input, WatchMask::CLOSE_WRITE)?;
    let mut buffer = [0; 1024];

    let mut window = Window::new("MeshView");
    let rainbow = include_bytes!("../res/rainbow.png");
    window.set_light(Light::StickToCamera);

    let scale = na::Vector3::new(1.0, 1.0, 1.0);
    let lines_width = 10.0;
    let lines_color = Some(na::Point3::new(0.0, 0.0, 0.0));

    let mesh = read_mesh(&opts.input)?;
    let mut object = window.add_mesh(mesh, scale);
    //object.set_color(1.0, 1.0, 0.0);
    object.set_lines_width(lines_width);
    object.set_texture_from_memory(rainbow, "rainbow");
    object.set_lines_color(lines_color);
    window.show();

    loop {
        if !window.render() {
            break;
        }

        if let Some(_) = inotify.read_events(&mut buffer)?.next() {
            // can cause race condition, so we sleep to try to avoid it
            std::thread::sleep(Duration::from_millis(100));
            println!("rerendering");
            window.remove_node(&mut object);
            let mesh = read_mesh(&opts.input)?;
            object = window.add_mesh(mesh, scale);
            //object.set_color(1.0, 1.0, 0.0);
            object.set_lines_width(lines_width);
            object.set_texture_from_memory(rainbow, "rainbow");
            object.set_lines_color(lines_color);
        }
    }

    Ok(())
}

fn read_mesh(path: &Path) -> Result<Rc<RefCell<Mesh>>, Error> {
    let mesh = exodus::read(path)?;
    let na_verts = mesh
        .vertices()
        .iter()
        .map(|p| na::Point3::new(p[0] as f32, p[1] as f32, p[2] as f32))
        .collect::<Vec<_>>();
    let mut na_faces = Vec::new();
    for block in mesh.blocks().iter() {
        for elem in block.elems().iter() {
            for side in elem.sides() {
                na_faces.push(na::Point3::new(
                    side.node(0) as u16,
                    side.node(1) as u16,
                    side.node(2) as u16,
                ));
            }
        }
    }
    let mut na_tex = Vec::with_capacity(na_verts.len());
    for _ in na_verts.iter() {
        na_tex.push(na::Point2::new(0.5, 0.5));
    }
    let mesh = Mesh::new(na_verts, na_faces, None, Some(na_tex), true);
    Ok(Rc::new(RefCell::new(mesh)))
}

// TODO clip plane
//let context = Context::get();
//context.enable(gl::CLIP_DISTANCE0);
//gl::ClipPlane(GL_CLIP_PLANE0, [1.0, 1.0, 1.0]);
