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
    //let material = Rc::new(RefCell::new(Box::new(NormalsMaterial::new()) as Box<dyn Material>));

    let mesh = read_mesh(&opts.input)?;
    let mut object = window.add_mesh(mesh, scale);
    //object.set_color(1.0, 1.0, 0.0);
    object.set_lines_width(lines_width);
    object.set_texture_from_memory(rainbow, "rainbow");
    object.set_lines_color(lines_color);
    //object.set_material(material.clone());
    window.show();

    loop {
        if !window.render() {
            break;
        }

        if inotify.read_events(&mut buffer)?.next().is_some() {
            // can cause race condition, so we sleep to try to avoid it
            std::thread::sleep(Duration::from_millis(100));
            // drain events to make sure we don't have duplicates
            while inotify.read_events(&mut buffer)?.next().is_some() {}
            println!("rerendering");
            window.remove_node(&mut object);
            let mesh = read_mesh(&opts.input)?;
            object = window.add_mesh(mesh, scale);
            //object.set_color(1.0, 1.0, 0.0);
            object.set_lines_width(lines_width);
            object.set_texture_from_memory(rainbow, "rainbow");
            object.set_lines_color(lines_color);
            //object.set_material(material.clone());
        }
    }

    Ok(())
}

/*pub struct Tile;

impl Tile {
    const LATTICE: [[u32; 3]; 8] = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ];
    const TETS: [mesh::Tet4; 5] = [
        mesh::Tet4::new([0, 1, 2, 5]),
        mesh::Tet4::new([0, 2, 7, 5]),
        mesh::Tet4::new([0, 7, 4, 5]),
        mesh::Tet4::new([2, 6, 7, 5]),
        mesh::Tet4::new([0, 2, 3, 7]),
    ];
}

pub struct Tris;

impl Tris {
    const LATTICE: [[u32; 3]; 4] = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 0, 1],
        [0, 0, 1],
        /*[0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],*/
    ];
    const TETS: [mesh::Tri3; 1] = [
        mesh::Tri3::new([0, 1, 2]),
        //mesh::Tri3::new([0, 2, 1]),
        //mesh::Tri3::new([0, 2, 3]),
    ];
}

pub struct Tets;

impl Tets {
    const LATTICE: [[i32; 3]; 5] = [
        [0, 0, 0],
        [1, 0, 1],
        [0, 1, 1],
        [-1, 0, 1],
        [0, 0, 2],
    ];
    const TETS: [mesh::Tet4; 2] = [
        mesh::Tet4::new([0, 1, 2, 3]),
        mesh::Tet4::new([1, 2, 3, 4]),
    ];
}*/

fn read_mesh(path: &Path) -> Result<Rc<RefCell<Mesh>>, Error> {
    /*use mesh::nalgebra::Point3;
    let mut mesh = mesh::Mesh::new();
    for c in Tile::LATTICE.iter() {
        mesh.add_vertex(Point3::new(-(c[0] as f32), c[2] as f32, c[1] as f32));
    }
    let mut block = mesh::Block::new();
    for t in Tile::TETS.iter() {
        block.add_elem(*t);
    }
    mesh.add_block(block);
    let mesh = mesh.boundary();*/

    let mesh = exodus::read(path)?.boundary();
    let na_verts = mesh
        .vertices()
        .iter()
        .map(|p| na::Point3::new(p[0] as f32, p[1] as f32, p[2] as f32))
        .collect::<Vec<_>>();
    let mut na_faces = Vec::new();
    let mut na_normals = Vec::new();
    for block in mesh.blocks().iter() {
        for elem in block.elems().iter() {
            //for elem in elem.sides() {
            na_faces.push(na::Point3::new(
                elem.node(1) as u16,
                elem.node(0) as u16,
                elem.node(2) as u16,
            ));
            let n = mesh.normal(&elem);
            let x = n[0] as f32;
            let y = n[1] as f32;
            let z = n[2] as f32;
            na_normals.push(na::Vector3::new(x, y, z));
            //}
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
