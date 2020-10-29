use anyhow::Error;
use clap::Clap;
use inotify::{Inotify, WatchMask};
use kiss3d::light::Light;
use kiss3d::resource::Mesh;
use kiss3d::window::Window;
use mesh::Element;
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

    let mesh = exodus::read(path)?.to_boundary();
    let mut na_verts = Vec::with_capacity(mesh.elems().len() * 3);
    let mut na_tex = Vec::with_capacity(mesh.elems().len() * 3);
    let mut na_faces = Vec::with_capacity(mesh.elems().len());
    let mut na_normals = Vec::with_capacity(mesh.elems().len());
    for (i, elem) in mesh.elems().iter().enumerate() {
        let elem_var = mesh.elem_var(0)[i];
        for node in elem.nodes() {
            let p = mesh.vertex(node);
            na_verts.push(na::Point3::new(p[0] as f32, p[1] as f32, p[2] as f32));
            na_tex.push(na::Point2::new(elem_var as f32, 0.5));
        }
        na_faces.push(na::Point3::new(
            i as u16 * 3 + 1,
            i as u16 * 3 + 0,
            i as u16 * 3 + 2,
        ));
        let n = mesh.normal(&elem);
        na_normals.push(na::Vector3::new(n[0] as f32, n[1] as f32, n[2] as f32));
    }
    let mesh = Mesh::new(na_verts, na_faces, Some(na_normals), Some(na_tex), true);
    Ok(Rc::new(RefCell::new(mesh)))
}

// TODO clip plane
//let context = Context::get();
//context.enable(gl::CLIP_DISTANCE0);
//gl::ClipPlane(GL_CLIP_PLANE0, [1.0, 1.0, 1.0]);
