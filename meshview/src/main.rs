use anyhow::Error;
use clap::Clap;
use inotify::{Inotify, WatchMask};
use kiss3d::light::Light;
use kiss3d::resource::Mesh;
use kiss3d::window::Window;
use mesh::Element1;
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

fn read_mesh(path: &Path) -> Result<Rc<RefCell<Mesh>>, Error> {
    //use mesh::BoundedElement;
    //use rand::Rng;
    //let mut rng = rand::thread_rng();

    /*let mut mesh = mesh::Mesh::new();
    mesh.add_elem_var("color");
    let offset = tessellate3d::Tile::LATTICE.len();
    let max = (2, 2, 2);

    fn transform(coord: u32, offset: usize, flipped: bool) -> f32 {
        (if flipped {
            1 - coord as usize + offset
        } else {
            coord as usize + offset
        }) as f32
    }

    for x in 0..max.0 {
        let fx = x % 2 == 1;
        for y in 0..max.1 {
            let fy = y % 2 == 1;
            for z in 0..max.2 {
                let fz = z % 2 == 1;
                for c in tessellate3d::Tile::LATTICE.iter() {
                    let cx = transform(c[0], x, fx);
                    let cy = transform(c[1], y, fy);
                    let cz = transform(c[2], z, fz);
                    mesh.add_vertex(mesh::nalgebra::Point3::new(cx, cy, cz));
                }
                let offset = offset * (x * max.1 * max.2 + y * max.2 + z);
                for t in tessellate3d::Tile::TETS.iter().skip(1).take(3) {
                    let color = rng.gen_range(0.0, 1.0);
                    let mut tet = mesh::Tet4::new([
                        t.node(0) + offset,
                        t.node(1) + offset,
                        t.node(2) + offset,
                        t.node(3) + offset,
                    ]);
                    if fx ^ fy ^ fz {
                        tet.flip();
                    }
                    if rng.gen_range(0.0, 1.0) > 0.2 {
                        mesh.add_elem(tet, &[color]);
                    }
                }
            }
        }
    }
    let mesh = mesh.to_boundary();*/

    let mesh = exodus::read(path)?;
    /*let var = mesh.add_elem_var("color");
    for var in mesh.elem_var_mut(var).iter_mut() {
        *var = rng.gen_range(0.0, 1.0);
    }*/
    println!("num tet {}", mesh.elems().len());
    let mesh = mesh.to_boundary();
    println!("num tri {}", mesh.elems().len());
    let mut na_verts = Vec::with_capacity(mesh.elems().len() * 3);
    let mut na_tex = Vec::with_capacity(mesh.elems().len() * 3);
    let mut na_faces = Vec::with_capacity(mesh.elems().len());
    let mut na_normals = Vec::with_capacity(mesh.elems().len());
    for (i, elem) in mesh.elems().iter().enumerate() {
        let elem_var = if mesh.elem_vars().is_empty() {
            0.5
        } else {
            mesh.elem_var(0)[i]
        };
        //let elem_var = mesh.elem_var(var)[i];

        //for elem in elem.sides() {
        let i = na_verts.len();
        for node in elem.nodes() {
            let p = mesh.vertex(node);
            na_verts.push(na::Point3::new(p[0] as f32, p[1] as f32, p[2] as f32));
            na_tex.push(na::Point2::new(elem_var as f32, 0.5));
        }
        na_faces.push(na::Point3::new(i as u16 + 0, i as u16 + 1, i as u16 + 2));
        let n = mesh.normal(&elem);
        na_normals.push(na::Vector3::new(n[0] as f32, n[1] as f32, n[2] as f32));
        //}
    }
    let mesh = Mesh::new(na_verts, na_faces, Some(na_normals), Some(na_tex), true);
    Ok(Rc::new(RefCell::new(mesh)))
}

// TODO clip plane
//let context = Context::get();
//context.enable(gl::CLIP_DISTANCE0);
//gl::ClipPlane(GL_CLIP_PLANE0, [1.0, 1.0, 1.0]);
