use kiss3d::context::Context;
use kiss3d::light::Light;
use kiss3d::resource::Mesh;
use kiss3d::window::Window;
use nalgebra as na;
use std::cell::RefCell;
use std::path::Path;
use std::rc::Rc;

fn main() {
    let coords = [
        [-0.5, -0.5, -0.5],
        [0.5, -0.5, -0.5],
        [0.5, 0.5, -0.5],
        [-0.5, 0.5, -0.5],
        [-0.5, -0.5, 0.5],
        [0.5, -0.5, 0.5],
        [0.5, 0.5, 0.5],
        [-0.5, 0.5, 0.5],
    ];
    let tets = [
        [0, 1, 2, 5],
        [0, 2, 7, 5],
        [0, 7, 5, 4],
        [2, 6, 7, 5],
        [0, 2, 3, 7],
    ];
    let tris = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]];
    let values = [1.0, 0.5, 0.8, 0.8, 0.8, 0.5, 0.8, 0.1];
    let mut na_verts = Vec::new();
    let mut na_tex = Vec::new();
    for coord in coords.iter() {
        na_verts.push(na::Point3::new(coord[0], coord[1], coord[2]));
    }
    for value in values.iter() {
        na_tex.push(na::Point2::new(*value, 0.5));
    }
    let mut na_faces = Vec::new();
    for tet in tets.iter() {
        for tri in tris.iter() {
            let face = na::Point3::new(tet[tri[0]], tet[tri[1]], tet[tri[2]]);
            na_faces.push(face);
        }
    }
    let mesh = Mesh::new(na_verts, na_faces, None, Some(na_tex), true);
    let rc_mesh = Rc::new(RefCell::new(mesh));

    let mut window = Window::new("MeshView");

    // TODO clip plane
    //let context = Context::get();
    //context.enable(gl::CLIP_DISTANCE0);
    //gl::ClipPlane(GL_CLIP_PLANE0, [1.0, 1.0, 1.0]);

    // Load colormap
    let rainbow = window.add_texture(&Path::new("res/rainbow.png"), "rainbow");

    // Add mesh
    let scale = na::Vector3::new(1.0, 1.0, 1.0);
    let mut object = window.add_mesh(rc_mesh, scale);
    //object.set_color(1.0, 1.0, 0.0);
    object.set_texture(rainbow);
    object.set_lines_width(10.0);
    object.set_lines_color(Some(na::Point3::new(0.0, 0.0, 0.0)));

    // Set light and render
    window.set_light(Light::StickToCamera);
    window.show();
    while window.render() {}
}
