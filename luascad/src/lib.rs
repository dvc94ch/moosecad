use hlua_derive::lua;
pub use implicit3d;
use std::collections::HashMap;

pub struct LuaSandbox<'lua> {
    lua: hlua::Lua<'lua>,
}

impl<'lua> LuaSandbox<'lua> {
    const USER_FUNCTION_NAME: &'static str = "__user_function__";
    const SANDBOX_ENV_NAME: &'static str = "__sandbox_env__";
    // This env is taken from:
    // http://stackoverflow.com/questions/1224708/how-can-i-create-a-secure-lua-sandbox
    const SANDBOX_ENV: &'static str = "{
        ipairs = ipairs,
        next = next,
        pairs = pairs,
        pcall = pcall,
        print = print,
        tonumber = tonumber,
        tostring = tostring,
        type = type,
        unpack = unpack,
        coroutine = { create = coroutine.create, resume = coroutine.resume,
            running = coroutine.running, status = coroutine.status,
            wrap = coroutine.wrap },
        string = { byte = string.byte, char = string.char, find = string.find,
            format = string.format, gmatch = string.gmatch, gsub = string.gsub,
            len = string.len, lower = string.lower, match = string.match,
            rep = string.rep, reverse = string.reverse, sub = string.sub,
            upper = string.upper },
        table = { insert = table.insert, maxn = table.maxn, remove = table.remove,
            sort = table.sort },
        math = { abs = math.abs, acos = math.acos, asin = math.asin,
            atan = math.atan, atan2 = math.atan2, ceil = math.ceil, cos = math.cos,
            cosh = math.cosh, deg = math.deg, exp = math.exp, floor = math.floor,
            fmod = math.fmod, frexp = math.frexp, huge = math.huge,
            ldexp = math.ldexp, log = math.log, log10 = math.log10, max = math.max,
            min = math.min, modf = math.modf, pi = math.pi, pow = math.pow,
            rad = math.rad, random = math.random, sin = math.sin, sinh = math.sinh,
            sqrt = math.sqrt, tan = math.tan, tanh = math.tanh },
        os = { clock = os.clock, difftime = os.difftime, time = os.time },
    }";

    pub fn new() -> Result<Self, hlua::LuaError> {
        let mut lua = hlua::Lua::new();
        lua.openlibs();
        lua.execute(&format!(
            "{env} = {val}",
            env = Self::SANDBOX_ENV_NAME,
            val = Self::SANDBOX_ENV,
        ))?;
        Ok(Self { lua })
    }

    pub fn load<'a, F>(&'a mut self, luaopen: &F)
    where
        F: Fn(hlua::LuaTable<hlua::PushGuard<&'a mut hlua::Lua<'lua>>>),
    {
        let env = self.lua.get(Self::SANDBOX_ENV_NAME).unwrap();
        luaopen(env);
    }

    pub fn execute<'a, T>(&'a mut self, script: &str) -> Result<T, hlua::LuaError>
    where
        T: for<'g> hlua::LuaRead<hlua::PushGuard<&'g mut hlua::PushGuard<&'a mut hlua::Lua<'lua>>>>,
    {
        self.lua
            .checked_set(Self::USER_FUNCTION_NAME, hlua::LuaCode(script))?;
        self.lua.execute(&format!(
            "debug.setupvalue({f}, 1, {env}); return {}();",
            f = Self::USER_FUNCTION_NAME,
            env = Self::SANDBOX_ENV_NAME,
        ))
    }
}

#[lua]
pub mod geom {
    use implicit3d::{self as sdf, Intersection, Object};
    use nalgebra as na;

    #[derive(Clone, Debug)]
    pub struct LObject {
        pub o: Box<dyn Object<f64>>,
    }

    impl LObject {
        pub fn plane_x(x: f64) -> Self {
            Self {
                o: Box::new(sdf::PlaneX::new(x)),
            }
        }

        pub fn plane_y(y: f64) -> Self {
            Self {
                o: Box::new(sdf::PlaneY::new(y)),
            }
        }

        pub fn plane_z(z: f64) -> Self {
            Self {
                o: Box::new(sdf::PlaneZ::new(z)),
            }
        }

        pub fn plane_neg_x(x: f64) -> Self {
            Self {
                o: Box::new(sdf::PlaneNegX::new(x)),
            }
        }

        pub fn plane_neg_y(y: f64) -> Self {
            Self {
                o: Box::new(sdf::PlaneNegY::new(y)),
            }
        }

        pub fn plane_neg_z(z: f64) -> Self {
            Self {
                o: Box::new(sdf::PlaneNegZ::new(z)),
            }
        }

        pub fn plane_hesian(nx: f64, ny: f64, nz: f64, p: f64) -> Self {
            Self {
                o: Box::new(sdf::NormalPlane::from_normal_and_p(
                    na::Vector3::new(nx, ny, nz),
                    p,
                )),
            }
        }

        pub fn plane_3_points(
            ax: f64,
            ay: f64,
            az: f64,
            bx: f64,
            by: f64,
            bz: f64,
            cx: f64,
            cy: f64,
            cz: f64,
        ) -> Self {
            Self {
                o: Box::new(sdf::NormalPlane::from_3_points(
                    &na::Point3::new(ax, ay, az),
                    &na::Point3::new(bx, by, bz),
                    &na::Point3::new(cx, cy, cz),
                )),
            }
        }

        pub fn sphere(radius: f64) -> Self {
            Self {
                o: Box::new(sdf::Sphere::new(radius)),
            }
        }

        pub fn i_cylinder(radius: f64) -> Self {
            Self {
                o: Box::new(sdf::Cylinder::new(radius)),
            }
        }

        pub fn i_cone(slope: f64) -> Self {
            Self {
                o: Box::new(sdf::Cone::new(slope, 0.0)),
            }
        }

        // TODO smooth optional
        pub fn cube(x: f64, y: f64, z: f64, smooth: f64) -> LObject {
            let b = Intersection::from_vec(
                vec![
                    Box::new(sdf::PlaneX::new(x / 2.0)),
                    Box::new(sdf::PlaneY::new(y / 2.0)),
                    Box::new(sdf::PlaneZ::new(z / 2.0)),
                    Box::new(sdf::PlaneNegX::new(x / 2.0)),
                    Box::new(sdf::PlaneNegY::new(y / 2.0)),
                    Box::new(sdf::PlaneNegZ::new(z / 2.0)),
                ],
                smooth,
            )
            .unwrap();
            Self { o: b }
        }

        // TODO smooth optional, r vs r1 && r2
        pub fn cylinder(length: f64, r1: f64, r2: f64, smooth: f64) -> Self {
            let cone = if (r1 - r2).abs() < f64::EPSILON {
                Box::new(sdf::Cylinder::new(r1)) as Box<dyn Object<f64>>
            } else {
                let slope = (r2 - r1).abs() / length;
                let offset = if r1 < r2 {
                    -r1 / slope - length * 0.5
                } else {
                    r2 / slope + length * 0.5
                };
                let mut cone = Box::new(sdf::Cone::new(slope, offset));
                let rmax = r1.max(r2);
                cone.set_bbox(&sdf::BoundingBox::new(
                    &na::Point3::new(-rmax, -rmax, f64::NEG_INFINITY),
                    &na::Point3::new(rmax, rmax, f64::INFINITY),
                ));
                cone
            };
            Self {
                o: Intersection::from_vec(
                    vec![
                        cone,
                        Box::new(sdf::PlaneZ::new(length / 2.0)),
                        Box::new(sdf::PlaneNegZ::new(length / 2.0)),
                    ],
                    smooth,
                )
                .unwrap(),
            }
        }

        pub fn translate(&self, x: f64, y: f64, z: f64) -> Self {
            Self {
                o: self.o.translate(&na::Vector3::new(x, y, z)),
            }
        }

        pub fn rotate(&self, x: f64, y: f64, z: f64) -> Self {
            Self {
                o: self.o.rotate(&na::Vector3::new(x, y, z)),
            }
        }

        pub fn scale(&self, x: f64, y: f64, z: f64) -> Self {
            Self {
                o: self.o.scale(&na::Vector3::new(x, y, z)),
            }
        }

        pub fn bend(&self, width: f64) -> Self {
            Self {
                o: Box::new(sdf::Bender::new(self.o.clone(), width)),
            }
        }

        pub fn twist(&self, height: f64) -> Self {
            Self {
                o: Box::new(sdf::Twister::new(self.o.clone(), height)),
            }
        }

        // TODO smooth optional and take a list of other
        pub fn r#union(&self, other: &Self, smooth: f64) -> Self {
            let v = vec![self.o.clone(), other.o.clone()];
            Self {
                o: sdf::Union::from_vec(v, smooth).unwrap(),
            }
        }

        // TODO smooth optional and take a list of other
        pub fn difference(&self, other: &Self, smooth: f64) -> Self {
            let v = vec![self.o.clone(), other.o.clone()];
            Self {
                o: sdf::Intersection::difference_from_vec(v, smooth).unwrap(),
            }
        }

        // TODO smooth optional and take a list of other
        pub fn intersection(&self, other: &Self, smooth: f64) -> Self {
            let v = vec![self.o.clone(), other.o.clone()];
            Self {
                o: sdf::Intersection::from_vec(v, smooth).unwrap(),
            }
        }

        pub fn volume(&self) -> Output {
            Output::Volume(self.o.clone())
        }

        pub fn boundary(&self) -> Output {
            Output::Boundary(self.o.clone())
        }

        #[lua(meta = "__tostring")]
        pub fn to_string(&self) -> String {
            format!("{:#?}", self.o)
        }
    }

    #[derive(Clone, Debug)]
    pub enum Output {
        Volume(Box<dyn Object<f64>>),
        Boundary(Box<dyn Object<f64>>),
    }

    impl Output {
        #[lua(meta = "__tostring")]
        pub fn to_string(&self) -> String {
            format!("{:#?}", self)
        }
    }
}

impl geom::Output {
    pub fn volume(&self) -> Option<&dyn implicit3d::Object<f64>> {
        if let Self::Volume(obj) = self {
            Some(&**obj)
        } else {
            None
        }
    }

    pub fn boundary(&self) -> Option<&dyn implicit3d::Object<f64>> {
        if let Self::Boundary(obj) = self {
            Some(&**obj)
        } else {
            None
        }
    }

    pub fn object_mut(&mut self) -> &mut dyn implicit3d::Object<f64> {
        match self {
            Self::Volume(obj) => &mut **obj,
            Self::Boundary(obj) => &mut **obj,
        }
    }
}

pub fn eval(script: &str) -> Result<HashMap<String, geom::Output>, hlua::LuaError> {
    let mut lua = LuaSandbox::new()?;
    lua.load(&geom::load);
    lua.execute(script)
}

#[cfg(test)]
mod tests {
    use super::eval;

    #[test]
    fn test_primitives() {
        let res =
            eval("return { cube = geom.cube(1, 1, 1, 0.0):translate(1.0, 1.0, 1.0):volume() }")
                .unwrap();
        assert!(res.contains_key("cube"));
        let res = eval("return { sphere = geom.sphere(0.5):volume() }").unwrap();
        assert!(res.contains_key("sphere"));
        let res = eval(
            r#"
            cube = geom.cube(1, 1, 1, 0.3)
            sphere = geom.sphere(0.5)
            diff = cube:difference(sphere, 0.3)
            union = cube:union(sphere, 0.3)
            return { xplicit = diff:scale(15, 15, 15):volume() }
        "#,
        )
        .unwrap();
        assert!(res.contains_key("xplicit"));
        let res = eval(
            r#"
            cube = geom.cube(1, 1, 1, 0.0)
            top = cube:intersection(geom.plane_x(0.5), 0.0)
            bottom = cube:intersection(geom.plane_neg_x(0.5), 0.0)
            return {
                cube = cube:volume(),
                top = top:boundary(),
                bottom = bottom:boundary(),
            }"#,
        )
        .unwrap();
        assert!(res.contains_key("cube"));
        assert!(res.contains_key("top"));
        assert!(res.contains_key("bottom"));
    }
}
