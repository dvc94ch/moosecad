const TEST_FILES: [&str; 3] = ["tests/hook-tri.e", "tests/hook-tetra.e", "tests/hook-out.e"];

fn print_file_info(path: &str) {
    let f = netcdf::open(path).unwrap();
    print_group(f.root().unwrap(), 0);
    //assert!(false);
}

fn print_group(g: netcdf::Group, i: usize) {
    let indent = "".repeat(i);
    println!("{}group.{}", indent, g.name());
    for var in g.variables() {
        println!(
            "{}var.{}: [{}; {:?}]",
            indent,
            var.name(),
            var.len(),
            var.vartype()
        );
        for attr in var.attributes() {
            println!("{}  attr.{} = {:?}", indent, attr.name(), attr.value().ok());
        }
        for dim in var.dimensions() {
            println!("{}  dim.{} = {}", indent, dim.name(), dim.len());
        }
    }
    for attr in g.attributes() {
        println!("{}attr.{} = {:?}", indent, attr.name(), attr.value().ok());
    }
    for dim in g.dimensions() {
        println!("{}dim.{} = {}", indent, dim.name(), dim.len());
    }
    for ty in g.types() {
        println!("{}ty.{:?}", indent, ty);
    }
    for group in g.groups() {
        print_group(group, i + 2);
    }
}

#[test]
fn netcdf_parse_tri_mesh() {
    print_file_info(TEST_FILES[0]);
}

#[test]
fn netcdf_parse_tetra_mesh() {
    print_file_info(TEST_FILES[1]);
}

#[test]
#[ignore] // requires enabling 64bit format
fn netcdf_parse_dataset() {
    print_file_info(TEST_FILES[2]);
}
