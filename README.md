# Moose CAD

- constructive solid geometry in lua
- tetrahedralize volumes using isosurface stuffing
- export volumes and boundaries to exodus2
- view meshes and simulation results with automatic reloading on file changes

## Example
```lua
cube = geom.cube(1, 1, 1, 0.0)
plane = geom.cube(1.0, 1.0, 0.01, 0.0)
top = cube:intersection(plane:translate(0.0, 0.0, 0.5), 0.0)
bottom = cube:intersection(plane:translate(0.0, 0.0, -0.5), 0.0)
return {
    cube = cube:volume(),
    top = top:boundary(),
    bottom = bottom:boundary(),
}
```

```sh
./moosecad -i examples/cube.lua -o cube.e --tessellation-resolution 1.0
#./moosesim -i cube.e -o cube.e -l bottom=0 -l top=100
./meshview -i cube.e
```

## License
MIT OR Apache-2.0
