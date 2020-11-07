cube = geom.cube(1, 1, 1, 0.0)
plane = geom.cube(1.0, 1.0, 0.01, 0.0)
top = cube:intersection(plane:translate(0.0, 0.0, 0.5), 0.0)
bottom = cube:intersection(plane:translate(0.0, 0.0, -0.5), 0.0)
return {
    cube = cube:volume(),
    top = top:boundary(),
    bottom = bottom:boundary(),
}
