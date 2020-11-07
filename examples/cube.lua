cube = geom.cube(1, 1, 1, 0.0)
top = cube:intersection(geom.plane_x(0.5), 0.0)
bottom = cube:intersection(geom.plane_neg_x(0.5), 0.0)
return {
    cube = cube:volume(),
    top = top:boundary(),
    bottom = bottom:boundary(),
}
