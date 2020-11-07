cube = geom.cube(1, 1, 1, 0.3)
sphere = geom.sphere(0.5)
diff = cube:difference(sphere, 0.3)
return diff:scale(15, 15, 15)
