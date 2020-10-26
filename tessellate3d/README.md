# tesselate3d
`tesselate3d` creates a tetrahedral volumetric mesh from an implicit function using the isosurface
stuffing algorithm. From a volumetric mesh a triangle surface mesh is trivially constructed by
removing faces that occur twice in the mesh. Typically rendering is performed using a surface mesh
while physics simulations require a volumetric mesh. This approach is much faster than first
creating a surface mesh and then computing a volumetric mesh from it.
