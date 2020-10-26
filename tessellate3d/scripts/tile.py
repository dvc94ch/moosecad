import numpy as np

cubic_lattice = [
    np.array([-.5, -.5, -.5]),
    np.array([ .5, -.5, -.5]),
    np.array([ .5,  .5, -.5]),
    np.array([-.5,  .5, -.5]),
    np.array([-.5, -.5,  .5]),
    np.array([ .5, -.5,  .5]),
    np.array([ .5,  .5,  .5]),
    np.array([-.5,  .5,  .5]),
]

cubic_tets = [
    [0, 1, 2, 5],
    [0, 2, 7, 5],
    [0, 7, 5, 4],
    [2, 6, 7, 5],
    [0, 2, 3, 7]
]

cubic_tets_coords = []
for tet in cubic_tets:
    coords = []
    for i in tet:
        coords.append(cubic_lattice[i])
    cubic_tets_coords.append(coords)

tet_tri_indices = [
    [0, 1, 2],
    [0, 3, 1],
    [1, 3, 2],
    [2, 3, 0],
]

def volume(tet):
    return 1/6 * abs(np.linalg.det(np.array([tet[0]-tet[3], tet[1]-tet[3], tet[2]-tet[3]])))

def normal(tri):
    return np.cross(tri[1]-tri[0], tri[2]-tri[0])

def centroid(tet):
    return sum(tet) * 1/4

if __name__ == '__main__':
    total_volume = 0;
    for tet in cubic_tets_coords:
        print('coords', tet)
        v = volume(tet)
        print('volume', v)
        c = centroid(tet)
        print('center', c)
        total_volume += v

        normals = []
        for tri in tet_tri_indices:
            n = normal(np.array([
                tet[tri[0]],
                tet[tri[1]],
                tet[tri[2]],
            ]))
            print('normal', n)
            normals.append(n)
        # should be zero
        print('sum_normals', sum(normals))

    # should be one
    print(total_volume)
