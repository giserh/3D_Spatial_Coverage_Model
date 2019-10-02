import pickle
import numpy as np
from scipy.spatial import KDTree
from utils3D import *


def algorithm(arg):
    type = arg['type']
    dataset = arg['dataset']
    angles = arg['angles']
    cells = arg['cells']
    iter = 100
    size = arg['size']
    Rv = arg['Rv']
    alpha = arg['alpha']
    path = arg['path']
    query= arg['query']

    pklpath = None

    if dataset == 2 and Rv == 0.05 and alpha == 45:
        pklpath = "data_all_2.pkl"
    if dataset == 2 and Rv == 0.1 and alpha == 45:
        pklpath = "data_all_2_R=0.1.pkl"
    if dataset == 2 and Rv == 0.15 and alpha == 45:
        pklpath = "data_all_2_R=0.15.pkl"

    if dataset == 1 and Rv == 0.05 and alpha == 45:
        pklpath = "data_all_1.pkl"

    fov_pts, fov_angles, x_range, y_range, z_range, scales = pickle.load(open(pklpath, "rb"))

    z_range[0] -= Rv

    pts = [pt[2] for pt in fov_pts]
    pts = np.array(pts)
    kdtree = KDTree(pts)

    np.random.seed(0)


    try:
        queries = pickle.load(open(query, "rb"))

        if abs(queries[0][4] - queries[0][1] - size) > 0.01:
            print("Resizing")
            queries2 = []
            for i in range(iter):
                print(i)
                id, min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons_index = queries[i]
                centerx, centery, centerz = (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2
                min_x = centerx - size / 2
                max_x = centerx + size / 2
                min_y = centery - size / 2
                max_y = centery + size / 2
                min_z = centerz - size / 2
                max_z = centerz + size / 2

                region_query_polygon = get_target_object_polygon(min_x, min_y, min_z, max_x, max_y, max_z)

                fov_polyhedrons_index2 = []
                for i in fov_polyhedrons_index:
                    pts = fov_pts[i]
                    fov_polyhedron = get_fov_polyhedron_from_pts(pts)
                    tmp = fov_polyhedron.intersect(region_query_polygon)
                    if len(tmp.toPolygons()) != 0:
                        fov_polyhedrons_index2.append(i)
                queries2.append([id, min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons_index2])
            queries = queries2

    except:
        queries = []
        it = 0
        while True:
            min_x, min_y, min_z, max_x, max_y, max_z = gen_cube(x_range, y_range, z_range, size)
            region_query_polygon = get_target_object_polygon(min_x, min_y, min_z, max_x, max_y, max_z)
            p = [(min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2]
            R = (max_x - min_x) * np.sqrt(3) / 2 + Rv
            indices = kdtree.query_ball_point(p, R)

            fov_polyhedrons_index = []
            for i in indices:
                pts = fov_pts[i]
                fov_polyhedron = get_fov_polyhedron_from_pts(pts)
                tmp = fov_polyhedron.intersect(region_query_polygon)
                if len(tmp.toPolygons()) != 0:
                    fov_polyhedrons_index.append(i)



            if len(fov_polyhedrons_index) != 0:
                queries.append([it, min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons_index])
                print(min_x, min_y, min_z, max_x, max_y, max_z)
                print("%d / %d" % (len(fov_polyhedrons_index), len(indices)))

                it += 1
                if it == iter:
                    break

        queries.sort(key=lambda x: len(x[-1]))
        pickle.dump(queries, open(query, "wb"))
    # print(queries)

    try:
        data = pickle.load(open(path, "rb"))
    except:
        data = {}

    if type == 1:

        for it in range(iter):
            id, min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons_index = queries[it]
            # print(len(fov_polyhedrons_index))
            print("%d/%d" % (it, iter))
            if data.__contains__(id):
                print("%d %d %f %f" % (id, len(fov_polyhedrons_index), data[id][0], data[id][1]))
                continue
            region_query_polygon = get_target_object_polygon(min_x, min_y, min_z, max_x, max_y, max_z)


            f = open("coverage.in", "w")
            f.write(str(min_x) + ' ' + str(max_x) + ' ' + str(min_y) + ' ' + str(max_y) + ' ' + str(min_z) + ' ' + str(max_z) + '\n')
            f.write(str(len(fov_polyhedrons_index)) + '\n')


            fov_polyhedrons = []
            for index in fov_polyhedrons_index:
                pts = fov_pts[index]
                for pt in pts:
                    f.write(str(pt[0]) + ' ' + str(pt[1]) + ' ' + str(pt[2]) + '\n')
                f.write('\n')
                fov_polyhedron = get_fov_polyhedron_from_pts(pts)
                fov_polyhedrons.append(fov_polyhedron)
            st = time.time()
# -----------------------------------
            if arg['method'] == 'Geo':

                ans1 = execute(region_query_polygon, fov_polyhedrons) # use this line to enable geometry calculation
# -----------------------------------
            if arg['method'] == 'MC':
                f.close()
                import subprocess
                subprocess.run("../build/alg1")
                # ans1 = execute_MC2(min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons)
                f = open("coverage.out", "r")
                ans1 = float(f.read())
                f.close()

            dur = time.time() - st
            data[id] = (ans1, dur)
            pickle.dump(data, open(path, "wb"))
# -----------------------------------
            print("%d %d %f %f" % (id, len(fov_polyhedrons), ans1, dur))

# geometry calculation is much too slow for ECM and WCM, so we ommit them in the public code.
    # if type == 2:
    #
    # ...


if __name__ == '__main__':


    arg = {'type': 1,
           'dataset': 2,
           'angles': None,
           'cells': None,
           'size': 0.1,
           'Rv': 0.05,
           'alpha': 45,
           'path': 'result1_MC.pkl',
           'query': 'query1.pkl',
           'method': 'MC'}
    algorithm(arg)

    arg = {'type': 1,
           'dataset': 2,
           'angles': None,
           'cells': None,
           'size': 0.1,
           'Rv': 0.05,
           'alpha': 45,
           'path': 'result1_Geo.pkl',
           'query': 'query1.pkl',
           'method': 'Geo'}
    algorithm(arg)
