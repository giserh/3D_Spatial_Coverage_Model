import sys, time
import pandas
from csg.core import CSG
from csg.geom import Vertex, Vector, Polygon
from scipy.spatial.transform import Rotation
import numpy as np
import pickle
import pandas as pd
from scipy.spatial import ConvexHull
# from meshReduction import reduce
from scipy.spatial import KDTree

def print_poly(obj):
    print("poly: ")
    polys = obj.toVerticesAndPolygons()
    for poly in polys:
        print(poly)


def get_fov_polyhedron_from_pts(pts):

    pts = [Vertex(p) for p in pts]

    face_index = [[0, 1, 2], [3, 0, 2], [4, 3, 2], [1, 4, 2], [1, 0, 4], [4, 0, 3]] # order is changed by pycsg
    faces = []
    for face in face_index:
        faces.append(Polygon([pts[i].clone() for i in face]))
    polyhedron = CSG.fromPolygons(faces)

    return polyhedron


def get_fov_polyhedron(pos, theta, alpha=45, R=0.1):

    # get polyhedron from pos and theta

    pts = [Vertex([0, 0, 0])]
    dx = [-1, 1]
    dz = [-1, 1]

    for i in range(2):
        for j in range(2):
            pts.append(Vertex([
                dz[j] * R * np.sin(alpha / 2 * np.pi / 180),

                dx[i] * R * np.sin(alpha / 2 * np.pi / 180),

                               R * np.cos(alpha / 2 * np.pi / 180),
                               ]))

    face_index = [[2, 1, 0],
             [4, 2,  0],
             [3, 4, 0],
             [1, 3, 0],
             [1, 2, 3],
             [3, 2, 4]]
    faces = []
    for face in face_index:
        faces.append(Polygon([pts[i].clone() for i in face]))
    polyhedron = CSG.fromPolygons(faces)

    # finish creating a polyhedron before rotation and translation

    r = Rotation.from_euler('ZYX', theta, degrees=False)
    rotvec = r.as_rotvec()
    angle = np.linalg.norm(rotvec)
    if angle == 0:
        rotvec = [0, 0, 1]
    else:
        rotvec /= angle

    # rotation and translation
    polyhedron.rotate(rotvec, angle * 180 / np.pi)
    polyhedron.rotate([0, 0, 1], 90)
    polyhedron.rotate([0, 1, 0], 180)

    polyhedron.translate(pos)



    # f = open("fovs.csv", "a")
    # polygons = polygon.toPolygons()
    # f.write(str(polygons[0].vertices[2].pos[0]) + ','
    #         + str(polygons[0].vertices[2].pos[1]) + ','
    #         + str(polygons[0].vertices[2].pos[2]) + ',')
    # f.write(str(polygons[0].vertices[1].pos[0]) + ','
    #         + str(polygons[0].vertices[1].pos[1]) + ','
    #         + str(polygons[0].vertices[1].pos[2]) + ',')
    # f.write(str(polygons[0].vertices[0].pos[0]) + ','
    #         + str(polygons[0].vertices[0].pos[1]) + ','
    #         + str(polygons[0].vertices[0].pos[2]) + ',')
    # f.write(str(polygons[2].vertices[0].pos[0]) + ','
    #         + str(polygons[2].vertices[0].pos[1]) + ','
    #         + str(polygons[2].vertices[0].pos[2]) + ',')
    # f.write(str(polygons[2].vertices[1].pos[0]) + ','
    #         + str(polygons[2].vertices[1].pos[1]) + ','
    #         + str(polygons[2].vertices[1].pos[2]) + '\n')

    return polyhedron

# Creates rectangular polygon for a range query
def get_target_object_polygon(min_x, min_y, min_z, max_x, max_y, max_z):

    res = CSG.cube([(min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2],
             [(max_x - min_x) / 2, (max_y - min_y) / 2, (max_z - min_z) / 2])
    return res

def volume(obj):
    # compute volume of polyhedron
    polygons = obj.toPolygons()
    v = 0
    for polygon in polygons:
        n = len(polygon.vertices)
        signs = []
        for i in range(1, n - 1):
            mat = np.matrix([polygon.vertices[0].pos,
                                    polygon.vertices[i].pos,
                                    polygon.vertices[i + 1].pos])
            # print(mat)
            det = np.linalg.det(mat)
            signs.append(np.sign(det))
            v += det
            # print(np.linalg.det(mat))
        if abs(np.mean(signs)) != 1:
            print("Wrong")
            print(signs)
    return v / 6

def sequence_union(fov_polyhedrons):

    n = len(fov_polyhedrons)
    if n == 0: return CSG.fromPolygons([])
    if n == 1: return fov_polyhedrons[0]

    unified_fov_polyhedron = fov_polyhedrons[0]
    for i in range(1, n):
        unified_fov_polyhedron = unified_fov_polyhedron.union(fov_polyhedrons[i])
    return unified_fov_polyhedron

#
# def sequence_union_reduction(fov_polyhedrons, thre, ratio):
#
#     n = len(fov_polyhedrons)
#     if n == 0: return CSG.fromPolygons([])
#     if n == 1: return fov_polyhedrons[0]
#
#     unified_fov_polyhedron = fov_polyhedrons[0]
#     for i in range(1, n):
#         if i % 10 == 0:
#             print("%d of %d" % (i, n))
#         unified_fov_polyhedron = unified_fov_polyhedron.union(fov_polyhedrons[i])
#         v, f, _ = unified_fov_polyhedron.toVerticesAndPolygons()
#         if len(v) > thre:
#             unified_fov_polyhedron = reduce(unified_fov_polyhedron, ratio)
#     return unified_fov_polyhedron


def cascaded_union(fov_polyhedrons):

    # using divide and conquer to get union

    n = len(fov_polyhedrons)
    if n == 0: return CSG.fromPolygons([])
    if n == 1: return fov_polyhedrons[0]
    left = cascaded_union(fov_polyhedrons[:n // 2])
    right = cascaded_union(fov_polyhedrons[n // 2:])
    unified_fov_polyhedron = left.union(right)

    return unified_fov_polyhedron

def execute(region_query_polygon, fov_polyhedrons):
    if len(fov_polyhedrons) == 0:
        return 0
    # Unify all FOV polygons
    unified_fov_polyhedron = cascaded_union(fov_polyhedrons)

    # Get coverage through intersection of target polygon with unified FOV polygons
    region_query_coverage_polygon = region_query_polygon.intersect(unified_fov_polyhedron)

    return (volume(region_query_coverage_polygon)/volume(region_query_polygon))


# def cascaded_union_reduction(fov_polyhedrons, thre, ratio):
#
#     # using divide and conquer to get union
#
#     n = len(fov_polyhedrons)
#     if n == 0: return CSG.fromPolygons([])
#     if n == 1: return fov_polyhedrons[0]
#     left = cascaded_union_reduction(fov_polyhedrons[:n // 2], thre, ratio)
#     right = cascaded_union_reduction(fov_polyhedrons[n // 2:], thre, ratio)
#     unified_fov_polyhedron = left.union(right)
#     v, f, _ = unified_fov_polyhedron.toVerticesAndPolygons()
#     if len(v) > thre:
#         unified_fov_polyhedron = reduce(unified_fov_polyhedron, ratio)
#     return unified_fov_polyhedron
#

# def execute_reduction(region_query_polygon, fov_polyhedrons, thre, ratio):
#     if len(fov_polyhedrons) == 0:
#         return 0
#     # Unify all FOV polygons
#     unified_fov_polyhedron = cascaded_union_reduction(fov_polyhedrons, thre, ratio)
#
#     # Get coverage through intersection of target polygon with unified FOV polygons
#     region_query_coverage_polygon = region_query_polygon.intersect(unified_fov_polyhedron)
#
#     return (volume(region_query_coverage_polygon)/volume(region_query_polygon))


def execute_MC(min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons):
    if len(fov_polyhedrons) == 0:
        return 0
    iter = 1000
    cnt = 0
    volumes = [volume(fov) for fov in fov_polyhedrons]

    for i in range(iter):
        pt = (np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y), np.random.uniform(min_z, max_z))

        for j, fov in enumerate(fov_polyhedrons):
            pts = fov.toVerticesAndPolygons()[0]
            tmp = ConvexHull([pts[0], pts[1], pts[2], pts[3], pts[4], pt])
            if abs(tmp.volume - volumes[j]) < 1e-10:
                cnt += 1
                break

    return cnt / iter


def execute_MC2(min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons, iter=1000):
    if len(fov_polyhedrons) == 0:
        return 0

    cnt = 0
    polys = []
    for fov in fov_polyhedrons:
        tmp = []
        for poly in fov.toPolygons():
            tmp.append(np.array([[pt.pos.x, pt.pos.y, pt.pos.z] for pt in poly.vertices]))
        polys.append(tmp)

    for i in range(iter):
        pt = np.array((np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y), np.random.uniform(min_z, max_z)))

        for fov in polys:
            flag = True
            ck = 0
            for poly in fov:
                tmp = np.linalg.det([poly[0] - pt, poly[1] - pt, poly[2] - pt])
                ck += tmp
                if tmp < 0:
                    flag = False
                    break

            if flag:
                cnt += 1
                break


    return cnt / iter



def execute_approx(min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons):

    # calculate an approximation of the coverage percentage as an upper bound

    div = 15
    xs = [min_x + i * (max_x - min_x) / div for i in range(div + 1)]
    ys = [min_y + i * (max_y - min_y) / div for i in range(div + 1)]
    zs = [min_z + i * (max_z - min_z) / div for i in range(div + 1)]

    whole = get_target_object_polygon(min_x, min_y, min_z, max_x, max_y, max_z)
    required_polys = []
    for polygon in fov_polyhedrons:
        if len(whole.intersect(polygon).toPolygons()) != 0:
            required_polys.append(polygon)
    # divede the whole cube into div * div * div cells

    cnt = 0
    for i in range(div):
        for j in range(div):
            for k in range(div):
                small = get_target_object_polygon(
                    xs[i], ys[j], zs[k], xs[i + 1], ys[j + 1], zs[k + 1])
                for polygon in required_polys:
                    if len(small.intersect(polygon).toPolygons()) != 0:
                        cnt += 1
                        break
                # a cell intersected by any fov is considered as covered

    return cnt * 1.0 / div ** 3



def read_data(st, ed, jump, R=0.1):
    # data = pandas.read_csv('extracted-data.csv')
    data = pandas.read_csv('../dataset_small.csv')


    data = data[st:ed:jump]
    data = data.reset_index(drop=True)

    N = len(data)
    print(data)
    # choose 1 fov every # jump

    # min_lat = np.min(data['GPS_Lat'])
    # max_lat = np.max(data['GPS_Lat'])
    # min_lng = np.min(data['GPS_Lon'])
    # max_lng = np.max(data['GPS_Lon'])
    # min_hgt = np.min(data['GPS_Alt'])
    # max_hgt = np.max(data['GPS_Alt'])
    min_lat = np.min(data['lat'])
    max_lat = np.max(data['lat'])
    min_lng = np.min(data['lng'])
    max_lng = np.max(data['lng'])
    min_hgt = np.min(data['hgt'])
    max_hgt = np.max(data['hgt'])

    # import matplotlib.pyplot as plt
    #
    # plt.subplot(3, 1, 1)
    # plt.hist(data['ATT_Yaw'])
    # plt.title('yaw')
    # plt.subplot(3, 1, 2)
    # plt.hist(data['ATT_Pitch'], bins=30)
    # plt.title('pitch')
    # plt.subplot(3, 1, 3)
    # plt.hist(data['ATT_Roll'], bins=30)
    # plt.title('roll')
    # plt.savefig('yaw_pitch_roll.jpg')
    # # compute the range of 3D pos

    lat_scale = 40075.017 / 360
    lng_scale = 40075.017 / 360 * np.cos((min_lat + max_lat) / 2 * np.pi / 180)
    hgt_scale = 0.001

    # print("range of latitude: %lf to %lf, %lf in total" % (min_lat, max_lat, max_lat - min_lat))
    # print("%lf in meters" % ((max_lat - min_lat) * lat_scale * 1000))
    # print("range of longitude: %lf to %lf, %lf in total" % (min_lng, max_lng, max_lng - min_lng))
    # print("%lf in meters" % ((max_lng - min_lng) * lng_scale * 1000))
    # print("range of altitude: %lf to %lf, %lf in total" % (min_hgt, max_hgt, max_hgt - min_hgt))
    #
    # print("range of pitch: %lf to %lf in radian, ranges %lf percent of pi"
    #       % (np.min(data['ATT_Pitch']), np.max(data['ATT_Pitch']),
    #          (np.max(data['ATT_Pitch']) - np.min(data['ATT_Pitch'])) / np.pi * 100))
    # print("range of roll: %lf to %lf in radian, ranges %lf percent of pi"
    #       % (np.min(data['ATT_Roll']), np.max(data['ATT_Roll']),
    #          (np.max(data['ATT_Roll']) - np.min(data['ATT_Roll'])) / np.pi * 100))


    # compute the scale factor into kilometer


    x_center = (min_lng + max_lng) / 2 * lng_scale
    y_center = (min_lat + max_lat) / 2 * lat_scale
    z_center = (min_hgt + max_hgt) / 2 * hgt_scale

    # center of the region

    x_range = [min_lng * lng_scale - x_center, max_lng * lng_scale - x_center]
    y_range = [min_lat * lat_scale - y_center, max_lat * lat_scale - y_center]
    z_range = [min_hgt * hgt_scale - z_center, max_hgt * hgt_scale - z_center]

    # range of the region
    # Note: set the origin to the center


    fov_pts = [None for i in range(N)]
    fov_angles = [tuple() for i in range(N)]


    for i in range(N):
        if i % 1000 == 0:
            print("Processing %d" % i)
        fov = get_fov_polyhedron([data['lng'][i] * lng_scale - x_center,
                                             data['lat'][i] * lat_scale - y_center,
                                             data['hgt'][i] * hgt_scale - z_center],
                                            [data['yaw'][i],
                                             data['pitch'][i],
                                             data['roll'][i]], 45, R)

        pts = fov.toVerticesAndPolygons()[0]
        # print(fov.toVerticesAndPolygons()[1])
        fov_pts[i] = pts

        fov_angles[i] = (data['yaw'][i], data['pitch'][i], data['roll'][i])

    # get a list of polygons


    scales = [lat_scale, lng_scale, hgt_scale]

    return fov_pts, fov_angles, x_range, y_range, z_range, scales


def gen_cube(x_range, y_range, z_range, size):

    # generate a random cube with the range

    min_x = np.random.uniform(x_range[0], x_range[1]) - size / 2
    max_x = min_x + size
    min_y = np.random.uniform(y_range[0], y_range[1]) - size / 2
    max_y = min_y + size
    min_z = np.random.uniform(z_range[0], z_range[1]) - size / 2
    max_z = min_z + size
    return min_x, min_y, min_z, max_x, max_y, max_z


def gen_cube_to_file(x_range, y_range, z_range, size, num, filename):
    cubes = np.array([gen_cube(x_range, y_range, z_range, size) for i in range(num)])
    data = {'min_x': cubes[:, 0], 'min_y': cubes[:, 1],
            'min_z': cubes[:, 2], 'max_x': cubes[:, 3],
            'max_y': cubes[:, 4], 'max_z': cubes[:, 5]}
    df = pd.DataFrame(data)
    df.to_csv(filename)

if __name__ == '__main__':
    fov_pts, fov_angles, x_range, y_range, z_range, scales = read_data(None, None, 1, R=0.05)
    pickle.dump((fov_pts, fov_angles, x_range, y_range, z_range, scales), open("data_all_2.pkl", "wb"))
    # fov_pts, fov_angles, x_range, y_range, z_range, scales = pickle.load(open("data_all_1.pkl", "rb"))
    # pts = [pt[0] for pt in fov_pts]
    # pts = np.array(pts)
    # print(pts.shape)
    # kdtree = KDTree(pts)
    # pickle.dump(kdtree, open("kdtree_all_1.pkl", "wb"))


# min_x, max_x = x_range
# min_y, max_y = y_range
# min_z, max_z = z_range


# print(execute(min_x, min_y, min_z, max_x, max_y, max_z, fov_polyhedrons))