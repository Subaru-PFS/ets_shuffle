import numpy as np
from coordinates.CoordTransp import CoordinateTransform as ctrans
from ets_shuffle import query_utils
import matplotlib.path as mppath


def flag_close_pairs(ra, dec, radius):
    ra, dec = ra*np.pi/180, dec*np.pi/180
    rsq = (radius*np.pi/180)**2
    assert ra.ndim == 1 and dec.ndim == 1 and ra.shape[0] == dec.shape[0],\
        "bad input"
    nval = ra.shape[0]
    xyz = np.zeros((nval, 3), dtype=np.float64)
    cdec = np.cos(dec)
    xyz[..., 0] = cdec*np.cos(ra)
    xyz[..., 1] = cdec*np.sin(ra)
    xyz[..., 2] = np.sin(dec)
    flags = np.zeros(nval, dtype=np.bool)

    for i in range(ra.shape[0]):
        distsq = np.sum((xyz[i]-xyz[i+1:])**2, axis=-1)
        tflags = distsq < rsq
        if np.any(tflags):
            flags[i] = True
        flags[i+1:] = np.logical_or(flags[i+1:], tflags)
    return flags


def update_coords_for_proper_motion(ra, dec, pmra, pmdec, epoch_cat,
                                    epoch_now):
    deltaRA = (epoch_now - epoch_cat)*pmra/3600000./np.cos(dec*np.pi/180.)
    deltaDE = (epoch_now - epoch_cat)*pmdec/3600000.
    return ra+deltaRA, dec+deltaDE


def plot_focal_plane(cameras, targets, gstargets):
    import matplotlib.pyplot as plt
    import matplotlib.path as mppath
    import matplotlib.patches as patches
    fig, ax = plt.subplots()
    for i in range(cameras.shape[0]):
        p = mppath.Path(cameras[i, (0, 1, 2, 3, 0), :], closed=True)
        patch = patches.PathPatch(p, fill=False)
        ax.add_patch(patch)
    ax.set_xlim(-300, 300)
    ax.set_ylim(-300, 300)
    ax.scatter(targets[:, 0], targets[:, 1], color='g', marker=".")
    ax.scatter(gstargets[:, 0], gstargets[:, 1], color='r')
    plt.show()


# returns a (6,4,2) float array with PFI coordinates of the 4 corners
# of the 6 guide star cams.
# This is a placeholder until we know how to obtain the real geometry
# rough values are taken from
# https://sumire.pbworks.com/w/file/58576912/agcamerastudy.pdf
def guidecam_geometry():
    # The AG camera is described by its 4 corners. Let's use a square with
    # 10 mm side length for now, at 25mm distance from the center
    agcoord0 = np.zeros((4, 2))
    dist = 250.
    sidelength = 15.
    agcoord0[0, :] = [-.5*sidelength, -.5*sidelength]
    agcoord0[1, :] = [+.5*sidelength, -.5*sidelength]
    agcoord0[2, :] = [+.5*sidelength, +.5*sidelength]
    agcoord0[3, :] = [-.5*sidelength, +.5*sidelength]
    agcoord0[:, 0] += dist
    # get the remaining cameras via 60-degree rotations
    agcoord = np.zeros((6, 4, 2))
    for i in range(0, 6):
        ang = (30. + 60.*i)*np.pi/180.
        ca, sa = np.cos(ang), np.sin(ang)
        for j in range(4):
            x, y = agcoord0[j, :]
            agcoord[i, j, :] = ca*x - sa*y, ca*y + sa*x
    return agcoord


# radius of the focal plane: roughly 220mm
# 1mm corresponds roughly to 10.2 arcseconds

def main():
    # Section with externally provided parameters

    obs_time = "2020-01-01 15:00:00"
    # telescope pointing
    raTel_deg, decTel_deg = 34., -3.
    # focal plane position angle
    pa_deg = 0.
    # maximum magnitude for guide stars
    guidestar_mag_max = 19
    # maximum magnitude for close neighours of guide stars
    guidestar_neighbor_mag_max = 21
    # minimum distance (in degrees) between guide star candidates
    guidestar_minsep_deg = 1./3600


    # guide star cam geometries
    agcoord = guidecam_geometry()

    # internal, technical parameters
    # set focal plane radius
    fp_rad_deg = 260. * 10.2/3600

    # Find guide star candidates

    # Query GAIA2 for a circular region containing all guide cam FOVs
    # Obtain all targets with g_mean_mag<=guidestar_neighbor_mag_max that have
    # proper motion information
    conn, table, coldict = query_utils.openGAIA2connection()
    racol, deccol = coldict["ra"], coldict["dec"]
    req_columns = [coldict["id"], racol, deccol, coldict["pmra"],
                   coldict["pmdec"], 'phot_g_mean_mag']
    constraints = [
        query_utils.build_circle_query(
            raTel_deg, decTel_deg, fp_rad_deg*1.2, coldict),
        query_utils.build_pm_query(coldict),
        query_utils.build_mag_query(guidestar_neighbor_mag_max, 0,
                                    'phot_g_mean_mag')]
    res = query_utils.run_query(conn, table, req_columns, constraints)
    # adjust for proper motion
    obs_year = float(obs_time[0:4])
    res[racol], res[deccol] = \
        update_coords_for_proper_motion(res[racol], res[deccol],
                                        res[coldict["pmra"]],
                                        res[coldict["pmdec"]], 2000., obs_year)

    # compute PFI coordinates
    tmp = np.array([res[racol], res[deccol]])
    tmp = ctrans(xyin=tmp,
                 za=0., mode="sky_pfi", inr=0., pa=pa_deg,
                 cent=np.array([raTel_deg, decTel_deg]),
                 time=obs_time)
    res["xypos"] = np.array([tmp[0, :], tmp[1, :]]).T

    # determine the subset of sources falling within the guide cam FOVs
    # For the moment I'm using matplotlib's path functionality for this task
    # Once the "pfi_sky" transformation direction is available in
    # pfs_utils.coordinates, we can do a direct polygon query for every camera,
    # which should be more efficient.
    targets = {}
    tgtcam = []
    for i in range(agcoord.shape[0]):
        p = mppath.Path(agcoord[i])
        # find all targets in the slighty enlarged FOV
        tmp = p.contains_points(res["xypos"], radius=1.)  # 1mm more
        tdict = {}
        for key, val in res.items():
            tdict[key] = val[tmp]
        # eliminate close neighbors
        flags = flag_close_pairs(tdict[racol], tdict[deccol],
                                 guidestar_minsep_deg)
        for key, val in tdict.items():
            tdict[key] = val[np.invert(flags)]
        # eliminate all targets which are not bright enough to be guide stars
        flags = tdict["phot_g_mean_mag"] < guidestar_mag_max
        for key, val in tdict.items():
            tdict[key] = val[flags]
        # eliminate all targets which are not really in the camera's FOV
        flags = p.contains_points(tdict["xypos"])  # 1mm more
        for key, val in tdict.items():
            tdict[key] = val[flags]
        # append the results for this camera to the full list
        tgtcam.append(tdict)
        for key, val in tdict.items():
            if key not in targets:
                targets[key] = val
            else:
                targets[key] = np.concatenate((targets[key], val))

    for i, d in enumerate(tgtcam):
        print("AG camera #{}".format(i))
        print(d[coldict["id"]])
    plot_focal_plane(agcoord, res["xypos"], targets["xypos"])


if __name__ == "__main__":
    main()
