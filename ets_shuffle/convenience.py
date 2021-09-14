import numpy as np


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
# data taken from email by Yuki Moritani, Apr 6 2020
def guidecam_geometry():
    import pfs.utils.coordinates.CoordTransp as ctrans
    agcoord = np.zeros((6, 4, 2))
    for i in range(0, 6):
        agcoord[i, 0, :] = ctrans.ag_pixel_to_pfimm(i, 0.5, 0.5)
        agcoord[i, 1, :] = ctrans.ag_pixel_to_pfimm(i, 1023.5, 0.5)
        agcoord[i, 2, :] = ctrans.ag_pixel_to_pfimm(i, 1023.5, 1023.5)
        agcoord[i, 3, :] = ctrans.ag_pixel_to_pfimm(i, 0.5, 1023.5)
    return agcoord
