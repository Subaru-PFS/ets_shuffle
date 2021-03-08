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
    agcoord0 = np.zeros((4, 2))
    dist = 241.292  # mm
    sidelength = 13.312  # mm
    # not yet used : npix = 1024  # pixels along each direction
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
