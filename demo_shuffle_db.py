import numpy as np
from pfs.utils.coordinates.CoordTransp import CoordinateTransform as ctrans
from ets_shuffle import query_utils
from ets_shuffle.convenience import (flag_close_pairs,
                                     update_coords_for_proper_motion,
                                     plot_focal_plane,
                                     guidecam_geometry)
import matplotlib.path as mppath


def get_db():
    from opdb import opdb
    hostname='localhost'
    port='5432'
    dbname='opdb'
    username='pfs'
    passwd='2394'

    db = opdb.OpDB(hostname, port, dbname, username, passwd)
    db.connect()
    return db


# 1mm on the focal plane corresponds roughly to 10.2 arcseconds

def main():
    # establish DB connection
    db = get_db()
    # Requested PFS design ID, will be passed as a module parameter later on
    pfs_design_id = -7188705301062151535
    # Get parameters related to this task
    df = db.fetch_by_id('pfs_design', pfs_design_id=pfs_design_id)
    print(df)
    obs_time = str(df['to_be_observed_at'][0])
    print(obs_time,type(obs_time))
    raTel_deg = df['ra_center_designed'][0]
    decTel_deg = df['dec_center_designed'][0]
    pa_deg = df['pa_designed'][0]
    print(raTel_deg, decTel_deg, pa_deg)

    # Section with externally provided parameters

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
