import numpy as np
import argparse
from pfs.utils.coordinates.CoordTransp import CoordinateTransform as ctrans
import pfs.datamodel
from ets_shuffle import query_utils
from ets_shuffle.convenience import (flag_close_pairs,
                                     update_coords_for_proper_motion,
                                     plot_focal_plane,
                                     guidecam_geometry)
import ets_fiber_assigner
import matplotlib.path as mppath

# radius of the focal plane: roughly 220mm
# 1mm corresponds roughly to 10.2 arcseconds

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--design_in", type=int, help="input pfsDesign ID")
    parser.add_argument("--design_in_dir", type=str, default=".", help="input pfsDesign directory")
    parser.add_argument("--design_out", type=int, help="output pfsDesign ID")
    parser.add_argument("--design_out_dir", type=str, default=".", help="output pfsDesign directory")
    parser.add_argument("--guidestar_mag_max", type=float, default=19., help="maximum magnitude for guide star candidates")
    parser.add_argument("--guidestar_neighbor_mag_min", type=float, default=21., help="minimum magnitude for objects in the vicinity of guide star candidates")
    parser.add_argument("--guidestar_minsep_deg", type=float, default=1./3600, help="radius of guide star candidate vicinity")
    parser.add_argument("--observation_time", type=str, default="2020-01-01 15:00:00", help="planned time of observation")
    args = parser.parse_args()

    input_design = pfs.datamodel.PfsDesign.read(args.design_in, args.design_in_dir)
    raTel_deg, decTel_deg, pa_deg = input_design.raBoresight, input_design.decBoresight, input_design.posAng

    # this should come from the pfsDesign as well, but is not yet in there
    # (DAMD-101)
    obs_time = args.observation_time

    guidestar_mag_max = args.guidestar_mag_max
    guidestar_neighbor_mag_min = args.guidestar_neighbor_mag_min
    guidestar_minsep_deg = args.guidestar_minsep_deg


    # guide star cam geometries
    agcoord = guidecam_geometry()

    # internal, technical parameters
    # set focal plane radius
    fp_rad_deg = 260. * 10.2/3600

    # Find guide star candidates

    # Query GAIA2 for a circular region containing all guide cam FOVs
    # Obtain all targets with g_mean_mag<=guidestar_neighbor_mag_min that have
    # proper motion information
    conn, table, coldict = query_utils.openGAIA2connection()
    racol, deccol = coldict["ra"], coldict["dec"]
    req_columns = [coldict["id"], racol, deccol, coldict["pmra"],
                   coldict["pmdec"], 'phot_g_mean_mag']
    constraints = [
        query_utils.build_circle_query(
            raTel_deg, decTel_deg, fp_rad_deg*1.2, coldict),
        query_utils.build_pm_query(coldict),
        query_utils.build_mag_query(guidestar_neighbor_mag_min, 0,
                                    'phot_g_mean_mag')]
    res = query_utils.run_query(conn, table, req_columns, constraints)
    # FIXME: run similar query, but without the PM requirement, to get a list of
    # potentially too-bright neighbours

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

    # Write the results to a new pfsDesign file. Data fields are according to
    # DAMD-101.
    # required data:
    # ra/dec of guide star candidates: in racol, deccol
    # PM information: in pmra, pmdec
    # parallax: currently N/A
    # flux: currently N/A
    # AgId: trivial to obtain from data structure
    # AgX, AgY (pixel coordinates): only computable with access to the full
    #   AG camera geometry
    output_design = input_design
    output_design.pfsDesignId = args.design_out
#    [...]
    output_design.write(dirName=args.design_out_dir)


if __name__ == "__main__":
    main()
