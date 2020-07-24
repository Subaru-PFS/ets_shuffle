import numpy as np
from ics.cobraOps.Bench import Bench
import ets_fiber_assigner.netflow as nf

# define locations of the input files
catalog_path = "data/"
fscience_targets = catalog_path+"pfs_preliminary_target_cosmology.dat"

def queryGAIA2(ra, dec, radius, minmag, maxmag, band):
    from astroquery.gaia import Gaia 
    if band != "g":
        raise NotImplementedError

    query_string = """SELECT source_id, ra, dec, phot_g_mean_mag, pmra, pmdec
        FROM gaiadr2.gaia_source
        WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{},{},{}))=1
        AND pmra IS NOT NULL
        AND pmdec IS NOT NULL
        AND phot_g_mean_mag<{}
        AND phot_g_mean_mag>{};""".format(ra, dec, radius, minmag, maxmag)
    job = Gaia.launch_job_async(query_string, dump_to_file=False)
    res = job.get_results()
    return res

def queryGAIA2_poly(corners, minmag, maxmag, band):
    from astroquery.gaia import Gaia 
    if band != "g":
        raise NotImplementedError
    print(*(corners.reshape((-1,))))
    query_string = """SELECT source_id, ra, dec, phot_g_mean_mag, pmra, pmdec
        FROM gaiadr2.gaia_source
        WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),
          POLYGON('ICRS',{},{},{},{},{},{},{},{}))=1
        AND pmra IS NOT NULL
        AND pmdec IS NOT NULL
        AND phot_g_mean_mag<{}
        AND phot_g_mean_mag>{};""".format(*(corners.reshape((-1,))), minmag, maxmag)
    print(query_string)
    job = Gaia.launch_job_async(query_string, dump_to_file=False)
    res = job.get_results()
    return res

if __name__ == "__main__":
    bench = Bench(layout="full")
    fp_rad_mm = bench.radius
    # TODO missing functionality: compute FP radius on the sky
#    fp_rad_deg = ???
    fp_rad_deg = 1.
    # determine a telescope pointing, for now just pointing at the center
    # of the prvided target catalog
    raTel_deg, decTel_deg = nf.telescopeRaDecFromFile(fscience_targets)
    # find bright stars within FOV
    res = queryGAIA2(raTel_deg, decTel_deg, fp_rad_deg, 15, -10, "g")
    print(res)
    # The AG camera is described by its 4 corners. Let's use a square with
    # 2.5 arcmin side length for now
    agcoord = np.zeros((4,2))
    sidelength = 2.5/60
    agcoord[0,:] = [raTel_deg, decTel_deg]
    agcoord[1,:] = [raTel_deg+sidelength, decTel_deg]
    agcoord[2,:] = [raTel_deg+sidelength, decTel_deg+sidelength]
    agcoord[3,:] = [raTel_deg, decTel_deg+sidelength]
    res = queryGAIA2_poly(agcoord, 20, -10, "g")
    print(res)
