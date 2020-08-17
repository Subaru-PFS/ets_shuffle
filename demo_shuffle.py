import numpy as np


# Translation tables for column names
coldicts = {
    "GAIA2" : {
        "ra" : "ra",
        "dec" : "dec",
        "pmra" : "pmra",
        "pmdec" : "pmdec",
        "id" : "source_id" },
    "SDSSdr12" : {
        "ra" : "RA_ICRS",
        "dec" : "DE_ICRS",
        "pmra" : "pmRA",
        "pmdec" : "pmDE",
        "id" : "objID" }
    }


def build_circle_query(ra, dec, radius, coldict):
    return """CONTAINS(POINT('ICRS',{},{}),CIRCLE('ICRS',{},{},{}))=1""".format(coldict["ra"], coldict["dec"], ra, dec, radius)


def build_polygon_query(coords, coldct):
    res = """CONTAINS(POINT('ICRS',{},{}),
          POLYGON('ICRS',""".format(coldict["ra"], coldict["dec"])
    for i in range(coords.shape[0]-1):
        res += str(coords[i,0]) + "," + str(coords[i,1]) + ","
    res += str(coords[-1,0]) + "," + str(coords[-1,1])
    res += "))=1"
    return res


def build_mag_query(minmag, maxmag, mag_col):
    return """{} between {} and {}""".format(mag_col, maxmag, minmag)


def build_pm_query(coldict):
    return """{} is not null and {} is not null""".format(coldict["pmra"], coldict["pmdec"])


def openGAIA2connection():
    from astroquery.gaia import Gaia
    return Gaia, "gaiadr2.gaia_source", coldicts["GAIA2"]


def openVizierSDSSdr12connection():
    from astroquery.utils.tap.core import TapPlus
    tap = TapPlus(url="http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap") 
    return tap, 'vizls."V/147/sdss12"', coldicts["SDSSdr12"]


def openVizierGAIA2connection():
    from astroquery.utils.tap.core import TapPlus
    tap = TapPlus(url="http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap") 
    return tap, 'vizls."I/345/gaia2"', coldicts["GAIA2"]


def run_query(conn, table, req_columns, constraints):
    colstring = ""
    for col in req_columns[:-1]:
        colstring += col+', '
    colstring += req_columns[-1]
    constr = ""
    for constraint in constraints[:-1]:
        constr += constraint+' AND '
    constr += constraints[-1]
    query_string = 'SELECT ' + colstring + ' FROM ' + table + ' WHERE ' + constr + ';'
    print(query_string)
    job = conn.launch_job_async(query_string, dump_to_file=False)
    res = job.get_results()
    return res


def flag_close_pairs(ra, dec, radius):
    ra, dec = ra*np.pi/180, dec*np.pi/180
    rsq = (radius*np.pi/180)**2
    assert(ra.ndim==1 and dec.ndim==1 and ra.shape[0]==dec.shape[0], "bad input")
    nval = ra.shape[0]
    xyz = np.zeros((nval,3), dtype=np.float64)
    cdec = np.cos(dec)
    xyz[..., 0] = cdec*np.cos(ra)
    xyz[..., 1] = cdec*np.sin(ra)
    xyz[..., 2] = np.sin(dec)
    flags = np.zeros(nval, dtype=np.bool)
    for i in range(ra.shape[0]):
        distsq = ((xyz[i,0]-xyz[i+1:,0])**2
                 +(xyz[i,1]-xyz[i+1:,1])**2
                 +(xyz[i,2]-xyz[i+1:,2])**2)
        tflags = distsq < rsq
        if np.any(tflags):
            flags[i] = True
        flags[i+1:] = flags[i+1:] | tflags
    return flags


def update_coords_for_proper_motion(ra, dec, pmra, pmdec, epoch_cat, epoch_now):
    deltaRA = (epoch_now - epoch_cat)*pmra / 3600000. / np.cos(dec*np.pi/180.)
    deltaDE = (epoch_now - epoch_cat)*pmdec / 3600000.
    return ra+deltaRA, dec + deltaDE


if __name__ == "__main__":
    # set focal plane radius
    fp_rad_deg = 1.
    # set some telescope pointing
    raTel_deg, decTel_deg = 34, -3


    # Find bright stars within FOV

    # Query GAIA via their own server
    conn, table, coldict = openGAIA2connection()
    req_columns = [coldict["id"], coldict["ra"], coldict["dec"], coldict["pmra"], coldict["pmdec"], 'phot_g_mean_mag']
    constraints = [build_circle_query(raTel_deg, decTel_deg, fp_rad_deg, coldict),
                   build_pm_query(coldict),
                   build_mag_query(15, 0, 'phot_g_mean_mag')]
    res = run_query(conn, table, req_columns, constraints)
    print(res)

    # Query Vizier's copy of GAIA
    conn, table, coldict = openVizierGAIA2connection()
    req_columns = [coldict["id"], coldict["ra"], coldict["dec"], coldict["pmra"], coldict["pmdec"], 'phot_g_mean_mag']
    constraints = [build_circle_query(raTel_deg, decTel_deg, fp_rad_deg, coldict),
                   build_pm_query(coldict),
                   build_mag_query(15, 0, 'phot_g_mean_mag')]
    res = run_query(conn, table, req_columns, constraints)
    print(res)

    # Query Vizier's copy of SDSS
    conn, table, coldict = openVizierSDSSdr12connection()
    req_columns = [coldict["id"], coldict["ra"], coldict["dec"], coldict["pmra"], coldict["pmdec"], 'gmag']
    constraints = [build_circle_query(raTel_deg, decTel_deg, fp_rad_deg, coldict),
                   build_pm_query(coldict),
                   build_mag_query(15, 0, 'gmag')]
    res = run_query(conn, table, req_columns, constraints)
    print(res)


    # Find guide star candidates

    # Query GAIA2 for a rectangular region
    conn, table, coldict = openGAIA2connection()
    req_columns = [coldict["id"], coldict["ra"], coldict["dec"], coldict["pmra"], coldict["pmdec"], 'phot_g_mean_mag']
    # The AG camera is described by its 4 corners. Let's use a square with
    # 2.5 arcmin side length for now
    # I'm also putting it into the center of the FP because it's just an example
    agcoord = np.zeros((4,2))
    sidelength = 2.5/60
    agcoord[0,:] = [raTel_deg, decTel_deg]
    agcoord[1,:] = [raTel_deg+sidelength, decTel_deg]
    agcoord[2,:] = [raTel_deg+sidelength, decTel_deg+sidelength]
    agcoord[3,:] = [raTel_deg, decTel_deg+sidelength]
    constraints = [build_polygon_query(agcoord, coldict),
                   build_pm_query(coldict),
                   build_mag_query(20, 0, 'phot_g_mean_mag')]
    res = run_query(conn, table, req_columns, constraints)
    flags = flag_close_pairs(res[coldict["ra"]], res[coldict["dec"]], 0.0000001)
    ranew, decnew = update_coords_for_proper_motion(res[coldict["ra"]], res[coldict["dec"]], res[coldict["pmra"]], res[coldict["pmdec"]], 2000., 2020.)
    print(decnew-res[coldict["dec"]])
    print(res)
    print(flags)
