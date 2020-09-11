# Translation tables for column names in different catalogues
coldicts = {
    "GAIA2": {
        "ra": "ra",
        "dec": "dec",
        "pmra": "pmra",
        "pmdec": "pmdec",
        "id": "source_id"},
    "SDSSdr12": {
        "ra": "RA_ICRS",
        "dec": "DE_ICRS",
        "pmra": "pmRA",
        "pmdec": "pmDE",
        "id": "objID"}
    }


def build_circle_query(ra, dec, radius, coldict):
    """Returns a query substring corresponding to a cone search of radius rad
    around a given ra, dec.

    Parameters
    ----------
    ra, dec, radius : float
        RA, DEC and cone radius in degrees
    coldict : dict(string, string)
        the translation table for the catalogue in use

    Returns
    -------
    string : the query substring
    """
    return """CONTAINS(POINT('ICRS',{},{}),CIRCLE('ICRS',{},{},{}))=1"""\
        .format(coldict["ra"], coldict["dec"], ra, dec, radius)


def build_polygon_query(coords, coldict):
    """Returns a query substring corresponding to a polygon query.

    Parameters
    ----------
    coords : numpy.ndarray(float) of shape (*,2)
        RA, DEC of the polygon vertices in deg
    coldict : dict(string, string)
        the translation table for the catalogue in use

    Returns
    -------
    string : the query substring
    """
    res = """CONTAINS(POINT('ICRS',{},{}),
          POLYGON('ICRS',""".format(coldict["ra"], coldict["dec"])
    for i in range(coords.shape[0]-1):
        res += str(coords[i, 0]) + "," + str(coords[i, 1]) + ","
    res += str(coords[-1, 0]) + "," + str(coords[-1, 1])
    res += "))=1"
    return res


def build_mag_query(minmag, maxmag, mag_col):
    """Returns a query substring corresponding to a magnitude selection for
    a given column.

    Parameters
    ----------
    minmag : float
        the magnitude corresponding to the faintest acceptable target
    maxmag : float
        the magnitude corresponding to the brigutest acceptable target
    mag_col : string
        the column name to use for the selection

    Returns
    -------
    string : the query substring
    """
    return """{} between {} and {}""".format(mag_col, maxmag, minmag)


def build_pm_query(coldict):
    """Returns a query substring selecting only targets that have proper motion
    information.

    Parameters
    ----------
    coldict : dict(string, string)
        the translation table for the catalogue in use

    Returns
    -------
    string : the query substring
    """
    return """{} is not null and {} is not null"""\
        .format(coldict["pmra"], coldict["pmdec"])


def openGAIA2connection():
    """Opens a connection to the GAIA2 server and returns relevant information

    Returns
    -------
    object : the database connection object
    string : the name of the table being queried
    dict(string, string) : the relevant translation table for column names
    """
    from astroquery.gaia import Gaia
    return Gaia, "gaiadr2.gaia_source", coldicts["GAIA2"]


def openVizierSDSSdr12connection():
    """Opens a connection to the SDSS DR1.2 table at VIZIER

    Returns
    -------
    object : the database connection object
    string : the name of the table being queried
    dict(string, string) : the relevant translation table for column names
    """
    from astroquery.utils.tap.core import TapPlus
    tap = TapPlus(url="http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap",
                  verbose=False)
    return tap, 'vizls."V/147/sdss12"', coldicts["SDSSdr12"]


def openVizierGAIA2connection():
    """Opens a connection to the GAIA2 table at VIZIER

    Returns
    -------
    object : the database connection object
    string : the name of the table being queried
    dict(string, string) : the relevant translation table for column names
    """
    from astroquery.utils.tap.core import TapPlus
    tap = TapPlus(url="http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap",
                  verbose=False)
    return tap, 'vizls."I/345/gaia2"', coldicts["GAIA2"]


def run_query(conn, table, req_columns, constraints):
    """Runs the provided query and returns the results

    Parameters
    ----------
    conn : object
        the database connection object
    table : string
        the table to be used
    req_columns : list[string]
        the required columns
    contraints : string
        the query string contraining geometry, magnitude, etc.

    Returns
    -------
    dict(string, numpy.ndarray) : the retrieved quantities
    """
    import numpy as np
    colstring = ""
    for col in req_columns[:-1]:
        colstring += col+', '
    colstring += req_columns[-1]
    constr = ""
    for constraint in constraints[:-1]:
        constr += constraint+' AND '
    constr += constraints[-1]
    query_string = ('SELECT ' + colstring + ' FROM ' + table
                    + ' WHERE ' + constr + ';')
#    print(query_string)
    job = conn.launch_job_async(query_string, dump_to_file=False, verbose=False)
    res = job.get_results()
    res2 = {}
    for key in res.keys():
        res2[key] = np.array(res[key])
    return res2
