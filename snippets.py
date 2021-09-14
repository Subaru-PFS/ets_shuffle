# Various potentially useful, unsorted information related to ets_shuffle

# Some handy functions for exploring a DB server
def getColumns(catalog):
    from astroquery.utils.tap.core import TapPlus
    tap = TapPlus(url="http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap") 
    query_string = """SELECT Top 1 * FROM {}""".format(catalog)
    job = tap.launch_job_async(query_string, dump_to_file=False)
    res = job.get_results()
    return res.columns


def VizierTables():
    from astroquery.utils.tap.core import TapPlus
    tap = TapPlus(url="http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap")
    tables = tap.load_tables(only_names=True)
    interestingNames = ("sdss", "2mass", "gaia", "tess", "usno")
    for table in tables:
        name = table.get_qualified_name().lower()
        for iname in interestingNames:
            if iname in name:
                print(table.get_qualified_name())

#snippets for coordinate conversion
# import astropy.units as u
# from astropy.coordinates import SkyCoord
#     crd = SkyCoord(ra*u.deg, dec*u.deg, frame="fk5")
#     test = crd.icrs
#     print(crd)
#     print(test)


# interesting catalogs
# viz7."J/A+A/625/L13/tessbpic"
# viz7."J/A+A/544/A81/snsdss"
# viz7."J/ApJ/780/92/2mass"
# viz7."J/A+A/625/A87/110gaia"
# viz7."J/A+A/551/A78/2mass"
# viz7."J/PAZh/44/124/sasdssbc"
# viz7."J/AJ/157/235/keptess"
# viz7."J/AJ/143/52/sdss_cl"
# viz7."J/ApJS/231/1/sdss"
# viz7."J/A+A/611/A68/SDSS2"
# viz7."J/A+A/611/A68/Gaia"
# viz7."J/ApJ/780/92/sdss"
# viz7."J/A+A/633/A98/gaia"
# viz7."J/A+A/597/A134/psb_sdss"
# viz7."J/MNRAS/406/1595/sdssh"
# viz7."J/A+A/608/A148/x2mass"
# viz7."J/ApJ/749/10/SDSS-obs"
# viz7."J/A+A/569/A124/guv_sdss"
# viz7."J/A+A/618/A20/TESSq"
# viz7."J/MNRAS/482/4570/gaia2wd"
# viz7."J/A+A/611/A68/SDSS1"
# viz7."J/A+A/608/A148/origaia"
# viz7."J/A+A/625/A87/129gaia"
# viz7."J/A+A/618/A20/TESSa"
# viz7."J/MNRAS/482/4570/gaiasdss"
# viz7."J/MNRAS/482/4570/gaiasdss"
# viz7."V/150/notessup"
# viz7."J/A+A/625/A87/122gaia"
# viz7."J/A+A/633/A30/tess"
# viz7."J/A+A/625/A87/116gaia"
# viz7."J/PAZh/44/124/sdss-bc"
# viz6."J/A+A/463/175/2mass"
# viz6."J/ApJS/177/103/sdss"
# vizA."J/A+A/532/A74/sdss"
# vizA."J/A+A/499/395/sdss1001"
# vizA."J/ApJ/703/L72/2mass"
# vizA."J/A+A/499/395/sdss0903"
# vizA."J/ApJ/699/800/sdss"
# vizA."J/A+A/499/395/sdss1335"
# vizA."J/A+A/499/395/sdss1353"
# vizA."J/A+A/499/395/sdss1206"
# vizls."V/139/sdss9"
# vizls."II/281/2mass6x"
# vizls."I/337/gaia"
# vizls."V/147/sdss12"
# vizls."II/306/sdss8"
# vizls."I/345/gaia2"
# vizls."I/347/gaia2dis"
# viz5."J/A+A/446/551/iso2mass"
