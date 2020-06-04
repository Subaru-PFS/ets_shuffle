"""HETDEX shot shuffling.
by Maximilian Fabricius 2011
"""
from datetime import datetime
import logging
import os
import sys
import time

import numpy

import pyhetdex.het.fplane as fplane_parser

from . import parang, visualize
from .handle_local_catalogue import LocalCatalogue
from . import findStars
from distutils.dir_util import mkpath


failure_code_dict = {-1: "Insufficient number of calibration stars.",
                     -2: "Insufficient number of guider or WFS stars.",
                     -3: "Bright star."}


class shotinfo(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []
        s.append("# "+str(datetime.now()))
        s.append("#ID       \t shot id")
        s.append("#RA       \t shot RA[deg]")
        s.append("#DEC      \t shot DEC[deg]")
        s.append("#PA       \t shot pa[deg]")
        s.append("#RA_new   \t new shot RA[deg]")
        s.append("#DEC_new  \t new shot DEC[deg]")
        s.append("#D        \t displacement distance")
        s.append("#N_IFU    \t number if calibration stars in IFUs")
        s.append("#FCODE    \t failure code")
        s.append("#RA_gp1   \t guide probe 1 star RA[deg]")
        s.append("#DEC_gp1  \t guide probe 1 star DEC[deg]")
        s.append("#RA_gp2   \t guide probe 2 star RA[deg]")
        s.append("#DEC_gp2  \t guide probe 2 star DEC[deg]")
        s.append("#RA_wfs1  \t wave front sensor 1 star RA[deg]")
        s.append("#DEC_wfs1 \t wave front sensor 1 star DEC[deg]")
        s.append("#RA_wfs2  \t wave front sensor 2 star RA[deg]")
        s.append("#DEC_wfs2 \t wave front sensor 2 star DEC[deg]")
        s.append("#COMMENT  \t comment")
        s.append("#")
        s.append("#%5s %12s %12s %12s %12s %12s %12s %6s %6s %12s %12s %12s"
                 " %12s %12s %12s %12s %12s %-12s" %
                 ("ID", "RA", "DEC", "PA", "RA_new", "DEC_new", "d", "N_IFU",
                  "FCODE", "RA_gp1", "DEC_gp1", "RA_gp2", "DEC_gp2", "RA_wfs1",
                  "DEC_wfs1", "RA_wfs2", "DEC_wfs2", "COMMENT"))
        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeShuffle(cls, f, sid, ra, dec, pa, shot, failure_code,
                     guideWFSSol):
        """Write something to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        dRa, dDec, d, n = shot
        s = (" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %6d %6d" %
             (sid, ra, dec, pa, ra+dRa, dec+dDec, d, n, failure_code))
        for star in guideWFSSol:
            if star is not None:
                s += " %12.6f %12.6f" % (star[2], star[3])
            else:
                s += " %12s %12s" % ("-", "-")
        f.write(s)
        f.write(" \"%s\"" % "" + "\n")
        f.flush()

    @classmethod
    def writeFailed(cls, f, sid, ra, dec, pa, failure_code):
        """Write failures to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = (" %5d %12.6f %12.6f %12.6f %12s %12s %12s %6s %6d" %
             (sid, ra, dec, pa, "-", "-", "-", "-", failure_code))
        for i in range(4):
            s += " %12s %12s" % ("-", "-")

        f.write(s)
        f.write(" \"%s\"" % failure_code_dict[failure_code] + "\n")

        f.flush()


class shotlist(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []

        s.append("# "+str(datetime.now()))
        s.append("#col  1: RA after shuffling (negative RA indicate failed"
                 " shots)")
        s.append("#col  2: DEC after shuffling")
        s.append("#col  3: original RA")
        s.append("#col  4: original DEC")
        s.append("#col  5: parallactic angle")
        s.append("#col  6: observation time (decimal days from survey start,"
                 " i.e 1st day of observing")
        s.append("#col  7: average pupil illumination [m^2]")
        s.append("#col  8: transparency on night of observation")
        s.append("#col  9: seeing on night of observation")
        s.append("#col 10: East (0) or West (1) track")
        s.append("#col 11: star ID for guide probe 1")
        s.append("#col 12: star ID for guide probe 2")
        s.append("#col 13: star ID for WFS 1")
        s.append("#col 14: star ID for WFS 2")
        # s.append("#col 11: extinction from Schlegel et al. / Peek & Graves"
        #           " 2010")

        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeShuffle(cls, f, sid, ra, dec, ra_org, dec_org, PA, obstime, ill,
                     transp, seeing, track, guideWFSSol):
        """Write something to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        sIDs = []
        for s in guideWFSSol:
            if s is None:
                sIDs.append("NA")
            else:
                sIDs.append("%09d%09d" % tuple(s[0:2]))

        s = ("%11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f"
             " %4d %s %s %s %s" % (ra, dec, ra_org, dec_org, PA, obstime, ill,
                                   transp, seeing, track, sIDs[0], sIDs[1],
                                   sIDs[2], sIDs[3]))
        f.write(s + "\n")
        f.flush()


class starCatalog(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []

        s.append("# "+str(datetime.now()))
        s.append("#col  1: SHOTID ")
        s.append("#col  2: ID ")
        s.append("#col  3: RA ")
        s.append("#col  4: DEC ")
        s.append("#col  5: u ")
        s.append("#col  6: g ")
        s.append("#col  7: r ")
        s.append("#col  8: i ")
        s.append("#col  9: z ")
        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeStars(cls, f, sid, stars):
        for star in stars:
            if star is not None:
                s = ("%06d %09d%09d %11.6f %11.6f %7.2f %7.2f %7.2f"
                     " %7.2f %7.2f" % ((sid,) + tuple(star[:9])))
                f.write(s + "\n")
        f.flush()


class ACAMCatalog(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []

        s.append("# "+str(datetime.now()))
        s.append("#col  1: SHOTID ")
        s.append("#col  2: ID ")
        s.append("#col  3: RA ")
        s.append("#col  4: DEC ")
        s.append("#col  4: x ")
        s.append("#col  4: y ")
        s.append("#col  5: u ")
        s.append("#col  6: g ")
        s.append("#col  7: r ")
        s.append("#col  8: i ")
        s.append("#col  9: z ")
        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeStars(cls, f, sid, stars):
        for star in stars:
            if star is not None:
                s = ("%06d %09d%09d %11.6f %11.6f %10.2f %10.2f %7.2f %7.2f"
                     " %7.2f ")
                f.write(s % ((sid,) + tuple(star[:9])) + "\n")
        f.flush()


def fastpickShuffle(args, oo, ra, dec, pa, shuffle_result):
    """Defines the criteria for picking a shuffle.

    Parameters
    ----------
    args : Namespace
        parsed command line arguments
    00 : ?
    ra, dec : ?
    pa : ?
    shuffle_result : ?

    Returns
    -------
    ?
    failure_code : int
        code defined in :data:`failure_code_dict` to indicate a failure or not
    """
    log = logging.getLogger('shuffle')
    config = args.config
    # now pick optimum shuffle based on maximum number of calibration stars and
    # minimum shuffling distance sort by number of ifu calibration stars and by
    # shuffling distance
    # nn = shuffle_result[:, 3]
    # dd = shuffle_result[:, 2]
    tmp_result = numpy.copy(shuffle_result)
    nn = tmp_result[:, 3]
    dd = tmp_result[:, 2]

    # Sungryong added ;; give equal n-star weight, if > numstarThreshold ;;
    # enough cal-star
    # print("nn before = ")
    # print(nn)
    numThres = 12
    ithreshold = numpy.where(nn > numThres)
    nn[ithreshold] = numThres
    # print("nn after = ")
    # print(nn)

    # distance mininum
    # print("dd = ")
    # print(dd)
    # ddmin = dd[ithreshold].min()
    # print("dd minimum for nn > numThres (1arcmin=0.0167, dThreshold=0.013) =
    # ")
    # print(ddmin)

    # ##### End : change the order

    # will sort by number of avail. IFU stars first and by shuffling distance
    # second
    ii = numpy.lexsort((dd, -nn))

    # print "Maximum number of IFU stars found is:", max(nn)
    # print "Searching for guider WF sensor solution ..."
    failure_code = 0
    dRa, dDec, d, n = shuffle_result[ii][0]

    # Sungryong added
    # print("In Pickshuffle 1st, 2nd, 3rd, lexsort result..")
    # print(shuffle_result[ii][0])
    # print(shuffle_result[ii][1])
    # print(shuffle_result[ii][2])
    # print("shuffle_result[ii].shape = ")
    # print(shuffle_result[ii].shape)
    # print("lexsort out >> ii.shape : ")
    # print(ii.shape)
    # print("ii.shape type : %d " % ii.shape)
    numavail = ii.shape[0]

    # ### Safe search
    log.info(">> Normal Search begins ... ")
    isort = 0
    for isort in range(0, numavail):
        dRa, dDec, d, n = shuffle_result[ii][isort]
        # check for bright star
        if findStars.hasBrightStar(ra+dRa, dec+dDec, pa, config, oo):
            failure_code = -3
            # return None, failure_code
        # find a solution for guide a wave front sensing stars
        # print "dRa,dDec:", dRa,dDec
        guideWFSSol = findStars.fastfindGuideWFSSol(ra+dRa, dec+dDec, pa, oo,
                                                    config)

        if guideWFSSol is not None:
            # if there is a solution pick this shuffle
            failure_code = 0
            break
            # return [[dRa, dDec, d, n], guideWFSSol],  failure_code
        else:
            failure_code = -2
            # return None, failure_code

    log.info(">> Normal pickShuffle : isort = %d,  failure_code = %d ",
             isort, failure_code)
    # print("loop pickShuffle = isort failure_code = ")
    # print(isort)
    # print(failure_code)

    # Sungryong added
    distThreshold = 0.0167
    if failure_code == 0 and d < distThreshold:
        log.info("pickShuffle : normal return")
        return [[dRa, dDec, d, n], guideWFSSol],  failure_code
    else:
        return None, failure_code


def shuffle(args, ra, dec, pa, ifu_centers, oo, ifu_size, edge,  radius,
            dstep):
    """Actual shuffle procedure.

    Walks a grid of allowed shuffle positions. At each positions checks how
    many stars fall into IFUs and finally picks a shuffle position based of a
    set of criteria that are defined in :func:`pickShuffle`.

    csc = calibration star candidates
    gsc = guide star candidates
    wsc = wavefront sensor star candidates

    Parameters
    ----------
    args : Namespace
        parsed command line arguments
    ??

    Returns
    -------
    ??
    """
    log = logging.getLogger('shuffle')
    config = args.config

    # find those stars which are suitable candidates for IFU calibration stars
    cal_star_cand = findStars.fastfindStarCand(oo, config, "ifu")
    log.info("Number of suitable IFU calibration in the FOV: %d",
             (len(cal_star_cand)))

    result = []
    rsq = radius**2.
    # walk grid of shuffle positions
    for x in (numpy.arange(-radius, radius+dstep, dstep) *
              numpy.cos(dec/180. * numpy.pi)):
        for y in (numpy.arange(-radius, radius+dstep, dstep)):
            # only consider shuffle distances of less than "radius"
            dsq = x*x + y*y
            if dsq > rsq:
                continue
            # figure out how many stars there are in IFUs
            insiders = findStars.inside(ra+x, dec+y, pa, ifu_centers,
                                        cal_star_cand, ifu_size, edge,
                                        PLOT=False)
            if len(insiders) > 0:
                n = len(insiders[insiders])
            else:
                n = 0
            # Sungryong added
            # print("searching grid : %f, %f, n = %d." % (ra+x, dec+y, n))

            result.append([x, y, numpy.sqrt(dsq), n])

    result = numpy.array(result)
    log.debug("shape of result = %s", result.shape)

    # picked_shuffle, failure_code = pickShuffle(args, oo, ra, dec, pa, result)
    picked_shuffle, failure_code = fastpickShuffle(args, oo, ra, dec, pa,
                                                   result)
    if picked_shuffle is not None:
        # look again which stars actually ended up inside IFUs
        [dRa, dDec, d, n], guideWFSSol = picked_shuffle
        insiders = findStars.inside(ra+dRa, dec+dDec, pa, ifu_centers,
                                    cal_star_cand, ifu_size, edge, PLOT=False)
        return ([dRa, dDec, d, n, guideWFSSol, cal_star_cand[insiders]],
                failure_code)
    else:
        return None, failure_code


def load_fplane_file(config):
    """Load the fplane file, get the x and y postions and convert them into
    degrees

    Parameters
    ----------
    config : :class:`ConfigParser.ConfigParser` instances
        configuration object

    Returns
    -------
    ifu_centers : 2d numpy array
        n * 2 array containing the centers of ifus in the focal plane
    """
    fplane_file = config.get("General", "fplane_file")
    fplane = fplane_parser.FPlane(fplane_file)
    ifu_centers = numpy.array([[ifu.x, ifu.y] for ifu in fplane.ifus])
    ifu_ihmpid = numpy.array([ifu.ifuslot for ifu in fplane.ifus])
    ifu_ids = numpy.array([[ifu.ifuslot] for ifu in fplane.ifus])

    return ifu_centers / 3600., ifu_ihmpid, ifu_ids


def do(args, data, radius=None, start=1, stop=numpy.Inf, orig_loc=None,
       targetID=None):
    """
    Main procedure. Loops over all shots and performs shuffling for them.

    Parameters
    ----------
    args : Namespace
        parsed command line arguments
    data : ndarray
        shot positions
    radius : int
        maximum shift of the shuffling
    start : int, optional
        ??
    stop : int, optional
        ??
    """
    config = args.config
    log = logging.getLogger('shuffle')

    ifu_centers, ifu_ihmpid, ifu_ids = load_fplane_file(config)

    result = []

    stop = min(max(data[:, 0]), stop)

    # timekeeping
    ntotal = stop-start
    count = 0
    mdt = [0]
    mmdt = 0

    localcat_obj = None
    localcat = config.get("General", "localcat")
    if localcat != "None":
        log.debug("A local catalogue will be used from %s", localcat)
        localcat_obj = LocalCatalogue()
        if not localcat_obj.read_catalogue(localcat):
            log.warning("Error: Reading local catalogue failed!")
            localcat_obj = None

    outinfo = config.get("General", "outinfo")
    log.debug("Information on sources will be stored in %s", (outinfo))
    outfile = config.get("General", "outfile")
    log.debug("Shuffled shots will be stored in %s", (outfile))
    outprobestars = config.get("General", "probestarcatalog")
    log.debug("Catalog of guide and WFS stars will be stored in %s",
              (outprobestars))

    outifustars = config.get("General", "ifustarcatalog")
    log.debug("Catalog of guide and WFS stars will be stored in %s",
              (outifustars))
    outacamstars = config.get("General", "acamstarcatalog")
    log.debug("Catalog of guide and WFS stars will be stored in %s",
              (outacamstars))

    with open(outinfo, 'a') as f_info, open(outfile, 'a') as f_out,\
      open(outprobestars, 'a') as f_probestarcat,\
      open(outacamstars, 'a') as f_acamstarcat,\
      open(outifustars, 'a') as f_ifustarcat:

        shotinfo.writeHeader(f_info)
        shotlist.writeHeader(f_out)
        starCatalog.writeHeader(f_probestarcat)
        ACAMCatalog.writeHeader(f_acamstarcat)
        starCatalog.writeHeader(f_ifustarcat)

        for sid in numpy.arange(start, stop+1):
            eta = (mmdt*(ntotal-count))/3600.
            log.debug("#########################################")
            log.debug("Running shot #%d", sid)
            if eta == 0:
                log.debug(" ETA = ? h")
            elif eta < 1.:
                log.debug(" ETA = %.1f min", (mmdt*(ntotal-count))/60.)
            else:
                log.debug(" ETA = %.1f h", (mmdt*(ntotal-count))/3600.)
            log.debug("#########################################")
            # timekeeping
            t1 = time.time()

            ######################################
            # here we start with the actual work
            ii = data[:, 0] == sid
            shot = data[ii][0]

            ra, dec = shot[1], shot[2]
            track = shot[7]
            if track == 0:
                log.debug("East Track")
                traj = 'EAST'
                tr = 'E'
            if track == 1:
                log.debug("West Track")
                traj = 'WEST'
                tr = 'W'

            pas = parang.parang(dec, az_input=args.az, verbose=False)
            if pas[0] == 0:
                log.error("No track is available for DEC=%f.", dec)
                sys.exit(0)
            if pas[0] == 1:
                log.warning("Warning: Only one track is available for DEC= %f,"
                            " two tracks are identical.", dec)
            pa = pas[int(track) + 1]
            az = pas[int(track) + 3]

            log.info("RA, DEC, PA, AZ, track = %s, %s, %s, %s, %s",
                     ra, dec, pa, az, track)

            # get sources from sdss or a local file
            if localcat_obj:
                source_name = localcat
            else:
                source_name = config.get("General", "catalog")

            log.info("Retrieving sources from %s ...", (source_name))
            if radius is None:
                radius = config.getfloat("General", "radius")
            boxsize = config.getfloat("General", "dfplane") + 2*radius
            oo, cattype = findStars.getStars(sid, ra, dec, pa, config,
                                             boxsize=boxsize,
                                             catalog=config.get("General",
                                                                "catalog"),
                                             local_catalogue=localcat_obj,
                                             cachedir=args.config.get('directories',
                                                                      'cache'))

            # Remove duplicates in the oo-catalog : 1.5 arcsec dup_radius
            log.debug("All dumped number of nearby stars: %d", len(oo))

            # Eliminate close pairs using "guide1_minsep"
            dd = config.getfloat("MagLimits", "guide1_minsep")
            if dd > 0:
                oo = findStars.removeDuplicates(oo, d=dd)
            log.debug("The final total number of nearby stars: %d", len(oo))
            log.info("Writing out probe star candidates")

            mkpath('gp_wfs')
            pstarlist = findStars.getGuideWFScand(ra, dec, pa, oo, config)
            fname_list = ['g1cand', 'g2cand', 'w1cand', 'w2cand']
            for pstar, fname in zip(pstarlist, fname_list):
                outname = '%s_%s_%s.%s' % (targetID, args.ifuslot, tr, fname)
                f_open = open(os.path.join('gp_wfs', outname), 'w')
                starCatalog.writeHeader(f_open)
                starCatalog.writeStars(f_open, sid, pstar)

            # shuffle
            log.info("Shuffling...")
            dstep = config.getfloat("General", "dstep")
            ifu_size = config.getfloat("General", "ifu_size")
            ifu_edge = config.getfloat("General", "ifu_edge")

            ifuslot_exclude = config.get("General", "ifuslot_exclude")
            ifuslot_exclude = [i.strip() for i in ifuslot_exclude.split(',')]
            shuffle_ifu_centers = []
            for (ifucen, ifuslot) in zip(ifu_centers, ifu_ihmpid):
                if ifuslot not in ifuslot_exclude:
                    shuffle_ifu_centers.append(ifucen)
            shuffle_ifu_centers = numpy.array(shuffle_ifu_centers)

            shuffle_result, failure_code = shuffle(args, ra, dec, pa,
                                                   shuffle_ifu_centers, oo,
                                                   ifu_size=ifu_size,
                                                   edge=ifu_edge,
                                                   radius=radius, dstep=dstep)
            if shuffle_result is None:
                if config.get("General", "catalog") == "SDSSDR9":
                    oo, cattype = findStars.getStars(sid, ra, dec, pa, config,
                                                     boxsize=boxsize,
                                                     catalog="USNOA2",
                                                     local_catalogue=localcat_obj,
                                                     cachedir=args.config.get('directories',
                                                                              'cache'))
                    # Remove duplicates in the oo-catalog :
                    # 1.5 arcsec dup_radius
                    log.info("SDSS Failed, so re-running with USNOA2")
                    log.debug("All dumped number of nearby stars: %d", len(oo))

                    # Eliminate close pairs using "guide1_minsep"
                    dd = config.getfloat("MagLimits", "guide1_minsep")
                    if dd > 0:
                        oo = findStars.removeDuplicates(oo, d=dd)
                    log.debug("The final total number of nearby stars: %d",
                              len(oo))
                    shuffle_result, failure_code = shuffle(args, ra, dec, pa,
                                                           shuffle_ifu_centers, oo,
                                                           ifu_size=ifu_size,
                                                           edge=ifu_edge,
                                                           radius=radius,
                                                           dstep=dstep)

            if shuffle_result is not None:
                dRa, dDec, d, n, guideWFSSol, insiders = shuffle_result
                log.info("dRa, dDec = %.5f, %.5f (D=%.1f\") resulted in %d"
                         " sources.", dRa, dDec, d*3600., n)
                log.info("Number of stars in IFUs: %d", n)
                log.info("IFU stars are:")
                for s in insiders:
                    log.info("   %09d%09d %11.6f %11.6f %7.2f %7.2f %7.2f"
                             " %7.2f %7.2f", *tuple(s[:9]))
                cal_star_cand = findStars.fastfindStarCand(oo, config, "ifu")
                acam_star_cand = findStars.fastfindStarCand(oo, config, "acam")
                result.append([ra+dRa, dec+dDec, pa, n, True])
                if config.getboolean("General", "interactive"):
                    guideWFSSol = visualize.visualizeDS9([ra, dec, pa],
                                                         result[-1],
                                                         ifu_centers,
                                                         ifu_ids, guideWFSSol,
                                                         cal_star_cand,
                                                         config, ifu_ihmpid,
                                                         orig_loc, targetID)
                shotinfo.writeShuffle(f_info, sid, ra, dec, pa,
                                      (dRa, dDec, d, n), failure_code,
                                      guideWFSSol)
                shotlist.writeShuffle(f_out, sid, ra+dRa, dec+dDec, ra, dec,
                                      pa, shot[3], shot[4], shot[5], shot[6],
                                      track, guideWFSSol)
                starCatalog.writeStars(f_probestarcat, sid, guideWFSSol)
                starCatalog.writeStars(f_ifustarcat, sid, insiders)
                # visualization of the solution

                ACAMCatalog.writeStars(f_acamstarcat, sid, acam_star_cand)
                if config.getboolean("General", "visualize"):
                    if not targetID:
                        image_name = ('image%04d_%0.3f_%0.3f_%s.jpg'
                                      % (sid, ra, dec, tr))
                    else:
                        image_name = ('%s_%s_%s.jpg'
                                      % (targetID, args.ifuslot, tr))
                    image_name = os.path.join(config.get('directories',
                                                         'images'),
                                              image_name)
                    visualize.visualize([ra, dec, pa], result[-1], ifu_centers,
                                        ifu_ids, guideWFSSol, cal_star_cand,
                                        config, image_name, ifu_ihmpid)
                if config.getboolean("General", "visualize_ACAM"):
                    if not targetID:
                        image_name = ('image%04d_%0.3f_%0.3f_%s.jpg'
                                      % (sid, ra, dec, tr))
                    else:
                        image_name = ('%s_%s_%s.jpg'
                                      % (targetID, args.ifuslot, tr))
                    image_name = os.path.join(config.get('directories',
                                                         'acam_images'),
                                              image_name)
                    acamstar_pos = visualize.visualize_acam([orig_loc[0],
                                                             orig_loc[1], pa],
                                                            result[-1],
                                                            ifu_centers,
                                                            ifu_ids,
                                                            guideWFSSol,
                                                            acam_star_cand,
                                                            config,
                                                            image_name,
                                                            targetID)
                else:
                    acamstar_pos = None

                if config.getboolean('visualisation', 'visualize_probestars'):

                    image_name = os.path.join(config.get('directories',
                                                         'gp_images'))
                    visualize.visualize_probestars(guideWFSSol, config,
                                                   image_name, targetID, tr)
            else:
                log.error("Shuffle failed: %s",
                          failure_code_dict[failure_code])
                result.append([ra, dec, pa, 0, False])
                shotinfo.writeFailed(f_info, sid, ra, dec, pa, failure_code)
                shotlist.writeShuffle(f_out, sid, ra, dec, ra, dec, pa,
                                      shot[3], shot[4], shot[5], shot[6],
                                      track, [None, None, None, None])
            # end of actual work
            ######################################
            # timekeeping
            t2 = time.time()
            dt = t2 - t1
            # limit to one query per s (SDSS complains otherwise
            while (dt < 1.):
                t2 = time.time()
                dt = t2 - t1

            mdt.append(dt)
            log.debug("time : %f s", dt)
            mmdt = numpy.median(mdt)
            log.debug("mtime : %f s", mmdt)
            # i += 1
            count += 1

    return numpy.array(result)
