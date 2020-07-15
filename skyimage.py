import numpy
from urllib.request import urlopen
from urllib.parse import urlencode
import matplotlib.pyplot as plt
from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO
from astroquery.skyview import SkyView
import astropy.units as aunit


__all__ = ["retrieve_image"]


def retrieve_image(ra, dec, size, res, logger, options=None):
    """Wrapper function for retrieving image from SDSS. If region outside SDSS
    coverage, it uses DSS image instead.

    Parameters
    ----------
    ra, dec : float
        Right ascension and declination in degrees
    size : float
        1D extent of the image in degrees
    res : int
        image resolution
        the returned array will normally have res*res pixels, but this may
        be adjusted if the desired size cannot be delivered for some reason
    logger : logging object
        will be used for log messages
    options : string (optional)
        additional options to pass to the image query

    Returns
    -------
    image array, scaling matrix, URL, data source
    """

    def SDSS_coverage(ra, dec, logger):

        url_sdssCoverage = 'http://www.sdss3.org/dr9/index.php'
        request_sdssCoverage = urlencode({'coverageRA': ra,
                                          'coverageDec': dec})

        try:
            for line in urlopen(url_sdssCoverage, request_sdssCoverage.encode()):
                if b'overlaps with the SDSS DR9 survey area.' in line:
                    return True
        except IOError:
            logger.error("Could not connect to %s", url_sdssCoverage)

        return False

    def blank_image(size_pix, scale):
        imarray = numpy.zeros((size_pix, size_pix))
        CD = numpy.matrix([[-1.*scale/3600., 0], [0, 1.*scale/3600.]])
        return imarray, CD, '', ''


    def retrieve_image_PS1(ra, dec, size, scale, logger):
        """
        Load an image from the PS1 image cut out service. Choose
        a scale that is a multiple of 0.25. Bin image up to match
        the input scale.
        """

        url_ps1_fn = "http://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url_ps1_fitscut = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
        request = urlencode({'ra': ra, 'dec': dec, 'filters': 'gri'})

        # image params
        scale_factor = int(scale/0.25)
        size_pix = int(size*3600./scale)*scale_factor

        if scale < 0.25:
            raise ValueError("Please choose a scale thats a multiple of 0.25")

        # first get the names of the images we need
        try:
            imlist = urlopen(url_ps1_fn, request.encode())
        except IOError:
            logger.warning("HTTP ERROR request of: %s %s", url_ps1_fn, request)
            return blank_image(size_pix/scale_factor, scale)

        # produce a dictionary of filters and filenames
        imdict = {}
        next(imlist)
        for line in imlist:
            els = line.split()
            if len(els) < 7:
                logger.warning("Query failed!")
                return blank_image(size_pix/scale_factor, scale)
            imdict[els[4]] = els[7]

        # make the image cutout
        try:
            # just one band as otherwise too much noise to see overplot
            opts = {'Green': imdict['r'], 'Red': imdict['i'], 'Blue': imdict['g'],
                    'Size': size_pix, 'Ra': ra, 'Dec': dec, 'AutoscaleMax': 99.5,
                    'AutoscaleMin': 50.0}
        except KeyError:
            logger.warning("Required PS1 images not found!")
            return blank_image(int(size_pix/scale_factor), scale)

        request_cutout = urlencode(opts)

        # get image
        try:
            imfile = urlopen(url_ps1_fitscut, request_cutout.encode())
        except IOError:
            logger.warning("HTTP ERROR request of: %s %s", url_ps1_fitscut,
                        request_cutout)
            return blank_image(size_pix/scale_factor, scale)

        try:
            imarray = plt.imread(BytesIO(imfile.read()), format='jpeg')
        except IOError:
            logger.warning("Image acquisition failed!")
            return blank_image(size_pix/scale_factor, scale)

        # bin it up
        rgb_arr = []
        nx, ny, nc = imarray.shape
        for c in range(nc):
            channel = imarray[:, :, c]
            nx, ny = channel.shape
            # nsize = int(nx/scale_factor)
            im_view = channel.reshape(int(nx/scale_factor), scale_factor,
                                      int(ny/scale_factor), scale_factor)
            binned = im_view.mean(axis=3).mean(axis=1)
            rgb_arr.append(numpy.array(255*(1.0 -
                                            binned/numpy.max(numpy.array(binned))),
                           dtype=int))

        # stick binned up channels into single RGB image
        rgb_binned = numpy.zeros((rgb_arr[0].shape[0], rgb_arr[0].shape[1], 3),
                                 dtype=int)
        rgb_binned[..., 0] = rgb_arr[0]
        rgb_binned[..., 1] = rgb_arr[1]
        rgb_binned[..., 2] = rgb_arr[2]

        # Get WCS
        opts['getWCS'] = 'true'
        request_wcs = urlencode(opts)
        try:
            wcsr = urlopen(url_ps1_fitscut, request_wcs.encode())
        except IOError:
            logger.warning("HTTP ERROR request of: %s %s", url_ps1_fitscut,
                        request_wcs)
            return blank_image(size_pix/scale_factor, scale)

        try:
            wcs_dict = json.load(wcsr)
            CD = numpy.matrix(scale_factor *
                              numpy.array(wcs_dict['cdmatrix']).reshape(2, 2))
        except KeyError:
            logger.warning("Could not get WCS info!")
            return blank_image(size_pix/scale_factor, scale)

        return rgb_binned, CD, url_ps1_fitscut+'?'+request_cutout, 'PS1'


    def retrieve_image_PANSTARRS(ra, dec, size, scale, logger):

        def getimages(ra, dec, size=240, filters="grizy"):

            """Query ps1filenames.py service to get a list of images

            ra, dec = position in degrees
            size = image size in pixels (0.25 arcsec/pixel)
            filters = string with filters to include
            Returns a table with the results
            """

            service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
            url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
                   "&filters={filters}").format(**locals())
            return Table.read(url, format='ascii')

        def geturl(ra, dec, size=240, output_size=None, filters="grizy",
                   format="jpg", color=False):

            """Get URL for images in the table

            ra, dec = position in degrees
            size = extracted image size in pixels (0.25 arcsec/pixel)
            output_size = output (display) image size in pixels (default = size).
                          output_size has no effect for fits format images.
            filters = string with filters to include
            format = data format (options are "jpg", "png" or "fits")
            color = if True, creates a color image (only for jpg or png format).
                    Default is return a list of URLs for single-filter
                    grayscale images.
            Returns a string with the URL
            """

            if color and format == "fits":
                raise ValueError("color images are available only "
                                 "for jpg or png formats")
            if format not in ("jpg", "png", "fits"):
                raise ValueError("format must be one of jpg, png, fits")
            table = getimages(ra, dec, size=size, filters=filters)
            url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
                   "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
            if output_size:
                url = url + "&output_size={}".format(output_size)
            # sort filters from red to blue
            flist = ["yzirg".find(x) for x in table['filter']]
            table = table[numpy.argsort(flist)]
            if color:
                if len(table) > 3:
                    # pick 3 filters
                    table = table[[0, len(table)//2, len(table)-1]]
                for i, param in enumerate(["red", "green", "blue"]):
                    url = url + "&{}={}".format(param, table['filename'][i])
            else:
                urlbase = url + "&red="
                url = []
                for filename in table['filename']:
                    url.append(urlbase+filename)
            return url

        def getgrayim(ra, dec, size=240, output_size=None, filter="g",
                      format="jpg"):

            """Get grayscale image at a sky position

            ra, dec = position in degrees
            size = extracted image size in pixels (0.25 arcsec/pixel)
            output_size = output (display) image size in pixels (default = size).
                          output_size has no effect for fits format images.
            filter = string with filter to extract (one of grizy)
            format = data format (options are "jpg", "png")
            Returns the image
            """

            if format not in ("jpg", "png"):
                raise ValueError("format must be jpg or png")
            if filter not in list("grizy"):
                raise ValueError("filter must be one of grizy")
            url = geturl(ra, dec, size=size, filters=filter,
                         output_size=output_size, format=format)
            r = requests.get(url[0])
            im = Image.open(BytesIO(r.content))
            return im


        size_pix = int(size*3600./scale)
        gim = getgrayim(ra, dec, size=size_pix, filter="g")
        imarray = numpy.array(gim)
        CD = numpy.matrix([[-1.*scale/3600., 0], [0, 1.*scale/3600.]])
        return imarray, CD, '', 'PS1'


    def retrieve_image_SDSS(ra, dec, size, options, scale, logger):
        """
        Load image from sdss-dr9 or dss server(jpeg) and return the image array
        and the url. Note that the transformation from world coordinate(ra,dec) to
        pixel position(x,y) is simple projection without rotation, i.e.
        x=-scale*(ra-ra0)*cos(dec)+x0; y=scale*(dec-dec0)+y0
        """
        url_sdssJPEG = 'http://skyservice.pha.jhu.edu/DR9/ImgCutout/getjpeg.aspx'
        size_pix = int(size*3600./scale)
        # options for the finding chart (see
        # http://skyserver.sdss3.org/dr9/en/tools/chart/chart.asp )

        if size_pix > 2048:
            size_pix = 2048
            scale = size*3600./size_pix

        request_sdss = urlencode({'ra': ra, 'dec': dec, 'scale': scale,
                                  'height': size_pix, 'width': size_pix,
                                  'opt': options})
        # url =
        # "http://skyservice.pha.jhu.edu/DR9/ImgCutout/getjpeg.aspx?"+request_sdss

        try:
            imfile = urlopen(url_sdssJPEG, request_sdss.encode())
        except IOError:
            logger.warning("HTTP ERROR request of: %s %s", url_sdssJPEG, request_sdss)
            return -99, -99, 'bogus', 'bogus'
        imarray = plt.imread(BytesIO(imfile.read()), format='jpeg')

        CD = numpy.matrix([[-1.*scale/3600., 0], [0, 1.*scale/3600.]])
        return imarray, CD, url_sdssJPEG+'?'+request_sdss, 'SDSS'


    def retrieve_image_SkyView(ra, dec, size, scale, logger):
        size_pix = min(512, int(size*3600./scale))
        scale = size*3600./size_pix

        paths = SkyView.get_images(position='%3.5f %2.5f' % (ra, dec),
                                   survey='DSS2 Blue', coordinates='J2000',
                                   height=size*aunit.deg, width=size*aunit.deg,
                                   pixels=str(size_pix), cache=False)
        if paths:
            hdu = paths[0][0].header
            imarray = paths[0][0].data
            imarray = imarray[::-1,:]
            CD = numpy.matrix([[hdu['CDELT1'], 0], [0, hdu['CDELT2']]])
        else:
            # Default to PS1
            logger.warning("Could not contact SDSS or DSS. Falling back to PS1!")
            return retrieve_image_PS1(ra, dec, size, scale, logger)

        return imarray, CD, '', 'DSS2'


    scale = size*3600./res
    return retrieve_image_SkyView(ra, dec, size, scale, logger)
    if SDSS_coverage(ra, dec, logger):
        return retrieve_image_SDSS(ra, dec, size, options, scale, logger)
    else:
        try:
            return retrieve_image_PS1(ra, dec, size, scale, logger)
            return retrieve_image_SkyView(ra, dec, size, scale, logger)
        except Exception:
            return retrieve_image_PANSTARRS(ra, dec, size, scale, logger)


if __name__ == "__main__":
    import logging

    level = logging.DEBUG
    fmt = '[%(levelname)s - %(filename)s - %(asctime)s] %(message)s'
    fmt = logging.Formatter(fmt)

    handler = logging.StreamHandler()
    handler.setFormatter(fmt)
    handler.setLevel(level)

    log = logging.getLogger('test')
    log.setLevel(logging.DEBUG)
    log.addHandler(handler)

    ra, dec, size, options, scale = 0.,0.,3.,"", 1.5
    img = retrieve_image(ra, dec, size, 2048, log, options)
    plt.imshow(img[0])
    plt.show()
