""">> sqlcl << command line query tool by Tamas Budavari <budavari@jhu.edu>"""


try:
    from urllib import urlopen, urlencode
except ImportError:
    from urllib.request import urlopen
    from urllib.parse import urlencode

formats = ['csv', 'xml', 'html']

astro_url = 'http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
public_url = astro_url

default_url = public_url
default_fmt = 'csv'


def filtercomment(sql):
    "Get rid of comments starting with --"
    import os
    fsql = ''
    for line in sql.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep
    return fsql


def query(sql, url=default_url, fmt=default_fmt):
    "Run query and return file object"
    fsql = filtercomment(sql)
    params = urlencode({'cmd': fsql, 'format': fmt})

    return urlopen(url+'?%s' % params)
