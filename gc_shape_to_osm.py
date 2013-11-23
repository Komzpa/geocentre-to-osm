#!/usr/bin/python
# -*- coding: utf-8 -*-
# based on the Massachusetts GIS script by christopher schmidt
# based on version 0.1 downloaded from http://boston.freemap.in/osm/files/mgis_to_osm.py
# still setup for MassGIS data.  Compare with shp files avalible at same source as above to see how it should work.
import sys
reload(sys)
sys.setdefaultencoding("utf-8")          # a hack to support UTF-8


def point_line_distance(point, startline, endline):
    """
    check if the "line" is actually a point
    if not use
    http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    Copypasted from lakewalker
    """

    if (startline == endline):
        return ((startline[0] - endline[0]) ** 2 +
               (startline[1] - endline[1]) ** 2) ** 0.5
    else:
        return abs((endline[0] - startline[0]) * (startline[1] - point[1]) -
                  (startline[0] - point[0]) * (endline[1] - startline[1])) / \
                  ((endline[0] - startline[0]) ** 2 + (endline[1] - startline[1]) ** 2) ** 0.5


def douglas_peucker(nodes, epsilon):
    """
    makes a linear curve smoother see also
    http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
    Copypasted from lakewalker
    """

    farthest_node = None
    farthest_dist = 0
    first = nodes[0]
    last = nodes[-1]

    for i in xrange(1, len(nodes) - 1):
        d = point_line_distance(nodes[i], first, last)
        if d > farthest_dist:
            farthest_dist = d
            farthest_node = i

    if farthest_dist > epsilon:
        seg_a = douglas_peucker(nodes[0:farthest_node + 1], epsilon)
        seg_b = douglas_peucker(nodes[farthest_node:], epsilon)
        nodes = seg_a[:-1] + seg_b
    else:
        return [nodes[0], nodes[-1]]
    return nodes

VERSION = "0.2"
# Version 0.2 changes to generalize the script some for use with other data by Dalep

iSource = "GcC"
iAttrib = "GeoCentre Consulting"

ignoreFcodes = []

# Files will be split when longer than this number of nodes
maxNodes = 300000000

# Set the maximum length of a way (in nodes) before it is split into
# shorter ways
Max_Waylength = 3000000

try:
    from osgeo import ogr
    from osgeo import osr
except:
    import ogr
    import osr


def parse_restrictions(filename):
    dr = ogr.GetDriverByName("ESRI Shapefile")
    poDS = dr.Open(filename)

    if poDS is None:
        raise "Open failed."

    poLayer = poDS.GetLayer(0)

    poLayer.ResetReading()

    ret = []

    poFeature = poLayer.GetNextFeature()
    while poFeature:

        tags = {}

        fff = poFeature.GetField("EDGE1ID")
        tags["from"] = fff
        fff = poFeature.GetField("EDGE2ID")
        tags["to"] = fff

        ret.append(tags)
        poFeature = poLayer.GetNextFeature()

    return ret


# ====================================
# Edit parse_shp_for_osm section to fit your data!
# change poFeature.GetField("    ") to contain only the shape column names for the data you want
# and  tags["   "] to match the osm tag names you wish to use for that data.
# some tags will require changing a number to a meaningful value like the Highway tag.  See the metadata for the meaning of these tags.
# For any measurements be sure to check the unit value of the original data, and convert if needed to the expected unit for osm.
# ====================================
def parse_shp_for_osm(filename):
    # ogr.RegisterAll()

    dr = ogr.GetDriverByName("ESRI Shapefile")
    poDS = dr.Open(filename)

    if poDS is None:
        raise "Open failed."

    poLayer = poDS.GetLayer(0)

    poLayer.ResetReading()

    ret = []

    poFeature = poLayer.GetNextFeature()
    while poFeature:

        tags = {}

        # WAY ID
        tags[iSource + ":way_id"] = int(poFeature.GetField("FEATID"))

        # FEATURE NAME
        if poFeature.GetField("NAME"):
            tags["name"] = poFeature.GetField("NAME").decode("cp1251", "ignore")

        fff = poFeature.GetField("TUNNEL")
        if fff != "N":
            tags["tunnel"] = "yes"

        fff = poFeature.GetField("BRIDGE")
        if fff != "N":
            tags["bridge"] = "yes"

        fff = poFeature.GetField("ONEWAY")
        if fff != "B":
            tags["oneway"] = "yes"
        fff = poFeature.GetField("LEVEL")
        tags["highway"] = {"1": "motorway",
                           "2": "trunk",
                           "3": "primary",
                           "4": "secondary",
                           "5": "unclassified",
                           "6": "residential",
                           "7": "service",
                           "8": "track"
                           }.get(str(fff), "path")

        geom = []

        rawgeom = poFeature.GetGeometryRef()
        for i in range(rawgeom.GetPointCount()):
            geom.append((rawgeom.GetX(i), rawgeom.GetY(i)))

        ret.append((geom, tags))
        poFeature = poLayer.GetNextFeature()

    return ret


# ====================================
# to do read .prj file for this data
# Change the Projcs_wkt to match your datas prj file.
# ====================================
projcs_wkt = \
    """GEOGCS["GCS_WGS_1984",
	DATUM["D_WGS_1984",
	SPHEROID["WGS_1984",6378137,298.257223563]],
	PRIMEM["Greenwich",0],
	UNIT["Degree",0.017453292519943295]]"""

from_proj = osr.SpatialReference()
from_proj.ImportFromWkt(projcs_wkt)

# output to WGS84
to_proj = osr.SpatialReference()
to_proj.SetWellKnownGeogCS("EPSG:4326")

tr = osr.CoordinateTransformation(from_proj, to_proj)


def unproject(point):
    pt = tr.TransformPoint(point[0], point[1])
    return (pt[1], pt[0])


def round_point(point, accuracy=8):
    return tuple([round(x, accuracy) for x in point])


def compile_nodelist(parsed_gisdata, first_id=1):
    nodelist = {}

    i = first_id
    for geom, tags in parsed_gisdata:
        if len(geom) == 0:
            continue

        for point in geom:
            r_point = round_point(point)
            if r_point not in nodelist:
                nodelist[r_point] = (i, unproject(point))
                i += 1

    return (i, nodelist)


def adjacent(left, right):
    left_left = round_point(left[0])
    left_right = round_point(left[-1])
    right_left = round_point(right[0])
    right_right = round_point(right[-1])

    return (left_left == right_left or
            left_left == right_right or
            left_right == right_left or
            left_right == right_right)


def glom(left, right):

    left = list(left)
    right = list(right)

    left_left = round_point(left[0])
    left_right = round_point(left[-1])
    right_left = round_point(right[0])
    right_right = round_point(right[-1])

    if left_left == right_left:
        left.reverse()
        return left[0:-1] + right

    if left_left == right_right:
        return right[0:-1] + left

    if left_right == right_left:
        return left[0:-1] + right

    if left_right == right_right:
        right.reverse()
        return left[0:-1] + right

    raise 'segments are not adjacent'


def glom_once(segments):
    if len(segments) == 0:
        return segments

    unsorted = list(segments)
    x = unsorted.pop(0)

    while len(unsorted) > 0:
        n = len(unsorted)

        for i in range(0, n):
            y = unsorted[i]
            if adjacent(x, y):
                y = unsorted.pop(i)
                x = glom(x, y)
                break

        # Sorted and unsorted lists have no adjacent segments
        if len(unsorted) == n:
            break

    return x, unsorted


def glom_all(segments):
    unsorted = segments
    chunks = []

    while unsorted != []:
        chunk, unsorted = glom_once(unsorted)
        chunk = douglas_peucker(chunk, 0.00005)
        chunks.append(chunk)

    return chunks


def compile_waylist(parsed_gisdata, blank_way_id):
    waylist = {}

    # Group by iSource:way_id
    for geom, tags in parsed_gisdata:
        way_key = tags.copy()
        way_key = (way_key[iSource + ':way_id'], tuple([(k, v) for k, v in way_key.iteritems()]))

        if way_key not in waylist:
            waylist[way_key] = []

        waylist[way_key].append(geom)

    ret = {}
    for (way_id, way_key), segments in waylist.iteritems():

        if way_id != blank_way_id:
            ret[way_key] = glom_all(segments)
        else:
            ret[way_key] = segments

    return ret


import time
from xml.sax.saxutils import escape


def shape_to_osm(shp_filename, restrictions_file, base_filename, blank_way_id):

    import_guid = time.strftime('%Y%m%d%H%M%S')

    print "parsing shpfile"
    parsed_features = parse_shp_for_osm(shp_filename)

    print "compiling nodelist"
    i, nodelist = compile_nodelist(parsed_features)

    print "compiling waylist"
    waylist = compile_waylist(parsed_features, blank_way_id)

    filenumber = 1
    objectcount = 0
    seen = set()
    GcC_to_OSM = {}

    print "constructing osm xml file"
    ret = []
    ret.append("<?xml version='1.0' encoding='UTF-8'?>")
    ret.append("<osm version='0.6' generator='shape_to_osm.py'>")

    for waykey, segments in waylist.iteritems():
        for segment in segments:
            # write the nodes
            for point in segment:
                id, (lat, lon) = nodelist[round_point(point)]
                if id not in seen:
                    seen.add(id)
                    # write node
                    ret.append("  <node id='-%d' action='create' visible='true' lat='%f' lon='%f' >" % (id, lat, lon))
#                    ret.append( "    <tag k=\"source\" v=\"%s_import_v%s_%s\" />" % (iSource, VERSION, import_guid) )
#                    ret.append( "    <tag k=\"attribution\" v=\"%s\" />" % (iAttrib) )
                    ret.append("  </node>")
                    objectcount += 1
                else:
                    pass
                    # print "Skipping node %d" %id

            # write the way
            ret.append("  <way id='-%d' action='create' visible='true'>" % i)

            ids = [nodelist[round_point(point)][0] for point in segment]

            count = 0
            for id in ids:
                ret.append("    <nd ref='-%d' />" % id)
            for k, v in waykey:
                ret.append("    <tag k=\"%s\" v=\"%s\" />" % (k, escape(str(v))))

                if k == (iSource + ":way_id"):
                    GcC_to_OSM[v] = i

            ret.append("    <tag k=\"attribution\" v=\"%s\" />" % (iAttrib))
            ret.append("    <tag k=\"source\" v=\"Commercial data. Delete immediately if gets to OSM!\" />")

            ret.append("  </way>")
            objectcount += 1

            i += 1

    # print GcC_to_OSM
#    print parse_restrictions(restrictions_file)
    restrictions = parse_restrictions(restrictions_file)
    i = 1
    for rest in restrictions:
        i += 1
        aa = ""
        aa += '<relation id="-%s">' % i
        aa += "    <tag k=\"type\" v=\"restriction\" />"
        aa += "    <tag k=\"restriction\" v=\"no_right_turn\" />"
        for k, v in rest.items():

            aa += "    <member type=\"way\" role=\"%s\" ref=\"-%s\" />" % (k, GcC_to_OSM[v])
        aa += '</relation>'
        ret.append(aa)

    ret.append("</osm>")

    osm_filename = "%s%d.osm" % (base_filename, filenumber)
    print "writing %s" % osm_filename
    fp = open(osm_filename, "w")
    fp.write("\n".join(ret))
    fp.close()

if __name__ == '__main__':
    import sys
    import os.path
    if len(sys.argv) < 2:
        print "%s filename.shp restrictions.shp filename.osm" % sys.argv[0]
        sys.exit()
    shape = sys.argv[1]
    restrictions = sys.argv[2]
    osm = sys.argv[3]

    id = "1.shp"
        # Left over from massGIS unknown usage, but works fine hardcoded to "1.shp" which was the valu on a test of the actual mass data,
        # id = os.path.basename(shape).split("_")[-1]
    shape_to_osm(shape, restrictions, osm, id)
