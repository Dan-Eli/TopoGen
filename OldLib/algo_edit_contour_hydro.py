#! /usr/bin/env python
# -*- coding: UTF-8 -*-
#####################################################################################################################################

"""
    This algorithm edits the contours (hypsographic features) in order to keep the coherence with the lakes, rivers and streams.
    The following constraint are implemented:   
        - A contour cannot cross a lake
        - A contour can cross a River or Stream only once
        - For the rivers, spot height (x,y,z) are used to keep the coherence with the contour
        - For the streams, the stream orientation (uphill to downhill) is used to keep the coherence with the contour   
    The following strategy is implemented
        - The edition is done in the following order
            - First clean the Lakes, after the Rivers at last the Streams
        - A contour that touch a lake can be moved on a River or a Stream
        - A contour that touch a River can be moved on a Stream only
        - A contour that touch a Stream cannot be moved on any feature
        - When a contour is moved it cannot touch to another contour
    
    Usage:
        import algo_edit_contour_hydro

    Limits and constraints
        The streams must be oriented
        Works better when the contours lines are closed.
        Works better when the permanent and intermittent part of a Lake or River are merges together
        Works better when the stream belonging to the major stream are not broken (merge stream based on Sthaler order)
          
"""
__revision__ = "--REVISION-- : $Id: algo_edit_contour_hydro.py 472 2011-09-29 13:39:22Z dpilon $"

#####################################################################################################################################

import sys, operator
import abc
from math import sqrt
from copy import deepcopy
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection
from lib_genmetal import MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, Parameters,SpatialContainer, Algorithm, Holder, InternalError
                         
########################################################

# ************** � D�truire
try:
   import pyfme
   logger = pyfme.FMELogfile()
except:
    pass


# Public constant
ALTI = "ALTI"
AREA_MINIMUM = "AREA MINIMUM"
CODE = "CODE"
HYPSO = "HYPSO"
LAKE = "LAKE"
RIVER = "RIVER"
STREAM = "STREAM"
Z = "Z"
Z_POINT = "Z_POINT"

_AFTER = "After"
_ALGO = "Algo"
_BEFORE = "Before"
_CLEAN_HYDRO = "Clean Hydrography"
_DELETED = "Hypso deleted"
_FIRST_COORDS = "First coords"
_NOT_PROCESSED = "Hypso not processed"
_INSIDE = "Inside"
_IS_ALL_EDITED = "Is all edited?"
_LAST_COORDS = "Last coords"
_NOT_SIMPLE = "Not simple"
_ORIENTATION = "Orientation"
_OUTSIDE = "Outside"
_SIMPLE = "Simple"
_COMPLEX = "Complex"
_TOUCH = "Touch"

class HydroUtil:
    """
    In this class we find many little utility routines.
    Some of these routine should be moved in the more general classs GenUtil as they are really general
    """
    
    @staticmethod
    def shortest_perimeter(p0_lr, p1_lr, line):
        """Extracts the line forming the shortest path between two points on a closed line

        *Parameters*
            - p0_lr_hydro: Start linear reference on the closed line
            - p1_lr_hydro: End linear reference on the closed line

        *Returns*:
            - LineString describing the shortest path on the closed line between start and end
            
        """
        
        sub_perimeter_a = HydroUtil.extract_sub_line(p0_lr, p1_lr, line)
        sub_perimeter_b = HydroUtil.extract_sub_line(p1_lr, p0_lr, line)
        if sub_perimeter_a.length < sub_perimeter_b.length:
            shortest_perimeter = sub_perimeter_a
        else:
            shortest_perimeter = sub_perimeter_b
             
        return shortest_perimeter
    
    @staticmethod
    def extract_sub_line (lr_a, lr_b, line):
        """Extracts a sub line of a line located between 2 linear references
        
        Note: Work on a closed or an open line. For an open line, if first linear reference is greater than the second,
              it will extact 2 lines  
        
        *Parameters*
            - lr_a: First linear reference
            - lr_b: second linear reference
            - line: LineString object to extract a sub line
            
        *Return*
            - LineString or MultiLineString sub line 
        
        """
         
        sub_line = LineString(line.coords)
         
        if lr_a < lr_b:
            # Manage some limit cases
            if lr_a < 0.0:
                lr_a = 0.0
            if lr_b > line.length:
                lr_b = line.length
            
            if lr_b >= line.length:
                # Nothing to do to extract pass the end of line
                pass
            else:
                # Cut and delete from lr_b to end
                tmp_lines = HydroUtil.cut_line_distance(sub_line, lr_b)
                sub_line = tmp_lines[0]
         
            if lr_a <= 0.0:
                # Nothing to 
                pass
            else:
                # Cut from start to lr_a
                tmp_lines = HydroUtil.cut_line_distance(sub_line, lr_a)
                sub_line = tmp_lines[1]
                
        else:             
            # Manage some limit cases
            if lr_a > line.length:
                lr_a = line.length
            if lr_b < 0.0:
                lr_b = 0.0
            
            sub_lines = []
            
            if lr_a >= line.length:
                # Nothing to extract pass the end of line
                pass
            else:
                tmp_lines = HydroUtil.cut_line_distance(sub_line, lr_a)
                sub_lines.append(tmp_lines[1])
            
            if lr_b <= 0.0:
                # Nothing to do to extract before the start
                pass
            else:
                tmp_lines = HydroUtil.cut_line_distance(sub_line, lr_b)
                sub_lines.append(tmp_lines[0])
                
            # Special case the extracted line is the same as the original line
            if len(sub_lines) == 0:
                sub_lines = [LineString(line.coords)]
                 
            #Merge the line if possible
            sub_line = HydroUtil.merge_lines(sub_lines)
                         
        return sub_line
    
    @staticmethod
    def multi_to_list (multi_feature):
        """Transform a simple or a multi geometry into a list of geometry
        
        *Parameters*
            - multi_feature: simple or multi part geometry
            
        *Return*
            - List of simple geometry (Point, LineString, Polygon)
            
        """
        if ( isinstance(multi_feature, MultiPoint) or \
             isinstance(multi_feature, MultiLineString) or \
             isinstance(multi_feature, MultiPolygon) or \
             isinstance(multi_feature, GeometryCollection) ):           
            lst_feature = [feature for feature in multi_feature]
        else:
            # It's not a multi feature but put it in a list anyway
            lst_feature = []
            lst_feature.append(multi_feature)
        
        return lst_feature
        
    @staticmethod
    def cut_line_distance(line, distance, copy_att=False):
        """Cut the line at a specific distance along the line
    
        *Parameters*:
            - line: LineString to cut
            - distance: distance in ground unit from the start of the line where to cut the line
            - copy_att: Flag to enable (True) or disable (False) the copy of the attributes in
                        in ma_properties for MA_LineString only
        
        *Returns*:
            - List of LineString or MA_LineString
    """
    
        if (isinstance(line, MA_LineString)):
            if line.is_dual():
                coords = line.coords_dual
            else:
                coords = list(line.coords)
        elif (isinstance(line, LineString)):
            coords = list(line.coords)
            if (copy_att):
                raise InternalError ("Cannot copy attributes on LineString only on MA_LineString")
        else:
            raise InternalError ("Cannot cut a %s feature" %(line.geom_type))
        
        if distance <= 0.0 or distance >= line.length:
            # Special case extract the same line
            lst_coords_out = [coords]
        else:
            lst_coords_out = [coords]  # This line will prevent errors on non simple lines
            last_i = len(coords)-1
            for i, p in enumerate(coords):
                if (i < last_i):
                    if i == 0:
                        pd = 0.
                    else:
                        pd += sqrt( (coords[i-1][1] - coords[i][1])**2 + (coords[i-1][0] - coords[i][0])**2 )
                    
                else:
                    pd = line.length
                if pd == distance:
                    # Create the list of coordinate
                    lst_coords_out =  [coords[:i+1], coords[i:] ]
                if pd > distance:
                    # Create the list of coordinate
                    cp = line.interpolate(distance)
                    lst_coords_out =  [list(coords[:i]) + [(cp.x, cp.y)], [(cp.x, cp.y)] + list(coords[i:])]
                    break
        
        lines_out = []
        # Create the lines from the ist of coordinate
        for coords_out in lst_coords_out:
            if isinstance(line, MA_LineString):
                line_out =  MA_LineString(coords_out)
            else:
                line_out =  LineString(coords_out)
            if (copy_att):
                # Copy the attributes
                line_out.ma_properties = deepcopy(line.ma_properties)
            lines_out.append(line_out)
                                    
        return lines_out
    
    @staticmethod
    def merge_lines (lines):
        """Merge a list of LineString together.  The start end can be within a tolerance
        
        *Parameters*
            - lines: List of LineString
            
        *Returns*
            - LineString or MultiLineString if the lines are not all merged
            
        Note: The lines that cannot be merged are left unmerged in the list of LineString
              There is a LineMerge method in Shapely but the coordinates must be exaclty the same in order for the lines to be merged.
              There method is good for small amount of lines otherwise it would need to be optimized 
        """
        coords = []
        for line in lines:
            coords.append(list(line.coords))
         
        # Loop over each line and try to find a line to merge with
        for i in range(1, len(coords)): 
            if len(lines) == 1: break 
            for j in range(len(coords)-1, 0,-1):
                merge = False
                if GenUtil.distance(coords[0][0], coords[j][0]) < GenUtil.ZERO:

                    # Add start first line with start current line
                    coords[0].reverse()
                    coords[0] += coords[j][1:]
                    merge = True
                elif GenUtil.distance(coords[0][0], coords[j][-1]) < GenUtil.ZERO:                
                    # Add start first line with end current line
                    coords[0] = coords[j] + coords[0][1:]
                    merge = True
                elif GenUtil.distance(coords[0][-1], coords[j][0]) < GenUtil.ZERO: 
                    # Add end first line with start current line
                    coords[0] += coords[j][1:]
                    merge = True
                elif GenUtil.distance(coords[0][-1], coords[j][-1]) < GenUtil.ZERO:
                    # Add end first line with end current line
                    coords[j].reverse()
                    coords[0] += coords[j][1:]
                    merge = True
                    
                if merge:
                    # Remove the line that was merged
                    del coords[j]
                    break
                
        if len(coords) == 1:
            if GenUtil.distance(coords[0][0], coords[0][-1]) < GenUtil.ZERO:
                coords[0][0] = coords[0][-1]
            line = LineString(coords[0])
        else:
            line = MultiLineString(coords)
             
        return line
    
    @staticmethod
    def are_lines_to_close (source_line, targets_line, distance_min):
        """???
        """

        buffer_distance = distance_min*2.
        colinear = False
        lr0 = buffer_distance
        lr1 = source_line.length - buffer_distance
        if lr1 < lr0:
            lr0 = source_line.length/2
            lr1 = lr0*1.01
        sub_source_line = HydroUtil.extract_sub_line(lr0, lr1, source_line)
        
        for target_line in targets_line:
            lr0 = buffer_distance
            lr1 = target_line.length - buffer_distance
            if lr1 < lr0:
                lr0 = target_line.length/2
                lr1 = lr0*1.01
            sub_target_line = HydroUtil.extract_sub_line(lr0, lr1, target_line)
#            print "area line to close ", sub_source_line.distance(sub_target_line), "min: ", distance_min
            if sub_source_line.distance(sub_target_line) <= distance_min*0.9:
                colinear = True
                break
        
        return colinear
                                
    @staticmethod
    def extend_coord(coord_a, coord_b, length):
        """Extend a line at the end of a line segment by a certain length
        
        *Parameters*
            - coord_a: First coordinate of the line segment
            - coord_b: Second coordinate of the line segment
            - length: Length of the extension
            
        *Return*
            Tuple (x,y) of the end of the extended line
        
        """
        
        ax = coord_a[0]
        ay = coord_a[1]
        bx = coord_b[0]
        by = coord_b[1]
        
        dist = GenUtil.distance(coord_a, coord_b)
        if dist < GenUtil.ZERO:
            dist = GenUtil.ZERO
            
        cx = bx + (bx-ax)/dist * length
        cy = by + (by-ay)/dist * length
        
        return (cx,cy)
        
    @staticmethod
    def get_orientation (lst_coords):
        """Test if a closed list of coordinates is forming a clockwise or anti clockwise loop
        
        The test is done by calculating the area of the closed line
        
        lst_coords: List of (x,y) tuple forming a closed polygon
        
        return: >0: Means clockwise
                <0: Means anti clockwise 
        """
        
        nbr_edges = len(lst_coords)-1
        area = 0.
        for i in range(nbr_edges):
            coords_a = lst_coords[i]
            coords_b = lst_coords[i+1]
            area += coords_a[0]*coords_b[1] - coords_a[1]*coords_b[0] 
        
        area /= 2.
            
        return area
    
    @staticmethod
    def mean_max_distance (line_a, line_b): # conflict.to_cut_hypso_line  conflict.to_paste_hypso_line  
        """Calculate the mean maximum distance between two LineString
        
        A series of distance samples are taken on each line, the mean of the fartest samples is calculated 
         
        Parameters:
            line_a, line_b: LineString object to calculate the mean of the maximum distance between them
             
        Return value:
            Mean distance between the 2 LineString
         
        """
        # Total number of samples to take
        samples = 30
        # Number of fartest sample that will be used to calculate the mean maximum distance 
        sub_samples = samples/6
             
        distances = []
        for i in range(samples):
            ratio = float(i)/float(samples)
            point_on_line_a = line_a.interpolate(ratio, normalized=True)
            point_on_line_b = line_b.interpolate(ratio, normalized=True)
            distances.append(point_on_line_a.distance(line_b))
            distances.append(point_on_line_b.distance(line_a))
         
        # We keep only a subset of the samples of the farthest point
        total = 0
        distances.sort(reverse=True)
        for i in range(sub_samples):
            total += distances[i]
         
        mean = total/sub_samples
             
        return mean

class Statistics(GenStatistics):
    """Class that contains the statistics for the talweg correction algorithm
    
    Attributes
        stat_names: Name of the statistics for the TalwegStatistics class. These name are
                    used by the Statistics class

    """
    
    def __init__(self):
        """Initialize the attributes of an object of the class"""
        
        GenStatistics.__init__(self)
        self.stats_names = ((_NOT_SIMPLE, _NOT_PROCESSED, _DELETED, _INSIDE, _OUTSIDE, GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS, ))
        
        self.add_iteration()
        
    def get_stats (self, dummy=None):
        """Extract the current statistics and build  a list of string that forms the statistical message"
        
        Parameters:
            dummy: Unused... only there for compatibility purpose
        
        """
        
        str_out = []
        str_out.append( "%s coherence algorithm Statistics" %(_CLEAN_HYDRO) )
        str_out.append( "-------------------------------------" )
        str_out.append("Hypsographic line not simple: " + str(self.get_stats_name_count_total( _NOT_SIMPLE)))
        str_out.append("Hypsographic line deleted (completely within buffer): " + str(self.get_stats_name_count_total( _DELETED)))
        str_out.append("Hypsographic line not processed (first/last point within buffer): " + str(self.get_stats_name_count_total( _NOT_PROCESSED)))
        str_out.append("Simple hypsography line edited: " + str(self.get_stats_name_count_total( _INSIDE)))
        str_out.append("Complex hypsography line edited: " + str(self.get_stats_name_count_total( _OUTSIDE)))
        str_out.append( "-------------------------------------" )
        
        return str_out
    
class Buffer(object):
    """ Class that manages the buffer
    """
    
    # The class variable _TERATION_RATIO controls the number of iteration and the ratio of the buffer value for each iteration
    # If the _ITERATION_RATIO = [1., 1/2, 1/5] and the buffer value is 20 there will be 3 iteration and the buffer value would be: 20, 10, 4
    _ITERATION_RATIO = [1., 7./10., 1./2., 1./3., 1./10.]
    _MAX_SUB_ITERATION = 10
    
    def __init__(self, max_buffer):
        """Constructor of the Buffer class
        
        Parameters:
            max_buffer: Maximum buffer value used to calculate the different iteration
        """
        
        self.max_buffer = float(max_buffer)
        self._iteration_value = []
        
        # Calculate the buffer value for each iteration 
        for ratio in Buffer._ITERATION_RATIO:
            value = ratio * self.max_buffer
            self._iteration_value.append(value)
                
    def next_iteration(self):
        """Extract the next buffer value for the next iteration
        
        Parameters:
            None
            
        Return value:
            Value (float) of the next buffer value of the next iteration
        """

        # Iterate over the buffer value list
        for self._i_iteration, dummy in enumerate(self._iteration_value):
            self._i_sub_iteration = -1
            value = self.get_current()
            yield value
        
    def next_sub_iteration(self):
        """Extract the next buffer value for the next sub iteration
        
        Parameters:
            None
            
        Return value:
            Value (float) of the next buffer value of the next sub iteration
            None if try to extract past the last value
        """
        
        # Return the buffer value for a sub iteration
        self._i_sub_iteration += 1
        if self._i_sub_iteration < Buffer._MAX_SUB_ITERATION:
            value = self.get_current()
        else:
            # No more sub iteration left (should not happen...)
            value = None
            
        return value
            
    def reset_sub_iteration(self):
        """Reset the sub iteration flag
        
        Parameters:
            None
            
        Return value:
            None
        """
        
        self._i_sub_iteration = -1
            
    def get_current(self):
        """Extract the current buffer value for the corresponding iteration and sub iteration
        
        Parameters:
            None
            
        Return value:
            Value (float) of the current iteration and sub iteration
        """
        
        iteration_value = self._iteration_value[self._i_iteration]
        nbr_steps = Buffer._MAX_SUB_ITERATION * 4
        if self._i_sub_iteration == -1:
            i_sub_iteration = 0
        else:
            i_sub_iteration = self._i_sub_iteration
        sub_iteration_value = iteration_value - (float(i_sub_iteration)/float(nbr_steps))

        return sub_iteration_value                
            
    def is_last_iteration(self):
        """Flag indicating if it is the last iteration
        
        Parameters:
            None
            
        Return value
            Boolean: True: It is the last iteration; False: Otherwise
        """
        
        if self._i_iteration+1 == len(Buffer._ITERATION_RATIO):
            last_iteration = True
        else:
            last_iteration = False
            
        return last_iteration

    def get_current_delta(self):
        """Extract the difference between the current buffer value and the max buffer value
        
        Parameters:
            None
            
        Return value
            Difference between the current buffer value and the max buffer value 
        """
                
        total_delta_buffer = self.max_buffer - self.get_current()
        if total_delta_buffer <= 0.:
            total_delta_buffer = GenUtil.ZERO
            
        return total_delta_buffer

class Conflicts(object):
    """Abstract class that manage the resolution of a conflicts between a contour and a lake, river or stream
    
    Attributes
        - hydro_feature to process: Hydrographic feature (LineString or Polygon)
        - buffer: Buffer object
        - max_displacment: distance maximum that a portion of an hypsographic line can be moved
        - max_displacement: Maximum distance a contour can be moved when edited
        - s_container_z_point: Spatial Container for the Z point
        - min_atrea: Area minimum to keep new contours 
        
    """
    
    # Define the class as abstract class
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, hydro_feature, buf, params, s_container_z_point = None):
        """Initialize the different attributes of the object Conflicts
        
        Parameters
        - hydro_feature: Hydrographic feature (LineString or Polygon)
        - hypso_line: Hypsographic feature (LineString)
        - s_container_z_point: Spatial container for the altimetric information
        - buf: Buffer object
        - params: Parameters of the AlgoEditContour_Hydro class
        """
        
        self.hydro_feature = hydro_feature
        self.buffer = buf
        self.max_displacement = params.max_displacement
        self.s_container_z_point = s_container_z_point
        self.resolution = params.resolution
        self.min_area = params.min_area    
    
    @abc.abstractmethod
    def is_edition_needed(self, hypso_line ):
        """Abstract method determining if some edition is needed for a contour
        
        Parameters:
            - hypso_line: LineString hypsographic feature
            
        Return value
            - Boolean: True: Editionis needed; False: No edition is needed
            
        """
        
    @abc.abstractmethod
    def _is_displacement_allowed(self, distance, max_displacement, hypso_line=None):
        """Abstract method determining if an edition is allowed"""

    @abc.abstractmethod
    def _is_all_edited(self, hypso_line):
        """Checks if all the conflicts are solved/edited"""
        
    @abc.abstractmethod    
    def _priorize_crossing(self):
        """Priorize the crossing for the edition"""
            
    @abc.abstractmethod
    def _micro_intersection_edition(self, hypso_line, s_container):
        """Abstract class for the micro edition needed for some intersection"""
    
    def _create_buffer(self):
        """Creates a buffer for the hypsographic feature
        
        Parameters:
            None
            
        Return value:
            - Buffered Polygon feature of the hydrographic feature
        """
        
        return self.hydro_feature.buffer(self.buffer.get_current(), resolution=3)
    
    def _creates_hydro_hypso_crossing(self, hypso_line):
        """Extract the list of intersections that the buffered hydrographic feature is making with the hypsographic feature
        
        *Parameters*
            - hypso_line: LineString object representing the hypsographic feature 
            
        *Return Value*
            Boolean: True: The crossing were created; Flase: the crossing were not created whatever the reason
        
        """
        
        valid = False
        while not valid:
            # Create a buffer around the hydrographic feature
            buf_hydro_pol = self._create_buffer()
            
            # Make a LineString of the polygon exterior
            self.buf_hydro_line = LineString(buf_hydro_pol.exterior.coords)
            # Strip the polygon of any polygon interiors
            self.buf_hydro_pol = Polygon(buf_hydro_pol.exterior.coords)
            # Detect intersection between hydrographic feature  and hypsographic feature            
            intersect_points = self.buf_hydro_line.intersection(hypso_line)
            intersect_points = GenUtil.make_iterable(intersect_points)
            lst_points = [1 for geometry in intersect_points if not isinstance(geometry, Point)]
            # check that all the intersections are Points (no LineString are accpeted)
            if len(lst_points) == 0:
                # All intersections are Point
                self._creates_crossing_pairs(intersect_points, hypso_line)
                valid = self._validate_inside_outside(hypso_line)
            else:
                valid = False
                
            self.buffer.next_sub_iteration()
            self.valid_crossing = valid

        return valid
    
    def _creates_crossing_pairs (self, intersect_points, hypso_line):
        """Determine where the hypsographic line is inside and outside of the buffered hydrographic feature
        
        From the ordered list of intersection points determine for all the pair of points if they are located inside or outside of 
        the buffered hydrographic feature
        
        *Parameters*
            intersect_points: List of Point where the hypsographic feature is intersecting the buffered hydrographic feature
            hypso_line: Hypsographic LineString to process
            
        *Return value*
            None
        """
                 
        self.crossing_pairs = []
        if len(intersect_points) >= 2:
            
            # Order the intersection point on the hypsographic feature from first to last point
            crossing_poins_lr = []
            for point in intersect_points:
                crossing_poins_lr.append(hypso_line.project(point))
            crossing_poins_lr.sort()
             
            # Determine the pair of crossing 
            nbr_crossing_pair = len(crossing_poins_lr)
            lst_first_second = []
            for i in range(nbr_crossing_pair-1):
                lst_first_second.append((i,i+1))
            if GenUtil.is_line_closed(hypso_line):
                # For a closed line the first and last crossing are joined
                lst_first_second.append((i+1,0))
           
            # Creates the list of pair of crossing
            for first_second in lst_first_second:
                first = first_second [0]
                second = first_second [1]
                p0 = hypso_line.interpolate(crossing_poins_lr[first])
                p1 = hypso_line.interpolate(crossing_poins_lr[second])
                crossing_pair = Holder()
                # set some attribute to the crossing pair
                self._set_crossing_pair(hypso_line, crossing_pair, p0, p1)
                self.crossing_pairs.append(crossing_pair)
         
            # Give unique ID to each crossing_pair
            for i, crossing_pair in enumerate(self.crossing_pairs):
                crossing_pair.id = i
        else:
            # should almost never arrive and there is nothing to do
            pass
        
        return
    
    def _validate_inside_outside(self, hypso_line):
        """Validate if the switch between inside and outside are valid
        
        For an open contour (not closed). The first and the last crossing must be located inside the hydrographic feature 
        and all the other must be a switch between inside and outside. If two consecutive crossing are inside or outside this
        mean that the hypsographic feature is not crossing the hydrographic but just touching it
        
        For a closed contour.  There must must be a perfect switch between inside and outside between each crossing any none
        alternation is a mistake. 
        """
            
        valid_crossing = True
        crossing_types = []
        
        for crossing_pair in self.crossing_pairs:
            # All the crossing must be Inside,  Outside otherwise it's not valid
            if crossing_pair.type in (_INSIDE, _OUTSIDE):
                crossing_types.append(crossing_pair.type)
            else:
                # It's a COMPLEX crossing and it's not valid 
                valid_crossing = False
        
        if valid_crossing and len(self.crossing_pairs) >= 1:
                                
            if GenUtil.is_line_closed(hypso_line):
                # For a closed lines there must be an even number of crossing
                if len(self.crossing_pairs) % 2 == 0:
                    state = crossing_types[0]
                    # Validate the switching between Inside and Outside
                    for i in xrange(len(crossing_types)):
                        if crossing_types[i] == state:
                            if state == _INSIDE:  
                                state = _OUTSIDE
                            else:
                                state = _INSIDE
                        else:
                            valid_crossing = False
                else:
                    valid_crossing = False
                            
            else:
                # For an open line first and last crossing must be Inside
                if crossing_types[0] == _INSIDE and crossing_types[-1] == _INSIDE:
                    state = _OUTSIDE
                    # Validate the switching between Inside and Outside
                    for i in xrange(1,len(crossing_types)-1):
                        if crossing_types[i] == state:
                            if state == _INSIDE:  
                                state = _OUTSIDE
                            else:
                                state = _INSIDE
                        else:
                            valid_crossing = False
                else:
                    valid_crossing = False
                    
        return valid_crossing
        
    def _set_crossing_pair (self, hypso_line, crossing_pair, p0, p1):
        """Sets various attribute for a crossing pair.
        
        *Parameters*
            - crossing_pair: Structure containing the information about the crossing pair
            - p0: First Point where the hypso is crossing the hydrophic feature
            - p1: Second Point where the hypso is crossing the hydrophic feature
            
        *Return value*
            - None
        """
        
        crossing_pair.p0_point = p0
        crossing_pair.p1_point = p1
        if ( GenUtil.distance (crossing_pair.p0_point.coords[0], crossing_pair.p1_point.coords[0]) < 10.*GenUtil.ZERO ):
            # The distance between the 2 crossing points is 2 small and it could bring calculation and/or rounding errors
            crossing_pair.type = _COMPLEX
        else:
            # Extract the sub hypsographic LineString between the 2 intersections points
            crossing_pair.lr0_hypso = hypso_line.project(p0)
            crossing_pair.lr1_hypso = hypso_line.project(p1)
            crossing_pair.p0_p1_line = HydroUtil.extract_sub_line(crossing_pair.lr0_hypso, crossing_pair.lr1_hypso, hypso_line)
            
            #Determine if Inside or Outside
            mid_point = crossing_pair.p0_p1_line.interpolate(crossing_pair.p0_p1_line.length/2.)
            if mid_point.within(self.buf_hydro_pol):
                # The sub hypsographic LineString is located inside the hydrographic feature
                crossing_pair.type = _INSIDE 
                lr0_hydro = self.buf_hydro_line.project(p0)
                lr1_hydro = self.buf_hydro_line.project(p1)
                
                # Calculates the 2 polygons formed by the sub hypsographic intersecting the hydrographic polygon 
                half_a_hydro_line = HydroUtil.extract_sub_line(lr0_hydro, lr1_hydro, self.buf_hydro_line)
                half_a_hydro_close_line = HydroUtil.merge_lines([half_a_hydro_line, crossing_pair.p0_p1_line])
                crossing_pair.half_a_hydro_pol = Polygon(half_a_hydro_close_line)
                 
                half_b_hydro_line = HydroUtil.extract_sub_line(lr1_hydro, lr0_hydro, self.buf_hydro_line)
                half_b_hydro_close_line = HydroUtil.merge_lines([half_b_hydro_line, crossing_pair.p0_p1_line])
                crossing_pair.half_b_hydro_pol = Polygon(half_b_hydro_close_line)
                
            else:
                # The sub hypsographic LineString is located inside the hydrographic feature
                crossing_pair.type = _OUTSIDE
             
        return
    
    def _get_sibbling(self, hypso_line, crossing_pair_current):
        """Find the adjacents right and left (sibbling) crossing for a specific crossing pair 
        
           For a closed a closed line the crossing pairs are circular and there is always a sibbling.
           For an open line, the is no sibbling for the first/last crossing pair
           
        *Parameters*
            - hypso_line: Hypsographic LineString
            - crossing_pair_current: Current crossing for which we want the sibbling crossing pair
            
        *Return*
            - tuple containing the sibbling crossing pair or None if there is no sibbling crossing pair
        """
        
        underflow = -1
        overflow = len(self.crossing_pairs)
        current = crossing_pair_current.id
        
        # Set the sibbling
        before = current - 1
        after = current + 1
        
        if before == underflow:
            if GenUtil.is_line_closed(hypso_line):
                #Fro a closed line, it's a circular array
                crossing_pair_before = self.crossing_pairs[-1]
            else:
                # No adjacent crossing pair
                crossing_pair_before = None
        else:
            crossing_pair_before = self.crossing_pairs[before]
            
        if after == overflow:
            if GenUtil.is_line_closed(hypso_line):
                # It's a circular array
                crossing_pair_after = self.crossing_pairs[0]
            else:
                # No adjacent crossing pair
                crossing_pair_after = None
        else:
            crossing_pair_after = self.crossing_pairs[after]
            
        return (crossing_pair_before, crossing_pair_after)
                
    def _get_nbr_crossing_inside(self):
        """Calculates the number of crossing pair of type INSIDE (located the buffered hydrographic feature)
        
        *Parameter*
            None
            
        *Return value*
            Number of crossing pair located inside the hypsographic feature
            
        """
        nbr_crossing_inside = 0
        for crossing_pair in self.crossing_pairs:
            if crossing_pair.type == _INSIDE:
                nbr_crossing_inside += 1
                
        return nbr_crossing_inside
    
    def _set_conflict_attribute (self, hypso_line, conflict):
        """Calculates different attributes for a conflict
        
        *Parameters*
            - hypso_line: Hypsographic MA_LineString
            - conflict: Conflict for which to calculate the attributes
            
        *Return*
            None
        """
        
        current_lr0_hypso = hypso_line.project(conflict.current.p0_point)
        current_lr1_hypso = hypso_line.project(conflict.current.p1_point)
        current_lr0_hydro = self.buf_hydro_line.project(conflict.current.p0_point)
        current_lr1_hydro = self.buf_hydro_line.project(conflict.current.p1_point)
        to_cut_hypso_line = HydroUtil.extract_sub_line(current_lr0_hypso, current_lr1_hypso, hypso_line)
        
        if conflict.type == _SIMPLE:
            # Extract some attributes for _SIMPLE case
            conflict.to_keep_hypso_line_a = HydroUtil.extract_sub_line(current_lr1_hypso, current_lr0_hypso, hypso_line)
            conflict.to_paste_hypso_line_a = HydroUtil.shortest_perimeter(current_lr0_hydro, current_lr1_hydro, self.buf_hydro_line)
            conflict.to_keep_hypso_line_a = HydroUtil.multi_to_list (conflict.to_keep_hypso_line_a)
            conflict.new_hypso_line_a = HydroUtil.merge_lines(conflict.to_keep_hypso_line_a +  [conflict.to_paste_hypso_line_a])
                                    
            # Check that the line has the same orientation has the original line
            conflict.new_hypso_line_a = self._orient_line (hypso_line, conflict.new_hypso_line_a)
            
            # For a simple case always create empty Polygon 
            conflict.new_hypso_line_b = None
            conflict.to_keep_hypso_line_b = None
            conflict.to_paste_hypso_line_b = None
            
        else:
            # Extract some attributes needed to solve _COMPLEX case
            before_lr0_hypso = hypso_line.project(conflict.before.p0_point)
#            before_lr1_hypso = hypso_line.project(conflict.before.p1_point)
            before_lr0_hydro = self.buf_hydro_line.project(conflict.before.p0_point)
#            before_lr1_hydro = self.buf_hydro_line.project(conflict.before.p1_point)
#            after_lr0_hypso = hypso_line.project(conflict.after.p0_point)
            after_lr1_hypso = hypso_line.project(conflict.after.p1_point)
#            after_lr0_hydro = self.buf_hydro_line.project(conflict.after.p0_point)
            after_lr1_hydro = self.buf_hydro_line.project(conflict.after.p1_point)
            
            # Extract information for the first new hypso line (line_a)
            conflict.to_keep_hypso_line_a = HydroUtil.extract_sub_line(after_lr1_hypso, before_lr0_hypso, hypso_line)    
            conflict.to_paste_hypso_line_a = HydroUtil.shortest_perimeter(before_lr0_hydro, after_lr1_hydro, self.buf_hydro_line)
            conflict.to_keep_hypso_line_a = HydroUtil.multi_to_list (conflict.to_keep_hypso_line_a)
            conflict.new_hypso_line_a = HydroUtil.merge_lines(conflict.to_keep_hypso_line_a +  [conflict.to_paste_hypso_line_a])
                                    
            # MAke sure that the new line has the same orientation has the original line even if it is edited
            conflict.new_hypso_line_a = self._orient_line (hypso_line, conflict.new_hypso_line_a)

            # Extract information for the second new hypso line (line_b)
            conflict.to_keep_hypso_line_b = HydroUtil.extract_sub_line(current_lr0_hypso, current_lr1_hypso, hypso_line)
            conflict.to_keep_hypso_line_b = HydroUtil.multi_to_list (conflict.to_keep_hypso_line_b)
            conflict.to_paste_hypso_line_b = HydroUtil.shortest_perimeter(current_lr0_hydro, current_lr1_hydro, self.buf_hydro_line )
            conflict.new_hypso_line_b = HydroUtil.merge_lines(conflict.to_keep_hypso_line_b +  [conflict.to_paste_hypso_line_b])
            if isinstance(conflict.new_hypso_line_b, LineString) and GenUtil.is_line_closed(conflict.new_hypso_line_b):
                conflict.new_hypso_pol_b = Polygon(conflict.new_hypso_line_b.coords)
                if conflict.new_hypso_pol_b.area >= self.min_area:
                    # The hypso line b is forming a valid polygone with the area over the minimum size
                    valid_hypso_line_b = True
                else:
                    # Minimum size is not respected
                    valid_hypso_line_b = False
            else:
                # Not a valid polygon
                valid_hypso_line_b = False
            
            if not valid_hypso_line_b:
                # Polygon not valid...
                conflict.new_hypso_pol_b = None
                conflict.new_hypso_line_b = None
                conflict.to_keep_hypso_line_b = None
                conflict.to_paste_hypso_line_b = None
                                                                        
        if conflict.type == _SIMPLE:
            to_cut_hypso_line = HydroUtil.extract_sub_line(current_lr0_hypso, current_lr1_hypso, hypso_line)
        else:
            to_cut_hypso_line = HydroUtil.extract_sub_line(before_lr0_hypso, after_lr1_hypso, hypso_line)
            
        conflict.distance = HydroUtil.mean_max_distance(to_cut_hypso_line, conflict.to_paste_hypso_line_a)    
        
        return
 
    
    def _manage_priorized_complex_conflict(self, hypso_line, min_max_crossing_pair):
        """Identify the complex conflict for a hydrographic feature that need priority conflict management
        
        The most upstream  crossing pair must not be included in a complex conflict. The most downstream crossing pair
        is included in the conflict that has the most chances to be solved.  A river need priorization when the number of
        inside crossing pair is odd.
        
        *Parameters*
            - hypso_line: Hypsographic MA_LineString
            - min_max_crossing_pair: List containing the inside crossing pair sorted by elevation (min to max)
            
        *Return value*
            None
        """
        
        # Simple inside is the counter a crossing that should be process as simple case and not complex
        simple_inside = 0
        
        min_z_crossing_pair_id = min_max_crossing_pair[0].id
        max_z_crossing_pair_id = min_max_crossing_pair[-1].id

        if GenUtil.is_line_closed(hypso_line):
        
            # From the crossing with the highest z (highest priority) process the other crossing located before and after in 2 different pass
            for direction in (_BEFORE, _AFTER):
                crossing_pair_ids = []
                # Make sure with a crossing pair INSIDE
                if direction == _BEFORE:
                    id_start = max_z_crossing_pair_id - 2
                else:
                    id_start = max_z_crossing_pair_id + 2
                # Manage the circular array
                id_start = id_start % len(self.crossing_pairs)
                # Loop over all the crossing pair until it reached the loweest z (lowest priority) crossing_pair
                for i in range(len(self.crossing_pairs)):
                    if direction == _BEFORE:
                        id_current = id_start - i
                    else:
                        id_current = id_start + i
                    # Manage the circular array
                    id_current = id_current % len(self.crossing_pairs)
                    if self.crossing_pairs[id_current].id == min_z_crossing_pair_id:
                        break
                    else:
                        crossing_pair_ids.append(id_current)
                    
                #Make sure to finish with a crossing pair INSIDE (so delete the last which is a OUTSIDE)
                if len(crossing_pair_ids) >= 1:
                    del crossing_pair_ids[-1]
                simple_inside += self._loop_complex_crossing_pairs (hypso_line, crossing_pair_ids)
                        
            # Special case to manage the conflict pair with the lowest priority conflict
            (before,after) = self._get_sibbling(hypso_line, min_max_crossing_pair[0])
            # The best result is to make a conflict with the sibbling crossing having the shortest crossing pair
            if before.p0_p1_line.length < after.p0_p1_line.length:
                current = before
            else:
                current = after
            (before,after) = self._get_sibbling(hypso_line, current)
            
            # Extract the ids of all the complex conflicts
            conflict_ids = []
            for conflict in self.conflicts:
                conflict_ids.append(conflict.before.id)
                conflict_ids.append(conflict.current.id)
                conflict_ids.append(conflict.after.id)
                
            # Check that no crossing pairs have been processed
            if ( before.id  not in conflict_ids and \
                 current.id not in conflict_ids and \
                 after.id   not in conflict_ids ):
                crossing_pair_ids = [before.id, current.id, after.id] 
                simple_inside += self._loop_complex_crossing_pairs (hypso_line, crossing_pair_ids)
        
        else:
            # Open contour
            crossing_pair_ids = []
            # Process from the highest crossing pair to the first crossing pair
            for id_cp in range(max_z_crossing_pair_id-2, -1, -1):
                crossing_pair_ids.append(id_cp)
            simple_inside += self._loop_complex_crossing_pairs (hypso_line, crossing_pair_ids)
                    
            crossing_pair_ids = []
            # Process from the highest crossing pair to the last crossing pair
            for id_cp in range(max_z_crossing_pair_id+2, len(self.crossing_pairs)):
                crossing_pair_ids.append(id_cp)
            simple_inside += self._loop_complex_crossing_pairs (hypso_line, crossing_pair_ids)
            
        # Check the number of crossing minus the number false complex case (crossing that shoul be condidered as simple case )
        if (self._get_nbr_crossing_inside()-simple_inside)%2 == 0:
            self.conflicts = []
            # If the real number of crossing is even reprocess the crossing as unpriorized conflict (like Lakes )
            self._manage_unpriorized_complex_conflict(hypso_line)
                                    
        return
    
    def _manage_unpriorized_complex_conflict(self, hypso_line):
        """Identify the complex conflict for a hydrographic feature that do not need priority conflict management
        
        
        *Parameters*
            - hypso_line: Hypsographic MA_LineString
            
        *Return value*
            None
        """   
        
        crossing_pair_ids = []
        if GenUtil.is_line_closed(hypso_line):
            # For a close line
            max_length = -1.0
            # Find the longest crossing pair of type OUTSIDE
            for crossing_pair in self.crossing_pairs:
                if crossing_pair.type == _OUTSIDE:
                    if crossing_pair.p0_p1_line.length > max_length:
                        max_length = crossing_pair.p0_p1_line.length
                        start_id = crossing_pair.id
            
            # Process all the crossing pair starting from the longest one
            for i in xrange(1, len(self.crossing_pairs)):
                # Process list as a circular list
                current_id = (start_id+i)%len(self.crossing_pairs)
                crossing_pair_ids.append(current_id)

        else:
            #For an open line process crossing pair from first to last
            for crossing_pair in self.crossing_pairs:
                crossing_pair_ids.append(crossing_pair.id)
            
            # Process the complex conflicts
        self._loop_complex_crossing_pairs (hypso_line, crossing_pair_ids)
                                
        return
        
    def _loop_complex_crossing_pairs (self, hypso_line, crossing_pair_ids):
        """Loop over the different crossing pair so that we can create the complex conflict
        
        *Parameters*
            - hypso_line: Hypsographic MA_LineString 
            - crossing_pair_ids: List of crossing pair id to from which to create complex conflict
            
        *Return value*
            - Number of crossing pair that should be considered as part of simple case and not complex case
        """
        
        simple_case = 0
        # We must have at least 3 conflicts (INSIDE-OUTSIDE-INSIDE) in order to process a complex conflict
        if len(crossing_pair_ids) >= 3:
            i = 0
            while i+2 < len(crossing_pair_ids):
                current = self.crossing_pairs[crossing_pair_ids[i+1]]
                (before, after) = self._get_sibbling(hypso_line, current)
                complex_case = True
                for inside in (before, after):
                    # Check if the crossing pair should be identified as part of a simple case istead of a complex case
                    inside_lr0_hypso = hypso_line.project(inside.p0_point)
                    inside_lr1_hypso = hypso_line.project(inside.p1_point)
                    to_cut_hypso_line = HydroUtil.extract_sub_line(inside_lr0_hypso, inside_lr1_hypso, hypso_line)
                    inside_lr0_hydro = self.buf_hydro_line.project(inside.p0_point)
                    inside_lr1_hydro = self.buf_hydro_line.project(inside.p1_point)
                    to_paste_hypso_line = HydroUtil.shortest_perimeter(inside_lr0_hydro, inside_lr1_hydro, self.buf_hydro_line)
                    distance = HydroUtil.mean_max_distance(to_cut_hypso_line, to_paste_hypso_line)
                    if distance < self.max_displacement:
                        # When the distance to edit the line is smaller than the tolerance; the crossing pair should be part of a simple case
                        simple_case += 1
                        complex_case = False
                    
                if complex_case:
                    self._create_complex_conflict(before, current, after)
                    i += 4
                else:
                    i += 2
                    
        return simple_case
           
    def _create_complex_conflict(self, before, current, after):
        """Creates and validate if a complex conflict can be created
        
        *Parameters*
            -before: Conflict preceding the current conflict 
            -current: Current conflict to process
            -after: Conflict following the current conflict
            
        *Return*
            None
        """
        
        conflict = Holder()
        conflict.type = _COMPLEX
        conflict.before = before
        conflict.current  = current
        conflict.after  = after
        if conflict.before.type != _INSIDE and conflict.current.type != _OUTSIDE and conflict.before.type != _INSIDE:
            # This configuration of the conflict is illegal. Should never happen
            print "Integrity error"
            0/0

        self.conflicts.append(conflict)
            
        return
    
    def _create_simple_conflict(self, dummy):
        """Create a simple conflict
        
        Process all the crossing inside and process them as simple case
        
        *Parameters*
            - dummy: No utility (for signature purpose only)
        """
        
        for crossing in self.crossing_pairs:
            if crossing.type == _INSIDE:
                conflict = Holder()
                conflict.current = crossing
                conflict.type = _SIMPLE
                self.conflicts.append(conflict)
        
        return
    
    def _build_conflicts(self, conflict_type, hypso_line):
        """Build the different type of conflict: SIMPLE and COMPLEX
        
        *Parameters*
            - conflict_type: Type of conflict: SIMPLE or COMPLEX
            
        *Return value*
            None
        """    
        
        # Creates the crossing pair
        if self._creates_hydro_hypso_crossing(hypso_line):
            self.conflicts = []
            
            if conflict_type == _SIMPLE:
                self._create_simple_conflict(hypso_line)
                                   
            if conflict_type == _COMPLEX:
                # Order crossing pairs from min to max
                min_max_crossing_pairs = self._priorize_crossing()
                if min_max_crossing_pairs is None:
                    self._manage_unpriorized_complex_conflict(hypso_line)
                else:
                    self._manage_priorized_complex_conflict(hypso_line, min_max_crossing_pairs)

        return
    
    def _next_to_edit(self, hypso_line):
        """Iterator over the conflict to solve
        
        *Parameters*
            - hypso_line: Hypsographic MA_LineString
        """
                
        # Creates and process simple case
        self._build_conflicts(_SIMPLE, hypso_line)
        for next_conflict in self.conflicts:
            self._set_conflict_attribute(hypso_line, next_conflict)
            yield next_conflict

        # Make the buffer just a little bit smaller
        self.buffer.next_sub_iteration()
        
        # Creates and process complex case
        self._build_conflicts(_COMPLEX, hypso_line)
        for next_conflict in self.conflicts:
            self._set_conflict_attribute(hypso_line, next_conflict)
            yield next_conflict
        
        return
                    
    def _set_conflict_status(self, conflict, solved):
        """Modify the type of crossing pair from Inside or Outside to touch
        
        *Parameters*
            - conflict: Conflict to process
            - solved: Boolean: True: conflict solved; False: conflict unsolved
        
        *Return value*
            - None
        
        """
        
        if solved:
            if conflict.type == _SIMPLE:
                conflict.current.type = _TOUCH
            elif conflict.type == _COMPLEX:
                conflict.before.type = _TOUCH
                conflict.after.type = _TOUCH
                
    
    def _add_intersection (self, s_container, intersection_to_smooth, point, buffer_value, to_paste_hypso_line):
        """Add one intersection into the list of intersection to smooth
        
        *Parameters*
            - s_container: SpatialContainer containing the featues
            - intersection_to_smooth: List of intersection to smooth
            - point : Point representing the position 
            - buffer_value: Value of the buffer when the intersection was created
            - to_paste_hypso_line: LineString containing the line to paste
            
        *Return value*
            None
        """
        
        if to_paste_hypso_line.length > self.buffer.max_buffer*.1:
            if to_paste_hypso_line.length > (self.buffer.max_buffer*2.05):
                buffer_value = self.buffer.max_buffer
            else:
                buffer_value =  to_paste_hypso_line.length/2.05
        
            # Create intersection objet
            intersection = Holder()
            intersection.point = point  
            intersection.buffer_value = buffer_value
            intersection.to_paste_hypso_line = to_paste_hypso_line
            info = self._extract_line_to_smooth (s_container, buffer_value, point)
            intersection.to_smooth_line = info.to_smooth_line 
            
            intersection_to_smooth.append(intersection)
        else:
            # The line to paste is to small no need for smoothing this intersection
            pass
        
        return
    
    def _add_intersections(self, s_container, conflict, intersections_to_smooth):
        """Add the diffrent intersection to smooth of a conflict into the list of intersection to smooth
        
        *Parameters*
            - s_container:  SpatialContainer containing the featues
            - conflict: Conflict resolved
            intersection_to_smooth:  List of intersection to smooth
        
        """
        
        if conflict.type == _SIMPLE:
            # Manage Simple conflict
            self._add_intersection(s_container, intersections_to_smooth, conflict.current.p0_point, self.buffer.get_current(), 
                                   conflict.to_paste_hypso_line_a )
            self._add_intersection(s_container, intersections_to_smooth, conflict.current.p1_point, self.buffer.get_current(), 
                                   conflict.to_paste_hypso_line_a)
                        
        else:
            # Manage Complex conflict
            self._add_intersection(s_container, intersections_to_smooth, conflict.after.p1_point, self.buffer.get_current(), 
                                   conflict.to_paste_hypso_line_a)
            self._add_intersection(s_container, intersections_to_smooth, conflict.before.p0_point, self.buffer.get_current(), 
                                   conflict.to_paste_hypso_line_a)

            if conflict.to_paste_hypso_line_b is not None:
                # It's not all the complex conflict whp have a second  line
                self._add_intersection(s_container, intersections_to_smooth, conflict.current.p0_point, self.buffer.get_current(), 
                                       conflict.to_paste_hypso_line_b)
                self._add_intersection(s_container, intersections_to_smooth, conflict.current.p1_point, self.buffer.get_current(), 
                                       conflict.to_paste_hypso_line_b)

                
        return

    def _smooth_line(self, line_to_smooth):
        """Smooth the line by creating a smooth curve between the first and last point of the line to smooth
        
        *Parameters*
            - line_to_smooth: LineString to smooth
            
        *Return value*
            - LineString smoothed line
        """
        
        # Extract the coordinate of the actual point to modify
        p0 = line_to_smooth.coords[0]
        p1 = line_to_smooth.interpolate(.25, normalized=True)
        p1 = p1.coords[0]
        p2 = line_to_smooth.interpolate(.5, normalized=True)
        p2 = p2.coords[0]
        p3 = line_to_smooth.interpolate(.75, normalized=True)
        p3 = p3.coords[0]
        p4 = line_to_smooth.coords[-1]
        
        
        # Calculate the new smoothed point
        p0_p4_mid = GenUtil.mid_point(p0, p4)
        new_coord1 = GenUtil.mid_point(p2, p0_p4_mid)
        
        p0_new_coord0_mid =  GenUtil.mid_point(p0, new_coord1)
        new_coord0 = GenUtil.mid_point(p1, p0_new_coord0_mid)
        
        new_coord1_p4_mid =  GenUtil.mid_point(new_coord1, p4)
        new_coord2 = GenUtil.mid_point(p3, new_coord1_p4_mid)
        
        # Create to new smoothed line
        smoothed_line = LineString((p0, new_coord0, new_coord1, new_coord2, p4))
        
        return smoothed_line
            
    def _extract_line_to_smooth (self, s_container, buffer_value, point):
        """Extract the sub hypso line to smotth
        """
        
        # Object to contain different information for the smoothing
        info = Holder()
        info.to_smooth_line = None
        info.hypso_line = None
        
        # Locate the hypso line located on or almost on the point
        search_pol = point.buffer(buffer_value/10.)
        features = s_container.get_features(bounds=search_pol.bounds)
        hypso_lines = [feature for feature in features if feature.ma_properties[CODE] == HYPSO]
        for info.hypso_line in hypso_lines:
            if point.distance(info.hypso_line) < GenUtil.ZERO:
                # We have found the hypso line wanted
                break

        if info.hypso_line is not None:
            lr_hypso = info.hypso_line.project(point)
            if GenUtil.is_line_closed(info.hypso_line):
                distance_avalaible = info.hypso_line.length/2.
            else:
                distance_avalaible = min((info.hypso_line.length-lr_hypso), lr_hypso)
            
            if distance_avalaible > buffer_value*.1:
                if distance_avalaible < buffer_value:
                    
                    buffer_value = distance_avalaible 
            
                # Determine the start/end of the sub line on the hypso line
                info.lr_hypso_plus_buffer = lr_hypso + buffer_value
                if lr_hypso + buffer_value <= info.hypso_line.length:
                    info.lr_hypso_plus_buffer = lr_hypso + buffer_value
                else:
                    # Manage closed line
                    info.lr_hypso_plus_buffer = (lr_hypso + buffer_value) - info.hypso_line.length
                    
                if  lr_hypso - buffer_value >= 0.0:
                    info.lr_hypso_minus_buffer = lr_hypso - buffer_value
                else:
                    # Manage closed line
                    info.lr_hypso_minus_buffer = info.hypso_line.length - (buffer_value-lr_hypso)
                    
                # Extract sub hypso line to smooth
                info.to_smooth_line = HydroUtil.extract_sub_line(info.lr_hypso_minus_buffer, info.lr_hypso_plus_buffer, 
                                                                 info.hypso_line)
            else:
                # To small nothing to smooth
                pass
                        
        return info

    
    def _smooth_intersections (self, s_container, intersections_to_smooth):
        
        for rec in intersections_to_smooth:
            buffer_value = rec.buffer_value
            #Find the hypso line on top of the intersection
            # Extract the hypso lines near the new sub hypso line
            info = self._extract_line_to_smooth(s_container, buffer_value, rec.point)
            
            if info.to_smooth_line is not None and \
               rec.to_smooth_line is not None and \
               info.to_smooth_line.almost_equals(rec.to_smooth_line, 4):
                smoothed_hypso_line = self._smooth_line(info.to_smooth_line)
                
                to_keep_hypso_line = HydroUtil.extract_sub_line(info.lr_hypso_plus_buffer, info.lr_hypso_minus_buffer, info.hypso_line)
                to_keep_hypso_line = HydroUtil.multi_to_list (to_keep_hypso_line)
                new_hypso_line = HydroUtil.merge_lines(to_keep_hypso_line +  [smoothed_hypso_line])                        
                info.hypso_line.update_coords(new_hypso_line.coords, s_container)
            else:
                # The intersection point is not on the line nothing to do
                print "***Probl�mes de ligne...."
                pass
                            
        return
                
    def _resolve (self, hypso_line, s_container, stats):
        """
        """
           
        intersections_to_smooth = []
        for conflict in self._next_to_edit(hypso_line):
            conflict_resolved = False
            if self._is_displacement_allowed (conflict.distance, self.max_displacement, hypso_line):
                # Check the constraints for the new_hypso_line_a
                if ( self._check_constraints (s_container, self.hydro_feature, hypso_line, conflict.to_paste_hypso_line_a, 
                                              conflict.to_keep_hypso_line_a, conflict.new_hypso_line_a, self.buffer) ):
                    # Check the constraints for the new_hypso_line_b
                    if self._check_new_hypso_line_b(hypso_line, conflict, s_container):
                        conflict_resolved = True
                
                if conflict_resolved:
                    
                    # Update the coordinate of the hypso line
                    hypso_line.update_coords(conflict.new_hypso_line_a.coords, s_container)
                    for new_hypso_line_b in conflict.new_hypso_lines_b:
                        # For complex case there is a new hypso line to create
                        # Create a new MA_LineString feature
                        new_hypso_line = MA_LineString(list(new_hypso_line_b.coords))
                        # Copy the attributes
                        new_hypso_line.ma_properties = hypso_line.ma_properties.copy()
                        new_hypso_line.ma_properties[_ORIENTATION] = HydroUtil.get_orientation(list(new_hypso_line.coords))
                        # Add in the spatial container
                        s_container.add_feature(new_hypso_line)
                    
                    # Update the statistics
                    if conflict.type == _SIMPLE:
                        stats.add_stats(_INSIDE)
                    else:
                        stats.add_stats(_OUTSIDE)
                         
                    self._add_intersections(s_container, conflict, intersections_to_smooth)
                         
            self._set_conflict_status(conflict, conflict_resolved) 
            
        self._micro_intersection_edition(hypso_line, s_container)
        is_all_edited = self._is_all_edited(hypso_line)
        self._smooth_intersections (s_container, intersections_to_smooth)       

        return is_all_edited
    
    
    def _detect_potential_conflict(self, buf_hydro_pol, s_container):
        """
        """
        
#        buf_hydro_pol = hydro_feature.buffer(buf_value, resolution=3)
        
        features = s_container.get_features(bounds=buf_hydro_pol.bounds)
        potential_hypso_lines = [feature for feature in features if feature.ma_properties[CODE] == HYPSO]
        potential_hypso_lines =  filter(buf_hydro_pol.intersects, potential_hypso_lines)
        
        hypso_lines = []
        for hypso_line in potential_hypso_lines:
            if GenUtil.is_line_closed(hypso_line):
                hypso_lines.append(hypso_line)
            else:
                p0 = Point(hypso_line.coords[0])
                p1 = Point(hypso_line.coords[-1])
                if p0.disjoint(buf_hydro_pol) and p1.disjoint(buf_hydro_pol):
                    hypso_lines.append(hypso_line)
        
        hypso_lines.sort(key=operator.attrgetter("length"))
        
        return hypso_lines
       
    def resolve_hypso_lines (self, s_container, stats):
        """
        """
            
        buf_hydro_pol = self._create_buffer()
        hypso_lines = self._detect_potential_conflict(buf_hydro_pol, s_container)
        
        all_edited = True
        for hypso_line in hypso_lines:
            self.buffer.reset_sub_iteration()
            if self.is_edition_needed(hypso_line):
                edited = self._resolve(hypso_line, s_container, stats)
            else:
                edited = True
                
            all_edited = all_edited and edited 
            
        return all_edited
            
                
    @abc.abstractmethod
    def _check_specific_constraints(self, topology):
        """
        """

    def _check_constraints (self, s_container, hydro_feature, ori_hypso_line, new_sub_hypso_line, old_sub_hypso_line, 
                            new_hypso_line, buf):
        """Abstract method to check the constraint  of an edited hypsographic feature
        """
        
        # Extract the hypso lines near the new sub hypso line 
        features = s_container.get_features(bounds=new_sub_hypso_line.bounds, remove_keys=[ori_hypso_line._sci_id])
        
        # Sort the feature by type
        hypso_features = [feature for feature in features if feature.ma_properties[CODE] == HYPSO]
        lake_features  = [feature for feature in features if feature.ma_properties[CODE] == LAKE ]
        river_features = [feature for feature in features if feature.ma_properties[CODE] == RIVER ]
        stream_features = [feature for feature in features if feature.ma_properties[CODE] == STREAM ]
        
        topology = Holder()
        
        topology.hydro_feature = hydro_feature
        topology.new_hypso_line = new_hypso_line
        topology.buffer = buf
        topology.constraints_ok = True
        topology.hypso_distance = [hypso for hypso in hypso_features if hypso.distance(new_sub_hypso_line) <= buf.get_current() *.05]
        topology.lake_intersect  = filter(new_sub_hypso_line.intersects, lake_features)
        topology.river_intersect = filter(new_sub_hypso_line.intersects, river_features)
        topology.stream_intersect = filter(new_sub_hypso_line.intersects, stream_features)
        
        #Check that the new sub hypso line is not to close from the remaining part of the hypso line
        # Check if the new hypso line is simple
        if not new_hypso_line.is_simple:
            topology.constraints_ok = False
            return False

        if old_sub_hypso_line is not None:
            if HydroUtil.are_lines_to_close(new_sub_hypso_line, old_sub_hypso_line, buf.get_current()):
                topology.constraints_ok = False
                return False
            
        # Check that there are no other hypso line very close from the new sub hypso line
        if len(topology.hypso_distance) != 0:
            topology.constraints_ok = False
            return False

        # Check the specific constraint for the Lake, River and Stream
        if topology.constraints_ok:
            self._check_specific_constraints(topology)
            
        if topology.constraints_ok:
            # Make sure to reset the flag of in order to edit the feature 
            for hydro_feature in topology.reset_edition:
                hydro_feature.ma_properties[_IS_ALL_EDITED]= False
                
        return topology.constraints_ok


    def _orient_line (self, old_hypso_line, new_hypso_line):
        """Reorient a LineString
        
        After being edited the orientation of an hysographic feature can be modified. This routine will
        replace the line in its original orientation
        
        *Parameters*
            - old_hypso_line: LineString of the original hypsographic feature
            - new_hypso_line: LineString of the edited hypsographic feature
            
        *Return*
            - lineString oriented according to the original orientation of the hypsographic feature
        
        """
        
        if GenUtil.is_line_closed(old_hypso_line):
            # For a closed loop we check if the sign of the result of the orientation function
            old_orientation = old_hypso_line.ma_properties[_ORIENTATION]
            new_hypso_coords = list(new_hypso_line.coords)
            new_orientation = HydroUtil.get_orientation(new_hypso_coords)
            if (old_orientation > 0. and new_orientation > 0.) or \
               (old_orientation < 0. and new_orientation < 0.):
                is_oriented = True
            else:
                is_oriented = False
                
        else:
            # for an open line we check if the first coordinate is still the first coordinate
            old_coord_first = old_hypso_line.ma_properties[_FIRST_COORDS]
            new_coord_first = new_hypso_line.coords[0]
            if GenUtil.distance(old_coord_first, new_coord_first) < GenUtil.ZERO:
                is_oriented = True
            else:
                is_oriented = False
                new_hypso_coords = list(new_hypso_line.coords)
        
        if is_oriented:
            line = new_hypso_line
        else:
            # Re orient the hypsographic feature
            new_hypso_coords.reverse()
            line = LineString(new_hypso_coords)
        
        return line                                
    
    def _check_new_hypso_line_b(self, hypso_line, conflict, s_container):
        """
        """
         
        conflict.new_hypso_lines_b = []
        constraints_ok = False
        if conflict.new_hypso_line_b is not None:
        # Process the loop line
            current_delta = self.buffer.get_current_delta()
            buf_loop_areas = conflict.new_hypso_pol_b.buffer(current_delta*(-1.05), resolution=3)
            lst_sub_loop_area = HydroUtil.multi_to_list(buf_loop_areas)
            if len(lst_sub_loop_area) == 1:
                # Let's put back the original Polygon before the buffering operations
                lst_sub_loop_area = [conflict.new_hypso_pol_b]
            else:
                # Let's put at alomost orginal position
                for i in range(len(lst_sub_loop_area)):
                    lst_sub_loop_area[i] = lst_sub_loop_area[i].buffer(current_delta*1.05, resolution=3)
                    
            # Now just keep the valid loop line above the minimum area size
            for sub_loop_area in lst_sub_loop_area:
                if ( isinstance(sub_loop_area, Polygon) and \
                     not sub_loop_area.is_empty  and \
                     sub_loop_area.area >= self.min_area ):
                    # Extract the exterior coords
                    sub_loop_line = LineString(sub_loop_area.exterior.coords)
                    conflict.new_hypso_lines_b.append(sub_loop_line)
                    
            # Check the constraints of each loop line
            # Now just keep the valid loop line above the minimum area size
            for new_hypso_line in conflict.new_hypso_lines_b:
                constraints_ok = True
                if not self._check_constraints (s_container, self.hydro_feature, hypso_line,  new_hypso_line, 
                                                None, new_hypso_line, self.buffer):
                    constraints_ok = False
                    break
                                
        else:
            # There is no hypso_line b we assume everything is OK
            constraints_ok = True
            
        return constraints_ok                
                
                
class LakeConflicts(Conflicts):
    
    def __init__(self,hydro_feature, buf, params):
        
        Conflicts.__init__(self, hydro_feature, buf, params)
    
    def is_edition_needed (self, hypso_line ):
        """
        """
        
        return True
        
    def _is_all_edited(self, hypso_line):
        """Checks if all the conflicts are solved
        
        Whatever the code of the hydrographic feature (LAKE, RIVER, STREAM) when the number of conflicts is 0
        the feature is all edited.
        
        When the number of conflicts remaining is not 0. It depends on the type of hydrographic feature:
                           
        For a Lake
            - If the number of conflicts with the buffered hydrographic feature is not 0 it is always considered as
              not all edited
        
        """
        
        if self.valid_crossing:
            # Count the number of conflicts left to edit
            if self._get_nbr_crossing_inside() == 0:
                all_edited = True
            else:
                all_edited = False
        else:
            # There is a problem and we consider that it's not all edited
            all_edited = False
                                
        return all_edited
    
    def _micro_intersection_edition(self, hypso_line, s_container):
        """No micro crossing needed for the Lake"""
        
        return    
    
    def _priorize_crossing(self):
        """There is no priorization to do on the Lakes"""
        
        return None
    
    def _check_specific_constraints(self, topology):
        """
        """
        
        
        if len(topology.lake_intersect) == 0:
            topology.constraints_ok = True
            topology.reset_edition = topology.river_intersect + topology.stream_intersect
        else:
            topology.constraints_ok = False
            topology.reset_edition = None
        
            
        return topology
    
    def _is_displacement_allowed(self, distance, max_displacement, dummy=None):
        """Method determining if the edition is allowed"""
        
        if distance <= max_displacement:
            edition_allowed = True
        else:
            edition_allowed = False
                
        return edition_allowed

                
class StreamConflicts(Conflicts):
    
    def __init__(self,hydro_feature, buf, params):
        
        Conflicts.__init__(self, hydro_feature, buf, params)
        
    
    def _create_buffer(self):
        """
        """
        
        buf_value = self.buffer.get_current()
        if self.hydro_feature.length >= buf_value*1.1:
            lr_a = buf_value*.99
            lr_b = self.hydro_feature.length - lr_a
            to_buf_feature = HydroUtil.extract_sub_line(lr_a, lr_b, self.hydro_feature)
        else:
            to_buf_feature = self.hydro_feature.interpolate(self.hydro_feature.length/2.)
            buf_value = self.hydro_feature.length/2.01
            
        buf_hydro_pol = to_buf_feature.buffer(buf_value, resolution = 3)
        
        return buf_hydro_pol
    
    def _micro_intersection_edition(self, hypso_line, s_container):
        """Micro edition needed for the Streams"""
        
        if self.buffer.is_last_iteration():
            if self._get_nbr_crossing_inside() == 1:
                intersections = self.hydro_feature.intersection(hypso_line)
                intersections = GenUtil.make_iterable(intersections)
                if len(intersections) >= 2:
                    for crossing_pair in self.crossing_pairs:
                        if crossing_pair.type == _INSIDE:
                            inside_crossing_pair = crossing_pair
                            
                    lr0_hypso = hypso_line.project(inside_crossing_pair.p0_point)
                    lr1_hypso = hypso_line.project(inside_crossing_pair.p1_point)
                    to_keep_hypso_line = HydroUtil.extract_sub_line(lr1_hypso, lr0_hypso, hypso_line)
                    to_keep_hypso_line = HydroUtil.multi_to_list(to_keep_hypso_line)
                    to_paste_hypso_line = LineString((inside_crossing_pair.p0_point.coords[0], inside_crossing_pair.p1_point.coords[0]))
                    new_hypso_line = HydroUtil.merge_lines(to_keep_hypso_line+[to_paste_hypso_line])
                    hypso_line.update_coords(new_hypso_line.coords, s_container)
                    
        return
    
    def _is_all_edited(self, hypso_line):
        """Checks if all the conflicts are solved
        
        Whatever the code of the hydrographic feature (LAKE, RIVER, STREAM) when the number of conflicts is 0
        the feature is all edited.
        
        When the number of conflicts remaining is not 0. It depends on the type of hydrographic feature:
              
        For a River
            - If the number of conflicts with the buffered hydrographic feature is 1 and 
                 the type of conflict is simple and
                 the distance is lower than a minimum treshold
            We consider the River as all edited otherwise it's not all edited
            
            Note: The treshold is used in the case when a contour intersect just a samll region of a River, 
                   we don't consider it as all edited
                           
        """
            
        if self.valid_crossing:
            # Count the number of conflicts left to solve
            if self._get_nbr_crossing_inside() == 0:
                # It's all finished
                all_edited = True
            else: 
                all_edited = False               
#                if not Conflicts.is_edition_needed(self.hydro_feature, hypso_line ):
#                    all_edited = True
#                else:
#                    all_edited = False                                
        else:
            # There is a problem and we consider that it's not all edited
            all_edited = False
            
        return all_edited
                    
    def _is_displacement_allowed(self, distance, max_displacement, hypso_line):
        """Method determining if the edition is allowed for a Stream"""
        
        
        if distance <= max_displacement:
            edition_allowed = True
        else:
            edition_allowed = False
                
        return edition_allowed

    def _identifiy_conflicts(self, hypso_line, min_crossing_pairs):
        """Identifi the Stream conflict
        
        The identification of conclist for the Stream is the same as for the River
        """
        
        self._identifiy_conflicts_stream_river(hypso_line, min_crossing_pairs)
        
        
    def _priorize_crossing(self):
        """Priorize the crossing of type inside on a stream
        
        The priorization of the crossing of the Stream is important because we start to resolve conflicts
        with the crossing have the least priority to terminate with the one the highest priority. In order
        to priorize the conflicts the Stream must be oriented. Uphill conflicts have greater priority than
        downhill conflicts. Each conflict of type Inside receive a priority level
        
        *Parameters*
            None
            
        *Return*
            - List conflicts ordered by priority 
        
        """
    
        if self._get_nbr_crossing_inside() >= 3:
            for crossing_pair in self.crossing_pairs:
                if crossing_pair.type == _INSIDE:
                    mid_point = crossing_pair.p0_p1_line.interpolate(crossing_pair.p0_p1_line.length/2.0)
                    crossing_pair.hydro_lr = self.hydro_feature.project(mid_point)
                    
            lst_hydro_lr = [ crossing_pair for crossing_pair in self.crossing_pairs if crossing_pair.type == _INSIDE]
            lst_hydro_lr.sort(key=operator.attrgetter("hydro_lr"), reverse=True)
            hydro_lr_min_max  = (lst_hydro_lr[0], lst_hydro_lr[-1])
        else:
            # There is only one or two crossing inside nothing to priorize
            hydro_lr_min_max = None
            
        return hydro_lr_min_max  
    
    
    def _loop_complex_crossing_pairs (self, hypso_line, crossing_pair_ids):
        """
        """
        
        # We must have at least 3 conflicts (INSIDE-OUTSIDE-INSIDE) in order to process a complex conflict
        if len(crossing_pair_ids) >= 3:
            i = 0
            while i+2 < len(crossing_pair_ids):
                current = self.crossing_pairs[crossing_pair_ids[i+1]]
                (before, after) = self._get_sibbling(hypso_line, current)
                self._create_complex_conflict(before, current, after)
                i += 4
                    
        return 0  
    
    def _create_simple_conflict(self, hypso_line):
        """
        """
               
        intersect_lrs = []
        hydro_hypso_intersect = self.hydro_feature.intersection(hypso_line)
        hydro_hypso_intersect = GenUtil.make_iterable(hydro_hypso_intersect)
        for intersect in hydro_hypso_intersect:
            if isinstance(intersect, Point):
                # It's a point
                p = intersect
            else:
                # It's a LineString so take one coordinate
                p = Point(intersect.coord[0])
            lr = self.hydro_feature.project(p)
            intersect_lrs.append(lr)
            
        if len(intersect_lrs)% 2 == 1:
            # There is an odd number of intersection so there is priorization to do
            intersect_lrs.sort()
            z_max_lr = intersect_lrs[0]
            z_max_point = self.hydro_feature.interpolate(z_max_lr)
        else:
            # There is an even number of intersection so no priorization to do
            z_max_point = None  
                    
        for crossing in self.crossing_pairs:
            if crossing.type == _INSIDE:
                if ( z_max_point is None or \
                     crossing.p0_p1_line.distance(z_max_point) > GenUtil.ZERO ):
                    conflict = Holder()
                    conflict.current = crossing
                    conflict.type = _SIMPLE
                    self.conflicts.append(conflict)
        
        return
    
    def _check_specific_constraints(self, topology):
        
        # Check if the new hypso line passes through the extremity of the straem
        # In order to check if the new hypso line passes by the extremity of the line we create two lines
        # at the extremity of the stream we a length of just a little over the work buffer.
        hydro_coords = list(topology.hydro_feature.coords)
        extended_first_coord = HydroUtil.extend_coord(hydro_coords[1], hydro_coords[0], topology.buffer.get_current()*1.05)
        extended_last_coord = HydroUtil.extend_coord(hydro_coords[-2], hydro_coords[-1], topology.buffer.get_current()*1.05)
        extremity_a_line = LineString((extended_first_coord, hydro_coords[0]))
        extremity_b_line = LineString((extended_last_coord, hydro_coords[-1]))
        disjoint = topology.new_hypso_line.disjoint(extremity_a_line) and topology.new_hypso_line.disjoint(extremity_b_line) 
        # For a stream feature the edited hypso line must not cross any lake, river and the hyso_line must not pass near the extremity of the stream
        if ( len(topology.lake_intersect) == 0 and \
             len(topology.river_intersect) == 0 and \
             len(topology.stream_intersect) <= 1 and \
             disjoint ):
            if len(topology.stream_intersect) == 0:
                topology.constraints_ok = True
                topology.reset_edition = []
            else:
                buf = Buffer(self.buffer.get_current())
                for dummy in buf.next_iteration():
                    params = Holder()
                    params.max_displacement = self.max_displacement
                    params.resolution = self.resolution
                    params.min_area = self.min_area
                    hydro_line = topology.stream_intersect[0]
                    conflicts = StreamConflicts(hydro_line, buf, params)
                    if conflicts.is_edition_needed(topology.new_hypso_line):
                        topology.constraints_ok = False
                        topology.reset_edition = None
                    else:
                        topology.constraints_ok = True
                        topology.reset_edition = []
                    del buf
                    del hydro_line
                    break
        else:
            topology.constraints_ok = False
            topology.reset_edition = None

        return 

    def is_edition_needed(self, hypso_line):
        """
        Check if this hypsographic and hydrographic features need edition
        
        Because there so many intersections of hydrographic and hypsographic features and because the majority of
        the intersections don't need edition, this routine is doing a fast evaluation to determine
        if  some edition is needed
        """
        
        #Assume some edition is needed
        edition_needed = True
        hydro_hypso_intersection = self.hydro_feature.intersection(hypso_line)
        hydro_hypso_intersection = GenUtil.make_iterable(hydro_hypso_intersection)
        if len(hydro_hypso_intersection) == 1:
            buf_hydro_pol = self.hydro_feature.buffer(self.buffer.get_current(), resolution=3)
            buf_hydro_line = LineString(buf_hydro_pol.exterior.coords)
            buf_hydro_hypso_intersection = buf_hydro_line.intersection(hypso_line)
            buf_hydro_hypso_intersection = GenUtil.make_iterable(buf_hydro_hypso_intersection)
            if ( len(buf_hydro_hypso_intersection) == 2 and \
                 isinstance(buf_hydro_hypso_intersection[0], Point) and \
                 isinstance(buf_hydro_hypso_intersection[1], Point) ):
                lr0_hydro = buf_hydro_line.project(buf_hydro_hypso_intersection[0])
                lr1_hydro = buf_hydro_line.project(buf_hydro_hypso_intersection[1])
                p0_p1_line = HydroUtil.shortest_perimeter(lr0_hydro, lr1_hydro, buf_hydro_line)
                if p0_p1_line.intersects(buf_hydro_line):                
                    edition_needed = False

        return edition_needed

class RiverConflicts(Conflicts):
    
    def __init__(self,hydro_feature, buf, params, s_container_z_point):
        
        self.s_container_z_point = s_container_z_point
        
        Conflicts.__init__(self, hydro_feature, buf, params, s_container_z_point)
        
    def _is_all_edited(self, hypso_line):
        """Checks if all the conflicts are solved
        
        Whatever the code of the hydrographic feature (LAKE, RIVER, STREAM) when the number of conflicts is 0
        the feature is all edited.
        
        When the number of conflicts remaining is not 0. It depends on the type of hydrographic feature:
        
        For the Stream:
            - If the number of conflicts with the non buffered hydrographic feature is 1 we consider the Stream as 
              all edited otherwiese it's not all edited
        
        """
            
        if self.valid_crossing:
            # Count the number of crossing left
            if self._get_nbr_crossing_inside() == 0:
                # There is no more crossing inside. there is no more edition to do
                all_edited = True
            else:
                if self._get_nbr_crossing_inside() == 1:
                    if len(self.conflicts) == 1:
                        self._set_conflict_attribute (hypso_line, self.conflicts[0])
                        if self._is_displacement_allowed(self.conflicts[0].distance, self.max_displacement) :
                            # We should edit this conflict because it is very small
                            all_edited = False
                        else:
                            all_edited = True
                    else:
                        # We never know might be editabe
                        all_edited = False
                else:
                    # Two or more crossing inside. The edition is not terminated
                    all_edited = False
        else:
            # There is a problem and we consider that it's not all edited
            all_edited = False
                    
        return all_edited
    
    def _micro_intersection_edition(self, hypso_line, s_container):
        """No micro edition needed for the River"""
        
        return
    
    def _is_displacement_allowed(self, distance, max_displacement, dummy=None):
        """Method determining if the edition is allowed for a River"""
        
        if self._get_nbr_crossing_inside() == 1:
            if distance <= max_displacement/5.0:
                edition_allowed = True
            else:
                edition_allowed = False
        else:
            if distance <= max_displacement:
                edition_allowed = True
            else:
                edition_allowed = False
                
        return edition_allowed
    
    def _identifiy_conflicts(self, hypso_line, min_max_crossing_pairs):
        """Identifi the Stream conflict
        
        The identification of conclist for the Stream is the same as for the River
        """
        
        self._identifiy_conflicts_stream_river(hypso_line, min_max_crossing_pairs)
        
    def _get_polygon_mean_z(self, pol):
        """Calculate for a polygon the mean Z value using the Z point
        
        *Parameters*
            - pol: Polygon for which to find a Z mean value
            
        *Return*
            - Real value of the Z mean value or None if there is no Z point located inside the polygon
        
        """
        
        # Extract all the Z point located spatialy within the polygon 
        z_points = self.s_container_z_point.get_features(bounds=pol.bounds)
        z_points = filter(pol.contains, z_points)
        if len(z_points) >= 1:
            # Calculate the mean value
            total = 0.
            for z_point in z_points:
                total += z_point.ma_properties[Z]
            mean_z = total / float(len(z_points))
        else:
            mean_z = None
        
        return mean_z
    
    def _priorize_crossing(self):
        """Priorize the crossing of type inside on a River
        
        To priorize the crossing on the river we try (not always easy and 100% satisfaction) to order
        the conflicts that are on the River from downhill to uphill.  Each conflicts on a River is forming
        2 polygons. We eavaluate the Z mean value of each polygon and assign to the conflicts the value
        of the maximum Z. the conflict are ordered according to the Z value.
        
        *Parameters*
            None
            
        *Return*
            - Tuple containing the crossing with the lowest and highest Z values
            
        """
    
        if self._get_nbr_crossing_inside() >= 3:
            for crossing_pair in self.crossing_pairs:
                if crossing_pair.type == _INSIDE:
                    # Only process conflict Inside
                    mean_a = self._get_polygon_mean_z(crossing_pair.half_a_hydro_pol)
                    mean_b = self._get_polygon_mean_z(crossing_pair.half_b_hydro_pol)
                    if mean_a is not None and mean_b is not None:
                        crossing_pair.mean_z = max(mean_a,mean_b)
                    else:
                        # When the conflicts is too small we cannot assign a Z value to the conflict
                        crossing_pair.mean_z = sys.maxint
                    
            lst_mean_z = [ crossing_pair for crossing_pair in self.crossing_pairs if crossing_pair.type == _INSIDE and
                                                                                   crossing_pair.mean_z != sys.maxint]
                               
            if len(lst_mean_z) >= 1:
                lst_mean_z.sort(key=operator.attrgetter("mean_z"))
                mean_min_max  = (lst_mean_z[0], lst_mean_z[-1])
            else:
                # In the event there is no z value to priorize the crossing;
                mean_min_max = None
        else:
            # There is only one or two crossing inside so nothing to priorize
            mean_min_max = None
                    
        return mean_min_max 
    
    def _check_specific_constraints(self, topology):
        """ Check the constraint for the river
        """
        
        # For a river the hypso line must not cross any lake must not cross any river but can cross a stream
        if len(topology.lake_intersect) == 0 and len(topology.river_intersect) == 0:
            topology.constraints_ok = True
            topology.reset_edition = topology.stream_intersect
        else:
            topology.constraints_ok = False
            topology.reset_edition = None
            
        return topology
    
    def is_edition_needed (self, hypso_line ):
        """
        """
        
        return True

class AlgoEditContourHydro(Algorithm):
    """This is the main class for the AlgoEditContourHydro algorithm 
    
    Attributes:
        - params: Parameters of the algorithm
        
    """

    def __init__(self, max_buffer,
                       max_displacement,
                       min_area,
                       resolution=16,
                       debug=False ):
        """Initialize the attributes of an object of the class

        Parameters:
            max_buffer: Maximum buffer to apply to the hydro feature
            max_displacement: Maximum distance that an hypsographic line can be edited
            min_area: During hypsographic edition, hypsographic features can be created, these features will be kept 
                      when the area size is over the min_area
            resolution: Number of decimal used by Shapely to buffer features
            debug: Flag to enable(TRUE)/disable(FLASE) print for debug output 

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        self.params.max_buffer = max_buffer
        self.params.min_area = min_area
        self.params.max_displacement = max_displacement
        self.params.test_crossing_line = True
        self.params.test_simple_line = True 
        self.params.resolution = resolution # Used for test and debug
        self.params.debug = debug
        
        self.stats = Statistics()      
    
    def _clean_up(self, hypso_lines):
        """Delete unnecessary attributes added on the features during the process
        
        Spatial container containing cyclic reference must be delete manually otherwise the
        garbage collector is unable to free the memory
        
        Parameters: None
        
        Return value: None
        
        """
        
        for line in hypso_lines:
            if _IS_ALL_EDITED in line.ma_properties: del line.ma_properties[_IS_ALL_EDITED]
            if _FIRST_COORDS in line.ma_properties: del line.ma_properties[_FIRST_COORDS]
            if _LAST_COORDS in line.ma_properties: del line.ma_properties[_LAST_COORDS]
            if _ORIENTATION in line.ma_properties: del line.ma_properties[_ORIENTATION]

    
    def _free_memory(self, s_container):
        """Delete variables in order to free some memory
         
        Spatial container containes cyclic reference that must be delete manually otherwise the
        garbage collector is unable to free the memory
         
        Parameters: None
         
        Return value: None
         
         """
         
        del s_container
       
    def _edit_hypso_line (self, s_container, hypso_line, new_hypso_line, conflicts, conflict):
        """Method to check the constraint of the new hypsographic line and edit the line
        
       
                     False: The line breaks some constraint
            
        """
       
        if conflict.type == _INSIDE:
            self.stats.add_stats(_INSIDE)
        else:
            self.stats.add_stats(_OUTSIDE)
        hypso_line.update_coords(new_hypso_line.coords, s_container)
                
        return
        
    def _set_conflicts (self, s_container_z_point, hydro_feature, buf):
        
           
        if hydro_feature.ma_properties[CODE] == LAKE:
            conflicts = LakeConflicts(hydro_feature, buf, self.params )
            
        if hydro_feature.ma_properties[CODE] == RIVER:
            conflicts = RiverConflicts(hydro_feature, buf, self.params, s_container_z_point )
            
        if hydro_feature.ma_properties[CODE] == STREAM:
            conflicts = StreamConflicts(hydro_feature, buf, self.params )
        
        return conflicts


    def _add_attributes(self, s_container):
         
        hypso_lines = [feature for feature in s_container.get_features() if feature.ma_properties[CODE] == HYPSO ]
        hydro_features = [feature for feature in s_container.get_features() if feature.ma_properties[CODE] == LAKE or 
                                                                               feature.ma_properties[CODE] == RIVER ]
        
        for hydro_feature in hydro_features:
            # Recode the hydro_feature
            if hydro_feature.ma_properties[CODE] == LAKE:
                pass # It's OK
            elif hydro_feature.ma_properties[CODE] == RIVER and isinstance(hydro_feature, MA_Polygon):
                hydro_feature.ma_properties[CODE] = RIVER
            else:
                hydro_feature.ma_properties[CODE] = STREAM
                 
        for hydro_feature in hydro_features:
            hydro_feature.ma_properties[_IS_ALL_EDITED] = False
         
        for hypso_line in hypso_lines:
            if GenUtil.is_line_closed(hypso_line):
                hypso_line.ma_properties[_ORIENTATION] = HydroUtil.get_orientation(list(hypso_line.coords))
            else:
                hypso_line.ma_properties[_FIRST_COORDS] = hypso_line.coords[0]
                hypso_line.ma_properties[_LAST_COORDS] = hypso_line.coords[-1]         
        
    def _pre_process (self, buffer_value, hydro_features, s_container):
        """
        Pre processing of the hydrographic features
        
        The following operation are done 
            - Delete the hypsographic features that are located completly inside the buffer zone
            - Flag as ALL_EDITED the hydrographic features STREAM that do not need edition; this operation is done for optimization purpose as 
              there are many STREAM and most of them don't need any edition.  It's much faster to identify the streams that do not need edition 
              than to try to edit them 
        """        
        

        for hydro_feature in hydro_features:
        
            # Buffer the hydro feature
#            buf_hydro_pol = hydro_feature.buffer(buffer_value, resolution=self.params.resolution)
            buf_hydro_pol = hydro_feature.buffer(buffer_value, resolution=3)
            buf_outer_ring_coords = list(buf_hydro_pol.exterior.coords)
            hypso_lines_to_del = []
            if isinstance(hydro_feature, MA_Polygon):
                # The hydro feature is a polygon so we have to make a hole in order to search only the "donut" part 
                inner_ring_coords = list(hydro_feature.exterior.coords)
                donut_hydro_pol = Polygon(buf_outer_ring_coords, [inner_ring_coords])
            else:
                # The hydro feature is a line string so no hole to make
                donut_hydro_pol = Polygon(buf_outer_ring_coords)
            
            # Find all the contours that are completly inside the donut
            potential_lines = s_container.get_features(bounds=donut_hydro_pol.bounds )
            hypso_lines = [potential_line for potential_line in potential_lines if potential_line.ma_properties[CODE] == HYPSO]
            hypso_lines_to_del =  filter(donut_hydro_pol.contains, hypso_lines)
            self.stats.add_stats(_DELETED, value=len(hypso_lines_to_del))
            # Delete the hypso line that lie completely in the donut region
            s_container.del_features(hypso_lines_to_del)
            
        return
                                  
    def _detect_potential_hypso(self, hydro_feature, buf_value, s_container):
        """
        """
        
        buf_hydro_pol = hydro_feature.buffer(buf_value, resolution=3)
        
        features = s_container.get_features(bounds=buf_hydro_pol.bounds)
        potential_hypso_lines = [feature for feature in features if feature.ma_properties[CODE] == HYPSO]
        potential_hypso_lines =  filter(buf_hydro_pol.intersects, potential_hypso_lines)
        
        hypso_lines = []
        for hypso_line in potential_hypso_lines:
            if GenUtil.is_line_closed(hypso_line):
                hypso_lines.append(hypso_line)
            else:
                p0 = Point(hypso_line.coords[0])
                p1 = Point(hypso_line.coords[-1])
                if p0.disjoint(buf_hydro_pol) and p1.disjoint(buf_hydro_pol):
                    hypso_lines.append(hypso_line)
        
        hypso_lines.sort(key=operator.attrgetter("length"))
        
        return hypso_lines

    
    def process(self):
        """Main routine for the edition of the contour that are not coherent with the hydrography
        
        
        Parameters: None
            
        """

        s_container = self.load_features()
        
        self._add_attributes (s_container)
        
        # Extract all the hydro features
        features = s_container.get_features()
        lake_features  = [feature for feature in features if feature.ma_properties[CODE] == LAKE ]
        river_features = [feature for feature in features if feature.ma_properties[CODE] == RIVER ]
        stream_features = [feature for feature in features if feature.ma_properties[CODE] == STREAM]
        
        z_point_features = [feature for feature in features if feature.ma_properties[CODE] == Z_POINT]
        s_container.del_features(z_point_features)
        s_container_z_point = SpatialContainer()
        s_container_z_point.add_features(z_point_features)
        
        self._pre_process (self.params.max_buffer, lake_features+river_features+stream_features, s_container)

        buf = Buffer(self.params.max_buffer)
        for iteration, buffer_value in enumerate(buf.next_iteration()):
            i = 0
            # Order by Lake, River and stream
            hydro_features = lake_features + river_features + stream_features
            
            # Only keep the hydro features that need edition
            to_edit_hydro_features = [hydro_feature for hydro_feature in hydro_features if not hydro_feature.ma_properties[_IS_ALL_EDITED] ]
              
            for hydro_feature in to_edit_hydro_features:
                i += 1
                try:
                    msg =  str(iteration) + " / " + str(i) +  "-" + str(len(to_edit_hydro_features))
                    logger.log(msg)
                except:
                    pass

                conflicts = self._set_conflicts (s_container_z_point, hydro_feature, buf)
                hydro_feature.ma_properties[_IS_ALL_EDITED] = conflicts.resolve_hypso_lines (s_container,  self.stats)

                    
        # Delete the hydrographic features
        s_container.del_features(lake_features+river_features+stream_features)
        
        # Clean up attributes
        self._clean_up(s_container.get_features())
        
        # Extract all the features for the output
        self.extract_features_out(s_container.get_features())
        
        self._free_memory(s_container)
        
        GenUtil.print_debug (self.params, "End of %s algorithm" %(_CLEAN_HYDRO))
        