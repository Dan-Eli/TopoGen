#! /usr/bin/env python
# -*- coding: UTF-8 -*-
#####################################################################################################################################

"""
    This algorithm change the position of the contours in order to to keep the consistency with the with the Lakes, Rivers and Streams.
    The following constraint are implemented:   
        - A contour cannot cross a lake
        - A contour can cross a River or Stream only once
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
        Works better when the contours lines are closed.
        Works better when the permanent and intermittent part of a Lake or River are merges together
        Works better when the stream belonging to the major stream are not broken
          
"""
__revision__ = "--REVISION-- : $Id: algo_talweg_correction.py 472 2011-09-29 13:39:22Z dpilon $"

#####################################################################################################################################

import sys, operator
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection
from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm, Holder, InternalError
                         
########################################################
# Public constant
ALTI = "ALTI"
AREA_MINIMUM = "AREA MINIMUM"
CODE = "CODE"
HYPSO = "HYPSO"   # Code for hypsographic features
LAKE = "LAKE"
RIVER = "RIVER"
STREAM = "STREAM"
Z = "Z"
Z_POINT = "Z_POINT"

#_BUF_HYDRO_POL = "Buf hydro pol" 
#_BUF_HYDRO_EXT_COORDS = "Buf hydro ext coords"

_ALGO = "Algo"
_CLEAN_HYDRO = "Clean Hydrography"
_COMPLEX = "Complex"
_DELETED = "Hypso deleted"
_FIRST = "First"
_LAST = "Last"
_FIRST_COORDS = "First coords"
_NOT_PROCESSED = "Hypso not processed"
_INSIDE = "Inside"
#_INFLEXION = "Inflexion"
_IS_CLOSED = "Is closed?"
_IS_ALL_EDITED = "Is all edited?"
_LAST_COORDS = "Last coords"
_AFTER = "After"
_NOT_SIMPLE = "Not simple"
_ORIENTATION = "Orientation"
_OUTSIDE = "Outside"
_BEFORE = "Before"
_PERIMETER = "Perimeter"
#######_PRIORITY = "Priorite"
#_SPATIAL_PRIORITY = "Spatial priorite"
_SIMPLE = "Simple"
_COMPLEX = "Complex"
_TOUCH = "Touch"

class HydroUtil:
    """
    In this class we find many little utility routines.
    Some of these routine should be moved in the more general classs GenUtil
    """
    
    @staticmethod
    def normalize_one_measure (measure, measure_min, measure_max):
        """ Normalize the list of value between 0 and 1
        
        Parameters:
            - measure: Measure to normalize
            - measure_min: Minimum value for this measure
            - measure_max: Maximum value for this measure
            
        Return value:
            Float [0..1]: Normalize value 
        
        """
        
        if (measure_max-measure_min !=  0.):
            norm_measure = (measure-measure_min) / (measure_max-measure_min)
        else:
            norm_measure = 0. 
        
        return norm_measure
    
    @staticmethod
    def shortest_perimeter(p0_lr, p1_lr, line):
        """Extracts the line forming the shortest path between two points on a closed line

        *Parameters*
            - p0_lr_hydro: First linear reference on the closed line
            - p1_lr_hydro: Second linear reference on the closed line

        *Returns*:
            - LineString describing the shortest path on the closed line
            
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
         
        sub_line = line
         
        if lr_a < lr_b:
             
            if lr_b >= line.length:
                None
            else:
                tmp_lines = HydroUtil.cut_line_distance(sub_line, lr_b)
                sub_line = tmp_lines[0]
         
            if lr_a <= 0.0:
                None
            else:
                tmp_lines = HydroUtil.cut_line_distance(sub_line, lr_a)
                sub_line = tmp_lines[1]
                 
        else:
             
            sub_lines = []
            if lr_b >= line.length:
                None
            else:
                tmp_lines = HydroUtil.cut_line_distance(sub_line, lr_a)
                sub_lines.append(tmp_lines[1])
                 
            if lr_a <= 0:
                None
            else:
                tmp_lines = HydroUtil.cut_line_distance(sub_line, lr_b)
                sub_lines.append(tmp_lines[0])
                 
             
            sub_line = HydroUtil.merge_lines(sub_lines)
                         
        return sub_line
    
    @staticmethod
    def multi_to_list (multi_feature):
        """Transform a simple or a multi geometry into a list of geometry
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
        """Cut the line at specific distance along the line
    
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
            lst_coords_out = [coords]
        else:
            lst_coords_out = [coords]  # This line will prevent errors on non simple lines
            last_i = len(coords)-1
            for i, p in enumerate(coords):
                if (i < last_i):
                    if i == 0:
                        pd = 0.
                    else:
                        pd += GenUtil.distance(coords[i-1], coords[i])    
                    
                else:
                    pd = line.length
                if pd == distance:
                    lst_coords_out =  [coords[:i+1], coords[i:] ]
                if pd > distance:
                    cp = line.interpolate(distance)
                    lst_coords_out =  [list(coords[:i]) + [(cp.x, cp.y)], [(cp.x, cp.y)] + list(coords[i:])]
                    break
        
        lines_out = []
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
        """
        """
        coords = []
        for line in lines:
            coords.append(list(line.coords))
         
        for i in range(1, len(coords)): 
            if len(lines) == 1: break 
            for j in range(len(coords)-1, 0,-1):
                merge = False
                if round(coords[0][0][0], 4) == round(coords[j][0][0], 4) and \
                   round(coords[0][0][1], 4) == round(coords[j][0][1], 4):
                    # Add start first line with start current line
                    coords[0].reverse()
                    coords[0] += coords[j][1:]
                    merge = True
                elif round(coords[0][0][0], 4) == round(coords[j][-1][0], 4) and \
                     round(coords[0][0][1], 4) == round(coords[j][-1][1], 4):
                    # Add start first line with end current line
                    coords[0] = coords[j] + coords[0][1:]
                    merge = True
                elif round(coords[0][-1][0], 4) == round(coords[j][0][0], 4) and \
                     round(coords[0][-1][1], 4) == round(coords[j][0][1], 4):
                    # Add end first line with start current line
                    coords[0] += coords[j][1:]
                    merge = True
                elif round(coords[0][-1][0], 4) == round(coords[j][-1][0], 4) and \
                     round(coords[0][-1][1], 4) == round(coords[j][-1][1], 4):
                    # Add end first line with end current line
                    coords[j].reverse()
                    coords[0] += coords[j][1:]
                    merge = True
                    
                if round(coords[0][0][0], 4) == round(coords[j][-1][0], 4) and \
                   round(coords[0][0][1], 4) == round(coords[j][-1][1], 4):
                    coords[0][0] = coords[0][-1]
                     
                if merge:
                    del coords[j]
                    break
                
        if len(coords) == 1:
            if round(coords[0][0][0], 4) == round(coords[0][-1][0], 4) and \
               round(coords[0][0][1], 4) == round(coords[0][-1][1], 4):
                coords[0][0] = coords[0][-1]
            line = LineString(coords[0])
        else:
            line = MultiLineString(coords)
             
        return line
    
    @staticmethod
    def are_lines_to_close (source, targets, distance_min):
        
        colinear = False
        if source.length > distance_min *6.00001:
            sub_sources = GenUtil.cut_line_distance(source, source.length - 3.*distance_min)
            sub_sources = GenUtil.cut_line_distance(sub_sources[0], 3.*distance_min)
            sub_source = sub_sources[1]
            
            for target in targets:
                if sub_source.distance(target) < distance_min/2.1:
                    colinear = True
                    break
                
        return colinear
    
    @staticmethod
    def extend_line(coord_a, coord_b, length):
        
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
        
    def get_stats (self, type=None):
        """Extract the current statistics and build  a list of string that forms the statistical message"
        
        Parameters:
            type: Unused... only there for compatibility purpose
        
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
    
    _MAX_SUB_ITERATION = 250.
    
    def __init__(self, max_buffer, max_iteration):
        
        self._i_iteration = None
        self._i_sub_iteration = None
        self._max_buffer = float(max_buffer)
        self._max_iteration = float(max_iteration)
        
        self._delta_buffer = self._max_buffer / self._max_iteration
        self._sub_delta_buffer = ((self._delta_buffer*.75) / Buffer._MAX_SUB_ITERATION)
        
    def next_iteration(self):
        
        if self._i_iteration is None:
            self._i_iteration = 0
        else:
            self._i_iteration += 1
        self._i_sub_iteration = 0
        if self._i_iteration > self._max_iteration:
            0/0
        
    def next_sub_iteration(self):
        
        self._i_sub_iteration += 1
        if self._i_sub_iteration > Buffer._MAX_SUB_ITERATION:
            0/0
            
    def reset_sub_iteration(self):
        
        self._i_sub_iteration = 0
            
    def get_number_iteration_left(self):
        
        return self._max_iteration - self._i_iteration
    
    def get_current(self):
        
        return self._max_buffer - (self._i_iteration*self._delta_buffer) - (self._i_sub_iteration*self._sub_delta_buffer) 
    
    def get_current_delta(self):
        
        total_delta_buffer = self._max_buffer - self.get_current()
        if total_delta_buffer == 0.:
            total_delta_buffer = GenUtil.ZERO
            
        return total_delta_buffer
        

class Conflicts(object):
    
    
    def __init__(self, hydro_feature, hypso_line, buffer, s_container_z_point, max_displacement, resolution):
        

        
        self.buffer = buffer
        self.s_container_z_point = s_container_z_point
        self.max_displacement = max_displacement
        self.resolution = resolution
        self.hypso_line = hypso_line
        self.hydro_feature = hydro_feature
        self.buffer.reset_sub_iteration()
    
    
    def _analyse_crossing (self):
                 
        self.crossing_pairs = []
        if len(self.crossing_points) >= 2:
            linear_refs = []
            for point in self.crossing_points:
                linear_refs.append(self.hypso_line.project(point))
            linear_refs.sort()
             
            nbr_crossing_pair = len(linear_refs)
     
            lst_first_second = []
            for i in range(nbr_crossing_pair-1):
                lst_first_second.append((i,i+1))
            if self.hypso_line.ma_properties[_IS_CLOSED]:
                lst_first_second.append((i+1,0))
           
            for first_second in lst_first_second:
                first = first_second [0]
                second = first_second [1]
                p0 = self.hypso_line.interpolate(linear_refs[first])
                p1 = self.hypso_line.interpolate(linear_refs[second])
                crossing_pair = Holder()
                self._set_crossing_pair(crossing_pair, p0, p1)
                      
                self.crossing_pairs.append(crossing_pair)
         
            # Give unique ID to each crossing_pair
            for i, crossing_pair in enumerate(self.crossing_pairs):
                crossing_pair.id = i
        
        return

    def _distance (self, conflict):
        """Calculate the mean distance between two LineString
         
        Parameters:
            line_a, line_b: LineString object to calculate the distace
             
        Return value:
            Mean distance between the 2 LineString
         
        """
         
        samples = 30
        sub_samples = samples/6
             
        distances = []
        for i in range(31):
            ratio = i/30.
            point_to_cut_line = conflict.to_cut_hypso_line.interpolate(ratio, normalized=True)
            point_to_paste_line = conflict.to_paste_hypso_line.interpolate(ratio, normalized=True)
            distances.append(point_to_cut_line.distance(conflict.to_paste_hypso_line))
            distances.append(point_to_paste_line.distance(conflict.to_cut_hypso_line))
         
        total = 0
        distances.sort(reverse=True)
        for i in range(sub_samples):
            total += distances[i]
         
        mean = total/sub_samples
             
        return mean
    
    def _validate_inside_outside(self):
            
        valid_crossing = True
        crossing_types = []
        
        # Check the type of crossing
        for crossing_pair in self.crossing_pairs:
            # All the crossing must be Inside,  Outside
            if crossing_pair.type in (_INSIDE, _OUTSIDE):
                crossing_types.append(crossing_pair.type)
            else:
                valid_crossing = False
        
        if valid_crossing and len(self.crossing_pairs) >= 1:
                                
            if self.hypso_line.ma_properties[_IS_CLOSED]:
                # We must always have an even number of crossing
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
                # For an open line first and last crossing are Inside
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
        
    def _set_crossing_pair (self, crossing_pair, p0, p1):
        
        crossing_pair.p0_point = p0
        crossing_pair.p1_point = p1
        if ( GenUtil.distance (crossing_pair.p0_point.coords[0], crossing_pair.p1_point.coords[0]) < 10.*GenUtil.ZERO ):
            crossing_pair.type = _COMPLEX
        else:
            crossing_pair.lr0_hypso = self.hypso_line.project(p0)
            crossing_pair.lr1_hypso = self.hypso_line.project(p1)
            crossing_pair.p0_p1_line = HydroUtil.extract_sub_line(crossing_pair.lr0_hypso, crossing_pair.lr1_hypso, self.hypso_line)
            mid_point = crossing_pair.p0_p1_line.interpolate(crossing_pair.p0_p1_line.length/2.)

            if mid_point.within(self.buf_hydro_pol):
                crossing_pair.type = _INSIDE
                # Add some measure
                lr0_hydro = self.buf_hydro_line.project(p0)
                lr1_hydro = self.buf_hydro_line.project(p1)
                half_a_hydro_line = HydroUtil.extract_sub_line(lr0_hydro, lr1_hydro, self.buf_hydro_line)
                half_a_hydro_close_line = HydroUtil.merge_lines([half_a_hydro_line, crossing_pair.p0_p1_line])
                crossing_pair.half_a_hydro_pol = Polygon(half_a_hydro_close_line)
                 
                half_b_hydro_line = HydroUtil.extract_sub_line(lr1_hydro, lr0_hydro, self.buf_hydro_line)
                half_b_hydro_close_line = HydroUtil.merge_lines([half_b_hydro_line, crossing_pair.p0_p1_line])
                crossing_pair.half_b_hydro_pol = Polygon(half_b_hydro_close_line)
                
            else:
                crossing_pair.type = _OUTSIDE
             
        return

    
    def _get_sibbling(self, crossing_pair_current):
    
        
        current = crossing_pair_current.id
        before = current - 1
        after = current + 1
        
        underflow = -1
        overflow = len(self.crossing_pairs)
        
        if before == underflow:
            if self.hypso_line.ma_properties[_IS_CLOSED]:
                crossing_pair_before = self.crossing_pairs[-1]
            else:
                crossing_pair_before = None
        else:
            crossing_pair_before = self.crossing_pairs[before]
            
        if after == overflow:
            if self.hypso_line.ma_properties[_IS_CLOSED]:
                crossing_pair_after = self.crossing_pairs[0]
            else:
                crossing_pair_after = None
        else:
            crossing_pair_after = self.crossing_pairs[after]
            
        return (crossing_pair_before, crossing_pair_after)
        
    
# 
#     def _is_priorisable(self):
#         """
#         """
#         
#         if self._get_nbr_crossing_left() <= 1:
#             priorisable = False
#         else:
#             if self.buffer.get_number_iteration_left() <= 2:
#                 priorisable = False
#             else:
#                 priorisable = True
#                 
#         return priorisable
        
    
    def _is_stream_editable(self):
        """
        """
        
        if self.hydro_feature.ma_properties[CODE] == STREAM:
            nbr_crossing_inside = self._get_nbr_crossing_inside()
            intersections = self.hydro_feature.intersection(self.hypso_line)
            intersections = GenUtil.make_iterable(intersections)
#            print " len(intersections), nbr_crossing_inside:", len(intersections), nbr_crossing_inside
            if len(intersections) == nbr_crossing_inside:
                is_editable = True
            else:
                is_editable = False
#            if len(intersections)%2 == 0:
#                is_editable = True
#            else:
#                if  len(intersections) == nbr_crossing_inside:
#                    is_editable = True
#                else:
#                    is_editable = False
        else:
            is_editable = True
            
#             intersections = self.hydro_feature.intersection(self.hypso_line)
#             intersections = GenUtil.make_iterable(intersections)
#             lst_geometry = [1 for geometry in intersections]
#             if len(lst_geometry) == 1:
#                 is_editable = False
                
        return is_editable
    
#     def is_all_edited(self):
#         
#         self._set_conflict_status()
#         
#         return self._all_edited
        
    def _get_nbr_crossing_inside(self):
        """
        """
        nbr_crossing_inside = 0
        for crossing_pair in self.crossing_pairs:
            if crossing_pair.type == _INSIDE:
                nbr_crossing_inside += 1
                
        return nbr_crossing_inside
    
    def is_all_edited(self):
        """
        """
        
        if self.valid_crossing:
            # Count the number of conflicts left to solve
            nbr_conflicts_left = self._get_nbr_crossing_inside()
            
            if nbr_conflicts_left == 0:
                # It's all finished
                all_edited = True
#                self._priorisable = False
            else:
                if self.hydro_feature.ma_properties[CODE] == STREAM:
                    # It's a stream so let's check with the original line how many crossing are left
                    intersections = self.hydro_feature.intersection(self.hypso_line)
                    intersections = GenUtil.make_iterable(intersections)
                    if len(intersections) == 1:
                        all_edited = True
                    else:
                        # There are 0 or more than 2 crossing we have to clean the line 
                        all_edited = False # There is more than one let's continue
                if self.hydro_feature.ma_properties[CODE] == RIVER:
                    if nbr_conflicts_left == 1:
                        if len(self.conflicts) == 1 and  \
                            self.conflicts[0].type == _SIMPLE and \
                            self.conflicts[0].distance < self.max_displacement/5.0 :
                            all_edited = False
                        else:
                            all_edited = True
                    else:
                        all_edited = False
                if self.hydro_feature.ma_properties[CODE] == LAKE:
                    # There must be no conflict left for lakes
                    all_edited = False                  
        else:
            # If the crossing are  not valid there is an error with the crossing so akk the crossing are not edited
            all_edited = False
                    
        return all_edited
                    
#     def _get_conflict_mid_coord(self, conflict):
#         
#         # Find a representative point for the conflict
#         if conflict.type == _SIMPLE:
#             mid_coord = GenUtil.mid_point(conflict.current.p0_point.coords[0], conflict.current.p1_point.coords[0])
#         else:
#             mid_coord_before = GenUtil.mid_point(conflict.before.p0_point.coords[0], conflict.before.p1_point.coords[0])
#             mid_coord_after     = GenUtil.mid_point(conflict.after.p0_point.coords[0], conflict.after.p1_point.coords[0])
#             mid_coord = GenUtil.mid_point(mid_coord_before, mid_coord_after)
#                 
#         return mid_coord                            
    

                         
        self.conflicts.sort(key=operator.attrgetter('priority'))

        return
    
    def _set_conflict_attribute (self, conflict):
        
        current_lr0_hypso = self.hypso_line.project(conflict.current.p0_point)
        current_lr1_hypso = self.hypso_line.project(conflict.current.p1_point)
        current_lr0_hydro = self.buf_hydro_line.project(conflict.current.p0_point)
        current_lr1_hydro = self.buf_hydro_line.project(conflict.current.p1_point)
        conflict.to_cut_hypso_line = HydroUtil.extract_sub_line(current_lr0_hypso, current_lr1_hypso, self.hypso_line)
        
        if conflict.type == _SIMPLE:
            # Extract some values for _SIMPLE case
            conflict.to_keep_hypso_line = HydroUtil.extract_sub_line(current_lr1_hypso, current_lr0_hypso, self.hypso_line)
            conflict.to_paste_hypso_line = HydroUtil.shortest_perimeter(current_lr0_hydro, current_lr1_hydro, self.buf_hydro_line)
            conflict.loop_area = Polygon()
        else:
            # Extract some values for _COMPLEX case
            
            before_lr0_hypso = self.hypso_line.project(conflict.before.p0_point)
            before_lr1_hypso = self.hypso_line.project(conflict.before.p1_point)
            before_lr0_hydro = self.buf_hydro_line.project(conflict.before.p0_point)
            before_lr1_hydro = self.buf_hydro_line.project(conflict.before.p1_point)
            after_lr0_hypso = self.hypso_line.project(conflict.after.p0_point)
            after_lr1_hypso = self.hypso_line.project(conflict.after.p1_point)
            after_lr0_hydro = self.buf_hydro_line.project(conflict.after.p0_point)
            after_lr1_hydro = self.buf_hydro_line.project(conflict.after.p1_point)
            
            conflict.to_keep_hypso_line = HydroUtil.extract_sub_line(after_lr1_hypso, before_lr0_hypso, self.hypso_line)    
            conflict.to_paste_hypso_line = HydroUtil.shortest_perimeter(before_lr0_hydro, after_lr1_hydro, self.buf_hydro_line)
            to_close_loop_line =  HydroUtil.shortest_perimeter(current_lr0_hydro, current_lr1_hydro, self.buf_hydro_line )
            loop_line = HydroUtil.merge_lines([conflict.to_cut_hypso_line, to_close_loop_line])
            if isinstance(loop_line, LineString) and GenUtil.is_line_closed(loop_line):
                conflict.loop_area = Polygon(loop_line.coords) 
            else:
                # Probably to small area we can assume an empty Polygon
                conflict.loop_area = Polygon()
                                                                        
        conflict.distance = self._distance(conflict)    
        conflict.to_keep_hypso_line = HydroUtil.multi_to_list (conflict.to_keep_hypso_line)
    
        return
    
    def _extract_crossing(self):
        
        valid = False
        while not valid:
            # Find the crossing
#            valid = self._find_crossing(self.buffer.get_current())
            # Buffer of the polygon
            
            buf_hydro_pol = self.hydro_feature.buffer(self.buffer.get_current(), resolution=self.resolution)
            # Make a LineString of the polygon exterior
            self.buf_hydro_line = LineString(buf_hydro_pol.exterior.coords)
            # Strip the polygon of any polygon interiors
            self.buf_hydro_pol = Polygon(buf_hydro_pol.exterior.coords)
            # Detect crossing between hydro and hypso feature
            self.crossing_points = self.buf_hydro_line.intersection(self.hypso_line)
            self.crossing_points = GenUtil.make_iterable(self.crossing_points)
            lst_points = [1 for geometry in self.crossing_points if not isinstance(geometry, Point)]
            if len(lst_points) == 0:
                self._analyse_crossing()
                valid = self._validate_inside_outside()
            else:
                valid = False
                
            self.buffer.next_sub_iteration()
            self.valid_crossing = valid
            
        return valid

    def _get_polygon_mean_z(self, pol):
        
        z_points = self.s_container_z_point.get_features(bounds=pol.bounds)
        z_points = filter(pol.contains, z_points)
        if len(z_points) >= 1:
            total = 0.
            for z_point in z_points:
                total += z_point.ma_properties[Z]
            mean_z = total / float(len(z_points))
        else:
            mean_z = None
        
        return mean_z
            
    
    def _priorize_crossing_stream(self):
    
        if len(self.crossing_pairs) >= 1:
            for crossing_pair in self.crossing_pairs:
                if crossing_pair.type == _INSIDE:
                    mid_point = crossing_pair.p0_p1_line.interpolate(crossing_pair.p0_p1_line.length/2.0)
                    crossing_pair.hydro_lr = self.hydro_feature.project(mid_point)
                    
            lst_hydro_lr = [ crossing_pair for crossing_pair in self.crossing_pairs if crossing_pair.type == _INSIDE]
            lst_hydro_lr.sort(key=operator.attrgetter("hydro_lr"), reverse=True)
            hydro_lr_min_max  = (lst_hydro_lr[0], lst_hydro_lr[-1])
        else:
            hydro_lr_min_max = None
            
        return hydro_lr_min_max
            
    def _priorize_crossing_river_opt1(self):
    
        for crossing_pair in self.crossing_pairs:
            if crossing_pair.type == _INSIDE:
                mean_a = self._get_polygon_mean_z(crossing_pair.half_a_hydro_pol)
                mean_b = self._get_polygon_mean_z(crossing_pair.half_b_hydro_pol)
                if mean_a is not None and mean_b is not None:
                    crossing_pair.mean_z = max(mean_a,mean_b)
                else:
                    crossing_pair.mean_z = sys.maxint
                
        lst_mean_z = [ crossing_pair for crossing_pair in self.crossing_pairs if crossing_pair.type == _INSIDE and
                                                                               crossing_pair.mean_z != sys.maxint]
                           
        if len(lst_mean_z) >= 1:
            lst_mean_z.sort(key=operator.attrgetter("mean_z"))
            mean_min_max  = (lst_mean_z[0], lst_mean_z[-1])
        else:
            mean_min_max = None
        
        return mean_min_max 
    
    def _priorize_crossing_river_opt2(self):
    
        for crossing_pair in self.crossing_pairs:
            if crossing_pair.type == _INSIDE:
                buffer_inside = crossing_pair.p0_p1_line.buffer(crossing_pair.p0_p1_line.length*2.5)
                z_points = self.s_container_z_point.get_features(bounds=buffer_inside.bounds)
                z_points = filter(buffer_inside.contains, z_points)
                z_points = filter(self.hydro_feature.contains, z_points)
                samples = []
                for z_point in z_points:
                    sample = Holder()
                    sample.distance = crossing_pair.p0_p1_line.distance(z_point)
                    sample.z = z_point.ma_properties[Z]
                    samples.append(sample)
                    
                if len(samples) >= 1:
                    samples.sort(key=operator.attrgetter("distance"))
                    samples = samples[0:29]
                    samples_z = [sample.z for sample in samples]
                    crossing_pair.mean_z = sum(samples_z)/len(samples_z)
                else:
                    crossing_pair.mean_z = sys.maxint
                    
        lst_sort = [ crossing_pair for crossing_pair in self.crossing_pairs if crossing_pair.type == _INSIDE and
                                                                               crossing_pair.mean_z != sys.maxint] 
        
        lst_sort.sort(key=operator.attrgetter("mean_z"))
        
        return lst_sort
    
    def _identifiy_conflicts_stream_river(self):
        
        
        if self.hydro_feature.ma_properties[CODE] == RIVER:
            min_max_mean = self._priorize_crossing_river_opt1()
#            min_max_mean = self._priorize_crossing_river_opt2()
        else:
            min_max_mean = self._priorize_crossing_stream()
        
        if min_max_mean is not None:
            
            min_z_crossing_pair = min_max_mean[0]
            max_z_crossing_pair = min_max_mean[-1]

            if self.hypso_line.ma_properties[_IS_CLOSED]:
                (before, after) = self._get_sibbling(min_z_crossing_pair)
                if before.p0_p1_line.length > after.p0_p1_line.length:
                    set_current = before
                else:
                    set_current = after
            
                current = set_current
                for i_dummy in xrange(self._get_nbr_crossing_inside()/2):
                    (dummy, current_plus_1) = self._get_sibbling(current)
                    (dummy, current_plus_2) = self._get_sibbling(current_plus_1)
                    (dummy, current_plus_3) = self._get_sibbling(current_plus_2)
                    if self._build_conflict(current_plus_1, current_plus_2, current_plus_3, max_z_crossing_pair.id):
                        (dummy, current) = self._get_sibbling(current_plus_3)
                    else:
                        break
                
                current = set_current        
                for i_dummy in xrange(self._get_nbr_crossing_inside()/2):
                    (current_minus_1, dymmy) = self._get_sibbling(current)
                    (current_minus_2, dymmy) = self._get_sibbling(current_minus_1)
                    (current_minus_3, dymmy) = self._get_sibbling(current_minus_2)
                    if self._build_conflict(current_minus_3, current_minus_2, current_minus_1, max_z_crossing_pair.id):
                        (current, dummy) = self._get_sibbling(current_minus_3)
                    else:
                        break
            
            else:
                
                i = 0
                while  i+3 <= len(self.crossing_pairs):
                    before = self.crossing_pairs[i]
                    current = self.crossing_pairs[i+1]
                    after = self.crossing_pairs[i+2]
                    if self._build_conflict(before, current, after, max_z_crossing_pair.id):
                        i += 4
                    else:
                        break
                        
                i = len(self.crossing_pairs)-1
                while  i-3 >= 0:
                    after = self.crossing_pairs[i]
                    current = self.crossing_pairs[i-1]
                    before = self.crossing_pairs[i-2]
                    if self._build_conflict(before, current, after, max_z_crossing_pair.id):
                        i -= 4
                    else:
                        break
                    
        else:
            
            self._identifiy_conflicts_lake()
                        
        return
    
    
    def _identifiy_conflicts_lake(self):
        
        
        if self._get_nbr_crossing_inside() >= 2:        
            
            length_max = -1
            current_max = None
            if self.hypso_line.ma_properties[_IS_CLOSED]:
                for crossing_pair in self.crossing_pairs:
                    if crossing_pair.type == _OUTSIDE:
                        if crossing_pair.p0_p1_line.length > length_max:
                            length_max = crossing_pair.p0_p1_line.length
                            current_max = crossing_pair
                
                current = current_max
                for i_dummy in xrange(self._get_nbr_crossing_inside()/2):
                    (dummy, current_plus_1) = self._get_sibbling(current)
                    (dummy, current_plus_2) = self._get_sibbling(current_plus_1)
                    (dummy, current_plus_3) = self._get_sibbling(current_plus_2)
                    if self._build_conflict(current_plus_1, current_plus_2, current_plus_3, sys.maxint):
                        (dummy, current) = self._get_sibbling(current_plus_3)
                    else:
                        break
            else:
                current_max_id = self.crossing_pairs[-1].id
                i = 0
                while  i+3 <= len(self.crossing_pairs):
                    before = self.crossing_pairs[i]
                    current = self.crossing_pairs[i+1]
                    after = self.crossing_pairs[i+2]
                    if self._build_conflict(before, current, after, sys.maxint):
                        i += 4
                    else:
                        break   
                        
        return
    
    def _build_conflict(self, before, current, after, stop_id):
        
        conflict = Holder()
        conflict.type = _COMPLEX
#        conflict.is_editable = True
        conflict.before = before
        conflict.current  = current
        conflict.after  = after
        if conflict.before.type != _INSIDE and conflict.current.type != _OUTSIDE and conflict.before.type != _INSIDE:
            0/0
        if conflict.before.id == stop_id or conflict.after.id == stop_id:
            build = False
        else:
            build = True
            self.conflicts.append(conflict)
            
        return build

                
#     def _list_sub_crossing(self, direction, min_z_crossing_pair, max_z_crossing_pair ):
#         
#         start_id = max_z_crossing_pair.id
#         stop_id = min_z_crossing_pair.id
#         lst_crossing_pair = []
#         
#         if start_id != stop_id:
#             current_id = start_id
#             current = max_z_crossing_pair
#             lst_crossing_pair.append(current)
#             while current_id != None and current_id != stop_id:
#                 (before, after) = self._get_sibbling(current)
#                 if direction == _BEFORE:
#                     current = before
#                 else:
#                     current = after
#                 if current != None:
#                     current_id =current.id
#                     lst_crossing_pair.append(current)
#                 else:
#                     current_id = None
#         
#             lst_crossing_pair.reverse()
#             
#         return lst_crossing_pair
    
    def _identify_conflicts(self, conflict_type):
    
        
        if self._extract_crossing():
            self.conflicts = []
            if self._is_stream_editable():
                if conflict_type == _SIMPLE:
                    for crossing in self.crossing_pairs:
                        if crossing.type == _INSIDE:
                            conflict = Holder()
                            conflict.current = crossing
                            conflict.type = _SIMPLE
        #                    conflict.is_editable = True
                            self.conflicts.append(conflict)
                            
                if conflict_type == _COMPLEX:
                    if self.hydro_feature.ma_properties[CODE] in (STREAM, RIVER):
                        if self._get_nbr_crossing_inside() >= 3:
                            self._identifiy_conflicts_stream_river()
                        else:
                            self._identifiy_conflicts_lake()
                    else:
                        self._identifiy_conflicts_lake()

        return
    
    def next_to_edit(self):
                
        self._identify_conflicts(_SIMPLE)
        for next_conflict in self.conflicts:
            if next_conflict.current.type == _INSIDE:
#                    if next_conflict.is_editable:
                self._set_conflict_attribute(next_conflict)
                if not self.is_all_edited():
                    yield next_conflict
            else:
                print "Internal corruption. Complex conflict. Check the result"

        self._identify_conflicts(_COMPLEX)
        for next_conflict in self.conflicts:
            if next_conflict.before.type == _INSIDE and \
               next_conflict.current.type == _OUTSIDE and \
               next_conflict.after.type == _INSIDE:
                self._set_conflict_attribute(next_conflict)
                if not self.is_all_edited():
                    yield next_conflict
            else:
                print "Internal corruption. Complex conflict. Check the result"
        
        return
                    
    def set_conflict_status(self, hypso_line, conflict, solved):
        
        if solved:
            self.hypso_line = hypso_line
            if conflict.type == _SIMPLE:
                conflict.current.type = _TOUCH
            elif conflict.type == _COMPLEX:
                conflict.before.type = _TOUCH
                conflict.after.type = _TOUCH
             


class AlgoEditContourHydro(Algorithm):
    """This is the main class for the talweg algorithm 
    
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
            test_crossing_line: Flag to enable(TRUE)/disable(FLASE) checking the line crossing constraint
            test_simple_line: Flag to enable(TRUE)/disable(FLASE) checking the simple line constraint
            test_sidedness:  to enable(TRUE)/disable(FLASE) for checking the sidedness constraint
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

    def _process_closed_line (self, lst_coords):
        """
        Special simplification for closed line as we want to have at least 4 vertices to have a closed line with an area > 0

        Parameters:
            line: The line to process
            
        Return value:
            Set containing the index of the vertices to keep
            
        """
        
        first = 0
        last = len(lst_coords)-1
                
        # Calculate the farthest index from the first/last
        mid_index = None
        if ((first+1) <= last):
            mid_index, dummy = self._find_farthest_point(lst_coords, first, last)
            
        # Calculate the farthest point from the lower half of the closed area (if a lower half exist)
        lower_index = None
        lower_dist = None
        if ( mid_index is not None and ( (first+1) <= mid_index )):
            lower_index, lower_dist = self._find_farthest_point(lst_coords, first, mid_index)
                
        # Calculate the farthest point from the upper half of the closed area (if a upper half exist)
        upper_index = None
        upper_dist = None
        if ( mid_index is not None and ( (mid_index+1) <= last )):
            upper_index, upper_dist = self._find_farthest_point(lst_coords, mid_index, last)
                
        # From here determine the points to keep on the line
        index = set()                # Set of points to keep
        
        # Always keep the first and the last point
        index.update([first, last])  

        if (mid_index is not None):
            # Add the mid point
            index.update([mid_index])
            
            fourth_vertice = None
            if (lower_index is not None and upper_index is not None):
                if (lower_dist > upper_dist):
                    fourth_vertice = lower_index
                else:
                    fourth_vertice = upper_index
            else:
                if (lower_index is not None): fourth_vertice = lower_index
                if (upper_index is not None): fourth_vertice = upper_index
                
            if fourth_vertice is not None: index.update([fourth_vertice]) 
                        
        return (list(index))
    



    def _clean_up(self, hypso_lines):
        """Delete unnecessary attributes added on the features during the process
        
        Spatial container containes cyclic reference that must be delete manually otherwise the
        garbage collector is unable to free the memory
        
        Parameters: None
        
        Return value: None
        
        """
        
        for line in hypso_lines:
            if _IS_ALL_EDITED in line.ma_properties: del line.ma_properties[_IS_ALL_EDITED]
            if _IS_CLOSED in line.ma_properties: del line.ma_properties[_IS_CLOSED]
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
   
    def _get_orientation (self, lst_coords):
        """
        Test if a closed list of coordinates is clockwise or anti clockwise
        
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
    
    def _orient_line (self, old_hypso_line, new_hypso_line):
        """
        Orient the line 
        
        hypso_line: LineString contour line to check for orientation correctness
        """
        
        if old_hypso_line.ma_properties[_IS_CLOSED]:
            # For a closed loop we check if the sign of the result of the orientation function
            old_orientation = old_hypso_line.ma_properties[_ORIENTATION]
            new_hypso_coords = list(new_hypso_line.coords)
            new_orientation = self._get_orientation(new_hypso_coords)
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
            new_hypso_coords.reverse()
            line = LineString(new_hypso_coords)
        
        return line
    
    def _check_constraints (self, s_container, hypso_line, hydro_feature, new_sub_hypso_line, new_hypso_line, buffer):
        """Method to check the constraint  of the new hypsographic line and edit the line
        """
        
        # check if the new hypso line is simple
        line_simple = new_hypso_line.is_simple
            
        features = s_container.get_features(bounds=new_sub_hypso_line.bounds, remove_keys=[hypso_line._sci_id])
        
        # Sort the feature by type
        hypso_features = [feature for feature in features if feature.ma_properties[CODE] == HYPSO]
        lake_features  = [feature for feature in features if feature.ma_properties[CODE] == LAKE ]
        river_features = [feature for feature in features if feature.ma_properties[CODE] == RIVER ]
        stream_features = [feature for feature in features if feature.ma_properties[CODE] == STREAM ]
        
        
        hypso_distance = [hypso for hypso in hypso_features if hypso.distance(new_sub_hypso_line) <= buffer.get_current() *.05]
        lake_intersect  = filter(new_sub_hypso_line.intersects, lake_features)
        river_intersect = filter(new_sub_hypso_line.intersects, river_features)
        stream_intersect = filter(new_sub_hypso_line.intersects, stream_features)
        
        if not line_simple or len(hypso_distance) != 0:
            # The hypso line must always be simple and never be to closed from another hypso line
            constraints_ok = False
        else:
            # There are different constraint for the lake, the river and the stream
            if hydro_feature.ma_properties[CODE] == LAKE:
                # For a lake the hypso line must not cross any lake but can cross a river or a stream
                if len(lake_intersect) == 0:
                    constraints_ok = True
                else:
                    constraints_ok = False
                reset_edition = river_intersect + stream_intersect
            if hydro_feature.ma_properties[CODE] == RIVER:
                # For a river the hypso line must not cross any lake must not cross any river but can cross a stream
                if len(lake_intersect) == 0 and len(river_intersect) == 0:
                    constraints_ok = True
                    reset_edition = stream_intersect
                else:
                    constraints_ok = False
                reset_edition = stream_intersect
            if hydro_feature.ma_properties[CODE] == STREAM:
                # Check if the new hypso line passes through the extremity of the straem
                # In order to check if the new hypso line passes by the extremity of the line we create two lines
                # at the extremity of the stream we a length of just a little over the work buffer.
                hydro_coords = list(hydro_feature.coords)
                extended_first_coord = HydroUtil.extend_line(hydro_coords[1], hydro_coords[0], buffer.get_current()*1.05)
                extended_last_coord = HydroUtil.extend_line(hydro_coords[-2], hydro_coords[-1], buffer.get_current()*1.05)
                multi_line_extremity = MultiLineString( ( ((extended_first_coord, hydro_coords[0])), ((extended_last_coord, hydro_coords[-1])) ) )
#                a.append(MA_LineString(((extended_first_coord, hydro_coords[0]))))
#                a.append(MA_LineString(((extended_last_coord, hydro_coords[-1]))))
#                if hydro_feature.length > (buffer.get_current() +1.0)*1.1:

#                    point = hydro_feature.interpolate(1.0)
#                    coords_1 = point.coords[0]
#                    point = hydro_feature.interpolate(hydro_feature.length-1.0)
#                    coords_n = point.coords[0]#
#                    coords_first = hydro_feature.coords[0]
#                    coords_last = hydro_feature.coords[-1]
#                    rescale_start = GenUtil.rescale_vector(coords_1, coords_first, (buffer.get_current()+1.)*1.05)
#                    rescale_last = GenUtil.rescale_vector(coords_n, coords_last, (buffer.get_current()+1.)*1.05)
#                    print coords_last, rescale_last, buffer
#                    multi_line_extremity = MultiLineString( ( ((rescale_start, coords_first)), ((coords_last, rescale_last)) ) )
                    # Check if the new hypso line is disjoint of the 2 lines created before
                disjoint = new_hypso_line.disjoint(multi_line_extremity)
                # For a stream the hypso line must not cross any lake, river and the hyso_line must not pass near the extremity of the stream
                if len(lake_intersect) == 0 and len(river_intersect) == 0 and len(stream_intersect) == 0 and disjoint:
                    constraints_ok = True 
                else:
                    constraints_ok = False
#                else:
#                    # It's to small we assume it's OK
#                    constraints_ok = True
                reset_edition = []

            if constraints_ok:
                # Meake sure to reset the flag of in order to edit the feature 
                for hydro_feature in reset_edition:
                    hydro_feature.ma_properties[_IS_ALL_EDITED]= False

        return constraints_ok
             
    
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

    def _resolve_loops(self, s_container, hypso_line, hydro_feature, new_sub_hypso_line, loop_area, buffer):
        """
        """
        
        lst_sub_loop_line = []
        constraint_loops = True
    
        if not loop_area.is_empty:
        # Process the loop line
            current_delta = buffer.get_current_delta()
            buf_loop_areas = loop_area.buffer(current_delta*(-1.05), self.params.resolution)
            lst_sub_loop_area = HydroUtil.multi_to_list(buf_loop_areas)
            
            for sub_loop_area in lst_sub_loop_area:
                sub_loop_area = sub_loop_area.buffer(current_delta*1.05, self.params.resolution)
                if ( isinstance(sub_loop_area, Polygon) and \
                     not sub_loop_area.is_empty ):
                    if (sub_loop_area.area >= self.params.min_area):
                        sub_loop_line= LineString(sub_loop_area.exterior.coords)
                        constraints_loops = self._check_constraints (s_container, hypso_line, hydro_feature, sub_loop_line, 
                                                                     sub_loop_line, buffer)
                        if constraints_loops:
                            lst_sub_loop_line.append(sub_loop_line)
                        else:
                            break
        else:
            # There is no loop line so we can assume everything is OK
            pass
           
        return (lst_sub_loop_line, constraint_loops)
        
    def _resolve_conflicts (self, s_container, s_container_z_point, hydro_feature, hypso_line, buffer):
        
           
        conflicts = Conflicts(hydro_feature, hypso_line, buffer, s_container_z_point, self.params.max_displacement, 
                              self.params.resolution )
        
        for conflict in conflicts.next_to_edit():

            conflict_resolved = False
            if conflict.distance < self.params.max_displacement:  
             
                if not HydroUtil.are_lines_to_close(conflict.to_paste_hypso_line, conflict.to_keep_hypso_line, buffer.get_current_delta()): 
                            
                    conflict.to_keep_hypso_line.append(conflict.to_paste_hypso_line)
                    new_hypso_line = HydroUtil.merge_lines(conflict.to_keep_hypso_line)
                                    
                    # Chack that the line has the same orientation has the original line
                    new_hypso_line = self._orient_line (hypso_line, new_hypso_line)
            
                    (lst_loop_line,constraint_loops)  = self._resolve_loops(s_container, hypso_line, hydro_feature, 
                                                                            conflict.to_paste_hypso_line, conflict.loop_area, buffer)
                    
                    constraint_inside =  self._check_constraints (s_container, hypso_line, hydro_feature, conflict.to_paste_hypso_line, 
                                                                  new_hypso_line, buffer)
                    
                    if constraint_loops and constraint_inside:
                        hypso_line.update_coords(new_hypso_line.coords, s_container)
                        if conflict.type == _SIMPLE:
                            self.stats.add_stats(_INSIDE)
                        else:
                            self.stats.add_stats(_OUTSIDE)
                        conflict_resolved = True
                        for loop_line in lst_loop_line:
                            new_hypso_line = hypso_line.cloner()
                            s_container.add_feature(new_hypso_line)
                            new_hypso_line.update_coords(loop_line.coords, s_container)
                    
            conflicts.set_conflict_status(hypso_line, conflict, conflict_resolved)

        return conflicts.is_all_edited()


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
            hypso_line.ma_properties[_IS_CLOSED] = GenUtil.is_line_closed(hypso_line)

            if hypso_line.ma_properties[_IS_CLOSED]:
                hypso_line.ma_properties[_ORIENTATION] = self._get_orientation(hypso_line.coords_dual)
            else:
                hypso_line.ma_properties[_FIRST_COORDS] = hypso_line.coords_dual[0]
                hypso_line.ma_properties[_LAST_COORDS] = hypso_line.coords_dual[-1]         
    
    def _pre_process_hypso_lines (self, buffer, hydro_feature, s_container):
        """
        Identifies all the hypso that could be edited for this hydro feature
        
        The following operation are done 
            - buffer the hydro feature
            - delete all the hypso features that are located completly in the donut are of the buffer
            - not retain (not delete!!!) open hypso lines where the first or last point are located in the donut area
        """
        
        
        # Buffer the hydro feature
        buf_hydro_pol = hydro_feature.buffer(buffer, resolution=self.params.resolution)

        buf_outer_ring_coords = list(buf_hydro_pol.exterior.coords)
        
        # Delete the contours that lie completely in the buffer zone
        if isinstance(hydro_feature, MA_Polygon):
            # The hydro feature is a polugon so make a hole
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

        # Find all feature crossing the outer ring of the buffered hydro feature
        buf_hydro_line = LineString(buf_outer_ring_coords)
        hypso_lines = s_container.get_features(bounds=buf_hydro_line.bounds )
        hypso_lines = [hypso_line for hypso_line in hypso_lines if hypso_line.ma_properties[CODE] == HYPSO]
        hypso_lines =  filter(buf_hydro_line.intersects, hypso_lines)
        
        if hydro_feature.ma_properties[CODE] == STREAM:
            # For the stream do not retain hypso line that crosses only one time the hydro feature. They do not need edition
            tmp_hypso_lines = []
            for hypso_line in hypso_lines:
                intersections = hydro_feature.intersection(hypso_line)
                intersections = GenUtil.make_iterable(intersections)
                if len(intersections) >= 2:
                    tmp_hypso_lines.append(hypso_line)
                hypso_lines = tmp_hypso_lines 

        
        # Don't retain open hypso line where the first or last vertice lie in the buffer zone to difficul to manage
        pre_process_hypso_lines = []
        for hypso_line in hypso_lines:
            if hypso_line.is_simple:
                if not hypso_line.ma_properties[_IS_CLOSED]:
                    hypso_coords = list(hypso_line.coords)
                    p0 = Point(hypso_coords[0])
                    p1 = Point(hypso_coords[-1])
                    if p0.disjoint(donut_hydro_pol) and p1.disjoint(donut_hydro_pol):
                        pre_process_hypso_lines.append(hypso_line)
                    else:
                        self.stats.add_stats(_NOT_PROCESSED)
                else:
                    pre_process_hypso_lines.append(hypso_line)
            else:
                self.stats.add_stats(_NOT_SIMPLE)                
                 
        return pre_process_hypso_lines
     
    def process(self):
        """Main routine for the Talweg coherence algorithm
        
        This algorithm will edit the contour line in order to fit with the Talweg line.
        It will prevent line crossing and sidedness errors.
        
        Parameters: None
            
        """
                
        s_container = self.load_features()
        max_iteration = 5
        buffer = Buffer(self.params.max_buffer, max_iteration)
        
        self._add_attributes (s_container)
        
        # Extract all the hydro features
        features = s_container.get_features()
        lake_features  = [feature for feature in features if feature.ma_properties[CODE] == LAKE ]
        river_features = [feature for feature in features if feature.ma_properties[CODE] == RIVER ]
        stream_features = [feature for feature in features if feature.ma_properties[CODE] == STREAM]
        # Order by Lake, River and stream
        hydro_features = lake_features + river_features + stream_features
        
        z_point_features = [feature for feature in features if feature.ma_properties[CODE] == Z_POINT]
        s_container.del_features(z_point_features)
        s_container_z_point = SpatialContainer()
        s_container_z_point.add_features(z_point_features)

#        try:
#           import pyfme
#           logger = pyfme.FMELogfile()
#        except:
#            pass
        
        for dummy in xrange(max_iteration):
#            print "iteration"
            buffer.next_iteration()
            hydro_features = [hydro_feature for hydro_feature in hydro_features if not hydro_feature.ma_properties[_IS_ALL_EDITED] ]
            for i, hydro_feature in enumerate(hydro_features):
#                try:
#                   msg =  "Iteratiom, feature: " + str(dummy) + " " + str(i)
#                   logger.log(msg)
#                except:
#                    pass
                hypso_lines = self._pre_process_hypso_lines (buffer.get_current(), hydro_feature, s_container)
                hypso_lines.sort(key=operator.attrgetter("length"))
                all_edited = True

                for j, hypso_line in enumerate(hypso_lines):
#                    try:
#                        msg =  "hypso: " + str(j)
#                        logger.log(msg)
#                    except:
#                        pass
                    edited = self._resolve_conflicts (s_container, s_container_z_point, hydro_feature, hypso_line, buffer)
                    all_edited = edited and all_edited

                if all_edited:
                    hydro_feature.ma_properties[_IS_ALL_EDITED] = True
                    
        # Clean up attributes
        features = s_container.get_features()
        lake_features    = [feature for feature in features if feature.ma_properties[CODE] == LAKE ]
        river_features   = [feature for feature in features if feature.ma_properties[CODE] == RIVER ]
        stream_features  = [feature for feature in features if feature.ma_properties[CODE] == STREAM]
        hypso_features   = [feature for feature in features if feature.ma_properties[CODE] == HYPSO ]
        s_container.del_features(lake_features+river_features+stream_features)
        self._clean_up(hypso_features)
        
        # Extract all the features for the output
        self.extract_features_out(s_container.get_features())
        
        self._free_memory(s_container)
        
        GenUtil.print_debug (self.params, "End of %s algorithm" %(_CLEAN_HYDRO))
        