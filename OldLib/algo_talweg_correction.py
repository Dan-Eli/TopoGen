#! /usr/bin/env python
# -*- coding: UTF-8 -*-

#####################################################################################################################################

"""
    This algorithm change the position of the talweg (contour reentrance) in order to graphically fit the talweg 
    with the position of the linear river/stream.
    The algorithm will only change the position of the hypsographic linear features all other features are not modified. 
    Hypsograp point features are only used to keep topological coherence (sidedness)  when a talweg is moved.
    Hydrographic linear features are used for talweg detection and are used for topological coherence (line crossing and sidedness) 
    when a talweg is moved.
    Hydrographic point and area features are only used to keep topological coherence (line crossing and sidedness) 
    when a talweg is moved.
    
    Usage:
        import algo_talweg_correction

    Limits and constraints
        Always works better when the line to process meet the OGC simple line.
        The rivers must be oriented from upstream (first vertice) to downstream (last vertice).
        Works better when the contours line are simplified and smoothed according to the 
        target scale 
          
"""
__revision__ = "--REVISION-- : $Id: algo_talweg_correction.py 472 2011-09-29 13:39:22Z dpilon $"

#####################################################################################################################################

import math, sys

from shapely.geometry import Point, LineString, Polygon, MultiPoint
from algo_bends import AlgoBends 
from shapely.prepared import prep

from lib_genmetal import MA_Point, MA_LineString, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm
                         
########################################################

# Public constant
HYDRO = "HYDRO"   # Code for hydrographic features
HYPSO = "HYPSO"   # Code for hypsographic features

#Internal key word constants
_TALWEG = "Talweg"

# Constant for the statistics
_STAT_EDITED_SIMPLE = "Talweg simple"   # Number of simple talweg edited
_STAT_EDITED_INVERTED = "Talweg inverted" # Number of inverted talweg edited
_STAT_EDITED_ERROR = "Talweg error"    # Number of error trying to edit a talweg (simple or inverted)

# Constant for the state of a talweg (a bend in the contour line)
_TALWEG_TYPE_EMPTY = "Talweg empty"        # A talweg is empty when no stream is found inside a bend
_TALWEG_TYPE_SIMPLE = "Talweg simple"      # A talweg is simple when the stream enters the contours and leave through the bend base line
_TALWEG_TYPE_INVERTED = "Talweg inverted"  # A talweg is inverted when the stream enters the bend base line and leave through the contours
_TALWEG_TYPE_OTHER = "Talweg other"        # A talweg other is a talweg that we cannot set as EMPTY, SIMPLE or INVERTED

_EDITED = "Edited"                         # The coordiantes of the bends are edited
_UNEDITED = "Unedited"                     # The coordiantes of the bend are unedited

_MAX_DISPLACEMENT = 'max_displacement'
_CODE = 'code'
   

class Statistics(GenStatistics):
    """Class that contains the statistics for the talweg correction algorithm
    
    Attributes
        stat_names: Name of the statistics for the TalwegStatistics class. These name are
                    used by the Statistics class

    """
    
    def __init__(self):
        """Initialize the attributes of an object of the class"""
        
        GenStatistics.__init__(self)
        self.stats_names = ((_STAT_EDITED_SIMPLE, _STAT_EDITED_INVERTED, _STAT_EDITED_ERROR, 
                             GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS, ))
        
        self.add_iteration()
        
    def get_stats (self, type=None):
        """Extract the current statistics and build  a list of string that forms the statistical message"
        
        Parameters:
            type: Unused... only there for compatibility purpose
        
        """
        
        str_out = []
        str_out.append( "%s coherence algorithm Statistics" %(_TALWEG) )
        str_out.append( "-------------------------------------" )
        str_out.append("Simple talweg edited: " + str(self.get_stats_name_count_total( _STAT_EDITED_SIMPLE)))
        str_out.append("Inverted talweg edited: " + str(self.get_stats_name_count_total( _STAT_EDITED_INVERTED)))
        str_out.append("Edition error (not edited): " + str(self.get_stats_name_count_total( _STAT_EDITED_ERROR)))
        
        return str_out
    
class AlgoTalwegCorrection(Algorithm):
    """This is the main class for the talweg algorithm 
    
    Attributes:
        - params: Parameters of the algorithm
        
    """

    def __init__(self, test_crossing_line=True,
                       test_simple_line=True, 
                       test_sidedness=True,
                       exageration=10., 
                       debug=False):
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
        self.params.test_crossing_line=test_crossing_line
        self.params.test_simple_line=test_simple_line 
        self.params.test_sidedness=test_sidedness
        self.params.debug=debug
        
        self.stats = Statistics()      
    
    def _detect_bend_location(self): 
        """Calculates the position of each individual bends in the line
        
        The position of the bends are calculated according to the definition of the bencds 
        in the orginal paper of Wang 1998.  
        """
        
        algo_bends = AlgoBends()
        for hypso_line in self.hypso_lines:
            algo_bends.features.append( hypso_line )
        
        algo_bends.process()
                    
    def _add_line_attributes(self): 
        """Add some attributes to the lines and bends
        
        Parameters: None
        
        Return value: None
          
        """
        
        # Add some attributes for hypso lines
        for line in self.hypso_lines:    
            line.vertice_edited = set()  # Tracks the vertice edited
            for bend in line.bends:
                # Add some attributes to the bend
                bend.type = _TALWEG_TYPE_EMPTY
                bend.status = _UNEDITED
                
        # Add some attributes for hydro line
        for line in self.hydro_lines:
            line.buffer_line = []
            line.buffer_weight = []
            for value in (.01, 1., 2., 5.):
                # Buffer are calculated and stored to avoid recalculation each time needed
                buffer = line.buffer(value, resolution=2)
                buffer_coords = list(buffer.exterior.coords)
                line.buffer_line.append(LineString(buffer_coords))
                line.buffer_weight.append(value)
            
    
    def _extract_centricity(self, talweg_hydro_point, talweg_line, hydro_line):
        """Calculate a measure of centricity of the hydro line whitin the talweg
        
        To determine the centricity we intersects the buffered hydro line with the
        talweg line. An acceptable result is when the intersection is 2 and only 2 points
        located just before and just after the intersection of the hydro line and the
        talweg line
        
        Parameters:
            - talweg_hydro_point: Point where the hydro crosses the talweg
            - talweg_line: LineString object defining the talweg
            - hydro_line: LineString object defining the hydro line
            
        Return value:
            - Float: Measure of the quality of the centricity
                0: Is the best measure
                max(hydro_line.buffer_weight): Is the worst value measure
        """
        
        centricity = 0.
        # Poistion of the hydro talweg crossing 
        dist_talweg_hydro = talweg_line.project(talweg_hydro_point)
        
        # Check for each buffer if it meets the criteria.  Starting with the biggest buffer
        for i in reversed(range(1, len(hydro_line.buffer_line))):
            buffer_line = hydro_line.buffer_line[i]
            buffer_weight = hydro_line.buffer_weight[i]
            talweg_crossing_points = buffer_line.intersection(talweg_line)
            # Check if the results 2 and only 2 points
            if ( isinstance(talweg_crossing_points, MultiPoint) and 
                 len(talweg_crossing_points) == 2):
                dist_0_talweg_hydro = talweg_line.project(talweg_crossing_points[0])
                dist_1_talweg_hydro = talweg_line.project(talweg_crossing_points[1])
                if (dist_0_talweg_hydro > dist_1_talweg_hydro):
                    dist_0_talweg_hydro,dist_1_talweg_hydro = dist_1_talweg_hydro, dist_0_talweg_hydro
                # Check if the hydro talweg crossing is located in the middle of the 2 points
                if ( dist_0_talweg_hydro < dist_talweg_hydro and 
                     dist_1_talweg_hydro > dist_talweg_hydro ):
                    # We have a good result and exit
                    centricity = buffer_weight
                    break 
                    
        return max(hydro_line.buffer_weight) - centricity
    
    def _extract_talweg_base_height_ratio(self, talweg_coords):
        """Extract the ratio between the height of the talweg and the base of the talweg
        
        Parameters: 
           - talweg_coords: List of tuple of coordinates (x,y)
        
        Return:
           - Float >0: Ratio between the height and the base of the talweg
           
        """
        bend_i = 0
        bend_j = len(talweg_coords)-1
        
        #Loop over each vertex of the bend in order to find its distance from the bend base
        p1 = talweg_coords[0]
        p2 = talweg_coords[-1]
        lst_height = []
        for i in range(bend_i+1, bend_j):
            dist = GenUtil.distance_line_point(p1, p2, talweg_coords[i])
            lst_height.append(dist)
        
        # Extract the index of the distance of the bend peak
        max_height = max(lst_height)
        base = GenUtil.distance(talweg_coords[0],talweg_coords[-1])

        # Avoid division by zero
        if max_height >= GenUtil.ZERO:
            ratio =  base/max_height
        else:
            ratio = 1.0e+10

        return ratio
        
    
    def _extract_bend_peak(self, bend_i, bend_j, line):
        """Locate the position of the bend peak along the talweg line
        
        The bend peak position is the place along the talweg which is the further from the bend base.
        But in order to create a more natural location of the bend peak.  The position of the bend peak
        can be modify to fit between 2 vertice. 
                  
        Parameters:
            bend_i: Index of the start of the talweg in the line
            bend_j: Index of the end of the talweg in the line
            line: MA_LineString containing the bend to process
            
        Return value: Tuple (k,l) containing the index of the bend peak. If k,l are the same the bend peak
                      is located exactly at this index along the line; if k,l are different the bend peak is located 
                      between index k and l along the line
        
        """
        
        lst_infos=[]
        p0 = line.coords_dual[0]
        p1 = line.coords_dual[-1]
        
        # Manage flat angle (180.) caused by the vertice densification
        lst_angles = GenUtil.compute_angles(line.coords_dual)
        for i in xrange(1, len(line.coords_dual)-1):
            angle = lst_angles[i]
            if (math.fabs(180.-angle) > GenUtil.ZERO):
                # Distance from the base to the vertex along the talweg
                dist = GenUtil.distance_line_point(p0, p1, line.coords_dual[i])
                lst_infos.append((dist,i,angle))
            else:
                # Do not include flat angle in thelist
                pass
                
        # Extract the three most far off vertex with the index and the angle
        lst_infos.sort(reverse=True)        
        d0,i0,a0 = lst_infos[0][0], lst_infos[0][1], lst_infos[0][2] 
        d1,i1,a1 = None, None, None
        d2,i2,a2 = None, None, None
        if len(lst_infos) >= 2:
            d1,i1,a1 = lst_infos[1][0], lst_infos[1][1], lst_infos[1][2]
        if len(lst_infos) >= 3:
            d2,i2,a2 = lst_infos[2][0], lst_infos[2][1], lst_infos[2][2]
            if (math.fabs(i2-i1) == 1 ):
                # Vertex i1 and i2 or not located before and after the bend peak
                # They are side by side
                d1,i1,a1 = None, None, None
                d2,i2,a2 = None, None, None
                
        if (a0 < 65.):
            # When the angle is small the peak vertex is always the farther
            peak_kl = (i0, i0)
        elif ( d1 is not None and 
               d2 is not None and
               math.fabs(d1-d2) / max(d1,d2) < 0.1 ):
            # If the distance between the 2 farther vertice is whithin 10% the peak vertex 
            # remains the peak vertex
            peak_kl = (i0, i0)
        elif ( d1 is not None and
               math.fabs(d0-d1) / d0 < 0.2):
            # The 2 farther vertice are within 20% distance the peak is between the 2 vertice
            peak_kl = (i0, i1)
        else:
            # Otherwise take the peak vertex
            peak_kl = (i0, i0)
            
        if peak_kl[0] > peak_kl[1]:
            peak_kl = (peak_kl[1], peak_kl[0])
        
        return peak_kl
              
    def _talweg_screening(self, max_displacement, hydro_talweg_point, talweg_line):
        """Check that the talweg conforms to certain criteria
        
        This method is only used to speed up the process by removing early in the process
        talweg that do not meet some minimal talweg critaria
        
        Parameters:
            - max_displacement: Maximum tolerance this talweg line can be moved
            - hydro_talweg_crossing_point: Point where the hydro crosses the talweg
            - talweg_line: MA_LineString defining the talweg
            
        Return value:
            - Boolean: Boolean indicating if the talweg meets (True) or do not meet (False)
                       the minimal criteria
                       
        """
        
        minimal_criteria = False
        
        # Validate that the base of the talweg is below a threshold otherwise the talweg is to big and the correction
        # is often unpredictable
        if (GenUtil.distance(talweg_line.coords_dual[0],talweg_line.coords_dual[-1]) < max_displacement * 4.):  
        
            # The base height ration of the talweg must below a threshold
            base_height_ratio = self._extract_talweg_base_height_ratio(talweg_line.coords_dual)
            if (base_height_ratio < 5. and base_height_ratio > .15):
            
                # Check that the hydro line is crossing the talweg line roughly in the middle of the talweg line
                hydro_talweg_ratio =  talweg_line.project(hydro_talweg_point, normalized=True)
                if hydro_talweg_ratio > .05 and hydro_talweg_ratio < .95:                 
                    # Check that the amplitude of the talweg is between an interval
                    amplitude =  self._extract_bend_amplitude(talweg_line.coords_dual)
                    if (amplitude > 60. and amplitude < 270.):
                        # All the minimal criteria are met
                        minimal_criteria = True
                        
        return minimal_criteria

    def _extract_talweg_type(self, talweg_line, hydro_line):
        """This method determines the type of talweg a talweg line and an hydro_line are forming
        
        Three results are possible _TALWEG_SIMPLE, _TALWEG_INVERTED and _TALWEG_UNKNOW. Here are the
        condition to meetr for each results
        _TALWEG_SIMPLE: - The bend base line must be greater than zero. Note: A bend base line of zero happens
                          with a closed line with only one bend.
                        - Inside the polygon formed by the talweg there must be one and only one hydro feature
                        - The bend base line must not cross the talweg line (line crossing itself).
                        - The hydro line and the base line must each one cross once and only once the 
                          bend base line. 
                        - The oriented hydro line must cross the talweg line before the bend base line
          
        _TALWEG_INVERTED: - The bend base line must be greater than zero. A bend base line of zero happens
                            with a closed line with only one bend.
                          - Inside the polygon formed by the talweg there must be one and only one hydro feature
                          - The bend base line must not cross the talweg line (line crossing itself).
                          - The hydro line and the base line must each one cross once and only once the 
                            bend base line. 
                          - The oriented hydro line must cross the talweg line before the bend base line.
        
        _TALWEG_OTHER: All the other cases that do not meet the conditions for _TALWEG_SIMPLE and
                         _TALWEG_INVERTED
        
        Parameters:
            - talweg_line: LineString representing the talweg line
            - hydro_line: oriented MA_LineString representing the hydrographic line
            
        Return value:
            - The type of talweg: _TALWEG_SIMPLE, _TALWEG_INVERTED, _TALWEG_COMPLEX 
        
        """
        
        talweg_type = _TALWEG_TYPE_OTHER
        lst_coords = list(talweg_line.coords)
                    
        # Extract the hydro line in the talweg polygon
        talweg_polygon = Polygon(talweg_line)
        features = self.s_container.get_features(bounds=talweg_polygon.bounds, filter="feature.ma_properties['code']=='%s'" %HYDRO)
        features = filter(talweg_polygon.intersects, features)
                
        # If there is more than one feature inside the talweg polygon it is to complex... 
        if (len(features) == 1 ):
           # If the hydro line is not the same as the one passed in parameters there is a problem
           if ( features[0] == hydro_line):            
                # Close the talweg line with the base line and check if the resulting line is simple
                base_line = LineString([lst_coords[0], lst_coords[-1]])
                lst_coords = lst_coords + [lst_coords[0]]
                talweg_line_closed = LineString(lst_coords)
                # Check that the bend base line is greater than zero to handle closed lines which are not processed.
                if (base_line.length > GenUtil.ZERO):
                    # Check that the base line do not cross the talweg line
                    if (talweg_line_closed.is_simple):
                        # Check that the hydro line intersects the bend base line in one and only one point
                        # Checks that the hydro line intersects the talweg line in one and only one point
                        base_hydro_point = hydro_line.intersection(base_line)
                        talweg_hydro_point = hydro_line.intersection(talweg_line)
                        if ( isinstance(base_hydro_point, Point) and
                             isinstance(talweg_hydro_point, Point) ):
                             # Check if a very small buffered hydro line will cross the talweg line in only 2 points
                             talweg_hydro_points = hydro_line.buffer_line[0].intersection(talweg_line)
                             if ( isinstance(talweg_hydro_points, MultiPoint) and 
                                  len(talweg_hydro_points) == 2):
                                 # Validate that the distance between the 2 intersection points is below 1m
                                 if GenUtil.distance(talweg_hydro_points[0].coords[0], talweg_hydro_points[1].coords[0]) < 1.0:
                                     # Check if the hydro line cross the talweg line first
                                     if ( hydro_line.project(talweg_hydro_point) < hydro_line.project(base_hydro_point) ):
                                          talweg_type = _TALWEG_TYPE_SIMPLE
                                     else:
                                        talweg_type = _TALWEG_TYPE_INVERTED

        return talweg_type
    
    def _extract_bend_amplitude(self, lst_coords):
        """Calculates the amplitude of a bend by adding the angle of each coordinate of the line
        
        Parameters:
          - lst_coords: List of (x,y) tuple

        Return value
          - Float > 0: Sum of the angles
        
        """
        
        lst_angles = GenUtil.compute_angles(lst_coords)
        angles = 0.
        del lst_angles[0]
        del lst_angles[-1]
        for angle in lst_angles:
            angles += (180.-angle)
            
        return angles

    def _detect_talweg (self):
        """Determine and flag within each bend of the contours the one that are real talwegs
        
        There are 2 types of talweg:
           - Simple: When a stream enter through the contour line delimiting the talweg and
                     leave the bend through the bend base line.
           - Inverted: When a stream enter through bend base line and
                     and leave through the bend base line.

        
        This is true because the streams are oriented.
        
        Parameters: None
            
        Return value:  None
            
        """
                     
        # Loop over all hypso lines
        for hypso_line in self.hypso_lines:
            
            # Loop over all bends            
            for bend in hypso_line.bends:
                i = bend.i
                j = bend.j
                
                # Filter first by bounds first to save some time
                bend_polygon = Polygon(hypso_line.coords_dual[i:j+1])
                features = self.s_container.get_features(bounds=bend_polygon.bounds, 
                                                         filter="feature.ma_properties['code']=='%s'" %HYDRO)
                if features:
                    features = filter(bend_polygon.intersects, features)
                
                # If there is more than one feature inside the bend polygon we don't handle that case to complex 
                if (len(features) == 1):
                    hydro_line = features[0]
                    talweg_line = MA_LineString(hypso_line.coords_dual[i:j+1])
                    talweg_type = self._extract_talweg_type(talweg_line, hydro_line)
                    if ( talweg_type == _TALWEG_TYPE_SIMPLE or 
                         talweg_type == _TALWEG_TYPE_INVERTED):
                        
                        bend.hydro_line = hydro_line # Keep a reference to which hydro line
                        bend.hydro_talweg_point = talweg_line.intersection(hydro_line ) # Keep a reference to the intersection 
                        bend.type = talweg_type
                    else:
                            # A bend with a big amplitude is difficult to modify
                        bend.type == _TALWEG_TYPE_OTHER
                else:
                    if (len(features) >= 2):
                        # Configuration too complex to be corrected; 2 hydro line in the same talweg
                        bend.type == _TALWEG_TYPE_OTHER
                                                                           
    def _extract_talweg_measure(self, ori_talweg_coords, new_talweg_coords, hydro_line):
        """Calculate different measures of the talweg line in order to quantify the quality the talweg
        
        The following measures are calculated
           - The relative position of the hydro talweg crossing point on the talweg line. The value is between 0. and 0.5
             A value near 0 means the hydro line  is crossing right in the middle of the talweg (best case)
           - Form coherence which is the difference between the circularity ration of the talweg before and after the talweg 
             displacement. The value is between 0 and 1. A samll value near 0 means the shape of the talweg 
             as not change a lot (best case) 
           - The amplitude of the talweg (sum of the angles) minus 180. The value is > 0.  A value near 0 means the talweg is
             forming almost a half circle (best case)
           - The centricity which is a measure if the line is too close from the talweg line. the value is > 0. A value of 0
             means that the hydro is not to close from the talweg line (best case)
             
        Note: All the measures are at their best when they are near 0
        
        Parameters:
            - ori_talweg_coords: List of x,y tuple forming the original talweg line
            - new_talweg_coords: List of x,y tuple forming the new talweg line
            - hydro_line: hydrographic MA_LineString to process
            
        Return value
            - A list of 4 the four measures
               - Relative position
               - Form coherence 
               - Talweg amplitude
               - Centricity
        """
        
        new_talweg_line = LineString(new_talweg_coords)
                
        # Check that the talweg is positionned in the middle
        talweg_hydro_point = hydro_line.intersection(new_talweg_line)
        if ( isinstance(talweg_hydro_point, Point)):
            hydro_talweg_crossing =  new_talweg_line.project(talweg_hydro_point, normalized=True)
            if ( hydro_talweg_crossing > 0.375 and  hydro_talweg_crossing < 0.625):
                # Position of the hydro crossing on the talweg
                mid_talweg_line = math.fabs(0.5-hydro_talweg_crossing)
                
                # Extract the circularity ratio of the talweg before and after
                new_talweg_pol = Polygon(new_talweg_coords)
                ori_talweg_pol = Polygon(ori_talweg_coords)
                new_circul_ratio = (4 * math.pi * new_talweg_pol.area) / new_talweg_pol.length**2.
                ori_circul_ratio = (4 * math.pi * ori_talweg_pol.area) / ori_talweg_pol.length**2.
                form_value = math.fabs(new_circul_ratio - ori_circul_ratio)
                        
                # Extract the amplitude of the talweg 
                amplitude = math.fabs(self._extract_bend_amplitude(new_talweg_coords) - 180.)
                
                # Extract the measure of centricity to determine if the hydro feature is in the middle of the talweg
                centricity = self._extract_centricity(talweg_hydro_point, new_talweg_line, hydro_line)
                
                measures = [mid_talweg_line, form_value, amplitude, centricity]
            else:
                # Talweg in not near the center of the line and it's not acceptable 
                measures = None
        else:
            # Hydro line does not intersect talweg line after line warpping 
            measures = None

        return measures
                
    def _normalize_one_measure (self, measure, measure_min, measure_max):
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
    
    def _normalize_talweg_measures (self, talweg_measures):
        """Normalize the talweg measures to values between [0..1]
        
        Parameters:
            bend_base_infos: Alist of tuple containing 5 information
                             - Length of the base
                             - Angle formed at the start of the bend base
                             - Angle formed at the end of the bend base
                             - Index of the start of the bend
                             - Index of the end of the bend
        
        Return value
            A tuple containing 3 values
                - The added normalized value [0..3]
                - The start index of the bend i on the line
                - The end index of the bend j on the line
        
        """
        
        norm_talweg_measures = []
        
        # Extract min max value for each measure
        lst_mid_length = [talweg_measure[1] for talweg_measure in talweg_measures]
        mid_length_min = min(lst_mid_length)
        mid_length_max = max(lst_mid_length)
        
        lst_form_value = [talweg_measure[2] for talweg_measure in talweg_measures]
        form_value_min = min(lst_form_value)
        form_value_max = max(lst_form_value)
        
        lst_amplitude = [talweg_measure[3] for talweg_measure in talweg_measures]
        amplitude_min = min(lst_amplitude)
        amplitude_max = max(lst_amplitude)
        
        lst_centricity = [talweg_measure[4] for talweg_measure in talweg_measures]
        centricity_min  = min(lst_centricity)
        centricity_max  = max(lst_centricity)
        
        for i in range(len(talweg_measures)):
            mid_length   = lst_mid_length[i]
            form_value  = lst_form_value[i] 
            amplitude    = lst_amplitude[i]
            centricity   = lst_centricity[i]
            bend_i       = talweg_measures[i][4]
            bend_j       = talweg_measures[i][5]
            peak         = talweg_measures[i][6]
            vertex_added = talweg_measures[i][7]
            talweg_ori_coord = talweg_measures[i][8]
            talweg_new_coord = talweg_measures[i][9]
            
            # Normalize the different values between [0..1]
            norm_mid_length = self._normalize_one_measure( mid_length, mid_length_min, mid_length_max)
            norm_form_value = self._normalize_one_measure( form_value, form_value_min, form_value_max)
            norm_amplitude = self._normalize_one_measure( amplitude, amplitude_min, amplitude_max)
            norm_centricity = self._normalize_one_measure( centricity, centricity_min, centricity_max)
            
            norm_measure = norm_mid_length + norm_amplitude + norm_form_value + norm_centricity
            
            norm_talweg_measures.append((norm_measure, bend_i, bend_j, peak, vertex_added, talweg_ori_coord, talweg_new_coord)) 
        
        return norm_talweg_measures
        
    def _calculate_warp_simple (self, peak, talweg_line, hydro_line, max_displacement):
        """Calculate the warp information to apply to the talweg in order to fit on the hydro line for a talweg simple.
        
        For a simple talweg the warpping information are calculated as follow:
           - Calculate the distance from the bend peak to the hydrographic line
           - The bend peak is offset by this distance (delta_x and delta_y)
           - The offset is propagated up to the beginning (bend i) and ending (bend j) of the bend
           
        Parameter:
           - peak: Index of the peak position along the talweg
           - talweg_line: MA_LineString contour line
           - hydro_line: LineString representing the hydrographic line that must fit with the talweg
           - max_displacement: Maximum displacement allowed
           
        Return
           - Warppping information needed for the GenUtil.warp_line utilities (see the doc of the method for more information)
                List of tuple containing the information to warp a LineString. 
                Each tuple contains three information (index, delta_x, delta_y)
                    - index: Vertex index number in the LineString (>=0 and < the number of vertex in the line string)
                    - delta_x: Shift to apply in the x coordinate 
                    - delta_y: Shift to apply in the y coordinate
           - None if the displacement is over the maximum displacement allowed
        
        """
            
        # Calculate the projection of the bend peak on the stream and this the max displacement of the talweg
        peak_coord = talweg_line.coords_dual[peak]
        new_talweg_lr = hydro_line.project(Point(peak_coord))     # Extract the linear reference
        new_talweg_point = hydro_line.interpolate(new_talweg_lr)  # Extract the (x,y) coordinate
        delta_x = new_talweg_point.x - peak_coord[0]
        delta_y = new_talweg_point.y - peak_coord[1]
        if ( (delta_x**2.0 + delta_y**2.0)**.5 <= max_displacement ):
            # Warpping is within the tolerance
            warp_infos = []                
            warp_infos.append((0, 0., 0.))
            warp_infos.append((peak, delta_x, delta_y))
            warp_infos.append((len(talweg_line.coords_dual)-1, 0., 0.))
        else:
            # Warpping is outside the tolerance
            warp_infos = None
                
        return warp_infos
            
    def _check_constraints (self, hypso_line, hydro_line, hypso_new_coords, i, j):
        """Check if the new talweg line respect the different constraints: self crossing, line-line crossing, sidedness
        
            Parameters:
                hypso_line: original MA_LineString hypso line
                hydro_line: MA_LineString of the hydrographic feature
                hypso_new_coords: The new coordinates of the hypso line (list of x,y tuple)
                i: Index of the start of the coordinate modification on the line
                j: Index of the end of the coordinate modification on the line
                
            Retrun:
                True: All the constraints are met
                False: One of the constraints is not met
        """
        
        constraints_ok = False
        
        if (hypso_new_coords != []):
            # Calculate the different geometry needed to check if the warp line violated some constraints
            line_simple_line = LineString(hypso_new_coords)
            line_crossing_line = LineString(hypso_new_coords[i:j+1])
            polygon_sidedness = Polygon(hypso_new_coords[i:j+1])
            
            # Checks for any constraints violation and do not log any error
            conflict_type = GenUtil.test_constraints (self, None, line_simple_line, line_crossing_line, 
                                                      polygon_sidedness, self.s_container, 
                                                      [hypso_line._sci_id, hydro_line._sci_id], 
                                                      False)
            if (conflict_type is None):
                # No topological errors in the moved talweg 
                constraints_ok = True
                    
        return constraints_ok

                                  
    def _move_vertex(self, ind, talweg_coords, hydro_line):
        """Move the vertex along the hydro_line in order to form a nice angle with the surrounding vertice  
        
        If the talweg peak is between 2 vertices an extra vertice is added and projected in order 
        to have a rounder talweg instead of a sqare talweg. An iteration process is done for finding
        the best position to add the vertice
        
        Parameters:
            - ind: Index of the vertex to move 
            - hypso_line: Hypsograpgraphic line to process
            - hydro_line: Hydrographic line to process
        
        Return values:
            - List of tuple (x,y) of the talweg line
            
        """
        
        # Some constants
        flat_angle = 1.0
        max_steps = 10
        
        # Extract the angle of the vertice after and before the vertice to move
        angle_a = GenUtil.compute_angle(talweg_coords[ind-2], talweg_coords[ind-1], talweg_coords[ind])
        angle_b = GenUtil.compute_angle(talweg_coords[ind], talweg_coords[ind+1], talweg_coords[ind+2])
        if (math.fabs(180. - angle_a) > flat_angle and
            math.fabs(180. - angle_b) > flat_angle):            
            max_angle = (max(angle_a,angle_b))
            complement = 180. - max_angle
            target_angle = max_angle + complement/2.
            max_distance = GenUtil.distance(talweg_coords[ind-1], talweg_coords[ind+1]) # Max distance to move the the vertex
            pace = max_distance / max_steps
            coord = talweg_coords[ind]
            distance_lr = hydro_line.project(Point(coord))
                    
            # Using an iterative process, loop to find the best position of the new vertex
            # Iteration stops when the it match the target angle or the line is extended at its maximum            
            for i in xrange(1, max_steps):
                tentative_distance = distance_lr - (i*pace)
                new_point = hydro_line.interpolate(tentative_distance)
                new_coord = new_point.coords[0]
                talweg_coords[ind] = new_coord
                
                # Computes new angles 
                angle_a = GenUtil.compute_angle (talweg_coords[ind-2], talweg_coords[ind-1], talweg_coords[ind])
                angle_b = GenUtil.compute_angle (talweg_coords[ind], talweg_coords[ind+1], talweg_coords[ind+2])
                if (angle_a > target_angle or angle_b > target_angle):
                    break
        else:
            # Angle between is too small trying to correct the angle will cause more problem...
            pass

        return talweg_coords
    
    def _add_vertex(self, i, j, lst_coords):
        """Add a coordinate in the midlle of 2 coordinates in a list 
        
        Parameters:
            - i: Start index position to add a coordinate 
            - j: End index position to add a coordinate
            lst_coords: List of coordinates
        
        Return values:
            - List of tuple (x,y) with the mid coordinate added
            
        """
        
        if (i != j):
            vertex_added = True
            mid_coord = GenUtil.mid_point(lst_coords[i], lst_coords[j])
            if (j-i == 1):
                peak = i+1
                lst_coords.insert(peak, mid_coord)
            else:
                line = LineString(lst_coords[i:j+1])
                mid_distance = line.project(Point(mid_coord))
                peak = j-1
                for k in range(i+1,j):
                    distance = line.project(Point(lst_coords[k]))
                    if distance > mid_distance:
                        if (math.fabs(distance-mid_distance) > 1.0):
                            peak = k+1
                            lst_coords.insert(peak, mid_coord)
                        else:
                            peak = k
                        break
                    else:
                        peak = k
        else:
            # Nothing to do
            vertex_added = False
            peak = i
                      
        return (peak,vertex_added)    
            
    def _are_vertice_editable(self, coord_i, coord_j, hypso_line):
        """Checks if the vertice between coord_i and coord_j are editable
        
        Vertice are not editable when they were moved before
        
        Parameters:
            - coord_i: Start of the search
            - coord_j: end of the search
            - hypso_line: Hypso line to search
        
        Return value:
            True: The coordinates are editable
            False: The coordinates are uneditable
            
        """
        
        to_edit = range(coord_i, coord_j+1)
        set_to_edit = set(to_edit)
        set_edited = hypso_line.vertice_edited & set_to_edit
        if (len(set_edited) == 0):
            editable = True
        else:
            editable = False
    
        return editable
    
    def _flag_vertice_as_edited(self, coord_i, coord_j, hypso_line):
        """Flags the vertice located between coord_i and coord_j as edited
        
        Parameters:
           - i: Start vertice index to flag
           - j: End vertice index to flag
           - hypso_line: LineString hypso line containing the vertice to flag 
        
        """
        
        for coord in range(coord_i, coord_j+1):
            # Insert the index in the set of edited vertice
            hypso_line.vertice_edited.add(coord)
        
    def _adjust_vertice_position(self, index, hypso_line):
        """Adjust the position of the vertices in the bends and in the set vertice_added
        
        When a vertex is added in a line once must adjust the bend position and the list of vertice edited
        
        Parameters:
            index: Index of the position of the vertice added
            hypso_line: MA_LineString with bend information
            
        Return value: None
        
        """

        # Adjust bend position 
        for bend in hypso_line.bends:
            if bend.i >= index:
                bend.i += 1
            if bend.j >= index:
                bend.j += 1
                
        # Adjust vertice edited
        lst_vertice = [vertice for vertice in hypso_line.vertice_edited if vertice >= index]
        lst_vertice.sort(reverse=True)
        for vertice in lst_vertice:
            hypso_line.vertice_edited.remove(vertice)
            hypso_line.vertice_edited.add(vertice+1)
            
    def _hydro_intersect_segment(self, bend_i, bend_j, hypso_line, hydro_line):
        """Determine between vertice the hydro line crosses the talweg line
        
        Parameters:
            - bend_i: Start index position in the hypso line
            - bend_j: End index position in the hypso line
            - hypso_line: MA_LineString hypsographic line
            - hydro_line: MA_LineString hydrographic line
        """
        
        talweg_coords = hypso_line.coords_dual[bend_i:bend_j+1]
        talweg_line = LineString(talweg_coords)
        hypso_hydro_point = talweg_line.intersection(hydro_line)
        distance_crossing = talweg_line.project(hypso_hydro_point)
        
        # Loop over each vertex to find where the hydro line crosses the talweg line
        for i in xrange(1, len(talweg_coords)):
            distance = talweg_line.project(Point(talweg_coords[i]))
            if distance == distance_crossing:
                #Special case where the hydro fits exactly on a vertex
                vertice_min = i-1
                vertice_max = i+1
                if i == 0:
                    vertice_min = i
                if i == len(talweg_coords)-1:
                    vertice_max = i
                break
            else:
                if distance > distance_crossing:
                    # Extract the vertice before and after the line crossing
                    vertice_min = i-1
                    vertice_max = i
                    break
                    
        return (vertice_min+bend_i,vertice_max+bend_i)
            
    def _move_talweg_simple (self, bend, hypso_line):
        """Displace the simple talweg position in order to fit with the hydropgraphic line
        
        For one simpele talweg this method is trying to displace the talweg to make it fit on the hydro line
                
        Parameters:
            - bend: Bend containing the talweg to edit
            - hypso_line: MA_LineString of the hypsographic line to edit
            
        Return value:
            - True: Talweg is moved
            - False: Talweg was not moved
            
        """
        
        talweg_measures = []
        talweg_edited = False
        max_coords = len(hypso_line.coords_dual)
        max_displacement = hypso_line.ma_properties[_MAX_DISPLACEMENT]
        
        # Number of vertice before and after the talweg definition we are trying to find the best talweg
        scope = 6
        
        # Extract the index interval on the talweg line where the hydro line intersects it
        (vertice_min ,vertice_max) = self._hydro_intersect_segment(bend.i, bend.j, hypso_line, bend.hydro_line)

        for i_delta in range(-scope,scope+1):
            for j_delta in range(-scope,scope+1):
                bend_i = bend.i + i_delta
                bend_j = bend.j + j_delta
                if (bend_j-bend_i >= 2      and  # A talweg is least 3 vertice
                    bend_i >= 0             and  # Don't overflow the list
                    bend_j <= max_coords-1  and  # Don't overflow the list
                    bend_i <= vertice_min   and  # Be sure to have hydro line inside the talweg
                    bend_j >= vertice_max):      # Be sure to have hydro line inside the talweg

                    # Check if we can modify the vertice
                    if (self._are_vertice_editable(bend_i+1, bend_j-1, hypso_line)):
                        talweg_ori_coords = list(hypso_line.coords_dual[bend_i:bend_j+1])
                        talweg_line = MA_LineString(talweg_ori_coords)
                        # Check if the poteltial talweg meets some minimal criteria's (for performance only...)
                        if (self._talweg_screening(max_displacement, bend.hydro_talweg_point, talweg_line) ):                                
                            peak_kl = self._extract_bend_peak(0, len(talweg_ori_coords)-1, talweg_line)
                            (peak,vertex_added)  = self._add_vertex(peak_kl[0], peak_kl[1], talweg_ori_coords)
                            talweg_ori_line = MA_LineString(talweg_ori_coords)
                    
                            warp_infos = self._calculate_warp_simple(peak, talweg_ori_line, bend.hydro_line, 
                                                                     hypso_line.ma_properties[_MAX_DISPLACEMENT])
                            
                            if (warp_infos is not None):
                                talweg_new_line = GenUtil.warp_line(talweg_ori_line,warp_infos)
                                talweg_new_coords = list(talweg_new_line.coords)
                                # Extract the measures on the potential talweg
                                measures = self._extract_talweg_measure(talweg_ori_coords, talweg_new_coords, 
                                                                        bend.hydro_line)
                                if measures is not None:
                                    # Accumulate the measure
                                    talweg_measures.append(measures+[bend_i, bend_j, peak, vertex_added, 
                                                                     talweg_ori_coords, talweg_new_coords])
                                else:
                                    # Do not count this measure
                                    pass
                            else:
                                # Warp info is None because outside of tolerance
                                pass
                    
        if (talweg_measures):
            # Normalize the measure to be able to compare them
            norm_talweg_measures = self._normalize_talweg_measures(talweg_measures)
            # Sort the measure the smallest is the best
            norm_talweg_measures.sort()
            for norm_talweg_measure in norm_talweg_measures:
                bend_i = norm_talweg_measure[1]
                bend_j = norm_talweg_measure[2]
                peak = norm_talweg_measure[3]
                vertex_added = norm_talweg_measure[4]   
                talweg_ori_coords = norm_talweg_measure[5]
                talweg_new_coords = norm_talweg_measure[6]
                if (vertex_added):
                    # If the peak is between 2 vertice move the vertice
                   talweg_new_coords =  self._move_vertex(peak, talweg_new_coords, bend.hydro_line)
                hypso_new_coords = list(hypso_line.coords_dual)
                hypso_new_coords[bend_i:bend_j+1] = talweg_new_coords

                talweg_ori_line = MA_LineString(talweg_ori_coords)
                talweg_new_line = MA_LineString(talweg_new_coords)
                # Check that the original potential talweg and the resulting potential talweg were valid talweg
                if (self._extract_talweg_type(talweg_ori_line, bend.hydro_line) == _TALWEG_TYPE_SIMPLE and
                    self._extract_talweg_type(talweg_new_line, bend.hydro_line) == _TALWEG_TYPE_SIMPLE):
                    # Check the modified talweg do not break any constraint         
                    if (self._check_constraints (hypso_line, bend.hydro_line, hypso_new_coords, bend_i, bend_j) ):
                        hypso_line.update_coords(hypso_new_coords, self.s_container)
                        self._flag_vertice_as_edited(bend_i+1, bend_j-1, hypso_line)
                        if (bend_j-bend_i != len(talweg_new_coords)-1):
                            # Adjust vertice position because a vertice was added
                            self._adjust_vertice_position(bend_i+peak_kl[1], hypso_line)
                        talweg_edited = True
                        break
                
        if (not talweg_edited):
            # Every tentative have failed log this talweg as an error
            talweg_line = LineString(hypso_line.coords_dual[bend.i:bend.j+1])
            GenUtil.add_err_position (self, None, talweg_line, _STAT_EDITED_ERROR)
                
        return talweg_edited
                        
    def _manage_talweg_simple (self):
        """Processes all the talwegs of all the lines
               
        Parameters: None
            
        Return value:  None
            
        """
        
        for hypso_line in self.hypso_lines:

            for bend in hypso_line.bends:
                
                if bend.type == _TALWEG_TYPE_SIMPLE:
                
                    if self._move_talweg_simple(bend, hypso_line):
                        self.stats.add_stats(_STAT_EDITED_SIMPLE)
                        bend.status = _EDITED
                        
    def _manage_talweg_inverted (self):
        """Manage the displacement of  the simple talweg position in order to fit with the hydropgraphic line
        
        The list normalized_measures contains the sorted optimized talweg position. This method
        is doing the following tasks until the first instance is meeting the following conditions
               - Warp the talweg if the warpping in within the tolerance
               - Check if the resulting talweg meet the criterias to be a valid talweg
               - Check that no topological constraints are violated 
                
        Parameters: None
            
        Return value:  None
            
        """
        
        infinite = 1.0e+100
        for hypso_line in self.hypso_lines:
            max_bends = len(hypso_line.bends)
            for bend_nbr, bend in enumerate(hypso_line.bends):
                if bend.type == _TALWEG_TYPE_INVERTED:
                    distance_before = infinite
                    distance_after = infinite
                    if (bend_nbr-1 >= 0):
                        if hypso_line.bends[bend_nbr-1].type == _TALWEG_TYPE_EMPTY:
                            bend_before = hypso_line.bends[bend_nbr-1]
                            bend_before.type = _TALWEG_TYPE_SIMPLE
                            bend_before.hydro_line = bend.hydro_line
                            peak_kl = self._extract_bend_peak(bend_before.i, bend_before.j, hypso_line)
                            peak_coord = self._extract_peak_coord ( peak_kl[0], peak_kl[1], hypso_line)
                            distance_before = bend.hydro_line.distance(Point(peak_coord))
                    if (bend_nbr+1 <= max_bends-1):
                        if hypso_line.bends[bend_nbr+1].type == _TALWEG_TYPE_EMPTY:
                            bend_after = hypso_line.bends[bend_nbr+1]
                            bend_after.type = _TALWEG_TYPE_SIMPLE
                            bend_after.hydro_line = bend.hydro_line
                            peak_kl = self._extract_bend_peak(bend_after.i, bend_after.j, hypso_line)
                            peak_coord = self._extract_peak_coord ( peak_kl[0], peak_kl[1], hypso_line)
                            distance_after = bend.hydro_line.distance(Point(peak_coord))
                            
                    target_bend = None
                    bend_edited = False
                    if (distance_before != infinite and distance_before <= distance_after):
                        target_bend = bend_before
                    if (distance_after != infinite and distance_after <= distance_before):
                        target_bend = bend_after
                    if target_bend is not None:
                        bend_edited = self._move_talweg_simple(target_bend, hypso_line)
                            
                    if bend_edited:
                        self.stats.add_stats(_STAT_EDITED_INVERTED)
                    
###    def _manage_vertex_added(self, bend_i, bend_j, talweg_new_coords, hypso_line):
###        """Remove the colinear vertex on the hypso feature
###        
###        Removing the colinear features is needed for adding and moving correctly the extra vertex
###        
###        Parameters: None
###        
###        Return value: None
###        
###        """
###        
###        ori_talweg_coords = hypso_line.coords_dual[bend_i:bend_j+1]
###        
###        lst_angles = GenUtil.compute_angles(ori_talweg_coords)
        

###    def _warp_line_inverted (self, hypso_line, hydro_line, i_bend, talweg_i, talweg_j, offset, max_displacement):
###        """Warp the talweg to solve the problem of a talweg inverted.
###        
###        For an inverted talweg, the line is warpped is 2 times:
###        Firt:
###           - Extract the bend peak and calculated the spatial position
###           - Calculate the distance from the bend peak to the hydrographic line
###           - Calculate the middle position of the bend base
###           - Warp the bend peak to fit of the hydro line (delta_x and delta_y)
###           - Warp the bend base but by a factor of .7 instead of 1 this help 
###             reduce the impact on the previous and nexct bend
###           - If offset is positive the warpping propagation must be done from bend-1 up to bend+2
###           - If offset is negative the warpping propagation must be done from bend-2 up to bend+1
###           - If the bend peak is between 2 vertex we take the vertex closer to the hydro line
###        Second:
###           - The twos bend preceding or following (depending of the type of offset) are 
###             reduced in height as they were shrink by the warping of the talweg
###           
###        Parameter:
###           - hypso_line: MA_LineString contour line
###           - hydro_line: LineString representing the hydrographic line that must fit with the talweg
###           - i_bend: Index of the bend to process in the hypso_line
###           - max_displacement: Maximum displacement allowed
###           - offset: If positive the talweg is moved in forward position;
###                     If negative the talweg is moved in backward position;
###           
###        Return
###           - Tuple of 3 value
###                 - List of (x,y) tuple of the new coordinates of the line or None if the line was not warped
###                 - Index of the start of the modification on the line
###                 - Index of the end of the modification on the line
###           - True: The line is warpped
###           - False: The line was not warpped mainly because it offsets the max_displacement tolerance
###        
###        """
###        
###        i_start = None
###        j_end = None
###        
###        if (offset < 0):
###            bend_talweg = i_bend - 1
###        else:
###            bend_talweg = i_bend + 1
###            
###        (peak_k, peak_l) = self._extract_bend_peak(talweg_i, talweg_j, hypso_line)
###        
###        # Take the shortest peak position depending on the offset
###        if (offset > 0):
###            peak_k,peak_l = peak_k, peak_k
###        else:
###            peak_k,peak_l = peak_l, peak_l
###            
###        # Extract the peak coord for the bend
###        peak_coord = self._extract_peak_coord (peak_k, peak_l, hypso_line)
###            
###        # Calculate the warping of the bend peak on the stream
###        talweg_peak_lr = hydro_line.project(Point(peak_coord))      # Extract the linear reference
###        talweg_peak_point = hydro_line.interpolate(talweg_peak_lr)  # Extract the (x,y) coordinate
###        peak_delta_x = talweg_peak_point.x - peak_coord[0]
###        peak_delta_y = talweg_peak_point.y - peak_coord[1]
###        peak_displacment = (peak_delta_x**2.0 + peak_delta_y**2.0)**.5
###        
###        # Calculate the warping of the base of the talweg
###        mid_base_coord = GenUtil.mid_point(hypso_line.coords_dual[bend_talweg], hypso_line.coords_dual[bend_talweg])
###        talweg_mid_base_lr = hydro_line.project(Point(peak_coord))      # Extract the linear reference
###        talweg_mid_base = hydro_line.interpolate(talweg_mid_base_lr)  # Extract the (x,y) coordinate
###        base_delta_x = (talweg_mid_base.x - peak_coord[0]) * 0.7
###        base_delta_y = (talweg_mid_base.y - peak_coord[1]) * 0.7
###        base_displacment = (base_delta_x**2.0 + base_delta_y**2.0)**.5
###        
###        # Check if the displacement is within the tolerance
###        if (peak_displacment <= hypso_line.max_displacement and
###            base_displacment <= hypso_line.max_displacement):
###        
###            # Build the warping information for the first displacement
###            warp_infos = []
###            i_start = hypso_line.bends[i_bend-2].i
###            j_end = hypso_line.bends[i_bend+2].j
###            warp_infos.append((hypso_line.bends[i_bend-2].i, 0., 0.))
###            if (hypso_line.bends[i_bend-2].i != talweg_i):
###                warp_infos.append((talweg_i, base_delta_x, base_delta_y))
###            warp_infos.append((peak_k, peak_delta_x, peak_delta_y))
###            if (talweg_j != hypso_line.bends[i_bend+2].j):
###                warp_infos.append((talweg_j, base_delta_x, base_delta_y))
###            warp_infos.append((hypso_line.bends[i_bend+2].j, 0., 0.))
###            
###            hypso_new_coords = GenUtil.warp_line(hypso_line,warp_infos)
###            
###            # Calcultae the reduction of the inverted talweg
###            i = hypso_line.bends[i_bend].i
###            j = hypso_line.bends[i_bend].j
###            (peak_k, peak_l) = self._extract_bend_peak(i, j, hypso_line)
###            peak_coord = self._extract_peak_coord (peak_k, peak_l, hypso_line)
###            base_line = LineString([hypso_line.coords_dual[i],hypso_line.coords_dual[i]] )
###            base_line_lr = base_line.project(Point(peak_coord))      # Extract the linear referenc
###            base_line_point = base_line.interpolate(base_line_lr)    # Extract the (x,y) coordinate
###            delta_x = (peak_coord[0] - base_line_point.x) * 0.3 
###            delta_y = (peak_coord[1] - base_line_point.y) * 0.3
###            invert_displacment = (delta_x**2.0 + delta_y**2.0)**.5
###            
###            # Check that the displacement is within the tolerance
###            if (invert_displacment <= hypso_line.max_displacement):           
###                # apply the reduction on the bend
###                warp_infos = []
###                warp_infos.append((i, 0., 0.))
###                warp_infos.append((peak_k, -delta_x, -delta_y))
###                warp_infos.append((j, 0., 0.))
###                hypso_new_coords = GenUtil.warp_line(LineString(hypso_new_coords),warp_infos)
###            else:
###                # Do not apply the disaplcement
###                pass
###            
###        else:
###            # There is no displacement
###            hypso_new_coords = []
###        
###        return (hypso_new_coords, i_start, j_end)
###
###    def _adjust_bend_base (self, hypso_line, hydro_line, bend):
###        """Modify the position of the bend base in order to have a better approximation of the bend.
###        
###        With a better approximation of the bend it is possible to calculate a better bend peak and create 
###        a more natural displacement for the talweg.  The bend base is adjusted by moving the bend base (i,j) 
###        by plus or minus two vertice.
###        A better position for the bend base is met when we find the maximum for the 3 following measures: 
###            - The length of the base line is the shortest
###            - We maximize the amplitude of the bend
###            - The bend base is forming a right angle (90 degrees) with the talweg
###                
###        Parameters: 
###            - hypso_line: MA_LineString hypsographic line to process
###            - hydro_line: MA_LineString hydrographic line to process
###            - bend: bend to process
###            
###        Return value:
###            - Sorted normalized talweg measures. It is a list of the values
###                    - Normalize talweg value
###                    - i: Index of the start of the bend
###                    - j: Index of the end of the bend
###            
###        """
###        
###        i = bend.i
###        j = bend.j
###        talweg_measures = []
###
###        for delta_i in range(-2,3):
###            for delta_j in range(-2,3):
###                talweg_measures.append(self._extract_talweg_measure (hypso_line, hydro_line, i+delta_i, j+delta_j) )
###
###        # Remove None values from the talweg_measures list
###        talweg_measures = [talweg_measure for talweg_measure in talweg_measures if talweg_measure is not None]
###        
###        # Normalize [0..1] the measures in order to be able to compare the measures
###        norm_talweg_measures = self._normalize_talweg_measures (talweg_measures)
###    
###        # Sort the normalized valued the smallest in the one that maximise the evaluation of a good talweg
###        norm_talweg_measures.sort()
###
###        return norm_talweg_measures
###
###    def _analyze_talweg_position (self):
###        """Analyze the talwegs position along a line and flag the talweg that do not meet criterias to be corrected  
###        
###        The bend analysis is done in 2 passes:
###        The first pass checks the inverted talweg: If an inverted talweg is just before or after a simple talweg we remove 
###                                                   the correction flag on the talweg inverted in order to correct only 
###                                                   the simple talweg 
###        The second pass checks the simple talweg: To be corrected a simple talweg must be between 2 empty talweg. 
###                                                  This will give some room/space for the correction of the talweg.
###        
###        
###        Parameters: None
###            
###        Return value:  None
###            
###        """
###        
###        for pass_nbr in range(2): # Manage the 2 passes
###            
###            for hypso_line in self.hypso_lines:    
###                max_bend_i = len(hypso_line.bends)-1
###                for i, bend in enumerate(hypso_line.bends):
###                    
###                    # Manage the condition for the inverted talweg
###                    if (pass_nbr == 0):
###                        if bend.type == _TALWEG_TYPE_INVERTED:
###                            # if the previous or the next bend is a bend simple we can convert the current bend type to talweg empty
###                            if (self._bend_type(i-1, hypso_line) == _TALWEG_TYPE_SIMPLE or
###                                self._bend_type(i+1, hypso_line) == _TALWEG_TYPE_SIMPLE):
###                                    bend.type = _TALWEG_TYPE_EMPTY
###                            else:
###                                # Check if we can transform the talweg inverted into a talweg simple by assign the talweg correction 
###                                # to one of its neighbour. A talweg simple is always easier to solve than a talweg inverted
###                                if (i != 0):
###                                    bend_before = hypso_line.bends[i-1]
###                                    talweg_coords = hypso_line.coords_dual[bend_before.i:bend.j+2]
###                                    talweg_line = LineString(talweg_coords)
###                                    talweg_type = self._extract_talweg_type(talweg_line, bend.hydro_line)
###                                    if (talweg_type == _TALWEG_TYPE_SIMPLE):
###                                        if self._bend_type(i+2, hypso_line) == _TALWEG_TYPE_EMPTY:
###                                            bend.type = _TALWEG_TYPE_EMPTY
###                                            bend_before.type = _TALWEG_TYPE_SIMPLE
###                                            bend_before.hydro_line = bend.hydro_line
###                                            bend_before.j += 1
###
###                                if (i != max_bend_i):
###                                    bend_after = hypso_line.bends[i+1]
###                                    talweg_coords = hypso_line.coords_dual[bend_after.i-1:bend_after.j+1]
###                                    talweg_line = LineString(talweg_coords)
###                                    talweg_type = self._extract_talweg_type(talweg_line, bend.hydro_line)
###                                    if (talweg_type == _TALWEG_TYPE_SIMPLE):
###                                        if self._bend_type(i+2, hypso_line) == _TALWEG_TYPE_EMPTY:
###                                            bend.type = _TALWEG_TYPE_EMPTY
###                                            bend_after.type = _TALWEG_TYPE_SIMPLE
###                                            bend_after.hydro_line = bend.hydro_line
###                                            bend_after.i -= 1
###                                                    
###                    # Manage the case for the simple talweg
###                    if (pass_nbr == 1):
###                        if bend.type == _TALWEG_TYPE_SIMPLE:
###                            # Check that the Previous and Next bend are empty in order to give some roome for the talweg correction
###                            if (self._bend_type(i-1, hypso_line) == _TALWEG_TYPE_EMPTY and
###                                self._bend_type(i+1, hypso_line) == _TALWEG_TYPE_EMPTY): 
###                                pass  # Nothing special to do; the talweg will be moves
###                            else:
###                                # Flag the current bend as to complex to be corrected
###                                # Improvment it would be nice when this case happens to choose the best talweg to keep  instead of the first one...
###                                bend.type = _TALWEG_TYPE_EMPTY
###
###    def _distance_hydro_talweg(self, hydro_line, hypso_line, i, j):
###        """Calculates the distance between the talweg and the hydro line
###        
###        To calculate the distance the method uses the distance between the 
###        talweg centroid (http://en.wikipedia.org/wiki/Centroid) and the hydro line.
###        The talweg is always almost a convex polygon and the centroid of a convex polygon is
###        always located inside the area
###        
###        Parameters:
###            - hydro_line: MA_LineString hydro line to process
###            - hypso_line: MA_LineString hypso line to process
###            - i: start index of the talweg on the hypso line
###            - j: end index of the talweg on the hypso line
###            
###        Return:
###            - Distance between the talweg and the hydro line
###          
###        """
###        
###        talweg = Polygon(hypso_line.coords_dual[i:j+1])
###        distance = talweg.centroid.distance(hydro_line)
###        
###        return distance
###    
###    def _adjust_bend_inverted (self, hypso_line, i_bend):
###        """Determines the bend to use in order to correct an inverted talweg
###        
###        Parameter:
###            - hypso_line: Ma_LineString to process
###            - i_bend: Index of the bend in thew hypso_line
###            
###        Return value:
###            - Sorted normalized talweg measures. It is a list of the values
###                    - Normalize talweg value
###                    - i: Index of the start of the bend
###                    - j: Index of the end of the bend
###            - Empty list []. An empty list means the previous or next talweg are not good talweg candidate.
###        
###        """
###        
###        norm_talweg_measures = []
###        offset = None
###        
###        if (hypso_line.bends[i_bend].type == _TALWEG_TYPE_INVERTED):
###            if (i_bend-2 >= 0) and (i_bend+1 <= len(hypso_line.bends)-2):
###                # Checks if the previous and next bend are empty 
###                # It is a too complex case to manage when previous or next bends are not empty 
###                if (hypso_line.bends[i_bend-2].type == _TALWEG_TYPE_EMPTY and
###                    hypso_line.bends[i_bend-1].type == _TALWEG_TYPE_EMPTY and
###                    hypso_line.bends[i_bend].type   == _TALWEG_TYPE_INVERTED and
###                    hypso_line.bends[i_bend+1].type == _TALWEG_TYPE_EMPTY and
###                    hypso_line.bends[i_bend+2].type == _TALWEG_TYPE_EMPTY):
###                    
###                    # For the previous and next talweg adjust the talweg to find the best talweg form
###                    hydro_line = hypso_line.bends[i_bend].hydro_line
###                    previous_norm_talweg_measures = self._adjust_bend_base (hypso_line, hydro_line, hypso_line.bends[i_bend-1])
###                    next_norm_talweg_measures = self._adjust_bend_base (hypso_line, hydro_line, hypso_line.bends[i_bend+1])
###
###                    # Distance from the previous bend to the hydro line
###                    if previous_norm_talweg_measures != []:
###                        i = previous_norm_talweg_measures[-1][1]
###                        j = previous_norm_talweg_measures[-1][2]
###                        talweg = Polygon(hypso_line.coords_dual[i:j+1])
###                        previous_distance_talweg_hydro = talweg.centroid.distance = talweg.centroid.distance(hydro_line)
###                    else:
###                        previous_distance_talweg_hydro = sys.float_info.max
###                                                                                 
###                    # Distance from the next bend to the hydro line
###                    if next_norm_talweg_measures != []:
###                        i = next_norm_talweg_measures[-1][1]
###                        j = next_norm_talweg_measures[-1][2]
###                        talweg = Polygon(hypso_line.coords_dual[i:j+1])
###                        next_distance_talweg_hydro = talweg.centroid.distance = talweg.centroid.distance(hydro_line)
###                    else:
###                        next_distance_talweg_hydro = sys.float_info.max
###                    
###                    # The best candidate talweg is the one closest to the hydro line
###                    if ( previous_distance_talweg_hydro < next_distance_talweg_hydro ):
###                        norm_talweg_measures = previous_norm_talweg_measures
###                        offset = -1
###                    else:
###                        norm_talweg_measures = next_norm_talweg_measures
###                        offset = +1
###                        
###                    # The talweg most also form an amplitude of 75 degrees to be considered a valid talweg
###                    i = norm_talweg_measures[0][1]
###                    j = norm_talweg_measures[0][2]
###                    lst_coords = hypso_line.coords_dual[i:j+1]
###                    amplitude = self._extract_bend_amplitude(lst_coords)
###                    if (amplitude < 75.):
###                        norm_talweg_measures = []
###                        
###        return (norm_talweg_measures, offset)
###    
###    def _move_talweg_inverted (self):
###        """Displace the simple talweg position in order to fit with the hydropgraphic line  
###                
###        Parameters: None
###            
###        Return value:  None
###            
###        """
###        
###        for hypso_line in self.hypso_lines:    
###            for i_bend, bend in enumerate(hypso_line.bends):
###                (norm_talweg_measures, offset) = self._adjust_bend_inverted (hypso_line, i_bend)
###                if (norm_talweg_measures != []):
###                    bend_edited = False
###                    # Loop over each measure 
###                    for norm_talweg_measure in norm_talweg_measures:
###                        talweg_i = norm_talweg_measure[1]
###                        talweg_j = norm_talweg_measure[2]
###                        (hypso_new_coords, i_start, j_end) = self._warp_line_inverted (hypso_line, bend.hydro_line, i_bend, 
###                                                                                       talweg_i, talweg_j, offset, 
###                                                                                       hypso_line.max_displacement)
###                        
###                        if (hypso_new_coords != []):
###                            constraints_ok = self._check_constraints (hypso_line, bend.hydro_line, hypso_new_coords, 
###                                                                      i_start, j_end)   
###                            if (constraints_ok):
###                                bend_edited = True
###                                break
###                    
###                    if (bend_edited):
###                        self.stats.add_stats(_STAT_EDITED_INVERTED)
###                        hypso_line.update_coords(hypso_new_coords)
###                    else:
###                        # Every tentative have failed log this talweg as an error
###                        talweg_line = LineString(hypso_line.coords_dual[bend.i:bend.j+1])
###                        GenUtil.add_err_position (self, hypso_line._sif_id, talweg_line, _STAT_EDITED_ERROR)

###    def _extract_peak_coord (self, index_k, index_l, line):
###        """Extract the (x,y) coordinate located add index_k, index_l
###        
###        If index_k and index_l are the same value we extract the coordinate located at index_k 
###        along the line if index_k and index_l are diffrent we extract the middlew point between 
###        inde_k and index l
###        
###        Parameters: 
###           - index_k: first index along the line
###           - index_l: Second index along the line
###           - line: MA_LineString to process
###           
###        Return value:
###            - tuple (x,y) of the requested value
###        """
###        
###        if (index_k == index_l):
###            coord = line.coords_dual[index_k]
###        else:
###            coord = GenUtil.mid_point(line.coords_dual[index_k], line.coords_dual[index_l])
###        
###        return coord
###    def _bend_type(self, i_bend, hypso_line):
###        """Extract the type of bend
###        
###        If the index is outside of the bound than it return _TALWEG_TYPE_EMPTY  
###        
###        Parameters:
###            i_bend: Index of the bend to extract
###            hypso_line: LineString with Bend definition
###            
###        Return:
###            The bend type
###        """
###        
###        try:
###            bend_type = hypso_line.bends[i_bend].type
###        except IndexError:
###            bend_type = _TALWEG_TYPE_EMPTY
###            
###        return bend_type
###
###    def _extract_angle_hydro_talweg(self, talweg_line, hydro_line):
###        """Calculate the angle at the intersection of the hydro line and the talweg
###        
###        Parameters:
###            - talweg_line: LineString of the talweg
###            - hydro_line: LineString of the hydro_line
###            
###        Return value:
###            - Tuple of float: Angles formed by the hydro and the talweg
###            
###        """
###
###        crossing_point =  hydro_line.intersection(talweg_line)
###        tuple_angle = (0,0)
###        if (isinstance(crossing_point, Point)):
###            distance_lr = talweg_line.project(crossing_point)
###            
###            # Determine a point on the talweg line just a little bit before and after the crossing point
###            # in order to dertermine an angle
###            delta = talweg_line.length/1000.
###            distance_lr = talweg_line.project(crossing_point)
###            talweg_before_crossing = talweg_line.interpolate(distance_lr-delta) 
###            talweg_after_crossing = talweg_line.interpolate(distance_lr+delta)
###            # Determine a point on the hydro line just a little bit before the crossing point
###            # in order to dertermine an angle
###            distance_lr = hydro_line.project(crossing_point)
###            hydro_before_crossing = hydro_line.interpolate(distance_lr+delta)
###            angle_hydro_a = GenUtil.compute_angle(talweg_before_crossing.coords[0], crossing_point.coords[0], 
###                                                  hydro_before_crossing.coords[0])
###            angle_hydro_b = GenUtil.compute_angle(talweg_after_crossing.coords[0], crossing_point.coords[0], 
###                                                  hydro_before_crossing.coords[0])
###    
###            tuple_angle = (angle_hydro_a, angle_hydro_b)
###            
###        return tuple_angle 
###    def _load_features_talweg(self):
###        """Load the features in the spatial container
###        
###        The features are loaded in the spatial container but the hydrographic area features
###        area splitted and loaded as there individual linear rings interior and exterior composing
###        the area.
###        
###        Parameters: None
###        
###        Return values: None
###        
###        """
###        
###        self.polygons_original = self.polygons
###        self.polygons = []
###        
###        # Load all the features except hydro polygon in the spatial container
###        self.s_container = self.load_features()
###        self.s_container.add_features(self.polygons_original)
###        
###        for polygon in self.polygons_original:
###            linear_rings = []
###            linear_rings.append(MA_LineString(list(polygon.exterior.coords)))
###            for interior in polygon.interiors:
###                linear_rings.append(MA_LineString(list(interior.coords)))
###                
###            # Sets some properties
###            for linear_ring in linear_rings:
###                linear_ring.id = -1
###                linear_ring.code = HYDRO
                
###            # Load the rings of the polygon in the spatial container
###            self.s_container.add_features(linear_rings)
###                
###                             
###                                                
###    def _extract_features_talweg(self):
###        """Prepare the features for the extraction
###        
###        When the features were loaded hydrographic features were transformed into linear ring features. This
###        method is removing the linear ring features and replaces them with the original polygon features
###        
###        Parameters: None
###        
###        Return value: None
###        
###        """
###        
###        # Extract the linear ring
###        linear_rings = self.s_container.get_features(filter="feature._sif_id==-1 and " + 
###                                                            "feature.code=='%s' and " %(HYDRO) + 
###                                                            "feature.feature_type==GenUtil.LINE_STRING" )
###        # Delete the linear rings in the spatial container
###        self.s_container.del_features(linear_rings)
###        
###        # Add the original polygon in the spatial container
###        self.s_container.add_features(self.polygons_original)
###        
###        # Extract all the features for the output
###        self.extract_features_out(self.s_container.get_features())
        
    def _clean_up(self):
        """Delete unnecessary attributes added on the features during the process
        
        Spatial container containes cyclic reference that must be delete manually otherwise the
        garbage collector is unable to free the memory
        
        Parameters: None
        
        Return value: None
        
        """
        
        for line in self.hypso_lines:
            if hasattr(line, "vertice_edited"): del line.vertice_edited
            if hasattr(line, "bends"): del line.bends
            
        for line in self.hydro_lines:
            if hasattr(line, "buffer_line"): del line.buffer_line
            if hasattr(line, "buffer_weight"): del line.buffer_weight    
    
    def _free_memory(self):
        """Delete variables in order to free some memory
        
        Spatial container containes cyclic reference that must be delete manually otherwise the
        garbage collector is unable to free the memory
        
        Parameters: None
        
        Return value: None
        
        """
        
        del self.hypso_lines
        del self.hydro_lines
        del self.s_container
        
    def process(self):
        """Main routine for the Talweg coherence algorithm
        
        This algorithm will edit the contour line in order to fit with the Talweg line.
        It will prevent line crossing and sidedness errors.
        
        Parameters: None
            
        """
        
        self.s_container = self.load_features()
     
        # Extract hypso and hydro lines
        self.hypso_lines = self.s_container.get_features(filter="feature.ma_properties['code']=='%s' and feature.feature_type==GenUtil.LINE_STRING" %(HYPSO) )
        self.hydro_lines = self.s_container.get_features(filter="feature.ma_properties['code']=='%s' and feature.feature_type==GenUtil.LINE_STRING" %(HYDRO) )

        # Detect the position of the bends on each line
        self._detect_bend_location()
        
        # Add extra information on each bend
        self._add_line_attributes()
        
        # Detect the talweg position using the bend position
        self._detect_talweg()
            
        # Move simple talweg
        self._manage_talweg_simple()
        
        # Move inverted talweg
#        self._manage_talweg_inverted()
                    
        # Clean up attributes
        self._clean_up()
        
        # Extract all the features for the output
        self.extract_features_out(self.s_container.get_features())
        
        self._free_memory()
        
        GenUtil.print_debug (self.params, "End of %s algorithm" %(_TALWEG))