#!/usr/local/bin/python
# -=- encoding: utf-8 -=-

#####################################################################################################################################

"""
This algorithm implements the FirstLastPivoter algorithm to move first/last vertices to a straighter portion of the line feature.

Usage:
    import algo_first_last_pivoter  OR
    from algo_first_last_pivoter import AlgoFirstLastPivoter

Limits and constraints
    None

__revision__ = "--REVISION-- : $Id: first_last_pivoter.py 128 2010-12-21 16:53:31Z langlois $"

"""

#####################################################################################################################################

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm

from itertools import count, izip

########################################################

#Internal key word constants
_ALGO = 'FirstLastPivoter'

# If psyco package is available, use it to speed up processing (2x)
try:
    import psyco
    psyco.full()
except ImportError:
    pass

class Statistics(GenStatistics):
    """Maintains statistics for the current algorithm

    Attributes
        stat_names: Dictionary of statistics counters names maintained by Statistics class. These counters
                    are then incremented by calls made to the Statistics class

    """

    def __init__(self):
        """Initialize attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO, 'TOTAL_PRE', 'TOTAL_POST', 'PIVOTED', 'UNPIVOTED',))


    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and create a list of strings that builds the statistical message"

        Parameters:
            type: Always set to a summary of statistics

        Return value:
            A ready to be output string containing the statisticals informations about the current
            algorithm processing.

        """

        title = "%s Algorithm Statistics" % (_ALGO)
        str_out = []
        str_out.append("=" * len(title))
        str_out.append("%s Algorithm Statistics" % (_ALGO))
        str_out.append("-" * len(title))

        str_out.append("Features IN  : " + str(self.get_stats_name_count_total('TOTAL_PRE')))
        str_out.append("Features OUT : " + str(self.get_stats_name_count_total('TOTAL_POST')))
        str_out.append("-" * len(title))
        str_out.append("   PIVOTED Features   : " + str(self.get_stats_name_count_total('PIVOTED')))
        str_out.append("   UNPIVOTED Features : " + str(self.get_stats_name_count_total('UNPIVOTED')))
        str_out.append("=" * len(title))

        return str_out


class AlgoFirstLastPivoter(Algorithm):
    """
    This is the main class for the algorithm.

    Attributes:
        - params: Parameters for the algorithm
    """

    def __init__(self,
                 keep_iter_results=False,
                 debug=False):
        """Initialize the attributes of an object of this class

        Parameters:
            keep_iter_results: Flag to enable(True)/disable(False) keeping of the iterative results
            debug: Flag to enable(True)/disable(False) debug output

        Return value:
            None

        """

        Algorithm.__init__(self)

        self.params = Parameters()
        self.params.keep_iter_results=keep_iter_results
        self.params.debug=debug

        self.stats = Statistics()

    def _set_line_attributes (self):
        """This method sets attributes of line features in the spatial container.

        The routine checks:
            - if first/last vertices are the same : if so, feature is closed otherwise it is open and should not be pivoted
            - if feature has less than 3 vertices it is at its simplest form and is not pivoted

        Parameters:
            None

        Return value:
            None

        """

        # Select all and only lines
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):

            # increment counter of line features loaded in the spatial container
            self.stats.add_stats('TOTAL_PRE')

            # initialize line feature attributes
            line.is_simplest = False
            line.is_closed = False
            line.is_pivoted = False

            # Verify line validity to be pivoted
            if (len(line.coords_dual) >= 3):
                # verify line closure
                if (GenUtil.distance(line.coords_dual[0], line.coords_dual[-1]) <= GenUtil.ZERO):
                    line.is_closed = True
                else:
                    line.is_simplest = True
            else:
                line.is_simplest = True

            # increment counter of unpivoted line features
            if (line.is_simplest):
                self.stats.add_stats('UNPIVOTED')

        return


    def Pivot_this_feature(self, line, userFunction='((f[i-1]*0.5)+f[i]+(f[i+1]*0.5))'):
        """
        This method actually implements algorithm for FirstLastPivoter. This is actually where, if valid, the line feature
        is pivoted using the userFunction. Line feature coordinates are either unchanged or pivoted in place. In both
        cases, line feature acquired the attribute *is_pivoted* indicating if line has been pivoted (True) or no (False).

        Parameter
            line : line feature to pivot
            userFunction : defines the actual formula used to select pivot vertex. This formula may be set by the developer. Variable
                           *i* refers to the current vertex while variable *f* refers to the flat angle result for the current vertex. It
                           is also possible to use variable *a* to instead refer to the angle result of the current vertex.

        Return value:
            None

        """

        uFunction = userFunction

        # Extract feature coordinates and copy in new list structure to work with
        featureCoords = line.coords_dual
        newCoords = list(featureCoords)

        coordsCount = len(featureCoords)
        nbVertex = coordsCount - 1

        angles = []
        angles_flats = []
        userResults = []

        for i in range(0, nbVertex):

            # p1 = sommet précédent au sommet courant
            # p2 = sommet courant dont on veut connaître l'angle
            # p3 = sommet suivant
            # Lorsque p2 est le sommet 0, le sommet précédent est l'avant-avant-dernier sommet de la ligne.
            # Ici, la fonction modulo (%) permet de circuler sur les sommets de la ligne tout en évitant le dernier.
            # Ainsi, si la ligne a 5 sommets et si i=0 alors (nbVertex+(i-1))%nbVertex = (5+(0-1))%5 => 4
            #                                   si i=1 => (5+(1-1))%5 => 0
            #                                   si i=4 => (5+(4-1))%5 => 3
            p1 = featureCoords[(nbVertex + (i - 1)) % nbVertex]
            p2 = featureCoords[i]
            p3 = featureCoords[i + 1]

            # Calculate the angle between both vertex segments and store the result
            angle = GenUtil.compute_angle(p1, p2, p3)
            angles.append(angle)

            # Calculate the flat angle for the current vertex and store the result
            angle_flat = abs(abs(angle) - 180)
            angles_flats.append(angle_flat)

        # Transform i-x strings in user function to the appropriate evaluation string
        idx = 0
        while (idx < 20):
            idx += 1
            lookupStr = 'i-%s' % idx
            if (uFunction.find(lookupStr) != -1):
                uFunction = uFunction.replace(lookupStr, '(nbVertex+(' + lookupStr + '))%nbVertex')

        # Transform i+x strings in user function to the appropriate evaluation string
        idx = 0
        while (idx < 20):
            idx += 1
            lookupStr = 'i+%s' % idx
            if (uFunction.find(lookupStr) != -1):
                uFunction = uFunction.replace(lookupStr, '(nbVertex+' + lookupStr + ')%nbVertex')

        # Copy angles and angles_flats lists in order to be able to apply user function.
        # User uses 'a' and 'f' variables to refer to either angles or angles_flats lists.
        a = list(angles)
        f = list(angles_flats)

        for i in range(0, nbVertex):
            userResults.append(eval(uFunction))

        # Find, in a list, the smallest value and its index
        # source: http://www.gossamer-threads.com/lists/python/python/616403
        minUserFunction,minUserFunctionIndex = min(izip(userResults, count()))
        maxUserFunction,maxUserFunctionIndex = max(izip(userResults, count()))

        value = minUserFunction
        index = minUserFunctionIndex

        # Transfer coordinates into a new list while reordering them based on the pivot vertex
        borneInf = index
        borneSup = index + coordsCount
        for i in range(borneInf, borneSup):
          idx = i%nbVertex
          idxNv = i - index
          newCoords[idxNv] = featureCoords[idx]

        line.is_pivoted = (newCoords[0] != featureCoords[0])
        line.update_coords(newCoords, self.s_container)
        line.is_simplest = True

        return


    def process(self):
        """
        This method processes each line feature in turn. It also maintains general statistics about the overall process.
        This method ends when all line features have been processed.

        Parameter
            None

        Return value
            None

        """

        GenUtil.print_debug (self.params, "Start of %s algorithm)" % (_ALGO))
        GenUtil.print_debug (self.params, "Parameter description:")
        GenUtil.print_debug (self.params, "  - Keep iterative results: %s" % (self.params.keep_iter_results))

        # Load the shapely features into the spatial container
        self.s_container = self.load_features()

        if (self.params.debug):
            nbr_lines = len(self.s_container.get_features(filter="feature.feature_type=GenUtil.LINE_STRING"))
            nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT"))
            GenUtil.print_debug (self.params, "Number of lines imported: %s" % (nbr_lines))
            GenUtil.print_debug (self.params, "Number of points imported: %s" % (nbr_points))

        self.stats.add_iteration()
        self._set_line_attributes()

        GenUtil.print_debug (self.params, 'Start of iteration # %s' % (self.stats.get_nbr_iteration()))

        # Iterate through all lines which are not at their simplest state
        line_simplified = True

        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING and not feature.is_simplest"):
            self.Pivot_this_feature(line)

            if (line.is_pivoted):
                self.stats.add_stats('PIVOTED')
            else:
                self.stats.add_stats('UNPIVOTED')

        nb_features = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"))
        self.stats.add_stats('TOTAL_POST', nb_features)

        GenUtil.print_debug (self.params, 'End of iteration # %s' % (self.stats.get_nbr_iteration()))

        self.features = [feature for feature in self.s_container.get_features()]
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))
        #self.stats.print_stats()

