
"""
Phase plot
"""
import numpy as np
import matplotlib.pyplot as plt
HOURS_TO_RADIANS = 2*np.pi/24


class PhaseDiagram(object):
    """
    Make phase plot. 
    """

    def __init__(self, refstd,
                 fig=None, rect=111, label='_', y_lim=(0, 1.5), radial_label_pos = 0):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.

        Parameters:

        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        """
        self.ax = plt.subplot(111, polar=True)
        
        # suppress the radial labels
#         plt.setp(self.ax.get_yticklabels(), visible=False)

        # set the circumference labels
        self.ax.set_xticks(np.linspace(0, 2*np.pi, 24, endpoint=False))
        self.ax.set_xticklabels(range(24))
#         self.ax.xaxis.grid(True,color='r',linestyle='-')
        
        
        self.ax.set_ylim(y_lim)
        self.ax.yaxis.grid(True,color='r',linestyle='--')

        # make the labels go clockwise
        self.ax.set_theta_direction(-1)

        # place 0 at the top
        self.ax.set_theta_offset(np.pi/2.0)    
        
        self.ax.set_rlabel_position(radial_label_pos)
        
        self.ax.grid()
        

        self.samplePoints = []
        
        
    def get_ax(self):
        return self.ax
    
    def add_sample(self, phase, ampl, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """

#         l, = self.ax.plot(NP.arccos(corrcoef), stddev,
#                           *args, **kwargs)  # (theta, radius)
        l, = self.ax.plot(HOURS_TO_RADIANS*phase, ampl,
                          *args, **kwargs)
        self.samplePoints.append(l)

        return l
    
    def add_text(self, phase, ampl, text, **kwargs):
        
        self.ax.annotate(text, (HOURS_TO_RADIANS*phase, ampl), **kwargs)
        return None

    def add_grid(self, *args, **kwargs):
        """Add a grid."""

        self.ax.grid(*args, **kwargs)


