# Authors: M. Jacobi, J. Kuske, M. Reichert
################################################################################
import os
import numpy              as np
import matplotlib.pyplot  as plt
import matplotlib         as mpl
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects
from matplotlib.collections import PatchCollection
from tqdm                   import tqdm
from matplotlib             import cm
from matplotlib.patches     import Arrow, FancyBboxPatch
from matplotlib.colors      import LogNorm, SymLogNorm
from matplotlib.colors      import ListedColormap, LinearSegmentedColormap
from matplotlib.animation   import FuncAnimation
from matplotlib.widgets     import Slider, Button
from wreader                import wreader
from h5py                   import File
from nucleus_multiple_class import nucleus_multiple
from ngamma_eq              import ngamma_eq
from winvn_class            import winvn

################################################################################

class FlowAnimation(object):
    """
       Class to create an animation of the abundances and their flow for the nuclear reaction network WinNet.
    """

    def __init__(
        self,
        path,
        fig,
        # Folder to save the frames in
        frame_dir         = None,
        # Flows
        plot_flow         = True,
        # flow_group       = 'flows',
        flow_min          = 1e-8,
        flow_max          = 1e1,
        separate_fission  = True,
        fission_minflow   = 1e-10,
        flow_cbar         = True,
        flow_adapt_prange = True,
        flow_prange       = 5,
        flow_adapt_width  = True,
        flow_maxArrowWidth= 2.0,
        flow_minArrowWidth= 0.3,
        cmapNameFlow      = 'viridis',
        # Mass fractions
        abun_cbar        = True,
        X_min            = 1e-8,
        X_max            = 1,
        cmapNameX        = 'inferno',
        # Mass bins
        plot_abar        = True,
        plotMassBins     = True,
        addMassBinLabels = True,
        alphaMassBins    = 0.5,
        cmapNameMassBins = 'jet',
        massBins         = [[1,1],[2,71],[72,93],[94,110],[111,144],[145,169],[170,187],[188,205],[206,252],[253,337]],
        massBinLabels    = ['','','1st peak','','2nd peak','rare earths', '', '3rd peak', '', 'fissioning'],
        # Magic numbers
        plot_magic       = True,
        # Additional plots
        additional_plot  = 'none',
        # Tracked nuclei
        trackedrange     = (1e-8, 1),
        # Additional mainout
        amainoutrange    = (5e-10, 1),
        # Energy
        energyrange      = (1e10, 1e20),
        # Timescales
        timescalerange   = (1e-12, 1e10),
        timerange        = (1e-5 , 1e5),
        # Mainout
        plot_mainout     = True,
        densityrange     = (1e-5, 1e12),
        temperaturerange = (0, 10),
        yerange          = (0.0, 0.55),
        plot_logo        = True,
        indicate_r_path  = False,
        winvn_path       = None,
        interactive      = False
       ):
        """
        Parameters
        ----------
        path : str
            Path to the WinNet data.
        fig : matplotlib.figure.Figure
            Figure to plot the animation on.
        frame_dir : str
            Folder to save the frames in, default is path/frames.
        plot_flow : bool
            Plot the flow of the abundances.
        flow_min : float
            Minimum value for the flow.
        flow_max : float
            Maximum value for the flow.
        separate_fission : bool
            Separate the fissioning region from the fission products in the flow plot.
        fission_minflow : float
            Minimum value to indicate a fission region.
        flow_cbar : bool
            Plot the colorbar for the flow.
        flow_adapt_prange : bool
            Adapt the color range of the flow to the data.
        flow_prange : float
            Range (in log10) of the flow.
        flow_adapt_width : bool
            Adapt the width of the flow arrows to the flow.
        flow_maxArrowWidth : float
            Maximum width of the flow arrows.
        flow_minArrowWidth : float
            Minimum width of the flow arrows.
        cmapNameFlow : str
            Name of the colormap for the flow.
        abun_cbar : bool
            Plot the colorbar for the mass fractions.
        X_min : float
            Minimum value for the mass fractions.
        X_max : float
            Maximum value for the mass fractions.
        cmapNameX : str
            Name of the colormap for the mass fractions.
        plot_abar : bool
            Plot the average mass number.
        plotMassBins : bool
            Plot the mass bins.
        addMassBinLabels : bool
            Add labels to the mass bins.
        alphaMassBins : float
            Transparency of the mass bins.
        cmapNameMassBins : str
            Name of the colormap for the mass bin color.
        massBins : list
            List of mass bins.
        massBinLabels : list
            List of labels for the mass bins.
        plot_magic : bool
            Plot the magic numbers in the abundance plot.
        additional_plot : str
            Additional plot to be made, possible values: 'None', 'timescales', 'energy', 'tracked'.
        trackedrange : tuple
            Range of the tracked nuclei plot.
        amainoutrange : tuple
            Range of the additional mainout plot.
        energyrange : tuple
            Range of the energy axis.
        timescalerange : tuple
            Range of the timescales.
        timerange : tuple
            Range of the time axis.
        plot_mainout : bool
            Plot the mainout data.
        densityrange : tuple
            Range of the density axis in the mainout plot.
        temperaturerange : tuple
            Range of the temperature axis in the mainout plot.
        yerange : tuple
            Range of the electron fraction axis in the mainout plot.
        plot_logo : bool
            Plot the WinNet logo.
        indicate_r_path : bool
            Indicate the r-process path.
        winvn_path : str
            Path to the winvn file in case r-process path should be indicated.
        interactive : bool
            Enable interactive mode.
        """


        # Make some sanity check and ensure that the user has the correct version of Matplotlib
        if (mpl.__version__ < '3.8.0'):
            print('Using old version of Matplotlib ('+str(mpl.__version__)+'), some features may not work.')
            print('Need 3.8 or higher.')

        # Set the paths
        # Data directory:
        # Remember where this file is located
        self.__script_path = os.path.dirname(os.path.abspath(__file__))
        self.__data_path = os.path.join(self.__script_path,"data")
        # Frame directory:
        if frame_dir is None:
            self.frame_dir = f'{path}/frames'
        else:
            self.frame_dir = frame_dir

        # WinNet run path
        self.path = path

        # Save the parameters in class variables
        self.X_min            = X_min             # Minimum value for the mass fractions
        self.X_max            = X_max             # Maximum value for the mass fractions
        self.plotMassBins     = plotMassBins      # Plot the mass bins
        self.plot_abar        = plot_abar         # Plot the average mass number
        self.addMassBinLabels = addMassBinLabels  # Add labels to the mass bins
        self.alphaMassBins    = alphaMassBins     # Transparency of the mass bins
        self.cmapNameMassBins = cmapNameMassBins  # Name of the colormap for the mass bin color
        self.massBins         = massBins          # List of mass bins
        self.massBinLabels    = massBinLabels     # List of labels for the mass bins
        self.additional_plot  = additional_plot.lower() # Additional plot to be made, possible values: 'None', 'timescales', 'energy', 'tracked'
        self.energyrange      = energyrange       # Range of the energy axis
        self.trackedrange     = trackedrange      # Range of the tracked nuclei plot
        self.timescalerange   = timescalerange    # Range of the timescales
        self.timerange        = timerange         # Range of the time axis
        self.plot_mainout     = plot_mainout      # Plot the mainout data
        self.plot_logo        = plot_logo         # Plot the WinNet logo
        self.interactive      = interactive       # Have it interactive
        self.flow_group       = 'flows'           # Name of the group in the HDF5 file that contains the flow data
        self.plot_flow        = plot_flow         # Plot the flow of the abundances
        self.fig              = fig               # Figure to plot the animation on
        self.separate_fission = separate_fission  # Separate the fissioning region from the fission products in the flow plot
        self.plot_magic       = plot_magic        # Plot the magic numbers in the abundance plot
        self.densityrange     = densityrange      # Range of the density axis in the mainout plot
        self.temperaturerange = temperaturerange  # Range of the temperature axis in the mainout plot
        self.yerange          = yerange           # Range of the electron fraction axis in the mainout plot
        self.cmapNameX        = cmapNameX         # Name of the colormap for the mass fractions
        self.cmapNameFlow     = cmapNameFlow      # Name of the colormap for the flow
        self.flow_adapt_width = flow_adapt_width  # Adapt the width of the flow arrows to the flow
        self.flow_maxArrowWidth = flow_maxArrowWidth # Maximum width of the flow arrows
        self.flow_minArrowWidth = flow_minArrowWidth # Minimum width of the flow arrows
        self.flow_min           = flow_min          # Minimum value for the flow
        self.flow_max           = flow_max          # Maximum value for the flow
        self.flow_adapt_prange  = flow_adapt_prange # Adapt the color range of the flow to the data
        self.fission_minflow    = fission_minflow   # Minimum value for the fission flow
        self.amainoutrange      = amainoutrange     # Range of the additional mainout plot
        self.indicate_r_path    = indicate_r_path   # Indicate the r-process path
        self.winvn_path         = winvn_path        # Path to the winvn file in case r-process path should be indicated

        if (self.flow_adapt_prange):
            self.flow_prange = flow_prange       # Range (in log10) of the flow
        else:
            self.flow_prange = np.log10(self.flow_max) - np.log10(self.flow_min)

        # Check which additional plot should be made
        if self.additional_plot == 'timescales':
            self.plot_timescales = True
            self.plot_energy = False
            self.plot_tracked = False
            self.plot_addmainout = False
        elif self.additional_plot == 'energy':
            self.plot_timescales = False
            self.plot_energy = True
            self.plot_tracked = False
            self.plot_addmainout = False
        elif self.additional_plot == 'tracked':
            self.plot_timescales = False
            self.plot_energy = False
            self.plot_tracked = True
            self.plot_addmainout = False
        elif self.additional_plot == 'mainout':
            self.plot_timescales = False
            self.plot_energy = False
            self.plot_tracked = False
            self.plot_addmainout = True
        elif self.additional_plot == 'none':
            self.plot_timescales = False
            self.plot_energy = False
            self.plot_tracked = False
            self.plot_addmainout = False
        else:
            raise ValueError(f"Additional plot {self.additional_plot} not recognized. Possible values: 'None', 'timescales', 'energy', 'tracked'")

        # Initialize a class to read WinNet data
        self.wreader = wreader(path)

        # Read the stable isotopes
        self.N_stab, self.Z_stab  = np.loadtxt(os.path.join(self.__data_path,'../../../class_files/data/stableiso.dat'),
                                                unpack=True, usecols=(1, 2), dtype=int)

        # Read the nuclear chart nuclei
        sunet_path = os.path.join(self.__data_path, "sunet_really_complete")
        nuclei_names = np.loadtxt(sunet_path,dtype=str)
        nm = nucleus_multiple(names=nuclei_names)
        self.__A_plot = nm.A
        self.__Z_plot = nm.Z
        self.__N_plot = nm.N
        self.__min_N, self.__max_N = np.min(self.__N_plot), np.max(self.__N_plot)
        self.__min_Z, self.__max_Z = np.min(self.__Z_plot), np.max(self.__Z_plot)


        # Other data that does not change like magic numbers
        self.nMagic = [8, 20, 28, 50, 82, 126, 184] # Magic numbers in neutrons
        self.zMagic = [8, 20, 28, 50, 82, 114]      # Magic numbers in protons
        self.magic_excess = 4                       # How long should the line stick out over the nuc. chart?

        # Timescale plotting parameters
        self.timescale_colors = ["C1","C2","C3","C4","C5","C6","C7","C8","C9","C0"]
        self.timescale_entries = [['ag','ga'],['ng','gn'],['an','na'],['np','pn'],['pg','gp'],['ap','pa'],["beta"],["bfiss"],["nfiss"],["sfiss"]]
        self.timescale_labels = [r"$\tau_{\alpha,\gamma}$",r"$\tau_{n,\gamma}$",r"$\tau_{\alpha,n}$",r"$\tau_{n,p}$",r"$\tau_{p,\gamma}$",
                                 r"$\tau_{\alpha,p}$",r"$\tau_{\beta}$",r"$\tau_{\rm{bfiss}}$",r"$\tau_{\rm{nfiss}}$",r"$\tau_{\rm{sfiss}}$"]

        # Energy plotting parameters
        self.energy_colors = ["k","C2","C3","C4","C5","C6","C7","C8","C9","C0"]
        self.energy_entries = ['tot', 'ag_ga', 'ng_gn', 'an_na', 'np_pn', 'pg_gp', 'ap_pa', 'beta', 'fiss']
        self.energy_labels = ['Total', r"$( \alpha,\gamma )$", r"$( n,\gamma )$", r"$( \alpha,n )$", r"$( n,p )$",
                              r"$( p,\gamma )$", r"$( \alpha,p )$", r"$\beta$", r"$\rm{fission}$"]
        self.energy_lw     = np.ones(len(self.energy_entries))
        self.energy_lw[0]  = 2

        # Set up the norm of the flow
        self.flow_norm = LogNorm(flow_min, flow_max, clip=True)
        # In case of adaptive flow ranges, keep track of the maximums to adapt the ranges with a rolling average
        self.flow_max_history = np.ones(5)*np.nan

        # If the flow is not plotted, then the fissioning region should not be separated
        # and the flow colorbar should not be plotted
        if (not self.plot_flow):
            self.separate_fission = False
            flow_cbar = False

        # Initialize the axes and figure
        self.init_axes()
        # Initialize the data
        self.init_data()
        # Initialize the plot
        self.init_plot()
        # Initialize the colorbars
        self.init_cbars(abun_cbar, flow_cbar)

        # Set up the limits of the plot
        self.limits_plot = self.ax.get_xlim(), self.ax.get_ylim()

        # For interactive flow range
        self.flow_max_offset = 0.0
        self.flow_min_offset = 0.0

        if self.indicate_r_path:
            self.__init_ngamma_eq()

        if self.interactive:
            self.__init_sunet_indicator()
            self.__interactive_box = None
            self.__interactive_textbox = None
            self.winvn = winvn(self.winvn_path)
            self.winvn.read_winvn()
            df = self.winvn.get_dataframe()
            # Set a tuple of N and Z as index
            df.set_index(['N','Z'], inplace=True)
            self.winvn.set_dataframe(df)


    def __init_ngamma_eq(self):
        """
           Initialize the ngamma_eq class.
        """
        self.ngamma_eq = ngamma_eq(self.winvn_path)
        self.ngamma_eq_plot   = self.ax.plot([],[],color='purple',lw=1.5,zorder=99)
        self.ngamma_eq_plot_o = self.ax.plot([],[],color='w',lw=2.5,zorder=98)



    def __init_sunet_indicator(self):
        sunet_path = self.wreader.template['net_source']
        nuclei_names = np.loadtxt(sunet_path,dtype=str)
        nm = nucleus_multiple(names=nuclei_names)

        self.__sunet_lines = []
        # Loop through Zs
        for Z in np.unique(nm.Z):
            # Find the Ns that are have larger than distance one to the next N
            mask = (nm.Z == Z)
            Ns = nm.N[mask]
            # Find the Ns that are have larger than distance one to the next N
            diff = np.diff(Ns)
            mask = np.where(diff > 1)[0]
            # Loop through the Ns
            # Plot left and right
            line = self.ax.plot([np.min(Ns)-0.5, np.min(Ns)-0.5], [Z-0.5, Z+0.5], color='red', zorder=1000, lw=1)
            self.__sunet_lines.append(line)
            line = self.ax.plot([np.max(Ns)+0.5, np.max(Ns)+0.5], [Z-0.5, Z+0.5], color='red', zorder=1000, lw=1)
            self.__sunet_lines.append(line)
            for i in mask:
                # Add a line to the plot
                line = self.ax.plot([Ns[i]+0.5, Ns[i]+0.5], [Z-0.5, Z+0.5], color='red', zorder=1000, lw=1)
                self.__sunet_lines.append(line)
                line = self.ax.plot([Ns[i+1]-0.5, Ns[i+1]-0.5], [Z-0.5, Z+0.5], color='red', zorder=1000, lw=1)
                self.__sunet_lines.append(line)

        # Same but for Z
        for N in np.unique(nm.N):
            # Find the Ns that are have larger than distance one to the next N
            mask = (nm.N == N)
            Zs = nm.Z[mask]
            # Find the Ns that are have larger than distance one to the next N
            diff = np.diff(Zs)
            mask = np.where(diff > 1)[0]
            # Loop through the Ns
            # Plot left and right
            line = self.ax.plot([N-0.5, N+0.5], [np.min(Zs)-0.5, np.min(Zs)-0.5], color='red', zorder=1000, lw=1)
            self.__sunet_lines.append(line)
            line = self.ax.plot([N-0.5, N+0.5], [np.max(Zs)+0.5, np.max(Zs)+0.5], color='red', zorder=1000, lw=1)
            self.__sunet_lines.append(line)
            for i in mask:
                # Add a line to the plot
                line = self.ax.plot([N-0.5, N+0.5], [Zs[i]+0.5, Zs[i]+0.5], color='red', zorder=1000, lw=1)
                self.__sunet_lines.append(line)
                line = self.ax.plot([N-0.5, N+0.5], [Zs[i+1]-0.5, Zs[i+1]-0.5], color='red', zorder=1000, lw=1)
                self.__sunet_lines.append(line)

        # Hide the lines
        for line in self.__sunet_lines:
            for l in line:
                l.set_visible(False)
        self.sunet_indication = False




    def init_axes(self):
        """
           Initialize the axes and everything figure related of the plot.
        """
        # Make the hatch linewidth smaller
        mpl.rcParams['hatch.linewidth'] = 0.8

        # Set up the figure
        self.ax = self.fig.gca()
        self.ax.set_aspect('equal')

        # Set up the axes for the nuclear chart
        self.__init_nucchart_ax()

        # Set up the axes for the mainout plot
        if self.plot_mainout:
            self.__init_axMainout()

        # Set up the axes for the mass fraction plot
        self.__init_axAbund()

        # Timescale stuff
        if self.plot_timescales:
            self.__init_axTimescales()

        # Energy stuff
        if self.plot_energy:
            self.__init_axEnergy()

        # Tracked nuclei
        if self.plot_tracked:
            self.__init_axTracked()

        # Additional mainout
        if self.plot_addmainout:
            self.__init_axAddMainout()

        # WinNet logo
        if self.plot_logo:
            self.__init_logo()

        # interactive stuff
        if self.interactive:
            self.__init_interactive()

        # Start with a running movie
        self.movie_paused = False

        # Make white behind the colorbars and plots
        self.ax.add_patch(patches.Rectangle((0.3, 0.84), 0.8, 0.15, fill=True, color='w', zorder=100, transform=self.fig.transFigure))
        self.ax.add_patch(patches.Rectangle((0.15, 0.745), 0.395, 0.2, fill=True, color='w', zorder=100, transform=self.fig.transFigure))



    def __init_nucchart_ax(self):
        """
           Initialize the axes and everything figure related of the nuclear chart.
        """
        # Set up the main figure of the nuclear chart, i.e., remove borders, ticks, etc.
        self.ax.xaxis.set_visible(False)
        self.ax.yaxis.set_visible(False)
        self.ax.spines[['right', 'top', "bottom", "left"]].set_visible(False)


    def __init_axAbund(self):
        """
           Initialize the axes and everything figure related of the mass fraction plot.
        """
        # Add summed mass fraction plot
        self.axAbund = plt.axes([0.15,0.78,0.35,0.15])
        self.axAbund.set_xlabel(r'Mass number $A$')
        self.axAbund.set_ylabel(r'X(A)')
        self.axAbund.set_xlim(0,300)
        self.axAbund.set_yscale('log')
        self.axAbund.set_ylim(self.X_min, self.X_max)
        if self.plot_abar:
            self.axAbund.axvline(np.nan, color='tab:red',label=r"$\bar{A}$")
            self.axAbund.legend(loc='upper right')
        # Add mass fraction bins if needed
        if self.plotMassBins:
            self.axMassFrac = self.axAbund.twinx()
            self.axMassFrac.set_ylabel(r'$X(A)$')
            self.axMassFrac.set_ylim(0,1)
            self.axMassFrac.set_yscale('linear')


    def __init_axTimescales(self):
        """
           Initialize the axes and everything figure related of the timescales plot.
        """
        self.axTimescales = plt.axes([0.15,0.55,0.20,0.17])
        self.axTimescales.set_xlabel('Time [s]')
        self.axTimescales.set_ylabel('Timescales [s]')
        self.axTimescales.set_yscale('log')
        self.axTimescales.set_ylim(self.timescalerange[0],self.timescalerange[1])
        self.axTimescales.set_xscale('log')
        self.axTimescales.set_xlim(self.timerange[0],self.timerange[1])

    def __init_axTracked(self):
        """
           Initialize the axes and everything figure related of the tracked nuclei plot.
        """
        self.axTracked = plt.axes([0.15,0.55,0.20,0.17])
        self.axTracked.set_xlabel('Time [s]')
        self.axTracked.set_ylabel('Mass fractions')
        self.axTracked.set_yscale('log')
        self.axTracked.set_ylim(self.trackedrange[0],self.trackedrange[1])
        self.axTracked.set_xscale('log')
        self.axTracked.set_xlim(self.timerange[0],self.timerange[1])

    def __init_axEnergy(self):
        """
           Initialize the axes and everything figure related of the energy plot.
        """
        self.axEnergy = plt.axes([0.15,0.55,0.20,0.17])
        self.axEnergy.set_xlabel('Time [s]')
        self.axEnergy.set_ylabel('Energy [erg/g/s]')
        self.axEnergy.set_yscale('log')
        self.axEnergy.set_ylim(self.energyrange[0],self.energyrange[1])
        self.axEnergy.set_xscale('log')
        self.axEnergy.set_xlim(self.timerange[0],self.timerange[1])

    def __init_axAddMainout(self):
        """
           Initialize the axes and everything figure related of the energy plot.
        """
        self.axAddMainout = plt.axes([0.15,0.55,0.20,0.17])
        self.axAddMainout.set_xlabel('Time [s]')
        self.axAddMainout.set_ylabel('Abundance')
        self.axAddMainout.set_yscale('log')
        self.axAddMainout.set_ylim(self.amainoutrange[0],self.amainoutrange[1])
        self.axAddMainout.set_xscale('log')
        self.axAddMainout.set_xlim(self.timerange[0],self.timerange[1])

    def __init_axMainout(self):
        """
           Initialize the axes and everything figure related of the mainout plot.
        """
        def make_patch_spines_invisible(ax):
            ax.set_frame_on(True)
            ax.patch.set_visible(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
        # Density
        self.axMainout = plt.axes([0.65,0.2,0.25,0.25])
        self.axMainout.set_xlabel('Time [s]')
        self.axMainout.set_ylabel(r'Density [g/cm$^3$]')
        self.axMainout.set_yscale('log')
        self.axMainout.set_ylim(self.densityrange[0],self.densityrange[1])
        self.axMainout.set_xscale('log')
        self.axMainout.set_xlim(self.timerange[0],self.timerange[1])
        # Temperature
        self.axMainout_temp = self.axMainout.twinx()
        self.axMainout_temp.set_ylabel(r'Temperature [GK]')
        self.axMainout_temp.set_ylim(self.temperaturerange[0],self.temperaturerange[1])
        self.axMainout_temp.spines["left"].set_position(("axes", -0.2))
        make_patch_spines_invisible(self.axMainout_temp)
        self.axMainout_temp.spines["left"].set_visible(True)
        self.axMainout_temp.yaxis.set_label_position('left')
        self.axMainout_temp.yaxis.set_ticks_position('left')
        self.axMainout_temp.yaxis.label.set_color("tab:red")
        self.axMainout_temp.yaxis.set_tick_params(colors="tab:red")
        # Ye
        self.axMainout_ye = self.axMainout.twinx()
        self.axMainout_ye.set_ylabel(r'Electron fraction')
        self.axMainout_ye.set_ylim(self.yerange[0],self.yerange[1])
        make_patch_spines_invisible(self.axMainout_ye)
        self.axMainout_ye.spines["right"].set_visible(True)
        self.axMainout_ye.yaxis.set_label_position('right')
        self.axMainout_ye.yaxis.set_ticks_position('right')
        self.axMainout_ye.yaxis.label.set_color("tab:blue")
        self.axMainout_ye.yaxis.set_tick_params(colors="tab:blue")


    def __init_logo(self):
        """
           Initialize the axes and everything figure related of the WinNet logo.
        """
        self.axLogo = plt.axes([0.75,0.45,0.15,0.15])
        self.axLogo.axis('off')
        self.axLogo.imshow(plt.imread(os.path.join(self.__data_path,'WinNet_logo.png')))

    def __init_interactive(self):

        self.ax_slider = plt.axes([0.18, 0.08, 0.72, 0.02], facecolor='lightgoldenrodyellow')
        self.ax_button = plt.axes([0.18-0.018, 0.08, 0.012, 0.022])  # Adjust `bottom` for centering

        self.play_button = Button(self.ax_button, "❚❚")
        self.play_button.label.set_color("red")
        self.play_button.label.set_fontsize(8)
        self.play_button.on_clicked(self.pause_movie)
        self.fig.canvas.mpl_connect('key_press_event', self.arrow_update)
        self.__interactive_ax = None

        # Add a bookmark at a certain time in the slider
        # Calculate neutron freeze-out time
        min_val = 1-self.wreader.mainout['yn']/self.wreader.mainout['yheavy']
        nfreezeout = np.argmin(abs(min_val))
        # Check if its really the freeze-out time by checking if it switches from larger one to smaller 1
        if (min_val[nfreezeout-1] < 0) and (min_val[nfreezeout+1] > 0):
            # self.ax_slider.plot([nfreezeout,nfreezeout],[0,0.25], color='tab:red', linestyle='-', linewidth=1)  # Bookmark indicators
            self.ax_slider.axvline(nfreezeout, color='tab:red', linestyle='-', linewidth=1)  # Bookmark indicators

        # Check which data could be shown
        self.available_data = []
        self.toggle_buttons = []
        if self.wreader.check_existence('tracked_nuclei') !=0:
            self.available_data.append('tracked_nuclei')
            self.toggle_buttons.append(Button(plt.axes([0.18-0.018, 0.05, 0.012, 0.022]), "⚛"))
            # Fontsize
            self.toggle_buttons[-1].label.set_fontsize(12)
            self.toggle_buttons[-1].label.set_color('k')
            self.toggle_buttons[-1].on_clicked(self.toggle_button_event)

        if self.wreader.check_existence('timescales') !=0:
            self.available_data.append('timescales')
            amount_buttons = len(self.toggle_buttons)
            self.toggle_buttons.append(Button(plt.axes([0.18-0.018+0.015*amount_buttons, 0.05, 0.012, 0.022]), r"$\tau$"))
            self.toggle_buttons[-1].label.set_fontsize(12)
            self.toggle_buttons[-1].label.set_color('tab:green')
            self.toggle_buttons[-1].on_clicked(self.toggle_button_event)

        if self.wreader.check_existence('energy') !=0:
            self.available_data.append('energy')
            amount_buttons = len(self.toggle_buttons)
            self.toggle_buttons.append(Button(plt.axes([0.18-0.018+0.015*amount_buttons, 0.05, 0.012, 0.022]), "⚡"))
            self.toggle_buttons[-1].label.set_fontsize(12)
            self.toggle_buttons[-1].label.set_color('tab:orange')
            self.toggle_buttons[-1].on_clicked(self.toggle_button_event)

        if self.wreader.check_existence('mainout') !=0:
            self.available_data.append('mainout')
            amount_buttons = len(self.toggle_buttons)
            self.toggle_buttons.append(Button(plt.axes([0.18-0.018+0.015*amount_buttons, 0.05, 0.012, 0.022]), "m"))
            self.toggle_buttons[-1].label.set_fontsize(12)
            self.toggle_buttons[-1].label.set_color('tab:blue')
            self.toggle_buttons[-1].on_clicked(self.toggle_button_event)
        self.active_button = None

        # Add zoom button
        self.zoom_button = Button(plt.axes([0.18-0.018+0.015*(len(self.toggle_buttons)+0.2), 0.05, 0.012, 0.022]), "+")
        self.zoom_button.label.set_fontsize(12)
        self.zoom_button.label.set_color('k')
        self.zoom_button.on_clicked(self.zoom_button_event)
        self.zoomed = False

        # Check if the sunet path exists
        supath = self.wreader.template['net_source']
        addbutton = 0
        if os.path.exists(supath):
            self.sunet_button = Button(plt.axes([0.18-0.018+0.015*(len(self.toggle_buttons)+1+0.2), 0.05, 0.012, 0.022]), "S")
            self.sunet_button.label.set_fontsize(12)
            self.sunet_button.label.set_color('k')
            self.sunet_button.on_clicked(self.sunet_button_event)
            addbutton += 1
        if self.wreader.check_existence('mainout') !=0:
            if not self.indicate_r_path:
                self.__init_ngamma_eq()
                self.ngamma_eq_plot[0].set_visible(self.indicate_r_path)
                self.ngamma_eq_plot_o[0].set_visible(self.indicate_r_path)
            self.r_path_button = Button(plt.axes([0.18-0.018+0.015*(len(self.toggle_buttons)+1+addbutton+0.2), 0.05, 0.012, 0.022]), "r")
            self.r_path_button.label.set_fontsize(12)
            self.r_path_button.label.set_color('k')
            self.r_path_button.on_clicked(self.r_path_button_event)
            addbutton += 1

        # Check if flow is plotted and add a button to change the flow range
        if self.plot_flow:
            self.flow_buttons = [Button(plt.axes([0.88, 0.905, 0.01, 0.02]), "+")]
            self.flow_buttons[-1].label.set_fontsize(12)
            self.flow_buttons[-1].label.set_color('k')
            self.flow_buttons[-1].on_clicked(self.flow_button_event)

            self.flow_buttons.append(Button(plt.axes([0.868, 0.905, 0.01, 0.02]), "-"))
            self.flow_buttons[-1].label.set_fontsize(12)
            self.flow_buttons[-1].label.set_color('k')
            self.flow_buttons[-1].on_clicked(self.flow_button_event)

            self.flow_buttons.append(Button(plt.axes([0.762, 0.905, 0.01, 0.02]),"+"))
            self.flow_buttons[-1].label.set_fontsize(12)
            self.flow_buttons[-1].label.set_color('k')
            self.flow_buttons[-1].on_clicked(self.flow_button_event)

            self.flow_buttons.append(Button(plt.axes([0.75, 0.905, 0.01, 0.02]), "-"))
            self.flow_buttons[-1].label.set_fontsize(12)
            self.flow_buttons[-1].label.set_color('k')
            self.flow_buttons[-1].on_clicked(self.flow_button_event)

            # Add a button to change reset
            self.flow_buttons.append(Button(plt.axes([0.775, 0.905, 0.01, 0.02]), "⟲"))
            self.flow_buttons[-1].label.set_fontsize(12)
            self.flow_buttons[-1].label.set_color('k')
            self.flow_buttons[-1].on_clicked(self.flow_button_event)

            # self.flow_button.on_clicked(self.flow_button_event)

    def flow_button_event(self, event):
        if event.inaxes == self.flow_buttons[0].ax:
            self.flow_max_offset += 0.5
        elif event.inaxes == self.flow_buttons[1].ax:
            self.flow_max_offset -= 0.5
        elif event.inaxes == self.flow_buttons[2].ax:
            self.flow_min_offset -= 0.5
        elif event.inaxes == self.flow_buttons[3].ax:
            self.flow_min_offset += 0.5
        elif event.inaxes == self.flow_buttons[4].ax:
            self.flow_max_offset = 0.0
            self.flow_min_offset = 0.0

        # Refresh the animation at current position
        if self.movie_paused:
            self.update_frame(self.slider_bar.val)

        self.fig.canvas.draw_idle()
        pass

    def r_path_button_event(self, event):
        self.indicate_r_path = not self.indicate_r_path
        self.ngamma_eq_plot[0].set_visible(self.indicate_r_path)
        self.ngamma_eq_plot_o[0].set_visible(self.indicate_r_path)
        self.fig.canvas.draw_idle()

    def sunet_button_event(self, event):
        self.sunet_indication = not self.sunet_indication
        for line in self.__sunet_lines:
            for l in line:
                l.set_visible(self.sunet_indication)
        self.fig.canvas.draw_idle()


    def zoom_button_event(self, event):
        self.__toggle_zoom()
        if self.zoomed:
            self.zoom_button.label.set_text("-")
        else:
            self.zoom_button.label.set_text("+")

    def toggle_button_event(self, event):
        # Find the button that was clicked
        for i, button in enumerate(self.toggle_buttons):
            if event.inaxes == button.ax:
                # Toggle the state: activate this button and deactivate others
                if self.active_button != button:
                    # Unpress the previously active button if there is one
                    if self.active_button:
                        for t in ['top','right','bottom','left']:
                            self.active_button.ax.spines[t].set_color('k')
                            self.active_button.ax.spines[t].set_linewidth(0.5)
                    # Set the new active button
                    # Ändere die Umrandung des Buttons
                    for t in ['top','right','bottom','left']:
                        button.ax.spines[t].set_color('red')
                        button.ax.spines[t].set_linewidth(2)

                    self.active_button = button
                    if not (self.__interactive_ax is None):
                        self.__interactive_ax.remove()
                    if self.available_data[i] == 'timescales':
                        self.__toggle_timescales()
                    elif self.available_data[i] == 'energy':
                        self.__toggle_energy()
                    elif self.available_data[i] == 'tracked_nuclei':
                        self.__toggle_tracked()
                    elif self.available_data[i] == 'mainout':
                        self.__toggle_addmainout()
                else:
                    # Unpress the button
                    for t in ['top','right','bottom','left']:
                        button.ax.spines[t].set_color('k')
                        button.ax.spines[t].set_linewidth(0.5)
                    self.active_button = None
                    if not (self.__interactive_ax is None):
                        self.__interactive_ax.remove()
                    self.__interactive_ax = None
                    self.plot_timescales = False
                    self.plot_energy = False
                    self.plot_tracked = False
                    self.plot_addmainout = False
                self.fig.canvas.draw_idle()
                return


    def init_data(self):
        """
           Initialize the data for the plot.
        """
        # Read all things related to the snapshots
        self.Z = self.wreader.Z
        self.N = self.wreader.N
        self.A = self.wreader.A
        self.X = self.wreader.X[0,:]
        # Get the number of timesteps, i.e., the number of snapshots
        self.n_timesteps = self.wreader.nr_of_snaps

        # Set up timescale data
        if self.plot_timescales:
            self.__init_data_timescales(-1)

        # Set up energy data
        if self.plot_energy:
            self.__init_data_energy(-1)

        # Set up tracked nuclei data
        if self.plot_tracked:
            self.__init_data_tracked(-1)

        # Set up additional mainout data
        if self.plot_addmainout:
            self.__init_data_addmainout(-1)

        # Set up mainout data
        if self.plot_mainout:
            self.mainout_time        = self.wreader.mainout['time']
            self.mainout_density     = self.wreader.mainout['dens']
            self.mainout_temperature = self.wreader.mainout['temp']
            self.mainout_ye          = self.wreader.mainout['ye']

        if self.indicate_r_path:
            self.__init_ngamma_eq()

        self.time = 0

        # Set up the Abar
        self.Abar = 1.0/np.sum(self.X/self.A)
        # Set up the sum of the mass fractions
        self.Asum, self.Xsum = self.sum_over_A(self.A, self.X)
        # Set up the mass fraction bins
        self.Xbins = np.zeros(len(self.massBins),dtype=float)

        # Create custom colormap, with colormap for abundances and also the background colors
        abucmap = mpl.colormaps[self.cmapNameX]
        abundance_colors = abucmap(np.linspace(0, 1, 256))
        massbin_colormap = mpl.colormaps[self.cmapNameMassBins]
        amount_mass_bins = len(self.massBins)
        massbin_colors   = massbin_colormap(np.linspace(0, 1, amount_mass_bins))
        massbin_colors[:,3] = self.alphaMassBins
        self.massbin_colors = massbin_colors
        newcolors = np.vstack((massbin_colors, abundance_colors))
        self.__abundance_colors = ListedColormap(newcolors)
        # Create the values for the colors
        dist = (np.log10(self.X_max)-np.log10(self.X_min))/256.0
        self.values = np.linspace(np.log10(self.X_min)-amount_mass_bins*dist, np.log10(self.X_max),num=256+amount_mass_bins+1,endpoint=True)

        # Create the background array that will contain numbers according to the background colors
        background_Y = np.empty((self.__max_N+1, self.__max_Z+1))
        background_Y[:,:] = np.nan
        for index, mbin in enumerate(self.massBins):
            mask = (self.__A_plot >= self.massBins[index][0]) & (self.__A_plot <= self.massBins[index][1])
            background_Y[self.__N_plot[mask],self.__Z_plot[mask]] = self.values[index]
        self.__background_Y = background_Y

        # Set up the axes for the nuclear chart
        self.n = np.arange(0, self.__max_N+1)
        self.z = np.arange(0, self.__max_Z+1)

        # Set up the abundance array
        self.abun = self.__background_Y

        # Set up the array for the fission region
        self.fis_region = np.zeros_like(self.abun)

        # Set up the flow arrays if necessary
        if (self.plot_flow):
            self.flow_N, self.flow_Z = np.array([0]), np.array([0])
            self.flow_dn, self.flow_dz = np.array([0]), np.array([0])
            self.flow = np.array([0])


    def init_plot(self):
        """
           Initialize the plots.
        """
        # Plot the nuclear chart
        self.abun_im = self.ax.pcolormesh(self.n,self.z,self.abun.T,
                       cmap = self.__abundance_colors,vmin=(min(self.values)),
                       vmax=(max(self.values)),linewidth=0.0,edgecolor="face")


        if (self.plot_flow):
            # Plot the flows as quiver
            self.quiver = self.ax.quiver(
                self.flow_N, self.flow_Z,
                self.flow_dn, self.flow_dz,
                self.flow,
                norm=self.flow_norm,
                cmap=self.cmapNameFlow, angles='xy', scale_units='xy', scale=1,
                units='xy', width=0.1, headwidth=3, headlength=4
                )

            if self.flow_adapt_width:
                # Create patchcollection of arrows
                width = (self.flow_maxArrowWidth-self.flow_minArrowWidth)/self.flow_prange
                with np.errstate(divide='ignore'):
                    arrowwidth = (np.log10(self.flow)-np.log10(self.flow_min))*width
                    arrowwidth = np.maximum(arrowwidth, self.flow_minArrowWidth)


                flow_arrows = [Arrow(self.flow_N[i],self.flow_Z[i],self.flow_dn[i],self.flow_dz[i],width=arrowwidth[i],color='k') for i in range(len(self.flow))]
                a = PatchCollection(flow_arrows, cmap=self.cmapNameFlow, norm=self.flow_norm)
                a.set_array(self.flow)
                self.flow_patch = self.ax.add_collection(a)

            # Plot the fission region of positive fission products
            fisspos = self.fis_region.T
            fisspos[:] = np.nan
            self.fis_im_pos = self.ax.pcolor(self.n,self.z,
                fisspos,hatch='//////', edgecolor='tab:red',facecolor='none',
                linewidth=0.0,zorder=1000
                )

            # Plot the fission region of negative fission products
            fissneg = self.fis_region.T
            fissneg[:] = np.nan
            self.fis_im_neg = self.ax.pcolor(self.n,self.z,
                fissneg,hatch='//////', edgecolor='tab:blue',facecolor='none',
                linewidth=0.0,zorder=1000
                )


        # Plot stable isotopes as black rectangles
        edgecolors = np.full((max(self.n+1),max(self.z+1)),"none")
        edgecolors[self.N_stab, self.Z_stab] = "k"
        edgecolors = edgecolors.T.ravel()
        self.stable_im = self.ax.pcolormesh(self.n,self.z,self.abun.T,
                          facecolor="none",linewidth=0.5,edgecolor=edgecolors)
        # Alternatively make a scatter
        # self.stable_im = self.ax.scatter(
        #     N_stab, Z_stab, c='k', marker='o', s=1)

        # Plot the sum of the mass fractions
        self.mafra_plot = self.axAbund.plot(self.Asum, self.Xsum, color='k', lw=1.2)

        if self.plotMassBins:
            # Plot the mass bins
            self.mafrabin_plot = [self.axMassFrac.bar(self.massBins[i][0], self.Xbins[i],
                                  width=self.massBins[i][1]+1-self.massBins[i][0], align='edge',
                                  color=self.massbin_colors[i], edgecolor='grey') for i in range(len(self.massBins))]

            for i,b in enumerate(self.massBins):
                px = (b[0])+1
                py = 1
                axis_to_data = self.axAbund.transAxes + self.axAbund.transData.inverted()
                data_to_axis = axis_to_data.inverted()
                trans = data_to_axis.transform((px,py))
                txt = self.axAbund.text(trans[0],1.1,self.massBinLabels[i],transform = self.axAbund.transAxes,
                                  ha='left',va='center',clip_on=False, fontsize=8,color=self.massbin_colors[i])
                txt.set_path_effects([PathEffects.withStroke(linewidth=0.5, foreground='k')])

        # Plot the average mass number as red vertical line
        if self.plot_abar:
            self.Abar_plot = self.axAbund.axvline(self.Abar, color='tab:red', lw=1)

        # Plot the mainout data
        if self.plot_mainout:
            self.mainout_dens_plot = self.axMainout.plot     (np.nan, np.nan, color='k'       , lw=2, label=r"$\rho$")
            self.mainout_temp_plot = self.axMainout_temp.plot(np.nan, np.nan, color='tab:red' , lw=2, label=r"T$_9$")
            self.mainout_ye_plot   = self.axMainout_ye.plot  (np.nan, np.nan, color='tab:blue', lw=2, label=r"Y$_e$")
            # Create the legend
            lines = [self.mainout_dens_plot[0], self.mainout_temp_plot[0], self.mainout_ye_plot[0]]
            self.axMainout.legend(lines, [l.get_label() for l in lines],loc='upper right')

            # Set the Textbox for the mainout data
            left_side = "$t$"+"\n"+rf"$\rho$"+"\n"+rf"$T_9$"+"\n"+rf"$Y_e$"
            right_side = f"= {self.format_time(self.mainout_time[-1])[0]}\n"+\
                         f"= {self.to_latex_exponent(self.mainout_density[-1])}\n"+\
                         f"= {self.to_latex_exponent(self.mainout_temperature[-1])}\n"+\
                         f"= {self.mainout_ye[-1]:.3f}"
            units = f"{self.format_time(self.mainout_time[-1])[1]}\n"+\
                    f"g/cm$^3$\n"+\
                    "GK\n"+\
                    ""

            # Make a background box
            y_pos = 0.07
            rect = patches.FancyBboxPatch((0.35, y_pos-0.01), 0.17, 0.135, transform=self.ax.transAxes, boxstyle="round,pad=0.01", ec="k", fc="lightgrey", zorder=1,alpha=0.5)
            self.ax.add_patch(rect)

            # Plot the text
            self.ax.text(0.35, y_pos, left_side, transform=self.ax.transAxes, fontsize=12,)
            self.Mainout_text = self.ax.text(0.37, y_pos, right_side, transform=self.ax.transAxes, fontsize=12,)
            self.Mainout_units = self.ax.text(0.48, y_pos, units, transform=self.ax.transAxes, fontsize=12,)


        # Set the arrow at the bottom left
        arrowLength = 20
        self.ax.arrow(-8, -8, arrowLength, 0, head_width=2, head_length=2, fc='k', ec='k')
        self.ax.arrow(-8, -8, 0, arrowLength, head_width=2, head_length=2, fc='k', ec='k')
        self.ax.text(arrowLength-8+3,0-8,'N',horizontalalignment='left',verticalalignment='center',fontsize=14,clip_on=True)
        self.ax.text(0-8,arrowLength-8+3,'Z',horizontalalignment='center',verticalalignment='bottom',fontsize=14,clip_on=True)

        # Plot the timescales
        if self.plot_timescales:
            self.__init_plot_timescales()
        # Plot the energy
        if self.plot_energy:
            self.__init_plot_energy()
        # Plot the additional mainout
        if self.plot_addmainout:
            self.__init_plot_addmainout()
        # Plot the tracked nuclei
        if self.plot_tracked:
            self.__init_plot_tracked()

        # Plot magic numbers
        if self.plot_magic:
            # First neutron numbers
            for n in self.nMagic:
                if (any(self.__N_plot == n)):
                    # Find minimum and maximum Z for this N
                    zmin = np.min(self.__Z_plot[self.__N_plot == n])
                    zmax = np.max(self.__Z_plot[self.__N_plot == n])
                    # Plot horizontal line from zmin - self.magic_excess to zmax + self.magic_excess
                    self.ax.plot([n-0.5, n-0.5],[zmin-self.magic_excess, zmax+self.magic_excess], color='k', ls="--", lw=0.5)
                    self.ax.plot([n+0.5, n+0.5],[zmin-self.magic_excess, zmax+self.magic_excess], color='k', ls="--", lw=0.5)
                    # write the number on the bottom, truncate text when out of range
                    self.ax.text(n,zmin-self.magic_excess-1,int(n),ha='center',va='top',clip_on=True)
            # Second proton numbers
            for z in self.zMagic:
                # Find minimum and maximum N for this Z
                if (any(self.__Z_plot == z)):
                    nmin = np.min(self.__N_plot[self.__Z_plot == z])
                    nmax = np.max(self.__N_plot[self.__Z_plot == z])
                    self.ax.plot([nmin-self.magic_excess, nmax+self.magic_excess],[z-0.5, z-0.5], color='k', ls="--", lw=0.5)
                    self.ax.plot([nmin-self.magic_excess, nmax+self.magic_excess],[z+0.5, z+0.5], color='k', ls="--", lw=0.5)
                    # write the number on the bottom
                    self.ax.text(nmin-self.magic_excess-1,z,int(z),ha='right',va='center',clip_on=True)

        if self.separate_fission:
            # Make a legend
            # Get the ratio of x and y of the figure
            ratio = self.fig.get_figwidth()/self.fig.get_figheight()
            self.ax.add_patch(mpl.patches.Rectangle((0.88, 0.70), 0.01, 0.01*ratio, fill=False, transform=self.ax.transAxes, hatch="////", edgecolor='tab:blue', lw=1))
            self.ax.text(0.9, 0.70, 'Fissioning region', transform=self.ax.transAxes, fontsize=8, verticalalignment='bottom', horizontalalignment='left')
            self.ax.add_patch(mpl.patches.Rectangle((0.88, 0.67), 0.01, 0.01*ratio, fill=False, transform=self.ax.transAxes, hatch="////", edgecolor='tab:red', lw=1))
            self.ax.text(0.9, 0.67, 'Fission products', transform=self.ax.transAxes, fontsize=8, verticalalignment='bottom', horizontalalignment='left')



    def __init_plot_timescales(self):
        ls = ["-","--"]
        self.ts_plot = [ [ self.axTimescales.plot(self.ts_time,self.ts_data[i][j], color=self.timescale_colors[i], ls=ls[j],
                            label=(self.timescale_labels[i] if (j == 0) else "")) for j in range(len(self.ts_data[i]))] for i in range(len(self.ts_data))]
        # Also make the background of the box non-transparent
        self.axTimescales.legend(loc='upper right', ncol=2, bbox_to_anchor=(1.3, 1.0), frameon=True, facecolor='white', edgecolor='black', framealpha=1.0, fontsize=8)

    def __init_plot_addmainout(self):
        self.addmainout_plot = [self.axAddMainout.plot(self.addmainout_time,self.addmainout_data[k],
                                                label=self.addmainout_label[i], lw=2) for i,k in enumerate(self.addmainout_data.keys())]
        # Also make the background of the box non-transparent
        self.axAddMainout.legend(loc='upper right', ncol=1, bbox_to_anchor=(1.15, 1.0), frameon=True, facecolor='white', edgecolor='black', framealpha=1.0, fontsize=8)


    def __init_plot_energy(self):
        self.energy_plot = [self.axEnergy.plot(self.energy_time,self.energy_data[i], color=self.energy_colors[i],
                                                label=self.energy_labels[i], lw=self.energy_lw[i]) for i in range(len(self.energy_data))]
        # Also make the background of the box non-transparent
        self.axEnergy.legend(loc='upper right', ncol=2, bbox_to_anchor=(1.3, 1.0), frameon=True, facecolor='white', edgecolor='black', framealpha=1.0, fontsize=8)

    def __init_plot_tracked(self):
        self.tracked_plot = [self.axTracked.plot(self.tracked_time,self.track_nuclei_data[i],
                                                label=self.track_nuclei_labels[i]) for i in range(len(self.track_nuclei_labels))]
        # Also make the background of the box non-transparent
        self.axTracked.legend(loc='upper right', ncol=2, bbox_to_anchor=(1.3, 1.0), frameon=True, facecolor='white', edgecolor='black', framealpha=1.0, fontsize=8)



    def init_cbars(self, abun_cbar, flow_cbar):
        """
           Initialize the colorbars.
        """
        # Set up the number of colorbars
        n_cbars = abun_cbar + flow_cbar
        if n_cbars == 0:
            self.cax = []
            return
        if n_cbars == 1:
            self.cax = [self.fig.add_axes([0.75, 0.88, 0.14, 0.02])]
        elif n_cbars == 2:
            self.cax = [self.fig.add_axes([0.58, 0.88, 0.14, 0.02]),
                   self.fig.add_axes([0.75, 0.88, 0.14, 0.02])]

        # Counter for the colorbars
        ii = 0
        # Plot the abundance colorbar
        if abun_cbar:
            # Make a custom colormap since the abundance one has the background colors in as well
            self.abun_cbar = self.fig.colorbar(mpl.cm.ScalarMappable(norm=LogNorm(vmin=self.X_min,vmax=self.X_max), cmap="inferno"),
                    cax=self.cax[ii], orientation='horizontal', label='')
            self.abun_cbar.ax.set_title('Mass fraction')
            ii += 1

        # Plot the flow colorbar
        if flow_cbar:
            self.flow_cbar=self.fig.colorbar(
                self.quiver,
                cax=self.cax[ii],
                orientation='horizontal',
                label='',
                )
            self.flow_cbar.ax.set_title('Flow')
            ii += 1


    def update_data(self, ii):
        """
           Update the data for the plot.
        """

        self.time = self.wreader.snapshot_time[ii]
        yy = self.wreader.Y[ii]
        xx = yy * self.A
        self.abun = self.__background_Y.copy()
        mask = xx > self.X_min
        self.abun[self.N[mask], self.Z[mask]] = np.log10(xx[mask])

        self.Asum, self.Xsum = self.sum_over_A(self.A, xx)
        self.Abar = 1.0/np.sum(yy)

        for i in range(len(self.massBins)):
            mask = (self.A >= self.massBins[i][0]) & (self.A <= self.massBins[i][1])
            self.Xbins[i] = np.sum(xx[mask])


        if self.plot_timescales:
            self.__init_data_timescales(ii)

        if self.plot_energy:
            self.__init_data_energy(ii)

        if self.plot_addmainout:
            self.__init_data_addmainout(ii)

        if self.plot_tracked:
            self.__init_data_tracked(ii)

        if self.plot_mainout:
            self.mainout_time        = self.wreader.mainout['time'][:ii]
            self.mainout_time        = self.mainout_time-self.mainout_time[0]
            self.mainout_density     = self.wreader.mainout['dens'][:ii]
            self.mainout_temperature = self.wreader.mainout['temp'][:ii]
            self.mainout_ye          = self.wreader.mainout['ye'][:ii]

        if self.indicate_r_path:
            self.ngamma_eq.calc_r_process_path(self.wreader.mainout['dens'][ii],self.wreader.mainout['temp'][ii],self.wreader.mainout['yn'][ii])


        if ii == 0: return # no flows in first timestep

        if self.plot_flow:
            fpath=f'{self.path}/WinNet_data.h5'
            flow_dict = self.wreader.flow_entry(ii, self.flow_group)
            nin = flow_dict['n_in']
            zin = flow_dict['p_in']
            nout = flow_dict['n_out']
            zout = flow_dict['p_out']
            flow = flow_dict['flow']
            dz = zout - zin
            dn = nout - nin
            mask = flow > min(self.flow_norm.vmin,self.fission_minflow)
            self.flow_N = nin[mask]
            self.flow_Z = zin[mask]
            self.flow_dn = dn[mask]
            self.flow_dz = dz[mask]
            self.flow = flow[mask]
            # Keep track of the maximum flows for the colorbar
            self.flow_max_history    = np.roll(self.flow_max_history,1)
            try:
                ## might raise error if flow is not above threshold
                self.flow_max_history[0] = np.max(self.flow)
            except:
                pass
        if self.separate_fission:
            self.handle_fission()

    def __init_data_timescales(self, ii):
        self.ts_time = self.wreader.tau['time'][:ii]
        self.ts_time = self.ts_time-self.ts_time[0]
        self.ts_data = [ [ self.wreader.tau['tau_'+str(self.timescale_entries[i][j])][:ii]
                            for j in range(len(self.timescale_entries[i]))] for i in range(len(self.timescale_entries)) ]

    def __init_data_addmainout(self, ii):
        self.addmainout_time = self.wreader.mainout['time'][:ii]
        self.addmainout_time = self.addmainout_time-self.addmainout_time[0]
        self.addmainout_label = [r"Y$_n$", r"Y$_p$", r"Y$_\alpha$", r"Y$_{\text{heavy}}$", r"Y$_{\text{light}}$"]
        self.addmainout_data = {}
        self.addmainout_data["yn"]  = self.wreader.mainout['yn'][:ii]
        self.addmainout_data["yp"]  = self.wreader.mainout['yp'][:ii]
        self.addmainout_data["ya"]  = self.wreader.mainout['ya'][:ii]
        self.addmainout_data["yheavy"]  = self.wreader.mainout['yheavy'][:ii]
        self.addmainout_data["ylight"]  = self.wreader.mainout['ylight'][:ii]

    def __init_data_energy(self, ii):
        self.energy_time = self.wreader.energy['time'][:ii]
        self.energy_time = self.energy_time-self.energy_time[0]
        self.energy_data = [ self.wreader.energy['engen_'+self.energy_entries[i]][:ii] for i in range(len(self.energy_entries)) ]

    def __init_data_tracked(self, ii, force_label_init=False):
        self.tracked_time = self.wreader.tracked_nuclei['time'][:ii]
        self.tracked_time = self.tracked_time-self.tracked_time[0]
        self.track_nuclei_data  = [ self.wreader.tracked_nuclei[n][:ii] for n in self.wreader.tracked_nuclei['names'] ]
        if ii == -1 or force_label_init:
            self.track_nuclei_labels= self.wreader.tracked_nuclei['latex_names']

    def handle_fission(self,):
        """
           Separate the fissioning region from the fission products in the flow plot.
        """

        fis_mask = np.abs(self.flow_dn + self.flow_dz) > 32
        fis_N = self.flow_N[fis_mask]
        fis_Z = self.flow_Z[fis_mask]
        fis_dn = self.flow_dn[fis_mask]
        fis_dz = self.flow_dz[fis_mask]
        fis_flow = self.flow[fis_mask]

        mask = self.flow>self.flow_norm.vmin

        self.flow_N = self.flow_N[((~fis_mask) & mask)]
        self.flow_Z = self.flow_Z[((~fis_mask) & mask)]
        self.flow_dn = self.flow_dn[((~fis_mask) & mask)]
        self.flow_dz = self.flow_dz[((~fis_mask) & mask)]
        self.flow = self.flow[((~fis_mask) & mask)]

        self.fis_region.fill(0)
        for nn, zz, dn, dz, ff in zip(fis_N, fis_Z, fis_dn, fis_dz, fis_flow):
            self.fis_region[nn, zz] -= ff
            self.fis_region[nn+dn, zz+dz] += ff
        self.fis_region[(self.fis_region == 0)] = np.nan


    def update_abun_plot(self,):
        self.abun_im.set_array(self.abun.T)
        self.mafra_plot[0].set_ydata(self.Xsum)
        if self.plot_abar:
            self.Abar_plot.set_xdata([self.Abar])
        if self.plotMassBins:
            [self.mafrabin_plot[i][0].set_height(self.Xbins[i]) for i in range(len(self.massBins))]
        if self.plot_timescales:
            [ [ self.ts_plot[i][j][0].set_xdata(self.ts_time) for j in range(len(self.ts_plot[i]))] for i in range(len(self.ts_plot)) ]
            [ [ self.ts_plot[i][j][0].set_ydata(self.ts_data[i][j]) for j in range(len(self.ts_plot[i]))] for i in range(len(self.ts_plot)) ]
        if self.plot_energy:
            [ self.energy_plot[i][0].set_xdata(self.energy_time) for i in range(len(self.energy_plot))]
            [ self.energy_plot[i][0].set_ydata(self.energy_data[i]) for i in range(len(self.energy_plot))]
        if self.plot_addmainout:
            [ self.addmainout_plot[i][0].set_xdata(self.addmainout_time) for i in range(len(self.addmainout_plot))]
            [ self.addmainout_plot[i][0].set_ydata(self.addmainout_data[k]) for i,k in enumerate(self.addmainout_data.keys())]
        if self.plot_tracked:
            [ self.tracked_plot[i][0].set_xdata(self.tracked_time) for i in range(len(self.tracked_plot))]
            [ self.tracked_plot[i][0].set_ydata(self.track_nuclei_data[i]) for i in range(len(self.tracked_plot))]

        if self.plot_mainout:
            self.mainout_dens_plot[0].set_xdata(self.mainout_time)
            self.mainout_dens_plot[0].set_ydata(self.mainout_density)
            self.mainout_temp_plot[0].set_xdata(self.mainout_time)
            self.mainout_temp_plot[0].set_ydata(self.mainout_temperature)
            self.mainout_ye_plot[0].set_xdata(self.mainout_time)
            self.mainout_ye_plot[0].set_ydata(self.mainout_ye)
            right_side = f"= {self.to_latex_exponent(self.format_time(self.mainout_time[-1])[0])}\n"+\
                         f"= {self.to_latex_exponent(self.mainout_density[-1])}\n"+\
                         f"= {self.to_latex_exponent(self.mainout_temperature[-1])}\n"+\
                         f"= {self.mainout_ye[-1]:.3f}"
            self.Mainout_text.set_text(right_side)
            units = f"{self.format_time(self.mainout_time[-1])[1]}\n"+\
                    f"g/cm$^3$\n"+\
                    "GK\n"+\
                    ""
            self.Mainout_units.set_text(units)


    def update_fission_plot(self,):
        if self.plot_flow:
            fisspos = self.fis_region.T.copy()
            fisspos[self.fis_region.T<0] = np.nan
            self.fis_im_pos.set_array(fisspos)
            fissneg = self.fis_region.T.copy()
            fissneg[self.fis_region.T>0] = np.nan
            self.fis_im_neg.set_array(fissneg)

    def update_flow_plot(self,):
        if self.plot_flow:
            # Adapt flowmin and flowmax
            if self.flow_adapt_prange:
                lmaxflow = (np.nanmean(np.log10(self.flow_max_history)))+0.5+self.flow_max_offset
                lminflow = lmaxflow-self.flow_prange-self.flow_min_offset
                self.flow_cbar.mappable.set_clim(vmin=10**lminflow, vmax=10**lmaxflow)
                self.flow_max = 10**lmaxflow
                self.flow_min = 10**lminflow

            # Plot with a quiver as the width does not have to be adapted
            if not self.flow_adapt_width:
                self.quiver.N=len(self.flow_N)
                self.quiver.set_offsets(np.array([self.flow_N, self.flow_Z]).T)
                self.quiver.set_UVC(self.flow_dn, self.flow_dz, self.flow)
            else:
                # Quiver does not allow for changing the width of the arrows
                # Therefore, we have to us a patchcollection and draw it everytime new.
                self.flow_patch.remove()
                width = (self.flow_maxArrowWidth-self.flow_minArrowWidth)/self.flow_prange
                with np.errstate(divide='ignore'):
                    arrowwidth = (np.log10(self.flow)-np.log10(self.flow_min))*width + self.flow_minArrowWidth
                flow_arrows = [Arrow(self.flow_N[i],self.flow_Z[i],self.flow_dn[i],self.flow_dz[i],width=arrowwidth[i],color='k') for i in range(len(self.flow))]
                a = PatchCollection(flow_arrows, cmap=self.cmapNameFlow, norm=self.flow_norm)
                a.set_array(self.flow)
                self.flow_patch = self.ax.add_collection(a)

    def update_ngamma_plot(self,):
        if self.indicate_r_path:
            self.ngamma_eq_plot[0].set_xdata(self.ngamma_eq.path_N)
            self.ngamma_eq_plot[0].set_ydata(self.ngamma_eq.path_Z)
            self.ngamma_eq_plot_o[0].set_xdata(self.ngamma_eq.path_N)
            self.ngamma_eq_plot_o[0].set_ydata(self.ngamma_eq.path_Z)


    def update_frame(self, ii):
        self.update_data(ii)

        self.update_abun_plot()
        self.update_fission_plot()
        self.update_flow_plot()
        self.update_ngamma_plot()
        self.update_slider(ii)
        self.update_interactive_text(ii)
        return ii


    def update_interactive_text(self, ii):
        if self.interactive:
            if self.__interactive_textbox is not None:
                # Get the text and only change the last line
                text = self.__interactive_textbox.get_text()
                text = text.split("\n")
                N = int(text[1].split(" = ")[1])
                Z = int(text[2].split(" = ")[1])
                X = 10**self.abun[N,Z]
                if X < self.X_min:
                    X = 0
                text[-1] = f"X = {X:.2e}"
                self.__interactive_textbox.set_text("\n".join(text))

    def save_frame(self, ii):
        self.update_frame(ii)
        # self.fig.canvas.draw()
        # self.fig.canvas.flush_events()
        plt.savefig(f'{self.frame_dir}/frame_{ii}.png',
                    dpi=300, bbox_inches='tight')
        return ii

    def get_funcanimation(self, frames=None, **kwargs):
        if frames is None:
            frames = range(self.n_timesteps)
        self.frames=frames
        self.animation = FuncAnimation(self.fig, self.update_frame,
            frames=frames, **kwargs)
        if self.interactive:
            self.time_slider()
        return self.animation

    def time_slider(self):
        self.slider_bar = Slider(self.ax_slider, '', 1, self.n_timesteps-1, valinit=1)
        self.slider_bar.valtext.set_text('')
        self.slider_bar.on_changed(self.on_slider_update)
        self.fig.canvas.mpl_connect('button_press_event', self.on_slider_click)
        self.fig.canvas.mpl_connect('button_release_event', self.on_slider_release)

    def on_slider_update(self, val):
        ii = int(self.slider_bar.val)
        self.update_frame(ii)

    def update_slider(self, ii):
        try:
            # avoid infinite recursion if misusing slider as timebar for animation
            self.slider_bar.eventson =False
            self.slider_bar.set_val(ii)
            self.slider_bar.eventson = True

            self.slider_bar.valtext.set_text('')
        except:
            pass

    def pause_movie(self, click_event):
        if self.movie_paused:
            # Only resume if the slider is not being dragged
            ii = int(self.slider_bar.val)
            new_seq = list(range(ii, self.frames[-1])) + list(range(self.frames[0], ii))
            self.animation._iter_gen = lambda: iter(new_seq)
            self.animation.frame_seq = self.animation.new_frame_seq()
            self.animation.resume()
            self.movie_paused = False
            self.play_button.label.set_text("❚❚")
            self.play_button.label.set_color("red")
            # Update the appearance of the button
        else:
            self.animation.pause()
            self.movie_paused = True
            self.play_button.label.set_text(" ▶")
            self.play_button.label.set_color("green")

            # Update the appearance of the button

    def on_slider_click(self, event):
        if event.inaxes == self.ax_slider:
            self.animation.pause()  # Pause the animation when clicking the slider
            self.movie_paused = True
            self.play_button.label.set_text(" ▶")
            self.play_button.label.set_color("green")
        elif event.inaxes == self.ax:
            toolbar = plt.get_current_fig_manager().toolbar
            active_tool = toolbar.mode  # Get the current active tool
            if active_tool == '':
                if event.dblclick:
                    self.__toggle_zoom()
                else:
                    # Get the coordinates of the click
                    x, y = event.xdata, event.ydata
                    nucl_n = int(np.round(x))
                    nucl_z = int(np.round(y))
                    # Check if abundances are nan there
                    if np.isnan(self.abun[nucl_n, nucl_z]):
                        if not (self.__interactive_box is None):
                            # Remove the rectangle if it exists
                            self.__interactive_box.remove()
                            self.__interactive_box = None
                        if not (self.__interactive_textbox is None):
                            # Remove the textbox if it exists
                            self.__interactive_textbox.remove()
                            self.__interactive_textbox = None
                    else:
                        if not (self.__interactive_box is None):
                            # Remove the rectangle if it exists
                            self.__interactive_box.remove()
                            self.__interactive_box = None
                        if not (self.__interactive_textbox is None):
                            # Remove the textbox if it exists
                            self.__interactive_textbox.remove()
                            self.__interactive_textbox = None

                        # Create a rectangle there
                        self.__interactive_box = self.ax.add_patch(
                            patches.Rectangle((nucl_n-0.5, nucl_z-0.5), 1, 1, linewidth=1, edgecolor='r', facecolor='none'))

                        # Create textbox with information
                        df = self.winvn.get_dataframe()
                        # Get the name of the nucleus
                        nucl_name = df.loc[(nucl_n, nucl_z), 'name']
                        text = nucl_name.capitalize() + "\n"
                        text+= f"N = {nucl_n}\nZ = {nucl_z}\nA = {nucl_n+nucl_z}"
                        # Add mass fraction
                        X = 10**self.abun[nucl_n, nucl_z]
                        # Set 0 if below limit
                        if X < self.X_min:
                            X = 0
                        text+= f"\nX = {X:.2e}"


                        # Also make a border around it
                        self.__interactive_textbox = self.ax.text(0.44, 0.23, text, transform=self.ax.transAxes, fontsize=10,
                                                                horizontalalignment='left',
                                                                verticalalignment='bottom', bbox=dict(facecolor='white',
                                                                edgecolor='black', boxstyle='round,pad=0.5'))
        self.fig.canvas.draw_idle()


    def on_slider_release(self, event):
        # Reset slider_dragging to False when the mouse is released
        if event.inaxes == self.ax_slider:
            self.on_slider_update(self.slider_bar.val)  # Final update after release

    def arrow_update(self, event):
        ii = int(self.slider_bar.val)
        if event.key in ['left', 'down']:
            self.movie_paused = False
            self.pause_movie(event)
            if ii !=self.slider_bar.valmin:
                self.slider_bar.set_val(ii-1)
        elif event.key in ['right', 'up']:
            self.movie_paused = False
            self.pause_movie(event)
            if ii !=self.slider_bar.valmax:
                self.slider_bar.set_val(ii+1)
        elif event.key == " ":
            self.pause_movie(event)
        elif ((event.key == "t") or (event.key == "e")
               or (event.key == 'n') or (event.key == 'd')
               or (event.key == 'm')):

            if not (self.__interactive_ax is None):
                self.__interactive_ax.remove()

            if event.key == "t":
                # Toggle timescales
                self.__toggle_timescales()
            elif event.key == 'e':
                # Toggle energy
                self.__toggle_energy()
            elif event.key == 'n':
                # Toggle tracked nuclei
                self.__toggle_tracked()
            elif event.key == 'm':
                # Toggle tracked nuclei
                self.__toggle_addmainout()
            elif event.key == "d":
                # Shut of additional plots
                self.plot_timescales = False
                self.plot_energy = False
                self.plot_tracked = False
                self.plot_addmainout = False
                self.__interactive_ax = None

            self.fig.canvas.draw_idle()

    def __toggle_timescales(self):
        if 'timescales' in self.available_data:
            # Toggle timescales
            if self.plot_timescales == True:
                self.plot_timescales = False
                self.plot_energy = False
                self.plot_tracked = False
                self.plot_addmainout = False
                self.__interactive_ax = None
            else:
                ii = int(self.slider_bar.val)
                self.__init_data_timescales(ii)
                self.__init_axTimescales()
                self.plot_timescales = True
                self.plot_energy = False
                self.plot_tracked = False
                self.plot_addmainout = False
                self.__init_plot_timescales()
                self.__interactive_ax = self.axTimescales

    def __toggle_energy(self):
        if 'energy' in self.available_data:
            # Toggle energy
            if self.plot_energy == True:
                self.plot_timescales = False
                self.plot_energy = False
                self.plot_tracked = False
                self.plot_addmainout = False
                self.__interactive_ax = None
            else:
                ii = int(self.slider_bar.val)
                self.__init_data_energy(ii)
                self.__init_axEnergy()
                self.plot_timescales = False
                self.plot_energy = True
                self.plot_tracked = False
                self.plot_addmainout = False
                self.__init_plot_energy()
                self.__interactive_ax = self.axEnergy

    def __toggle_addmainout(self):
        if 'mainout' in self.available_data:
            # Toggle additional mainout information
            if self.plot_addmainout == True:
                self.plot_timescales = False
                self.plot_energy = False
                self.plot_tracked = False
                self.plot_addmainout = False
                self.__interactive_ax = None
            else:
                ii = int(self.slider_bar.val)
                self.__init_data_addmainout(ii)
                self.__init_axAddMainout()
                self.plot_timescales = False
                self.plot_energy = False
                self.plot_tracked = False
                self.plot_addmainout = True
                self.__init_plot_addmainout()
                self.__interactive_ax = self.axAddMainout

    def __toggle_tracked(self):
        if 'tracked_nuclei' in self.available_data:
            # Toggle energy
            if self.plot_tracked == True:
                self.plot_timescales = False
                self.plot_energy = False
                self.plot_tracked = False
                self.plot_addmainout = False
                self.__interactive_ax = None
            else:
                ii = int(self.slider_bar.val)
                self.__init_data_tracked(ii, force_label_init=True)
                self.__init_axTracked()
                self.plot_timescales = False
                self.plot_energy = False
                self.plot_tracked = True
                self.plot_addmainout = False
                self.__init_plot_tracked()
                self.__interactive_ax = self.axTracked

    def __toggle_zoom(self):
        if self.zoomed:
            self.ax.set_xlim(self.limits_plot[0])
            self.ax.set_ylim(self.limits_plot[1])
            self.zoomed = False
        else:
            self.ax.set_xlim(-5.5, 100)
            self.ax.set_ylim(-6.5, 50)
            self.zoomed = True
        self.fig.canvas.draw_idle()



    @staticmethod
    def to_latex_exponent(x):
        """
          Convert a number to a latex string with exponent.
        """
        if x == 0:
            return '0'
        exponent = np.floor(np.log10(x))
        mantissa = x / 10**exponent
        return f'{mantissa:.2f}'+r'$\times 10^{'+f'{int(exponent)}'+r'}$'


    @staticmethod
    def format_time(time):
        """
          Format a time in seconds to a more human readable format.
        """
        if time < 1:
            return time*1000, 'ms'
        if time < 60:
            return time, 's'
        if time < 3600:
            return time/60, 'min'
        if time < 3600*24:
            return time/3600, 'h'
        if time < 3600*24*365:
            return time/86400, 'd'
        if time < 3600*24*365*1000:
            return time/31536000, 'y'
        if time < 3600*24*365*1000*1000:
            return time/31536000/1000, 'ky'
        else:
            return time/31536000_000_000, 'My'

    def sum_over_A(self,A,X):
        A_unique = np.arange(max(A)+1)
        X_sum = np.zeros_like(A_unique,dtype=float)
        for a, x in zip(A,X):
            X_sum[int(a)] += x
        return A_unique, X_sum

################################################################################

# Funcanimation

if __name__ == "__main__":
    run_path = "../Example_NSM_dyn_ejecta_rosswog_hdf5"
    fig = plt.figure(figsize=(15, 8))

    # z_drip , n_drip = np.loadtxt('/home/mjacobi/frdm/dripline.dat', unpack=True)
    # plt.plot(n_drip, z_drip, 'k--', lw=1)

    anim = FlowAnimation(run_path, fig, flow_group='flows')
    ani = anim.get_funcanimation(interval=10, frames=range(1, anim.n_timesteps))
    plt.show()
