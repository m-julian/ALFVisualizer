
import numpy as np
import pyvista
from pyvista import _vtk

IPYGANY_MAP = {'reds': 'Reds',
               'spectral': 'Spectral'}


# shamelessly copied from matplotlib.colors
hexcolors = {
    'aliceblue':            '#F0F8FF',
    'antiquewhite':         '#FAEBD7',
    'aqua':                 '#00FFFF',
    'aquamarine':           '#7FFFD4',
    'azure':                '#F0FFFF',
    'beige':                '#F5F5DC',
    'bisque':               '#FFE4C4',
    'black':                '#000000',
    'blanchedalmond':       '#FFEBCD',
    'blue':                 '#0000FF',
    'blueviolet':           '#8A2BE2',
    'brown':                '#654321',
    'burlywood':            '#DEB887',
    'cadetblue':            '#5F9EA0',
    'chartreuse':           '#7FFF00',
    'chocolate':            '#D2691E',
    'coral':                '#FF7F50',
    'cornflowerblue':       '#6495ED',
    'cornsilk':             '#FFF8DC',
    'crimson':              '#DC143C',
    'cyan':                 '#00FFFF',
    'darkblue':             '#00008B',
    'darkcyan':             '#008B8B',
    'darkgoldenrod':        '#B8860B',
    'darkgray':             '#A9A9A9',
    'darkgreen':            '#006400',
    'darkgrey':             '#A9A9A9',
    'darkkhaki':            '#BDB76B',
    'darkmagenta':          '#8B008B',
    'darkolivegreen':       '#556B2F',
    'darkorange':           '#FF8C00',
    'darkorchid':           '#9932CC',
    'darkred':              '#8B0000',
    'darksalmon':           '#E9967A',
    'darkseagreen':         '#8FBC8F',
    'darkslateblue':        '#483D8B',
    'darkslategray':        '#2F4F4F',
    'darkslategrey':        '#2F4F4F',
    'darkturquoise':        '#00CED1',
    'darkviolet':           '#9400D3',
    'deeppink':             '#FF1493',
    'deepskyblue':          '#00BFFF',
    'dimgray':              '#696969',
    'dimgrey':              '#696969',
    'dodgerblue':           '#1E90FF',
    'firebrick':            '#B22222',
    'floralwhite':          '#FFFAF0',
    'forestgreen':          '#228B22',
    'fuchsia':              '#FF00FF',
    'gainsboro':            '#DCDCDC',
    'ghostwhite':           '#F8F8FF',
    'gold':                 '#FFD700',
    'goldenrod':            '#DAA520',
    'gray':                 '#808080',
    'green':                '#008000',
    'greenyellow':          '#ADFF2F',
    'grey':                 '#808080',
    'honeydew':             '#F0FFF0',
    'hotpink':              '#FF69B4',
    'indianred':            '#CD5C5C',
    'indigo':               '#4B0082',
    'ivory':                '#FFFFF0',
    'khaki':                '#F0E68C',
    'lavender':             '#E6E6FA',
    'lavenderblush':        '#FFF0F5',
    'lawngreen':            '#7CFC00',
    'lemonchiffon':         '#FFFACD',
    'lightblue':            '#ADD8E6',
    'lightcoral':           '#F08080',
    'lightcyan':            '#E0FFFF',
    'lightgoldenrodyellow': '#FAFAD2',
    'lightgray':            '#D3D3D3',
    'lightgreen':           '#90EE90',
    'lightgrey':            '#D3D3D3',
    'lightpink':            '#FFB6C1',
    'lightsalmon':          '#FFA07A',
    'lightseagreen':        '#20B2AA',
    'lightskyblue':         '#87CEFA',
    'lightslategray':       '#778899',
    'lightslategrey':       '#778899',
    'lightsteelblue':       '#B0C4DE',
    'lightyellow':          '#FFFFE0',
    'lime':                 '#00FF00',
    'limegreen':            '#32CD32',
    'linen':                '#FAF0E6',
    'magenta':              '#FF00FF',
    'maroon':               '#800000',
    'mediumaquamarine':     '#66CDAA',
    'mediumblue':           '#0000CD',
    'mediumorchid':         '#BA55D3',
    'mediumpurple':         '#9370DB',
    'mediumseagreen':       '#3CB371',
    'mediumslateblue':      '#7B68EE',
    'mediumspringgreen':    '#00FA9A',
    'mediumturquoise':      '#48D1CC',
    'mediumvioletred':      '#C71585',
    'midnightblue':         '#191970',
    'mintcream':            '#F5FFFA',
    'mistyrose':            '#FFE4E1',
    'moccasin':             '#FFE4B5',
    'navajowhite':          '#FFDEAD',
    'navy':                 '#000080',
    'oldlace':              '#FDF5E6',
    'olive':                '#808000',
    'olivedrab':            '#6B8E23',
    'orange':               '#FFA500',
    'orangered':            '#FF4500',
    'orchid':               '#DA70D6',
    'palegoldenrod':        '#EEE8AA',
    'palegreen':            '#98FB98',
    'paleturquoise':        '#AFEEEE',
    'palevioletred':        '#DB7093',
    'papayawhip':           '#FFEFD5',
    'peachpuff':            '#FFDAB9',
    'peru':                 '#CD853F',
    'pink':                 '#FFC0CB',
    'plum':                 '#DDA0DD',
    'powderblue':           '#B0E0E6',
    'purple':               '#800080',
    'raw_sienna':           '#965434',
    'rebeccapurple':        '#663399',
    'red':                  '#FF0000',
    'rosybrown':            '#BC8F8F',
    'royalblue':            '#4169E1',
    'saddlebrown':          '#8B4513',
    'salmon':               '#FA8072',
    'sandybrown':           '#F4A460',
    'seagreen':             '#2E8B57',
    'seashell':             '#FFF5EE',
    'sienna':               '#A0522D',
    'silver':               '#C0C0C0',
    'skyblue':              '#87CEEB',
    'slateblue':            '#6A5ACD',
    'slategray':            '#708090',
    'slategrey':            '#708090',
    'snow':                 '#FFFAFA',
    'springgreen':          '#00FF7F',
    'steelblue':            '#4682B4',
    'tan':                  '#D2B48C',
    'teal':                 '#008080',
    'thistle':              '#D8BFD8',
    'tomato':               '#FF6347',
    'turquoise':            '#40E0D0',
    'violet':               '#EE82EE',
    'wheat':                '#F5DEB3',
    'white':                '#FFFFFF',
    'whitesmoke':           '#F5F5F5',
    'yellow':               '#FFFF00',
    'yellowgreen':          '#9ACD32',
    'tab:blue':             '#1f77b4',
    'tab:orange':           '#ff7f0e',
    'tab:green':            '#2ca02c',
    'tab:red':              '#d62728',
    'tab:purple':           '#9467bd',
    'tab:brown':            '#8c564b',
    'tab:pink':             '#e377c2',
    'tab:gray':             '#7f7f7f',
    'tab:olive':            '#bcbd22',
    'tab:cyan':             '#17becf'}

color_char_to_word = {
        'b': 'blue',
        'g': 'green',
        'r': 'red',
        'c': 'cyan',
        'm': 'magenta',
        'y': 'yellow',
        'k': 'black',
        'w': 'white'}


PARAVIEW_BACKGROUND = [82/255., 87/255., 110/255.]


def hex_to_rgb(h):
    """Return 0 to 1 rgb from a hex list or tuple."""
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16)/255. for i in (0, 2, 4))


def string_to_rgb(string):
    """Convert a literal color string (i.e. white) to a color rgb.

    Also accepts hex strings or single characters from the following list.

        b: blue
        g: green
        r: red
        c: cyan
        m: magenta
        y: yellow
        k: black
        w: white

    """
    # if a single character
    if len(string) == 1:

        # Convert from single character to full hex
        if string.lower() not in color_char_to_word:
            raise ValueError('Single character string must be one of the following:'
                             f'\n{str(color_char_to_word.keys())}')

        colorhex = hexcolors[color_char_to_word[string.lower()]]

    # if a color name
    elif string.lower() in hexcolors:
        colorhex = hexcolors[string.lower()]

    elif string.lower() in 'paraview' or string.lower() in 'pv':
        # Use the default ParaView background color
        return PARAVIEW_BACKGROUND

    # try to convert to hex
    else:
        try:
            return hex_to_rgb(string)
        except:
            raise ValueError('Invalid color string or hex string.')

    return hex_to_rgb(colorhex)


def get_cmap_safe(cmap):
    """Fetch a colormap by name from matplotlib, colorcet, or cmocean."""
    try:
        from matplotlib.cm import get_cmap
    except ImportError:
        raise ImportError('cmap requires matplotlib')
    if isinstance(cmap, str):
        # check if this colormap has been mapped between ipygany
        if cmap in IPYGANY_MAP:
            cmap = IPYGANY_MAP[cmap]

        # Try colorcet first
        try:
            import colorcet
            cmap = colorcet.cm[cmap]
        except (ImportError, KeyError):
            pass
        else:
            return cmap
        # Try cmocean second
        try:
            import cmocean
            cmap = getattr(cmocean.cm, cmap)
        except (ImportError, AttributeError):
            pass
        else:
            return cmap
        # Else use Matplotlib
        cmap = get_cmap(cmap)
    elif isinstance(cmap, list):
        for item in cmap:
            if not isinstance(item, str):
                raise TypeError('When inputting a list as a cmap, each item should be a string.')
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(cmap)

    return cmap


def parse_color(color, opacity=None, default_color=None):
    """Parse color into a vtk friendly rgb list.

    Values returned will be between 0 and 1.

    """
    if color is None:
        if default_color is None:
            color = pyvista.global_theme.color
        else:
            color = default_color
    if isinstance(color, str):
        color = string_to_rgb(color)
    elif len(color) == 3:
        pass
    elif len(color) == 4:
        color = color[:3]
    else:
        raise ValueError(f"""
    Invalid color input: ({color})
    Must be string, rgb list, or hex color string.  For example:
        color='white'
        color='w'
        color=[1, 1, 1]
        color='#FFFFFF'""")
    if opacity is not None and isinstance(opacity, (float, int)):
        color = [color[0], color[1], color[2], opacity]
    return color


"""An internal module for wrapping the use of mappers."""


def make_mapper(mapper_class):
    """Wrap a mapper.

    This makes a mapper wrapped with a few convenient tools for managing
    mappers with scalar bars in a consistent way since not all mapper classes
    have scalar ranges and lookup tables.
    """

    class MapperHelper(mapper_class):
        """A helper that dynamically inherits the mapper's class."""

        def __init__(self, *args, **kwargs):
            self._scalar_range = None
            self._lut = None

        @property
        def scalar_range(self):
            if hasattr(self, 'GetScalarRange'):
                self._scalar_range = self.GetScalarRange()
            return self._scalar_range

        @scalar_range.setter
        def scalar_range(self, clim):
            if hasattr(self, 'SetScalarRange'):
                self.SetScalarRange(*clim)
            if self.lookup_table is not None:
                self.lookup_table.SetRange(*clim)
            self._scalar_range = clim

        @property
        def lookup_table(self):
            if hasattr(self, 'GetLookupTable'):
                self._lut = self.GetLookupTable()
            return self._lut

        @lookup_table.setter
        def lookup_table(self, lut):
            if hasattr(self, 'SetLookupTable'):
                self.SetLookupTable(lut)
            self._lut = lut

    return MapperHelper()


def change_actor_color(self, actor_name, plotter, color=None, style=None, scalars=None,
                clim=None, show_edges=None, edge_color=None,
                point_size=5.0, line_width=None, opacity=1.0,
                flip_scalars=False, lighting=None, n_colors=256,
                interpolate_before_map=True, cmap=None, label=None,
                reset_camera=None, scalar_bar_args=None, show_scalar_bar=None,
                multi_colors=False, name=None, texture=None,
                render_points_as_spheres=None, render_lines_as_tubes=False,
                smooth_shading=None, ambient=0.0, diffuse=1.0, specular=0.0,
                specular_power=100.0, nan_color=None, nan_opacity=1.0,
                culling=None, rgb=None, categories=False, silhouette=False,
                use_transparency=False, below_color=None, above_color=None,
                annotations=None, pickable=True, preference="point",
                log_scale=False, pbr=False, metallic=0.0, roughness=0.5,
                render=True, component=None, **kwargs):
    """Add any PyVista/VTK mesh or dataset that PyVista can wrap to the scene.

    This method is using a mesh representation to view the surfaces
    and/or geometry of datasets. For volume rendering, see
    :func:`pyvista.BasePlotter.add_volume`.

    Parameters
    ----------
    mesh : pyvista.DataSet or pyvista.MultiBlock
        Any PyVista or VTK mesh is supported. Also, any dataset
        that :func:`pyvista.wrap` can handle including NumPy
        arrays of XYZ points.

    color : str or 3 item list, optional, defaults to white
        Use to make the entire mesh have a single solid color.
        Either a string, RGB list, or hex color string.  For example:
        ``color='white'``, ``color='w'``, ``color=[1, 1, 1]``, or
        ``color='#FFFFFF'``. Color will be overridden if scalars are
        specified.

    style : str, optional
        Visualization style of the mesh.  One of the following:
        ``style='surface'``, ``style='wireframe'``, ``style='points'``.
        Defaults to ``'surface'``. Note that ``'wireframe'`` only shows a
        wireframe of the outer geometry.

    scalars : str or numpy.ndarray, optional
        Scalars used to "color" the mesh.  Accepts a string name
        of an array that is present on the mesh or an array equal
        to the number of cells or the number of points in the
        mesh.  Array should be sized as a single vector. If both
        ``color`` and ``scalars`` are ``None``, then the active
        scalars are used.

    clim : 2 item list, optional
        Color bar range for scalars.  Defaults to minimum and
        maximum of scalars array.  Example: ``[-1, 2]``. ``rng``
        is also an accepted alias for this.

    show_edges : bool, optional
        Shows the edges of a mesh.  Does not apply to a wireframe
        representation.

    edge_color : str or 3 item list, optional, defaults to black
        The solid color to give the edges when ``show_edges=True``.
        Either a string, RGB list, or hex color string.

    point_size : float, optional
        Point size of any nodes in the dataset plotted. Also
        applicable when style='points'. Default ``5.0``.

    line_width : float, optional
        Thickness of lines.  Only valid for wireframe and surface
        representations.  Default None.

    opacity : float, str, array-like
        Opacity of the mesh. If a single float value is given, it
        will be the global opacity of the mesh and uniformly
        applied everywhere - should be between 0 and 1. A string
        can also be specified to map the scalars range to a
        predefined opacity transfer function (options include:
        'linear', 'linear_r', 'geom', 'geom_r').  A string could
        also be used to map a scalars array from the mesh to the
        opacity (must have same number of elements as the
        ``scalars`` argument). Or you can pass a custom made
        transfer function that is an array either ``n_colors`` in
        length or shorter.

    flip_scalars : bool, optional
        Flip direction of cmap. Most colormaps allow ``*_r``
        suffix to do this as well.

    lighting : bool, optional
        Enable or disable view direction lighting. Default ``False``.

    n_colors : int, optional
        Number of colors to use when displaying scalars. Defaults to 256.
        The scalar bar will also have this many colors.

    interpolate_before_map : bool, optional
        Enabling makes for a smoother scalars display.  Default is
        ``True``.  When ``False``, OpenGL will interpolate the
        mapped colors which can result is showing colors that are
        not present in the color map.

    cmap : str, list, optional
        Name of the Matplotlib colormap to use when mapping the
        ``scalars``.  See available Matplotlib colormaps.  Only
        applicable for when displaying ``scalars``. Requires
        Matplotlib to be installed.  ``colormap`` is also an
        accepted alias for this. If ``colorcet`` or ``cmocean``
        are installed, their colormaps can be specified by name.

        You can also specify a list of colors to override an
        existing colormap with a custom one.  For example, to
        create a three color colormap you might specify
        ``['green', 'red', 'blue']``.

    label : str, optional
        String label to use when adding a legend to the scene with
        :func:`pyvista.BasePlotter.add_legend`.

    reset_camera : bool, optional
        Reset the camera after adding this mesh to the scene.

    scalar_bar_args : dict, optional
        Dictionary of keyword arguments to pass when adding the
        scalar bar to the scene. For options, see
        :func:`pyvista.BasePlotter.add_scalar_bar`.

    show_scalar_bar : bool
        If ``False``, a scalar bar will not be added to the
        scene. Defaults to ``True``.

    multi_colors : bool, optional
        If a ``MultiBlock`` dataset is given this will color each
        block by a solid color using matplotlib's color cycler.

    name : str, optional
        The name for the added mesh/actor so that it can be easily
        updated.  If an actor of this name already exists in the
        rendering window, it will be replaced by the new actor.

    texture : vtk.vtkTexture or np.ndarray or bool, optional
        A texture to apply if the input mesh has texture
        coordinates.  This will not work with MultiBlock
        datasets. If set to ``True``, the first available texture
        on the object will be used. If a string name is given, it
        will pull a texture with that name associated to the input
        mesh.

    render_points_as_spheres : bool, optional
        Render points as spheres rather than dots.

    render_lines_as_tubes : bool, optional
        Show lines as thick tubes rather than flat lines.  Control
        the width with ``line_width``.

    smooth_shading : bool, optional
        Enable smooth shading when ``True`` using either the 
        Gouraud or Phong shading algorithm.  When ``False``, use
        flat shading.
        Automatically enabled when ``pbr=True``.

    ambient : float, optional
        When lighting is enabled, this is the amount of light in
        the range of 0 to 1 (default 0.0) that reaches the actor
        when not directed at the light source emitted from the
        viewer.

    diffuse : float, optional
        The diffuse lighting coefficient. Default 1.0.

    specular : float, optional
        The specular lighting coefficient. Default 0.0.

    specular_power : float, optional
        The specular power. Between 0.0 and 128.0.

    nan_color : str or 3 item list, optional, defaults to gray
        The color to use for all ``NaN`` values in the plotted
        scalar array.

    nan_opacity : float, optional
        Opacity of ``NaN`` values.  Should be between 0 and 1.
        Default 1.0.

    culling : str, optional
        Does not render faces that are culled. Options are
        ``'front'`` or ``'back'``. This can be helpful for dense
        surface meshes, especially when edges are visible, but can
        cause flat meshes to be partially displayed.  Defaults to
        ``False``.

    rgb : bool, optional
        If an 2 dimensional array is passed as the scalars, plot
        those values as RGB(A) colors. ``rgba`` is also an
        accepted alias for this.  Opacity (the A) is optional.  If
        a scalars array ending with ``"_rgba"`` is passed, the default
        becomes ``True``.  This can be overridden by setting this
        parameter to ``False``.

    categories : bool, optional
        If set to ``True``, then the number of unique values in
        the scalar array will be used as the ``n_colors``
        argument.

    silhouette : dict, bool, optional
        If set to ``True``, plot a silhouette highlight for the
        mesh. This feature is only available for a triangulated
        ``PolyData``.  As a ``dict``, it contains the properties
        of the silhouette to display:

            * ``color``: ``str`` or 3-item ``list``, color of the silhouette
            * ``line_width``: ``float``, edge width
            * ``opacity``: ``float`` between 0 and 1, edge transparency
            * ``feature_angle``: If a ``float``, display sharp edges
                exceeding that angle in degrees.
            * ``decimate``: ``float`` between 0 and 1, level of decimation

    use_transparency : bool, optional
        Invert the opacity mappings and make the values correspond
        to transparency.

    below_color : str or 3 item list, optional
        Solid color for values below the scalars range
        (``clim``). This will automatically set the scalar bar
        ``below_label`` to ``'Below'``.

    above_color : str or 3 item list, optional
        Solid color for values below the scalars range
        (``clim``). This will automatically set the scalar bar
        ``above_label`` to ``'Above'``.

    annotations : dict, optional
        Pass a dictionary of annotations. Keys are the float
        values in the scalars range to annotate on the scalar bar
        and the values are the the string annotations.

    pickable : bool, optional
        Set whether this actor is pickable.

    preference : str, optional
        When ``mesh.n_points == mesh.n_cells`` and setting
        scalars, this parameter sets how the scalars will be
        mapped to the mesh.  Default ``'points'``, causes the
        scalars will be associated with the mesh points.  Can be
        either ``'points'`` or ``'cells'``.

    log_scale : bool, optional
        Use log scale when mapping data to colors. Scalars less
        than zero are mapped to the smallest representable
        positive float. Default: ``True``.

    pbr : bool, optional
        Enable physics based rendering (PBR) if the mesh is
        ``PolyData``.  Use the ``color`` argument to set the base
        color. This is only available in VTK>=9.

    metallic : float, optional
        Usually this value is either 0 or 1 for a real material
        but any value in between is valid. This parameter is only
        used by PBR interpolation. Default value is 0.0.

    roughness : float, optional
        This value has to be between 0 (glossy) and 1 (rough). A
        glossy material has reflections and a high specular
        part. This parameter is only used by PBR
        interpolation. Default value is 0.5.

    render : bool, optional
        Force a render when ``True``.  Default ``True``.

    component :  int, optional
        Set component of vector valued scalars to plot.  Must be
        nonnegative, if supplied. If ``None``, the magnitude of
        the vector is plotted.

    **kwargs : dict, optional
        Optional developer keyword arguments.

    Returns
    -------
    vtk.vtkActor
        VTK actor of the mesh.

    Examples
    --------
    Add a sphere to the plotter and show it with a custom scalar
    bar title.

    >>> import pyvista
    >>> sphere = pyvista.Sphere()
    >>> sphere['Data'] = sphere.points[:, 2]
    >>> plotter = pyvista.Plotter()
    >>> _ = plotter.add_mesh(sphere,
    ...                      scalar_bar_args={'title': 'Z Position'})
    >>> plotter.show()

    Plot using RGB on a single cell.  Note that since the number of
    points and the number of cells are identical, we have to pass
    ``preference='cell'``.

    >>> import pyvista
    >>> import numpy as np
    >>> vertices = np.array([[0, 0, 0], [1, 0, 0], [.5, .667, 0], [0.5, .33, 0.667]])
    >>> faces = np.hstack([[3, 0, 1, 2], [3, 0, 3, 2], [3, 0, 1, 3], [3, 1, 2, 3]])
    >>> mesh = pyvista.PolyData(vertices, faces)
    >>> mesh.cell_data['colors'] = [[255, 255, 255],
    ...                               [0, 255, 0],
    ...                               [0, 0, 255],
    ...                               [255, 0, 0]]
    >>> plotter = pyvista.Plotter()
    >>> _ = plotter.add_mesh(mesh, scalars='colors', lighting=False,
    ...                      rgb=True, preference='cell')
    >>> plotter.camera_position='xy'
    >>> plotter.show()

    Note how this varies from ``preference=='point'``.  This is
    because each point is now being individually colored, versus
    in ``preference=='point'``, each cell face is individually
    colored.

    >>> plotter = pyvista.Plotter()
    >>> _ = plotter.add_mesh(mesh, scalars='colors', lighting=False,
    ...                      rgb=True, preference='point')
    >>> plotter.camera_position='xy'
    >>> plotter.show()

    """

    # Avoid mutating input
    if scalar_bar_args is None:
        scalar_bar_args = {'n_colors': n_colors}
    else:
        scalar_bar_args = scalar_bar_args.copy()

    if show_edges is None:
        show_edges = self.plotter._theme.show_edges

    if edge_color is None:
        edge_color = self.plotter._theme.edge_color

    if show_scalar_bar is None:
        show_scalar_bar = self.plotter._theme.show_scalar_bar

    if lighting is None:
        lighting = self.plotter._theme.lighting

    # supported aliases
    clim = kwargs.pop('rng', clim)
    cmap = kwargs.pop('colormap', cmap)
    culling = kwargs.pop("backface_culling", culling)

    if not render_points_as_spheres:
        render_points_as_spheres = self.plotter._theme.render_points_as_spheres

    if not nan_color:
        nan_color = self.plotter._theme.nan_color
    nan_color = list(parse_color(nan_color))
    nan_color.append(nan_opacity)

    if color:
        color = self.plotter._theme.color

    if culling is True:
        culling = 'backface'

    rgb = kwargs.pop('rgba', rgb)

    self.mapper = make_mapper(_vtk.vtkDataSetMapper)
    self.mapper.SetInputData(self.mesh)
    self.mapper.GetLookupTable().SetNumberOfTableValues(n_colors)
    if interpolate_before_map:
        self.mapper.InterpolateScalarsBeforeMappingOn()

    actor = _vtk.vtkActor()
    prop = _vtk.vtkProperty()
    actor.SetMapper(self.mapper)
    actor.SetProperty(prop)

    # Handle making opacity array =========================================

    if use_transparency and np.max(opacity) <= 1.0:
        opacity = 1 - opacity
    elif use_transparency and isinstance(opacity, np.ndarray):
        opacity = 255 - opacity

    # Set actor properties ================================================

    # select view style
    if not style:
        style = 'surface'
    style = style.lower()
    if style == 'wireframe':
        prop.SetRepresentationToWireframe()
        if color is None:
            color = self.plotter._theme.outline_color
    elif style == 'points':
        prop.SetRepresentationToPoints()
    elif style == 'surface':
        prop.SetRepresentationToSurface()
    else:
        raise ValueError('Invalid style.  Must be one of the following:\n'
                            '\t"surface"\n'
                            '\t"wireframe"\n'
                            '\t"points"\n')

    prop.SetPointSize(point_size)
    prop.SetAmbient(ambient)
    prop.SetDiffuse(diffuse)
    prop.SetSpecular(specular)
    prop.SetSpecularPower(specular_power)

    # interpolation, no pbr or smooth shading
    prop.SetInterpolationToFlat()

    # edge display style
    if show_edges:
        prop.EdgeVisibilityOn()

    rgb_color = parse_color(color, default_color=self.plotter._theme.color)
    prop.SetColor(rgb_color)
    if isinstance(opacity, (float, int)):
        prop.SetOpacity(opacity)
    prop.SetEdgeColor(parse_color(edge_color))

    if render_points_as_spheres:
        prop.SetRenderPointsAsSpheres(render_points_as_spheres)
    if render_lines_as_tubes:
        prop.SetRenderLinesAsTubes(render_lines_as_tubes)

    # legend label
    if label:
        if not isinstance(label, str):
            raise TypeError('Label must be a string')
        geom = pyvista.Triangle()
        if scalars is not None:
            geom = pyvista.Box()
            rgb_color = parse_color('black')
        geom.points -= geom.center
        self._labels.append([geom, label, rgb_color])

    # lighting display style
    if not lighting:
        prop.LightingOff()

    # set line thickness
    if line_width:
        prop.SetLineWidth(line_width)

    self.Modified
    self.renderer.Modified()
    return actor
