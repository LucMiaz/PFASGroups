from PIL import Image
import svgutils.transform as sg
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
import svgutils.transform as sg


def merge_raster(imgs,buffer,ncols):
    """Merge multiple raster images into a grid layout.
    
    Parameters
    ----------
    imgs : list of PIL.Image
        List of PIL Image objects to arrange in grid
    buffer : int
        Spacing (in pixels) between images
    ncols : int
        Number of columns in the grid
    
    Returns
    -------
    tuple of (PIL.Image, int, int)
        - Merged image containing all input images in grid layout
        - Total width of the merged image
        - Total height of the merged image
    
    Notes
    -----
    Images are arranged left-to-right, top-to-bottom. Row height is determined
    by tallest image in that row.
    """
    coords = []
    tot_height = 0
    max_width = 0
    for imgs in [imgs[i:i+ncols] for i in range(0, len(imgs), ncols)]:
        tot_width = 0
        for img in imgs:
            coords.append((img,(tot_width,tot_height)))
            tot_width += img.width + buffer
        max_width = max(max_width, tot_width)
        tot_height = max(*[img.height for img in imgs],0)+ tot_height
    fig = Image.new("RGBA",(max_width, tot_height))
    for img, coord in coords:
        fig.paste(img,coord)
    return fig,tot_width,tot_height

def merge_svg(imgs,buffer, ncols):
    """Merge multiple SVG images into a single SVG with grid layout.
    
    Parameters
    ----------
    imgs : list of svgutils.transform.SVGFigure
        List of SVG figure objects to arrange in grid
    buffer : int
        Spacing (in pixels) between images
    ncols : int
        Number of columns in the grid
    
    Returns
    -------
    tuple of (svgutils.transform.SVGFigure, int, int)
        - Merged SVG figure containing all input SVGs
        - Total width of the merged figure
        - Total height of the merged figure
    
    Notes
    -----
    SVG dimensions are extracted from width/height attributes (supports 'px' and 'pt' units).
    Each input SVG is repositioned using moveto() to its grid location.
    """
    coords = []
    tot_height = 0
    max_width = 0
    for imgs in [imgs[i:i+ncols] for i in range(0, len(imgs), ncols)]:
        tot_width = 0
        for img in imgs:
            coords.append((img,(tot_width,tot_height)))
            tot_width += int(img.width.replace('px','').replace('pt','')) + buffer
        max_width = max(max_width, tot_width)
        tot_height = max(*[int(img.height.replace('px','').replace('pt','')) for img in imgs],0)+ tot_height
    fig = sg.SVGFigure(f"{max_width}px",f"{tot_height}px")
    for img, coord in coords:
        # get the plot object
        plot1 = img.getroot()
        plot1.moveto(coord[0],coord[1])
        # append plot to figure
        fig.append([plot1])
    return fig,tot_width,tot_height


def draw_images(imgs, buffer = 1 ,ncols = 2, svg = False):
    """Merge images into a grid layout (SVG or raster).
    
    Dispatcher function that calls either merge_svg or merge_raster based on image type.
    
    Parameters
    ----------
    imgs : list
        List of images (SVG figures or PIL Images)
    buffer : int, default=1
        Spacing between images in pixels
    ncols : int, default=2
        Number of columns in grid layout
    svg : bool, default=False
        If True, treats images as SVG figures. If False, treats as PIL Images.
    
    Returns
    -------
    tuple of (figure, int, int)
        - Merged figure (SVGFigure or PIL.Image)
        - Total width
        - Total height
    """
    if svg is True:
        fig,tot_width,tot_height = merge_svg(imgs,buffer, ncols)
    else:
        fig,tot_width,tot_height = merge_raster(imgs,buffer,ncols)
    return fig,tot_width,tot_height

def plot_mol(mol,**kwargs):
    """Plot a single RDKit molecule with customizable drawing options.
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object to draw
    **kwargs : dict
        Drawing options:
        - svg (bool): If True, generate SVG. If False, generate raster image
        - subwidth (int): Image width in pixels, default=300
        - subheight (int): Image height in pixels, default=300
        - fixedBondLength (int): Bond length in pixels, default=20
        - addAtomIndices (bool): Show atom indices, default=True
        - addBondIndices (bool): Show bond indices, default=False
        - bondLineWidth (float): Bond line thickness, default=1.5
        - maxFontSize (int): Maximum font size, default=16
        - minFontSize (int): Minimum font size, default=13
        - buffer (int): Spacing for grid (single molecule, so unused), default=2
        - ncols (int): Grid columns (single molecule, so unused), default=1
    
    Returns
    -------
    tuple of (figure, int, int)
        - Figure object (SVGFigure or PIL.Image)
        - Total width
        - Total height
    
    Notes
    -----
    Uses black-and-white atom palette (useBWAtomPalette). Returns result from
    draw_images which allows consistent interface with plot_mols.
    """
    imgs = []
    if kwargs.get('svg',False) is True:
        d2d = Draw.MolDraw2DSVG(kwargs.get('subwidth',300),kwargs.get('subheight',300))
    else:
        d2d = Draw.MolDraw2DCairo(kwargs.get('subwidth',300),kwargs.get('subheight',300))
    dopts = d2d.drawOptions()
    dopts.useBWAtomPalette()
    dopts.fixedBondLength = kwargs.get('fixedBondLength',20)
    dopts.addAtomIndices = kwargs.get('addAtomIndices',True)
    dopts.addBondIndices = kwargs.get('addBondIndices',False)
    dopts.bondLineWidth = kwargs.get('bondLineWidth',1.5)
    dopts.maxFontSize = kwargs.get('maxFontSize',16)
    dopts.minFontSize = kwargs.get('minFontSize',13)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    imgs.append(d2d.GetDrawingText())
    if kwargs.get('svg',False) is True:
        imgs = [sg.fromstring(img) for img in imgs]
    else:
        imgs = [Image.open(BytesIO(img)) for img in imgs]
    return draw_images(imgs, buffer = kwargs.get('buffer',2), ncols = kwargs.get('ncols',1), svg = kwargs.get('svg',False))

def plot_mols(mols,**kwargs):
    """Plot multiple RDKit molecules in a grid with customizable drawing options.
    
    Parameters
    ----------
    mols : list of rdkit.Chem.Mol
        List of RDKit molecule objects to draw
    **kwargs : dict
        Drawing options (same as plot_mol):
        - svg (bool): If True, generate SVG. If False, generate raster images
        - subwidth (int): Width per molecule in pixels, default=300
        - subheight (int): Height per molecule in pixels, default=300
        - fixedBondLength (int): Bond length in pixels, default=20
        - addAtomIndices (bool): Show atom indices, default=True
        - addBondIndices (bool): Show bond indices, default=False
        - bondLineWidth (float): Bond line thickness, default=1.5
        - maxFontSize (int): Maximum font size, default=16
        - minFontSize (int): Minimum font size, default=13
        - buffer (int): Spacing between molecules in pixels, default=2
        - ncols (int): Number of columns in grid, default=1
    
    Returns
    -------
    tuple of (figure, int, int)
        - Merged figure (SVGFigure or PIL.Image) with all molecules
        - Total width
        - Total height
    
    Examples
    --------
    >>> from rdkit import Chem
    >>> mols = [Chem.MolFromSmiles("CCF"), Chem.MolFromSmiles("CC(F)F")]
    >>> fig, width, height = plot_mols(mols, ncols=2, svg=True)
    
    Notes
    -----
    All molecules use black-and-white atom palette. Each molecule is drawn
    independently then combined into a grid layout.
    """
    imgs = []
    for mol in mols:
        if kwargs.get('svg',False) is True:
            d2d = Draw.MolDraw2DSVG(kwargs.get('subwidth',300),kwargs.get('subheight',300))
        else:
            d2d = Draw.MolDraw2DCairo(kwargs.get('subwidth',300),kwargs.get('subheight',300))
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.fixedBondLength = kwargs.get('fixedBondLength',20)
        dopts.addAtomIndices = kwargs.get('addAtomIndices',True)
        dopts.addBondIndices = kwargs.get('addBondIndices',False)
        dopts.bondLineWidth = kwargs.get('bondLineWidth',1.5)
        dopts.maxFontSize = kwargs.get('maxFontSize',16)
        dopts.minFontSize = kwargs.get('minFontSize',13)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        imgs.append(d2d.GetDrawingText())
    if kwargs.get('svg',False) is True:
        imgs = [sg.fromstring(img) for img in imgs]
    else:
        imgs = [Image.open(BytesIO(img)) for img in imgs]
    return draw_images(imgs, buffer = kwargs.get('buffer',2), ncols = kwargs.get('ncols',1), svg = kwargs.get('svg',False))