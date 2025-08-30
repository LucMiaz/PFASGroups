from PIL import Image
import svgutils.transform as sg
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
import svgutils.transform as sg


def merge_raster(imgs,buffer,ncols):
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
    if svg is True:
        fig,tot_width,tot_height = merge_svg(imgs,buffer, ncols)
    else:
        fig,tot_width,tot_height = merge_raster(imgs,buffer,ncols)
    return fig,tot_width,tot_height

def plot_mol(mol,**kwargs):
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