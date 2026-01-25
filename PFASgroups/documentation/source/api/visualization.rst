Visualization Module
====================

This module provides functions for drawing and visualizing molecules.

.. module:: PFASgroups.draw_mols
   :synopsis: Molecular visualization utilities

Overview
--------

The visualization module provides utilities for:

- Drawing individual molecules with annotations
- Creating grid layouts of multiple molecules
- Generating both raster (PNG) and vector (SVG) output
- Highlighting specific atoms and bonds

Main Functions
--------------

plot_mol
^^^^^^^^

.. py:function:: plot_mol(mol, **kwargs)

   Draw a single molecule with customizable options.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param svg: Generate SVG output instead of PNG
   :type svg: bool, optional
   :param subwidth: Image width in pixels (default: 300)
   :type subwidth: int, optional
   :param subheight: Image height in pixels (default: 300)
   :type subheight: int, optional
   :param addAtomIndices: Add atom index labels (default: True)
   :type addAtomIndices: bool, optional
   :param addBondIndices: Add bond index labels (default: False)
   :type addBondIndices: bool, optional
   :param fixedBondLength: Fixed bond length in pixels (default: 20)
   :type fixedBondLength: int, optional
   :param bondLineWidth: Bond line width (default: 1.5)
   :type bondLineWidth: float, optional
   :param maxFontSize: Maximum font size (default: 16)
   :type maxFontSize: int, optional
   :param minFontSize: Minimum font size (default: 13)
   :type minFontSize: int, optional
   :param buffer: Pixel buffer between images (default: 2)
   :type buffer: int, optional
   :param ncols: Number of columns in grid (default: 1)
   :type ncols: int, optional
   :returns: Tuple of (figure, width, height)
   :rtype: tuple

   **Example:**

   .. code-block:: python

      from PFASgroups import plot_mol
      from rdkit import Chem

      mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
      
      # Basic plot
      fig, w, h = plot_mol(mol)

      # SVG output with larger size
      fig, w, h = plot_mol(mol, svg=True, subwidth=400, subheight=400)

      # Without atom indices
      fig, w, h = plot_mol(mol, addAtomIndices=False)

plot_mols
^^^^^^^^^

.. py:function:: plot_mols(mols, **kwargs)

   Draw multiple molecules in a grid layout.

   :param mols: List of RDKit molecule objects
   :type mols: list[rdkit.Chem.Mol]
   :param kwargs: Same options as plot_mol
   :returns: Tuple of (figure, width, height)
   :rtype: tuple

   **Example:**

   .. code-block:: python

      from PFASgroups import plot_mols
      from rdkit import Chem

      mols = [
          Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O"),
          Chem.MolFromSmiles("FC(F)(F)C(F)(F)S(=O)(=O)O"),
          Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)O")
      ]
      
      # Grid layout
      fig, w, h = plot_mols(mols, ncols=2)

      # SVG output
      fig, w, h = plot_mols(mols, svg=True, ncols=3)

draw_images
^^^^^^^^^^^

.. py:function:: draw_images(imgs, buffer=1, ncols=2, svg=False)

   Combine multiple images into a grid layout.

   :param imgs: List of images (PIL Image or SVG objects)
   :type imgs: list
   :param buffer: Pixel buffer between images
   :type buffer: int, optional
   :param ncols: Number of columns in grid
   :type ncols: int, optional
   :param svg: Whether images are SVG format
   :type svg: bool, optional
   :returns: Tuple of (combined_figure, total_width, total_height)
   :rtype: tuple

Helper Functions
----------------

merge_raster
^^^^^^^^^^^^

.. py:function:: merge_raster(imgs, buffer, ncols)

   Merge raster images (PIL) into a grid.

   :param imgs: List of PIL Image objects
   :type imgs: list[PIL.Image.Image]
   :param buffer: Pixel buffer between images
   :type buffer: int
   :param ncols: Number of columns
   :type ncols: int
   :returns: Tuple of (merged_image, width, height)
   :rtype: tuple

merge_svg
^^^^^^^^^

.. py:function:: merge_svg(imgs, buffer, ncols)

   Merge SVG images into a grid.

   :param imgs: List of SVG figure objects
   :type imgs: list
   :param buffer: Pixel buffer between images
   :type buffer: int
   :param ncols: Number of columns
   :type ncols: int
   :returns: Tuple of (merged_svg, width, height)
   :rtype: tuple

PFAS-Specific Visualization
---------------------------

plot_pfasgroups
^^^^^^^^^^^^^^^

For PFAS-specific visualization with group highlighting, use the main
``plot_pfasgroups`` function from the core module:

.. code-block:: python

   from PFASgroups import plot_pfasgroups

   # Visualize PFAS groups
   plot_pfasgroups("FC(F)(F)C(F)(F)C(=O)O")

   # Multiple molecules
   plot_pfasgroups([
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O"
   ], ncols=2)

   # Save to file
   plot_pfasgroups(smiles, path="output.png")

   # SVG output
   plot_pfasgroups(smiles, svg=True, path="output.svg")

See :doc:`core` for full documentation of ``plot_pfasgroups``.

Output Formats
--------------

Raster (PNG)
^^^^^^^^^^^^

Default output format using RDKit's Cairo renderer:

.. code-block:: python

   fig, w, h = plot_mol(mol, svg=False)
   
   # fig is a PIL Image object
   fig.save("molecule.png")

Vector (SVG)
^^^^^^^^^^^^

SVG output for publication-quality graphics:

.. code-block:: python

   fig, w, h = plot_mol(mol, svg=True)
   
   # fig is an svgutils SVGFigure object
   fig.save("molecule.svg")

Jupyter Notebook
^^^^^^^^^^^^^^^^

For Jupyter notebook display:

.. code-block:: python

   from PFASgroups import plot_mol
   from rdkit import Chem
   from IPython.display import display, SVG

   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   
   # Display PNG
   fig, w, h = plot_mol(mol)
   display(fig)

   # Display SVG
   fig, w, h = plot_mol(mol, svg=True)
   display(SVG(fig.to_str()))

Customization Examples
----------------------

Custom Colors
^^^^^^^^^^^^^

.. code-block:: python

   from rdkit.Chem import Draw
   from rdkit import Chem

   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   
   # Use RDKit's drawing options directly
   d2d = Draw.MolDraw2DCairo(400, 400)
   dopts = d2d.drawOptions()
   
   # Custom atom colors
   dopts.setAtomPalette({
       6: (0.2, 0.2, 0.2),   # Carbon: dark gray
       9: (0.0, 0.8, 0.0),   # Fluorine: green
       8: (1.0, 0.0, 0.0),   # Oxygen: red
   })
   
   d2d.DrawMolecule(mol)
   d2d.FinishDrawing()

Highlighting Atoms
^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from rdkit.Chem import Draw
   from rdkit import Chem

   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   
   # Highlight specific atoms (e.g., fluorine atoms)
   fluorine_indices = [i for i, a in enumerate(mol.GetAtoms()) 
                       if a.GetSymbol() == 'F']
   
   img = Draw.MolToImage(mol, highlightAtoms=fluorine_indices)
   img.save("highlighted.png")
