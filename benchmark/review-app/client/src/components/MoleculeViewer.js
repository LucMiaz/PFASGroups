import React, { useEffect, useRef, useState } from 'react';

function MoleculeViewer({ smiles, width = 300, height = 200 }) {
  const canvasRef = useRef(null);
  const [error, setError] = useState(null);
  const [rdkit, setRdkit] = useState(null);

  useEffect(() => {
    const loadRDKit = async () => {
      try {
        // Try to load RDKit.js
        const RDKitModule = await import('@rdkit/rdkit');
        const rdkitInstance = await RDKitModule.initRDKitModule();
        setRdkit(rdkitInstance);
      } catch (err) {
        console.error('Error loading RDKit:', err);
        setError('RDKit.js not available. Showing SMILES string instead.');
      }
    };

    loadRDKit();
  }, []);

  useEffect(() => {
    if (rdkit && smiles && canvasRef.current) {
      try {
        // Clear previous drawing
        const canvas = canvasRef.current;
        const ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Create molecule from SMILES
        const mol = rdkit.get_mol(smiles);
        if (mol) {
          // Generate 2D coordinates
          mol.normalize_depiction();
          
          // Draw molecule to canvas
          const svg = mol.get_svg(width, height);
          
          // Convert SVG to canvas
          const img = new Image();
          const svgBlob = new Blob([svg], { type: 'image/svg+xml;charset=utf-8' });
          const url = URL.createObjectURL(svgBlob);
          
          img.onload = () => {
            ctx.drawImage(img, 0, 0);
            URL.revokeObjectURL(url);
          };
          
          img.src = url;
          
          // Clean up
          mol.delete();
          setError(null);
        } else {
          setError('Invalid SMILES string');
        }
      } catch (err) {
        console.error('Error drawing molecule:', err);
        setError('Error rendering molecule');
      }
    }
  }, [rdkit, smiles, width, height]);

  if (error) {
    return (
      <div 
        style={{ 
          width, 
          height, 
          display: 'flex', 
          alignItems: 'center', 
          justifyContent: 'center',
          border: '1px solid #ddd',
          borderRadius: '4px',
          backgroundColor: '#f8f9fa',
          flexDirection: 'column',
          padding: '10px'
        }}
      >
        <div style={{ fontSize: '12px', color: '#6c757d', marginBottom: '10px' }}>
          {error}
        </div>
        <div style={{ 
          fontFamily: 'monospace', 
          fontSize: '11px', 
          wordBreak: 'break-all',
          textAlign: 'center'
        }}>
          {smiles}
        </div>
      </div>
    );
  }

  if (!rdkit) {
    return (
      <div 
        style={{ 
          width, 
          height, 
          display: 'flex', 
          alignItems: 'center', 
          justifyContent: 'center',
          border: '1px solid #ddd',
          borderRadius: '4px',
          backgroundColor: '#f8f9fa'
        }}
      >
        Loading molecular viewer...
      </div>
    );
  }

  return (
    <div style={{ textAlign: 'center' }}>
      <canvas
        ref={canvasRef}
        width={width}
        height={height}
        style={{
          border: '1px solid #ddd',
          borderRadius: '4px',
          backgroundColor: 'white'
        }}
      />
      <div style={{ 
        fontSize: '10px', 
        color: '#6c757d', 
        marginTop: '5px',
        fontFamily: 'monospace'
      }}>
        {smiles}
      </div>
    </div>
  );
}

export default MoleculeViewer;