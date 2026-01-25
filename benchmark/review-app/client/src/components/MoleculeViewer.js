import React, { useEffect, useRef, useState } from 'react';

function MoleculeViewer({ smiles, width = 300, height = 200 }) {
  const containerRef = useRef(null);
  const [error, setError] = useState(null);
  const [rdkit, setRdkit] = useState(null);

  useEffect(() => {
    const loadRDKit = async () => {
      try {
        // Load RDKit from npm package
        const RDKit = await import('@rdkit/rdkit');
        // The package exports initRDKitModule in different ways
        const initRDKitModule = RDKit.default || RDKit.initRDKitModule || RDKit;
        
        // Call the initialization function
        let rdkitInstance;
        if (typeof initRDKitModule === 'function') {
          rdkitInstance = await initRDKitModule();
        } else if (typeof initRDKitModule.initRDKitModule === 'function') {
          rdkitInstance = await initRDKitModule.initRDKitModule();
        } else {
          throw new Error('Could not find initRDKitModule function');
        }
        
        setRdkit(rdkitInstance);
      } catch (err) {
        console.error('Error loading RDKit:', err);
        setError('RDKit.js not available. Showing SMILES string instead.');
      }
    };

    loadRDKit();
  }, []);

  useEffect(() => {
    if (rdkit && smiles && containerRef.current) {
      try {
        // Clear previous content
        containerRef.current.innerHTML = '';

        // Create molecule from SMILES
        const mol = rdkit.get_mol(smiles);
        if (mol && mol.is_valid()) {
          // Generate SVG directly
          const svg = mol.get_svg_with_highlights(JSON.stringify({
            width: width,
            height: height
          }));
          
          // Inject SVG directly into container
          containerRef.current.innerHTML = svg;
          
          // Clean up
          mol.delete();
          setError(null);
        } else {
          setError('Invalid SMILES string');
          if (mol) mol.delete();
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
      <div
        ref={containerRef}
        style={{
          border: '1px solid #ddd',
          borderRadius: '4px',
          backgroundColor: 'white',
          display: 'inline-block',
          minWidth: width,
          minHeight: height
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