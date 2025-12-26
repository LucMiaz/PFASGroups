import React, { useState } from 'react';

function MoleculeViewer({ smiles, width = 300, height = 200 }) {
  const [imageError, setImageError] = useState(false);
  
  // Use PubChem's structure rendering service - reliable and fast
  const imageUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/PNG?image_size=${width}x${height}`;
  
  if (imageError) {
    return (
      <div style={{
        width: width,
        height: height,
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        justifyContent: 'center',
        backgroundColor: '#f8f9fa',
        border: '1px solid #dee2e6',
        borderRadius: '8px',
        padding: '20px',
        fontFamily: 'monospace',
        fontSize: '14px',
        wordBreak: 'break-all',
        textAlign: 'center'
      }}>
        <div style={{ marginBottom: '10px', fontWeight: 'bold', color: '#6c757d' }}>
          SMILES:
        </div>
        <div style={{ color: '#212529' }}>
          {smiles}
        </div>
      </div>
    );
  }
  
  return (
    <div style={{
      width: width,
      height: height,
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      backgroundColor: 'white',
      border: '1px solid #dee2e6',
      borderRadius: '8px',
      overflow: 'hidden',
      padding: '5px'
    }}>
      <img
        src={imageUrl}
        alt={`Molecule: ${smiles}`}
        style={{
          maxWidth: '100%',
          maxHeight: '100%',
          objectFit: 'contain'
        }}
        onError={() => setImageError(true)}
      />
    </div>
  );
}

export default MoleculeViewer;
