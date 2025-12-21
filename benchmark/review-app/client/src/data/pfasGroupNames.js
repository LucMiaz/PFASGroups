// PFAS Group Names mapping (ID to Name)
export const PFAS_GROUP_NAMES = {
  1: "Perfluoroalkyl carboxylic acids (PFCAs)",
  2: "Polyfluoroalkyl carboxylic acid (PolyFCAs)",
  3: "Perfluoroalkyl dicarboxylic acids (PFdiCAs)",
  4: "Perfluoroalkylether carboxylic acids (PFECAs)",
  5: "Polyfluoroalkylether carboxylic acid (PolyFECAs)",
  6: "Perfluoroalkyl sulfonic acids (PFSAs)",
  7: "Polyfluoroalkyl sulfonic acid (PolyFSAs)",
  8: "Perfluoroalkyl disulfonic acids (PFdiSAs)",
  9: "Perfluoroalkylether sulfonic acids (PFESAs)",
  10: "Polyfluoroalkylether sulfonic acid (PolyFESAs)",
  11: "Perfluoroalkyl phosphonic acids (PFPAs)",
  12: "Polyfluoroalkyl phosphonic acid (PolyFPAs)",
  13: "Perfluoroalkane sulfonamides (FASAs)",
  14: "Polyfluoroalkane sulfonamides (PolyFASAs)",
  15: "Fluorotelomer alcohols (FTOHs)",
  16: "Fluorotelomer iodides",
  17: "Fluorotelomer acrylates",
  18: "Fluorotelomer methacrylates",
  19: "Per- and polyfluoroalkyl aldehydes (PFALs)",
  20: "Perfluoroalkane sulfonyl fluorides (PASFs)",
  21: "Perfluoroalkyl sulfonyl chlorides",
  22: "Perfluoroalkyl acid fluorides",
  23: "Per- and polyfluoroalkyl ketones",
  24: "Semifluorinated n-alkanes",
  25: "Fluorinated cycloalkanes",
  26: "Fluorinated aromatics",
  27: "Side-chain fluorinated aromatic polymers (SCFPs)",
  28: "Per- and polyfluoroalkyl amines",
  29: "Per- and polyfluoroalkyl amides",
  30: "Per- and polyfluoroalkyl imines",
  31: "Perfluoropolyethers (PFPEs)",
  32: "Polyfluoropolyethers (PolyFPEs)",
  33: "Branched perfluoroalkyl carboxylic acids",
  34: "Branched perfluoroalkyl sulfonic acids",
  35: "Cyclic perfluoroalkyl carboxylic acids",
  36: "Cyclic perfluoroalkyl sulfonic acids",
  37: "Perfluoroalkyl thiols",
  38: "Perfluoroalkyl disulfides",
  39: "Perfluoroalkyl nitro compounds",
  40: "Perfluoroalkyl nitriles",
  41: "Fluorinated silanes",
  42: "Fluorinated siloxanes",
  43: "Fluorinated stannanes",
  44: "Perfluorinated ethers (cyclic)",
  45: "Polyfluorinated ethers (cyclic)",
  46: "Fluorotelomer sulfonates",
  47: "Fluorotelomer sulfonic acids",
  48: "Fluorotelomer carboxylates",
  49: "Per- and polyfluoroalkyl betaines",
  50: "Per- and polyfluoroalkyl ureas"
};

export function getGroupName(groupId) {
  return PFAS_GROUP_NAMES[groupId] || `Group ${groupId}`;
}

export function getShortGroupName(groupId) {
  const fullName = PFAS_GROUP_NAMES[groupId];
  if (!fullName) return `G${groupId}`;
  
  // Extract alias if present (text in parentheses)
  const aliasMatch = fullName.match(/\(([^)]+)\)$/);
  if (aliasMatch) {
    return aliasMatch[1];
  }
  
  // Otherwise truncate if too long
  if (fullName.length > 35) {
    return fullName.substring(0, 32) + '...';
  }
  
  return fullName;
}
