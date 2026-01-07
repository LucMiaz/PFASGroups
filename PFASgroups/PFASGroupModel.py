from rdkit import Chem
import numpy as np
class PFASGroup():
    def __init__(self, id, name, smarts1, smarts2, smartsPath, constraints,**kwargs):
        self.id = id
        self.name = name
        self.smarts1 = smarts1
        self.smarts2 = smarts2
        self.smartsPath = smartsPath
        self.max_dist_from_CF = kwargs.get('max_dist_from_CF', 0)
        if self.smarts1 !="" and self.smarts1 is not None:
            self.smarts1 = Chem.MolFromSmarts(self.smarts1)
            self.smarts1.UpdatePropertyCache()
            Chem.GetSymmSSSR(self.smarts1)
            self.smarts1.GetRingInfo().NumRings()
            if self.smarts2 !="" and self.smarts2 is not None:
                self.smarts2 = Chem.MolFromSmarts(self.smarts2)
                self.smarts2.UpdatePropertyCache()
                Chem.GetSymmSSSR(self.smarts2)
                self.smarts2.GetRingInfo().NumRings()
            else:
                self.smarts2 = None
        else:
            self.smarts1 = None
            self.smarts2 = None
        self.constraints = constraints
    def __str__(self):
        return self.name
    def constraint_gte(self, formula_dict):
        success = True
        for e,n in self.constraints.get('gte',{}).items():
            success = success and formula_dict.get(e,0)>=n
        return success
    def constraint_lte(self, formula_dict):
        success = True
        for e,n in self.constraints.get('lte',{}).items():
            success = success and formula_dict.get(e,0)<=n
        return success
    def constraint_eq(self, formula_dict):
        success = True
        for e,n in self.constraints.get('eq',{}).items():
            success = success and formula_dict.get(e,0)==n
        return success
    def constraint_only(self, formula_dict):
        success = True
        if 'only' in self.constraints.keys():
            tot = sum(formula_dict.values())
            nn = 0
            for e in self.constraints['only']:
                nn += formula_dict.get(e,0)
            success = success and tot == nn
        return success
    def constraint_rel(self, formula_dict):
        success = True
        for e,v in self.constraints.get('rel',{}).items():
            n = sum([formula_dict.get(x,0) for x in v.get('atoms',[])])
            success = success and formula_dict.get(e,0)==n/v.get('div',1)+v.get('add',0)+sum([formula_dict.get(x,0) for x in v.get('add_atoms',[])])
        return success
    def formula_dict_satisfies_constraints(self,formula_dict):
        if len(self.constraints.keys())==0:
            return True
        success = True
        process = [None,self.constraint_rel,self.constraint_only,self.constraint_eq, self.constraint_lte, self.constraint_gte]
        k = process.pop()
        while success and k is not None:
            success = k(formula_dict)
            k = process.pop()
        return success