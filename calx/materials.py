import xrayutilities as xu
import xraylib_np
from . import consts as _consts
from scipy.constants import torr,bar,k,N_A,R
import numpy as np

# This module holds relevant materials of the
# xrayutilities materials class,

class MaterialCollection:
    """ Dummy class collections of materials (dict-like)."""
    def __init__(self, **entries):
        self.__dict__.update(entries)
    def __setitem__(self,key,value):
        self.__dict__.update({key:value})

_amorphous = dict()
_crystal = dict()
_gas = dict()

amorphous = MaterialCollection()
crystal = MaterialCollection()
gas = MaterialCollection()

def _get_transmission(self,d,E='config'):
    """ calculate the transmittion after thickness d (in m) of material at energy E (in eV)."""
    return np.exp(-d*1e6/self.absorption_length(E))

xu.materials.Material.transmission = _get_transmission

def _get_CS_Photo(self,E=None):
  """Returns the total photoabsorption cross section in m^2
     ID is the element symbol  E is the photon energy in eV
     NOTE: This is per molecule if chemical formula given
  """
  E = np.atleast_1d(E)
  CS=0
  if hasattr(self,'base') and self.base:
      base = self.base
  elif hasattr(self,'lattice'):
      base = []
      for atom,position,amount,dontknow in self.lattice.base():
          base.append((atom,amount))
  else:
      raise Exception('Did not find a base object!')

  for atom,amount in base:
    CS+=xraylib_np.CS_Photo(np.atleast_1d(atom.num),E)*amount*atom.weight
  return CS

xu.materials.Material.cs_photo = _get_CS_Photo


crystal['Si'] = xu.materials.Si
crystal['Ge'] = xu.materials.Ge
crystal['GaAs'] = xu.materials.GaAs
crystal['Al'] = xu.materials.Al
crystal['Diamond'] = xu.materials.C
crystal['Be'] = xu.materials.material.Crystal("Be", \
        xu.materials.spacegrouplattice.SGLattice(\
        194, 2.2858, 3.5843, atoms=[xu.materials.elements.Be, ], \
        pos=['2c', ]))
# crystal['Al2O3'] = xu.materials.material.Crystal("Al2O3", \
#     xu.materials.spacegrouplattice.SGLattice(167, 0.4759026, 1.299084, atoms=[xu.materials.elements.Al, xu.materials.elements.O],pos=['18e','12c'])

amorphous['B4C'] = xu.materials.material.Amorphous('B4C',2520,[('B',4),('C',1)])
amorphous['Mo'] = xu.materials.material.Amorphous('Mo',10220,[('Mo',1)])
amorphous['polyimide'] = xu.materials.material.Amorphous('polyimide',1430,[('C',22),('H',10),('N',2),('O',5)])
amorphous['mylar'] = xu.materials.material.Amorphous('mylar',1400,[('C',10),('H',8),('O',4)])
amorphous['polycarbonate'] = xu.materials.material.Amorphous('polycarbonate',1200,[('C',16),('H',14),('O',3)])
amorphous['Si3N4'] = xu.materials.material.Amorphous('Silicon nitride',3440,[('Si',3),('N',4)])
# amorphous['air'] = xu.materials.material.Amorphous('air',1000,[('N',1.562),('O',.42),('C',.0003),('Ar',.0094)])
amorphous['SiC'] = xu.materials.material.Amorphous('Silicon carbide',3210,[('Si',1),('C',1)])


# more useful values and constants
#elementName = DummyClassDict(_consts.elementName)
#meltPoint = DummyClassDict(_consts.meltPoint)
#density = DummyClassDict(_consts.Density)

class Gas(xu.materials.material.Amorphous):
    def __init__(self,name, pressure=bar, temperature=295, molecule_size=1, atoms=None, cij=None):
        """pressure in Pascal, temperature in Kelvin"""
        self.pressure = pressure
        self.temperature  = temperature
        self.molecule_size = molecule_size
        super(Gas,self).__init__(name,0,atoms=atoms,cij=cij)


    def _getdensity(self):
        """
        calculates the mass density of an material from the atomic composition and the average molecule size (ideal gas).

        Returns
        -------
        mass density in kg/m^3
        """
        num_dens = self.pressure/k/self.temperature
        return self._get_composition_mass()*num_dens*self.molecule_size
    density = property(_getdensity)

    def _get_composition_mass(self):
        w = 0
        for atom,occ in self.base:
            w += atom.weight * occ
        return w

gas['air'] = Gas('air',molecule_size=1.9917,atoms=[('N',1.562),('O',.42),('C',.0003),('Ar',.0094)])
gas['He'] = Gas('He',molecule_size=1,atoms=[('He',1)])
gas['N'] = Gas('N',molecule_size=2,atoms=[('N',1)])



def get_chem_formula(material_base):
    """creates string representation from a chemical formula inside a Material class instance."""
    string_rep = ''
    for atom,amount in material_base:
        string_rep += atom.basename + str(amount)
    return string_rep


def get_CS_Photo(ID,E=None):
  """Returns the total photoabsorption cross section in m^2
     ID is the element symbol
     E is the photon energy (default is current LCLS value)
     NOTE: This is per molecule if chemical formula given
  """
  ID=checkID(ID)
  E = getE(E)/1000.
  form=periodictable.formulas.parse_formula(ID)
  CS=0
  for id,atoms in form.atoms.items():
    CS+=xraylib.CS_Photo(elementZ[id.symbol],E)*atoms*AtomicMass[id.symbol]/c['NA']/u['cm']**2
  return CS

def Dose(ID,FWHM,J=None,E=None):
    """ Computes the absorbed dose in a solid in eV/atom (photoabsoption cross sections)
        ID is chemical fomula : 'Si'
        FWHM is the FEL spot size
        J is the LCLS pulse energy in Joules (default is LCLS value)
        E is photon energy in eV or keV (default is LCLS value)
        If no density is specified will use default value
    """
    ID=checkID(ID)
    E = getE(E)/1000.
#    if J==None:
#      J=pypsepics.get("SIOC:SYS0:ML00:AO627")/1000.
    wo=FWHM/1.18
    dose=2*J*CS_Photo(ID,E)*u['eV']/n.pi/wo**2/nAtoms(ID)
    return dose
