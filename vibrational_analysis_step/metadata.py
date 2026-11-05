# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from VibrationalAnalysis"""

metadata = {}

"""Description of the computational models for VibrationalAnalysis.

Hamiltonians, approximations, and basis set or parameterizations,
only if appropriate for this code. For example::

    metadata["computational models"] = {
        "Hartree-Fock": {
            "models": {
                "PM7": {
                    "parameterizations": {
                        "PM7": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": True,
                            "code": "mopac",
                        },
                        "PM7-TS": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": False,
                            "code": "mopac",
                        },
                    },
                },
            },
        },
    }
"""
# metadata["computational models"] = {
# }

"""Description of the VibrationalAnalysis keywords.

(Only needed if this code uses keywords)

Fields
------
description : str
    A human readable description of the keyword.
takes values : int (optional)
    Number of values the keyword takes. If missing the keyword takes no values.
default : str (optional)
    The default value(s) if the keyword takes values.
format : str (optional)
    How the keyword is formatted in the MOPAC input.

For example::
    metadata["keywords"] = {
        "0SCF": {
            "description": "Read in data, then stop",
        },
        "ALT_A": {
            "description": "In PDB files with alternative atoms, select atoms A",
            "takes values": 1,
            "default": "A",
            "format": "{}={}",
        },
    }
"""
# metadata["keywords"] = {
# }

"""Properties that VibrationalAnalysis produces.
`metadata["results"]` describes the results that this step can produce. It is a
dictionary where the keys are the internal names of the results within this step, and
the values are a dictionary describing the result. For example::

    metadata["results"] = {
        "total_energy": {
            "calculation": [
                "energy",
                "optimization",
            ],
            "description": "The total energy",
            "dimensionality": "scalar",
            "methods": [
                "ccsd",
                "ccsd(t)",
                "dft",
                "hf",
            ],
            "property": "total energy#Psi4#{model}",
            "type": "float",
            "units": "E_h",
        },
    }

Fields
______

calculation : [str]
    Optional metadata describing what subtype of the step produces this result.
    The subtypes are completely arbitrary, but often they are types of calculations
    which is why this is name `calculation`. To use this, the step or a substep
    define `self._calculation` as a value. That value is used to select only the
    results with that value in this field.

description : str
    A human-readable description of the result.

dimensionality : str
    The dimensions of the data. The value can be "scalar" or an array definition
    of the form "[dim1, dim2,...]". Symmetric tringular matrices are denoted
    "triangular[n,n]". The dimensions can be integers, other scalar
    results, or standard parameters such as `n_atoms`. For example, '[3]',
    [3, n_atoms], or "triangular[n_aos, n_aos]".

methods : str
    Optional metadata like the `calculation` data. `methods` provides a second
    level of filtering, often used for the Hamiltionian for *ab initio* calculations
    where some properties may or may not be calculated depending on the type of
    theory.

property : str
    An optional definition of the property for storing this result. Must be one of
    the standard properties defined either in SEAMM or in this steps property
    metadata in `data/properties.csv`.

type : str
    The type of the data: string, integer, or float.

units : str
    Optional units for the result. If present, the value should be in these units.
"""
metadata["results"] = {
    "energy": {
        "description": "The energy at the structure analyzed",
        "dimensionality": "scalar",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "ZPE": {
        "description": "The anharmonic (VPT2) zero-point energy, with its constant G0",
        "dimensionality": "scalar",
        "property": "ZPE#VPT2#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".4f",
    },
    "harmonic ZPE": {
        "description": "The harmonic zero-point energy from the same Hessian",
        "dimensionality": "scalar",
        "property": "harmonic ZPE#VPT2#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".4f",
    },
    "ZPE correction": {
        "description": "The anharmonic correction to the zero-point energy",
        "dimensionality": "scalar",
        "property": "ZPE correction#VPT2#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".4f",
    },
    "G0": {
        "description": "The constant term G0 of the VPT2 zero-point energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": "1/cm",
        "format": ".3f",
    },
    "harmonic frequencies": {
        "description": "The harmonic frequencies",
        "dimensionality": "[n_modes]",
        "property": "harmonic frequencies#VPT2#{model}",
        "type": "json",
        "units": "1/cm",
        "format": ".2f",
    },
    "fundamentals": {
        "description": "The anharmonic fundamental frequencies",
        "dimensionality": "[n_modes]",
        "property": "anharmonic fundamentals#VPT2#{model}",
        "type": "json",
        "units": "1/cm",
        "format": ".2f",
    },
    "overtones": {
        "description": "The first overtones, 2 nu_i",
        "dimensionality": "[n_modes]",
        "property": "overtones#VPT2#{model}",
        "type": "json",
        "units": "1/cm",
        "format": ".2f",
    },
    "combination bands": {
        "description": "The combination bands nu_i + nu_j",
        "dimensionality": "[n_modes, n_modes]",
        "property": "combination bands#VPT2#{model}",
        "type": "json",
        "units": "1/cm",
        "format": ".2f",
    },
    "chi": {
        "description": "The anharmonicity constants",
        "dimensionality": "[n_modes, n_modes]",
        "property": "anharmonicity constants#VPT2#{model}",
        "type": "json",
        "units": "1/cm",
        "format": ".3f",
    },
    "resonances": {
        "description": "The Fermi resonances treated, as text",
        "dimensionality": "scalar",
        "type": "string",
    },
    "flagged modes": {
        "description": "The modes VPT2 should not be trusted for (1-based)",
        "dimensionality": "[n_flagged]",
        "type": "json",
    },
    "n_modes": {
        "description": "The number of vibrational modes",
        "dimensionality": "scalar",
        "type": "integer",
    },
}
