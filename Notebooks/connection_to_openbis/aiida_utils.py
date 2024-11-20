from aiida.common.links import LinkType
import io
import tempfile
import os
from aiida import orm
import subprocess
from pathlib import Path
from typing import Optional, Tuple
import numpy as np
from ase import Atoms
from ase.io import read,write
import json
from ase.io.jsonio import encode
from aiida.common.exceptions import NotExistentAttributeError
import random

def find_bandgap(bandsdata, number_electrons=None, fermi_energy=None):
    """
    Tries to guess whether the bandsdata represent an insulator.
    This method is meant to be used only for electronic bands (not phonons)
    By default, it will try to use the occupations to guess the number of
    electrons and find the Fermi Energy, otherwise, it can be provided
    explicitely.
    Also, there is an implicit assumption that the kpoints grid is
    "sufficiently" dense, so that the bandsdata are not missing the
    intersection between valence and conduction band if present.
    Use this function with care!

    :param (float) number_electrons: (optional) number of electrons in the unit cell
    :param (float) fermi_energy: (optional) value of the fermi energy.

    :note: By default, the algorithm uses the occupations array
      to guess the number of electrons and the occupied bands. This is to be
      used with care, because the occupations could be smeared so at a
      non-zero temperature, with the unwanted effect that the conduction bands
      might be occupied in an insulator.
      Prefer to pass the number_of_electrons explicitly

    :note: Only one between number_electrons and fermi_energy can be specified at the
      same time.

    :return: (is_insulator, gap), where is_insulator is a boolean, and gap a
             float. The gap is None in case of a metal, zero when the homo is
             equal to the lumo (e.g. in semi-metals).
    """

    def nint(num):
        """
        Stable rounding function
        """
        if num > 0:
            return int(num + 0.5)
        else:
            return int(num - 0.5)

    if fermi_energy and number_electrons:
        raise EitherNumberOfElectronsOrFermiEnergyError()

    assert bandsdata.units == "eV"
    stored_bands = bandsdata.get_bands()

    if len(stored_bands.shape) == 3:
        # I write the algorithm for the generic case of having both the
        # spin up and spin down array

        # put all spins on one band per kpoint
        bands = np.concatenate(list(stored_bands), axis=1)
    else:
        bands = stored_bands

    # analysis on occupations:
    if fermi_energy is None:
        num_kpoints = len(bands)

        if number_electrons is None:
            try:
                _, stored_occupations = bandsdata.get_bands(also_occupations=True)
            except KeyError as exc:
                raise FermiEnergyOrOccupationsNotPresentError() from exc

            # put the occupations in the same order of bands, also in case of multiple bands
            if len(stored_occupations.shape) == 3:
                # I write the algorithm for the generic case of having both the
                # spin up and spin down array

                # put all spins on one band per kpoint
                occupations = np.concatenate(list(stored_occupations), axis=1)
            else:
                occupations = stored_occupations

            # now sort the bands by energy
            # Note: I am sort of assuming that I have an electronic ground state

            # sort the bands by energy, and reorder the occupations accordingly
            # since after joining the two spins, I might have unsorted stuff
            bands, occupations = (
                np.array(y)
                for y in zip(
                    *[
                        zip(*j)
                        for j in [
                            sorted(
                                zip(i[0].tolist(), i[1].tolist()),
                                key=lambda x: x[0],
                            )
                            for i in zip(bands, occupations)
                        ]
                    ]
                )
            )
            number_electrons = int(
                round(sum([sum(i) for i in occupations]) / num_kpoints)
            )

            homo_indexes = [
                np.where(np.array([nint(_) for _ in x]) > 0)[0][-1]
                for x in occupations
            ]
            if (
                len(set(homo_indexes)) > 1
            ):  # there must be intersections of valence and conduction bands
                return False, None, None, None
            else:
                homo = [_[0][_[1]] for _ in zip(bands, homo_indexes)]
                try:
                    lumo = [_[0][_[1] + 1] for _ in zip(bands, homo_indexes)]
                except IndexError as exc:
                    raise NeedMoreBandsError() from exc

        else:
            bands = np.sort(bands)
            number_electrons = int(number_electrons)

            # find the zero-temperature occupation per band (1 for spin-polarized
            # calculation, 2 otherwise)
            number_electrons_per_band = 4 - len(stored_bands.shape)  # 1 or 2
            # gather the energies of the homo band, for every kpoint
            homo = [
                i[int(number_electrons / number_electrons_per_band) - 1] for i in bands
            ]  # take the nth level
            try:
                # gather the energies of the lumo band, for every kpoint
                lumo = [
                    i[int(number_electrons / number_electrons_per_band)] for i in bands
                ]  # take the n+1th level
            except IndexError as exc:
                raise NeedMoreBandsError() from exc

        if number_electrons % 2 == 1 and len(stored_bands.shape) == 2:
            # if #electrons is odd and we have a non spin polarized calculation
            # it must be a metal and I don't need further checks
            return False, None, None, None

        # if the nth band crosses the (n+1)th, it is an insulator
        gap = min(lumo) - max(homo)
        if gap == 0.0:
            return False, 0.0, None, None
        elif gap < 0.0:
            return False, gap, None, None
        else:
            return True, gap, max(homo), min(lumo)

    # analysis on the fermi energy
    else:
        # reorganize the bands, rather than per kpoint, per energy level

        # I need the bands sorted by energy
        bands.sort()

        levels = bands.transpose()
        max_mins = [(max(i), min(i)) for i in levels]

        if fermi_energy > bands.max():
            raise FermiEnergyAndBandsEnergiesError(where="above")
        if fermi_energy < bands.min():
            raise FermiEnergyAndBandsEnergiesError(where="below")

        # one band is crossed by the fermi energy
        if any(i[1] < fermi_energy and fermi_energy < i[0] for i in max_mins):
            return False, 0.0, None, None

        # case of semimetals, fermi energy at the crossing of two bands
        # this will only work if the dirac point is computed!
        elif any(i[0] == fermi_energy for i in max_mins) and any(
            i[1] == fermi_energy for i in max_mins
        ):
            return False, 0.0, None, None
        # insulating case
        else:
            # Take the max of the band maxima below the fermi energy.
            homo = max([i[0] for i in max_mins if i[0] < fermi_energy])
            # Take the min of the band minima above the fermi energy.x
            lumo = min([i[1] for i in max_mins if i[1] > fermi_energy])

            gap = lumo - homo
            if gap <= 0.0:
                raise WrongCodeError()
            return True, gap, homo, lumo

# not used
def get_preceding_workchains(node):
    preceding_workchains = set()

    def recursive_trace(node):
        if isinstance(node, orm.WorkChainNode):
            preceding_workchains.add(node)
            # Recursively check incoming nodes for workchains
        for link in node.base.links.get_incoming().all():
            recursive_trace(link.node)

    # Start tracing from the provided node
    recursive_trace(node)

    return preceding_workchains

def get_all_preceding_main_workchains(node):
    """By MAIN workchain it is meant a workchain that is not called by a workchain"""
    main_workchains = set()
    visited = set()  # Track visited nodes to avoid infinite recursion

    # Helper function to check if a node is a MAIN workchain
    def is_main_workchain(n):
        return isinstance(n, orm.WorkChainNode) and not n.base.links.get_incoming(link_type=LinkType.CALL_WORK).all()

    # Recursive function to trace back and find all preceding MAIN workchains
    def trace_back_main_workchains(n):
        if n.pk in visited:
            return  # Avoid processing the same node again
        visited.add(n.pk)

        # Check if the current node is a MAIN workchain
        if is_main_workchain(n):
            main_workchains.add(n)

        # Regardless of whether it's a MAIN workchain, continue tracing backward
        for input_link in n.base.links.get_incoming().all():
            trace_back_main_workchains(input_link.node)

    # Start tracing back from the given node
    trace_back_main_workchains(node)

    # Return the PKs of all unique MAIN workchains found
    return [wc for wc in main_workchains]

def get_qe_output_parameters(outputs):
    parameters=[
        'energy',
        'volume',
        'fft_grid',
        'energy_xc',
        'occupations',
        'total_force',
        'energy_ewald',
        'energy_units',
        'fermi_energy',
        'forces_units',
        'stress_units',
        'energy_hartree',
        'energy_accuracy',
        'energy_smearing',
        'energy_xc_units',
        'number_of_bands',
        'smooth_fft_grid',
        'symmetries',
        'symmetries_units',
        'total_force_units',
        'energy_ewald_units',
        'fermi_energy_units',
        'inversion_symmetry',
        'lattice_symmetries',
        'number_of_k_points',
        'energy_one_electron',
        'number_of_electrons',
        'energy_hartree_units',
        'magnetization_angle1',
        'magnetization_angle2',
        'number_of_atomic_wfc',
        'number_of_symmetries',
        'energy_accuracy_units',
        'no_time_rev_operations',
        'spin_orbit_calculation',
        'non_colinear_calculation',
        'energy_one_electron_units',
        'number_of_spin_components',
        'number_of_bravais_symmetries',
    ]
    return {parameter: outputs[parameter] for parameter in parameters if parameter in outputs}

def get_qe_input_parameters(outputs):
    parameters=[
        'lsda',
         'degauss',
         'rho_cutoff',
         'wfc_cutoff',
         'smearing_type',
         'constraint_mag',
         'number_of_atoms',
         'do_magnetization',
         'energy_threshold',
         'rho_cutoff_units',
         'spin_orbit_domag',
         'wfc_cutoff_units',
         'number_of_species',
         'has_electric_field',
         'time_reversal_flag',
         'monkhorst_pack_grid',
         'energy_accuracy_units',
         'energy_smearing_units',
         'has_dipole_correction',
         'monkhorst_pack_offset',
         'lda_plus_u_calculation',
         'spin_orbit_calculation',
         'starting_magnetization',
         'dft_exchange_correlation',
         'do_not_use_time_reversal',
         'non_colinear_calculation',
         'number_of_spin_components',
    ]
    return {parameter: outputs[parameter] for parameter in parameters if parameter in outputs}

def get_dft_parameters_qe(inputs,outputs):
    """Retrieves from QE workchains the parameters needed to create the DFT object
    in input the inputs of QE workchain and the output_parameters. Will be simplified when QeAppWorkchain
    bugs for not exposing some of the outputs will be fixed
    """
    parameters={}
    parameters['code']=inputs.pw.code.description
    parameters['functional']=outputs['dft_exchange_correlation']
    parameters['plus_u'] = outputs['lda_plus_u_calculation']
    parameters['spin_orbit_couplig'] = outputs['spin_orbit_calculation']
    parameters['uks'] = outputs['lsda']
    parameters['charge']= inputs.pw.parameters['SYSTEM']['tot_charge']
    parameters['vdw_corr']=inputs.pw.parameters['SYSTEM'].get('vdw_corr','none')

    return parameters

def get_dft_parameters_cp2k(code,dft_para):
    """Retrieves from CP2K workchains teh parameters to define the DFT object. Very preliminary"""
    parameters={}
    parameters['code']=code
    parameters['functional']='PBE'
    parameters['plus_u'] = False
    parameters['spin_orbit_couplig'] = False
    parameters['uks'] = dft_para.get('uks',False)
    parameters['charge']= dft_para.get('charge',False)
    parameters['vdw_corr']=dft_para.get('vdw',False)
    
    return parameters

def geo_to_png(ase_geo):
    ase_geo.write('ase_geo.png')
    return 'ase_geo.png'
    
def guess_dimensionality(ase_geo: Optional[Atoms] = None, thr_vacuum: float = 5) -> Optional[Tuple[int, Tuple[bool, bool, bool]]]:
    """Guess the dimensionality of a structure. thr_vacuum in Å.
    returns:
    -int dimensionality
    -(Bool,Bool,Bool) PBC
    """
    if ase_geo is None or not hasattr(ase_geo, "positions") or not hasattr(ase_geo, "cell"):
        return None

    sys_size = np.ptp(ase_geo.positions, axis=0)
    cell_lengths = np.diagonal(ase_geo.cell)

    # Check if the structure has no meaningful cell (cell lengths close to 0)
    if np.any(cell_lengths < 0.1):
        return 0, (False, False, False)

    # Determine vacuum presence in each direction
    has_vacuum = sys_size + thr_vacuum < cell_lengths
    # Invert vacuum to indicate bulk presence
    has_bulk = ~has_vacuum  # Logical NOT to invert (True where no vacuum)
    dimensionality = 3 - np.sum(has_vacuum)  # 3D: no vacuum, 2D: 1 vacuum, etc.

    return dimensionality, tuple(has_bulk)

def creator_of_structure(struc):
    the_creator = struc.creator
    previous_wc = the_creator
    while previous_wc is not None:
        the_creator=previous_wc
        previous_wc = previous_wc.caller
    return the_creator

def is_structure_optimized(structure):
    creator = creator_of_structure(structure)
    geo_opt = False
    cell_opt = False
    cell_free = ''
    if creator is None:
        return geo_opt,cell_opt,cell_free
    if creator.process_label == 'Cp2kGeoOptWorkChain':
        geo_opt = True
        cell_opt = creator.label == 'CP2K_CellOpt'
        if cell_opt:
            cell_free = creator.inputs.sys_params.get_dict()['cell_opt_constraint']
        return geo_opt,cell_opt,cell_free
    if creator.process_label == 'QeAppWorkChain':
        for wc in creator.called_descendants:
            if wc.process_label == 'PwRelaxWorkChain':
                geo_opt = True
                cell_free = wc.inputs.base.pw.parameters.get_dict().get('CELL', {}).get('cell_dofree', '')
                if cell_free != '':
                    cell_opt = True
                break            
        return geo_opt,cell_opt,cell_free

# Assuming 'data' is an AiiDA Data object
def aiida_data_to_json(data):
    """Exports AiiDA xxData object as .json. Does not work for StructureData"""
    # Create a temporary file path
    temp_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.json', delete=False)
    try:
        temp_file.close()  # Close it to allow export to overwrite the file
        
        # Export the Data object to the temporary file in JSON format
        data.export(temp_file.name, fileformat='json', overwrite=True)
        
        # Read the contents of the temporary file
        with open(temp_file.name, 'r') as f:
            json_string = f.read()
    finally:
        # Clean up the temporary file
        os.remove(temp_file.name)
    
    return json_string

def create_obis_object(obtype=None,parameters=None):
    """Function to create oBIS object given parameters"""
    obisuuid = ''.join(random.choice('abcdefghijklmnopqrstuvwxyz1234567890') for _ in range(8))
    print(f'oBIS object {obtype} {obisuuid} created with parameters: {parameters.keys()}')
    return obisuuid

def get_uuids_from_oBIS():
    """TBD get from openbis list of already exported uuids we should make sure this is fast"""
    return {'wc_uuids':[],'structure_uuids':[]}

def structure_to_atomistic_model(structure,uuids):
    """Check if this atomistic model is already in OBIS otherwise create.
    Output: uuid of the oBIS object for linking
    """
    uuid = structure.uuid
    #check if th eatomistic model is already in oBIS
    if uuid in uuids:
        obisuuid=1234
        print(f'atomistic model already in oBIS with uuid {obisuuid} ')
        return obisuuid
    #check if the geometry is optimized
    ase_geo = structure.get_ase()
    dimensionality = guess_dimensionality(ase_geo)
    dictionary={
        'uuid':structure.uuid,
        'dimensionality':dimensionality[0],
        'PBC':dimensionality[1],
        'optimized':is_structure_optimized(structure),
        'structure': encode(ase_geo),
        'eln_preview':geo_to_png(ase_geo) # writes teh file ase_geo.png and returns teh filename
    }
    
    obobject = create_obis_object(obtype='atomistic_model',parameters=dictionary)
    return obobject

def create_and_export_AiiDA_archive(uuid):
    """Create archive.aiida as temporary file, with nodes only from a MAIN workchain. 
    To be sent to AiiDA_nodes object"""
    
    # Define the output file path
    output_file = "archive.aiida"

    # Define the command, converting the integer `pk` to a string
    command = [
        "verdi", "archive", "create", output_file,
        "--no-call-calc-backward", "--no-call-work-backward",
        "--no-create-backward", "-N", uuid
    ]

    try:
        # Execute the command
        result = subprocess.run(command, capture_output=True, text=True)

        # Check for errors
        if result.returncode != 0:
            print(f"An error occurred: {result.stderr}")
        else:
            # Capture the absolute path to the created file
            created_file_path = Path(output_file).resolve()
            
            # Create the AiiDA_nodes object in openBIS
            obobject = create_obis_object(obtype='aiida_nodes',parameters={'uuid':uuid,'aiida.archive':output_file,'comments':''})
                       
    finally:
        # Ensure the file is deleted
        file_to_delete = Path(output_file)
        if file_to_delete.exists():
            file_to_delete.unlink()  # Delete the file
            print(f"File {file_to_delete} has been deleted.")
        else:
            print(f"File {file_to_delete} does not exist, nothing to delete.")
    return obobject 

def PwRelaxWorkChain_export(workchain,uuids): # is SUB of QeAppWorkChain
    pw_input_parameters = workchain.inputs.base.pw.parameters.get_dict()
    output_parameters_dict = workchain.outputs.output_parameters.get_dict()

    dft_object_parameters=get_dft_parameters_qe(workchain.inputs.base,output_parameters_dict)
    # Create the DFT object in openBIS
    dft_object = create_obis_object(obtype='dft',parameters=dft_object_parameters)
    
    dictionary = {
        'uuid': workchain.uuid,
        'DFT':dft_object, #link/incorporate DFT object
        'cell_opt' : 'CELL' in pw_input_parameters,
        'force_conv_thr' : pw_input_parameters['CONTROL']['forc_conv_thr'] ,
        'constraints' : False,
        'output_parameters' : get_qe_output_parameters(workchain.outputs.output_parameters.get_dict()) ,
        'input_parameters' : get_qe_input_parameters(workchain.outputs.output_parameters.get_dict())
    }
    

    
    # create oBIS GEO_OPT object
    obobject = create_obis_object(obtype='geo_opt',parameters=dictionary)
    
    input_structure = workchain.inputs.structure
    # if missing create oBIS object and obtain uuid
    input_structure = structure_to_atomistic_model(input_structure,uuids)
    
    # TBD link ipnut_structure, that is a oBIS uuid, as parent
    print(f'atomistic model {input_structure} is parent of oBIS GEO_OPT {obobject} ')
    
    output_structure = workchain.outputs.output_structure
    # if missing create oBIS object and obtain uuid
    output_structure = structure_to_atomistic_model(output_structure,uuids)
    
    # TBD link ipnut_structure, taht is a oBIS uuid, as parent
    print(f'atomistic model {output_structure} is child of oBIS GEO_OPT {obobject} ')    
    
    return obobject

def BandsWorkChain_export(workchain,uuids): # is SUB of QeAppWorkChain
    try:
        root_in = workchain.inputs.bands
        root_out = workchain.outputs.bands
    except NotExistentAttributeError:
        root_in=workchain.inputs.bands_projwfc
        root_out = workchain.outputs.bands_projwfc
    
    dft_object_parameters=get_dft_parameters_qe(root_in.bands,root_out.scf_parameters.get_dict())
    # Create the DFT object in openBIS
    dft_object = create_obis_object(obtype='dft',parameters=dft_object_parameters)
    
    output_parameters=get_qe_output_parameters(root_out.scf_parameters.get_dict())
    input_parameters = get_qe_input_parameters(root_out.scf_parameters.get_dict())
    dictionary = {
        'uuid':workchain.uuid,
        'DFT':dft_object,
        'input_structure':workchain.inputs.structure,
        'outputs':{
            'pdos':aiida_data_to_json(root_out.projwfc.Dos),
            'pbands':aiida_data_to_json(root_out.projwfc.bands),
            'projections':aiida_data_to_json(root_out.projwfc.projections),
            'band_structure':aiida_data_to_json(root_out.band_structure)
        },
        'output_parameters' : output_parameters ,
        'input_parameters' : input_parameters,
        'bandgap':find_bandgap(root_out.band_structure, number_electrons=output_parameters['number_of_electrons'])[1]
    }
    
    # Create BANDSTRUCURE object
    obobject = create_obis_object(obtype='bands',parameters=dictionary)
    
    input_structure = workchain.inputs.structure
    # if missing create oBIS object and obtain uuid
    input_structure = structure_to_atomistic_model(input_structure,uuids)
    
    # TBD link ipnut_structure, that is a oBIS uuid, as parent
    print(f'atomistic model {input_structure} is parent of oBIS BANDS {obobject} ')
    
    return obobject

def PdosWorkChain_export(workchain,uuids): # is SUB of QeAppWorkChain    
    root_in = workchain.inputs
    root_out = workchain.outputs    

    dft_object_parameters=get_dft_parameters_qe(root_in.scf,root_out.nscf.output_parameters.get_dict())
    # Create the DFT object in openBIS
    dft_object = create_obis_object(obtype='dft',parameters=dft_object_parameters)    
    
    dictionary = {
        'uuid':workchain.uuid,
        'DFT':dft_object,
        'input_structure':workchain.inputs.structure,
        'outputs':{
            'dos': aiida_data_to_json(workchain.outputs.dos.output_dos),
            'pdos':aiida_data_to_json(workchain.outputs.projwfc.Dos)
        },
        'output_parameters' : get_qe_output_parameters(root_out.nscf.output_parameters.get_dict()) ,
        'input_parameters' : get_qe_input_parameters(root_out.nscf.output_parameters.get_dict())        
    }  
    
    # create oBIS PDOS object
    obobject = create_obis_object(obtype='pdos',parameters=dictionary)
    
    input_structure = workchain.inputs.structure
    # if missing create oBIS object and obtain uuid
    input_structure = structure_to_atomistic_model(input_structure,uuids)
    # TBD link ipnut_structure, that is a oBIS uuid, as parent
    print(f'atomistic model {input_structure} is parent of oBIS PDOS {obobject} ')
    
    return obobject

def VibroWorkChain_export(workchain,uuids): # is SUB of QeAppWorkChain
    # outputs are not exported so we look for a PwBaseWorkChain
    for wkc in workchain.called_descendants:
        if wkc.process_label == 'PwBaseWorkChain':
            root_in = wkc.inputs
            root_out = wkc.outputs
            break
    
    dft_object_parameters=get_dft_parameters_qe(root_in,root_out.output_parameters.get_dict())
    # Create the DFT object in openBIS
    dft_object = create_obis_object(obtype='dft',parameters=dft_object_parameters)    
        
    dictionary = {
        'uuid':workchain.uuid,
        'DFT':dft_object,
        'input_structure':workchain.inputs.structure,
        'outputs':{
            'phonon_bands':aiida_data_to_json(workchain.outputs.phonon_bands),
            'phonon_pdos':aiida_data_to_json(workchain.outputs.phonon_pdos),
            'phonon_thermo': aiida_data_to_json(workchain.outputs.phonon_thermo),
        }        
    } 
    
    # Create VIBSPEC object
    obobject = create_obis_object(obtype='vibspec',parameters=dictionary)  
    
    input_structure = workchain.inputs.structure
    # if missing create oBIS object and obtain uuid
    input_structure = structure_to_atomistic_model(input_structure,uuids)    
    # TBD link ipnut_structure, that is a oBIS uuid, as parent
    print(f'atomistic model {input_structure} is parent of oBIS VIBSPEC {obobject} ')    
    
    return obobject

def Cp2kGeoOptWorkChain_export(workchain,uuids): # Can be both MAIN and SUB. Do not export in case is SUB
    sys_params = workchain.inputs.sys_params.get_dict()
    dft_params = workchain.inputs.dft_params.get_dict()
    code=workchain.inputs.code.description
    output_parameters = workchain.outputs.output_parameters.get_dict()
    input_parameters = {}
    try:
        input_parameters = workchain.outputs.final_input_parameters.get_dict()
    except NotExistentAttributeError:
        pass
    input_structure = workchain.inputs.structure
    
    
    dft_object_parameters=get_dft_parameters_cp2k(code,dft_params)
    # Create the DFT object in openBIS
    dft_object = create_obis_object(obtype='dft',parameters=dft_object_parameters)    
    
    dictionary = {
        'uuid': workchain.uuid,
        'DFT':dft_object, #link/incorporate DFT object
        'cell_opt' : 'cell' in sys_params,
        'force_conv_thr' : '' ,
        'constraints' : sys_params['constraints'] != "",
        'output_parameters' : output_parameters ,
        'input_parameters' : input_parameters
    }
    
    # create oBIS GEO_OPT object
    obobject = create_obis_object(obtype='geo_opt',parameters=dictionary)   
    
    input_structure = workchain.inputs.structure
    # if missing create oBIS object and obtain uuid
    input_structure = structure_to_atomistic_model(input_structure,uuids)
    
    # TBD link ipnut_structure, that is a oBIS uuid, as parent
    print(f'atomistic model {input_structure} is parent of oBIS GEO_OPT {obobject} ')
    
    output_structure = workchain.outputs.output_structure
    # if missing create oBIS object and obtain uuid
    output_structure = structure_to_atomistic_model(output_structure,uuids)
    
    # TBD link ipnut_structure, taht is a oBIS uuid, as parent
    print(f'atomistic model {output_structure} is child of oBIS GEO_OPT {obobject} ')    
        
    # TBD output trajectory....
    #output_trajectory = aiida_data_to_json(workchain.outputs.output_trajectory)    
    
    return obobject

def Cp2kStmWorkChain_export(workchain,uuids):
    dft_params = workchain.inputs.dft_params.get_dict()
    spm_params = workchain.inputs.spm_params.get_dict()
    cp2k_code = workchain.inputs.cp2k_code.description
    spm_code=workchain.inputs.spm_code.description
    output_parameters = {}
    try:
        output_parameters = workchain.outputs.dft_output_parameters.get_dict()
    except NotExistentAttributeError:
        pass
    input_parameters = {}
    
    dft_object_parameters=get_dft_parameters_cp2k(spm_code,dft_params)
    # Create the DFT object in openBIS
    dft_object = create_obis_object(obtype='dft',parameters=dft_object_parameters)     
        
    dictionary = {
        'uuid': workchain.uuid,
        'DFT':dft_object, #link/incorporate DFT object
        'Bias':spm_params['--energy_range'],
        'Iso':spm_params['--isovalues'],
        'Height':spm_params['--heights'],
        'p-tip':spm_params['--p_tip_ratios'],
        'output_parameters' : output_parameters ,
        'input_parameters' : input_parameters,
    }
    
    # create oBIS GEO_OPT object
    obobject = create_obis_object(obtype='spm',parameters=dictionary)
    
    input_structure = workchain.inputs.structure
    # if missing create oBIS object and obtain uuid
    input_structure = structure_to_atomistic_model(input_structure,uuids)
    
    # TBD link ipnut_structure, that is a oBIS uuid, as parent
    print(f'atomistic model {input_structure} is parent of oBIS GEO_OPT {obobject} ')   
    return obobject

workchain_exporters={
    'PwRelaxWorkChain':PwRelaxWorkChain_export,
    'BandsWorkChain':BandsWorkChain_export,
    'PdosWorkChain':PdosWorkChain_export,
    'VibroWorkChain':VibroWorkChain_export,
    'Cp2kGeoOptWorkChain':Cp2kGeoOptWorkChain_export,
    'Cp2kStmWorkChain':Cp2kStmWorkChain_export,
}

def export_workchain(workchain):
    workchains_to_export = get_all_preceding_main_workchains(workchain)
    
    # check available uuids in oBIS form AiiDA_nodes objects
    list_uuids_oBIS = get_uuids_from_oBIS()
    
    #create individual oBIS objects 
    for main_wc in workchains_to_export:
        if main_wc.is_finished_ok: # if not we do not parse it but it will still be in the AiiDA archive
            if main_wc.uuid not in list_uuids_oBIS['wc_uuids']:
                print(f'dealing with main WC {main_wc.pk}')
                
                #create global .aiida for main_wc and AiiDA_nodes openBIS object with the archive as dataset
                AiiDA_archive=create_and_export_AiiDA_archive(main_wc.uuid)
                if main_wc.process_label in workchain_exporters:
                        #check if wc.uuid already in openBIS
                        #if not in openbis create pertinent oBIS object and populate it
                        export=workchain_exporters[main_wc.process_label](main_wc,list_uuids_oBIS['structure_uuids'])
                        # TBD create parent link
                        print(f'AiiDA_archive {AiiDA_archive} is parent of {export}')
                else:
                    print(main_wc.process_label,' checking sub_workchains')        
                    #all workchains called by teh main workchain
                    wc_tree = main_wc.called_descendants
                    for wc in wc_tree:
                        #we export only properties related workchains: geo_opt, bands,...
                        if wc.process_label in workchain_exporters:
                            #check if wc.uuid already in openBIS
                            #if not in openbis create pertinent oBIS object and populate it
                            export=workchain_exporters[wc.process_label](wc,list_uuids_oBIS['structure_uuids'])
                            # TBD create parent link
                            print(f'AiiDA_archive {AiiDA_archive} is parent of {export}')
                        #else:
                        #    print(wc.process_label,' should not be exported')            
        
    return