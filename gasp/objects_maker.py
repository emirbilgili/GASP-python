# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Objects Maker module:

This module contains functions for creating the (singleton) objects
used by the genetic algorithm during the search.

"""

from gasp import general
from gasp import population
from gasp import geometry as geo
from gasp import variations
from gasp import energy_calculators
from gasp import organism_creators
from gasp import development

from pymatgen.core.structure import Structure

import os
import math


def make_objects(parameters):
    """
    Constructs the needed objects for the genetic algorithm search.

    Returns a dictionary containing the objects.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file
    """

    # to hold all the objects
    objects_dict = {}

    # make the geometry object
    if 'Geometry' not in parameters:
        geometry = geo.Bulk()
    elif parameters['Geometry'] in (None, 'default'):
        geometry = geo.Bulk()
    elif 'shape' not in parameters['Geometry']:
        geometry = geo.Bulk()
    elif parameters['Geometry']['shape'] == 'cluster':
        geometry = geo.Cluster(parameters['Geometry'])
    elif parameters['Geometry']['shape'] == 'wire':
        geometry = geo.Wire(parameters['Geometry'])
    elif parameters['Geometry']['shape'] == 'sheet':
        geometry = geo.Sheet(parameters['Geometry'])
    elif parameters['Geometry']['shape'] == 'interface':
        geometry = geo.Substrate_2D(parameters['Geometry'])
    # TODO: add any other non-bulk geometries here
    else:
        geometry = geo.Bulk()

    objects_dict['geometry'] = geometry
    substrate_search = False
    if geometry.shape == 'interface':
        substrate_search = True

    # make the composition space object
    if 'CompositionSpace' in parameters:
        composition_space = general.CompositionSpace(
            parameters['CompositionSpace'], sub_search=substrate_search)
    else:
        print('Input file must contain a "CompositionSpace" block.')
        print("Quitting...")
        quit()

    objects_dict['composition_space'] = composition_space

    # make the constraints object
    if 'Constraints' in parameters:
        cons_min_num = 2
        constraints_parameters = parameters['Constraints']
        if substrate_search:
            cons_min_num = 1
            if 'sub_search' not in constraints_parameters:
                constraints_parameters['sub_search'] = True
        if 'min_num_atoms' in constraints_parameters:
            if constraints_parameters['min_num_atoms'] < cons_min_num:
                print('The value passed to the "min_num_atoms" keyword in the '
                      'Constraints block must greater than or equal to 2.')
                print('Quitting...')
                quit()
        constraints = development.Constraints(constraints_parameters,
                                              composition_space)
    else:
        constraints = development.Constraints('default', composition_space)

    objects_dict['constraints'] = constraints

    # make the development object
    if 'Development' in parameters:
        developer = development.Developer(parameters['Development'], geometry)
    else:
        developer = development.Developer('default', geometry)

    objects_dict['developer'] = developer

    # make the redundancy guard object
    if 'RedundancyGuard' in parameters:
        redundancy_guard = development.RedundancyGuard(
            parameters['RedundancyGuard'], geometry)
    else:
        redundancy_guard = development.RedundancyGuard('default', geometry)

    objects_dict['redundancy_guard'] = redundancy_guard

    # make the id generator
    id_generator = general.IDGenerator()
    objects_dict['id_generator'] = id_generator

    # make the organism creators
    initial_organism_creators = make_organism_creators(
        parameters, composition_space, constraints)

    # if more than one organism creator, sort them so that the attempts-based
    # ones are at the front and the successes-based ones are at the back
    if len(initial_organism_creators) > 1:
        initial_organism_creators.sort(key=lambda x: x.is_successes_based)

    objects_dict['organism_creators'] = initial_organism_creators

    # the number of energy calculations to run at a time
    if 'NumCalcsAtOnce' not in parameters:
        num_calcs_at_once = 1
    elif parameters['NumCalcsAtOnce'] in (None, 'default'):
        num_calcs_at_once = 1
    else:
        num_calcs_at_once = parameters['NumCalcsAtOnce']

    objects_dict['num_calcs_at_once'] = num_calcs_at_once

    # get the run title
    if 'RunTitle' not in parameters:
        run_dir_name = 'garun'
    elif parameters['RunTitle'] in (None, 'default'):
        run_dir_name = 'garun'
    else:
        run_dir_name = 'garun_' + str(parameters['RunTitle'])

    objects_dict['run_dir_name'] = run_dir_name

    # make the energy calculator
    energy_calculator = make_energy_calculator(parameters, geometry,
                                               composition_space)
    objects_dict['energy_calculator'] = energy_calculator

    # make the stopping criteria
    stopping_criteria = make_stopping_criteria(parameters, composition_space)
    objects_dict['stopping_criteria'] = stopping_criteria

    # determine which variations should have non-zero default fractions
    do_permutation = False
    if len(composition_space.get_all_swappable_pairs()) > 0:
        do_permutation = True

    # get the number of atoms per composition
    if len(composition_space.endpoints) > 1:
        atoms_per_comp = len(composition_space.endpoints)
    else:
        atoms_per_comp = \
            composition_space.endpoints[0].reduced_composition.num_atoms

    # see if numstoichmut can be done w/o violating the min or max number of
    # atoms constraints
    do_atomsmut = False
    if composition_space.objective_function == 'pd':
        if constraints.min_num_atoms != constraints.max_num_atoms:
            do_atomsmut = True
    elif composition_space.objective_function == 'epa':
        bottom = int(math.ceil(constraints.min_num_atoms/atoms_per_comp))
        top = int(math.floor(constraints.max_num_atoms/atoms_per_comp))
        if top > bottom:
            do_atomsmut = True

    # set default fractions for the variations
    default_variation_fractions = {}
    if do_permutation and do_atomsmut:
        default_variation_fractions['permutation'] = 0.1
        default_variation_fractions['num_atoms_mut'] = 0.1
        default_variation_fractions['structure_mut'] = 0.1
        default_variation_fractions['mating'] = 0.7
    elif not do_permutation and do_atomsmut:
        default_variation_fractions['permutation'] = 0.0
        default_variation_fractions['num_atoms_mut'] = 0.1
        default_variation_fractions['structure_mut'] = 0.1
        default_variation_fractions['mating'] = 0.8
    elif do_permutation and not do_atomsmut:
        default_variation_fractions['permutation'] = 0.1
        default_variation_fractions['num_atoms_mut'] = 0.0
        default_variation_fractions['structure_mut'] = 0.1
        default_variation_fractions['mating'] = 0.8
    elif not do_permutation and not do_atomsmut:
        default_variation_fractions['permutation'] = 0.0
        default_variation_fractions['num_atoms_mut'] = 0.0
        default_variation_fractions['structure_mut'] = 0.2
        default_variation_fractions['mating'] = 0.8

    # make the variations
    variations_list = make_variations(parameters, default_variation_fractions,
                                      composition_space)

    # check that at least one variation has been used
    if len(variations_list) == 0:
        print('At least one variation must be used. Either leave entire '
              '"Variations" block blank to use default variations, or specify '
              'at least one variation within the "Variations" block.')
        print('Quitting...')
        quit()

    # check that the variations' fraction variables sum to 1
    frac_sum = 0.0
    for variation in variations_list:
        frac_sum = frac_sum + variation.fraction
    if frac_sum < 0.999 or frac_sum > 1.001:
        print("The Variations' fraction values must sum to 1.")
        print('Quitting...')
        quit()

    objects_dict['variations'] = variations_list

    # make the pool, selection, and composition fitness weight
    if 'Pool' not in parameters:
        pool = population.Pool(None, composition_space, run_dir_name)
    else:
        if 'num_promoted' in parameters['Pool']:
            if parameters['Pool']['num_promoted'] < 1:
                print('At least one organism must be promoted in the Pool.')
                print('Quitting...')
                quit()
            else:
                pool = population.Pool(parameters['Pool'], composition_space,
                                       run_dir_name)
        else:
            pool = population.Pool(parameters['Pool'], composition_space,
                                   run_dir_name)

    if 'Selection' not in parameters:
        selection = general.SelectionProbDist(None, pool.size)
    else:
        selection = general.SelectionProbDist(parameters['Selection'],
                                              pool.size)

    if 'CompositionFitnessWeight' not in parameters:
        comp_fitness_weight = general.CompositionFitnessWeight(None)
    else:
        if 'max_weight' in parameters['CompositionFitnessWeight']:
            if parameters['CompositionFitnessWeight']['max_weight'] < 0 or \
                    parameters['CompositionFitnessWeight']['max_weight'] > 1:
                print('The maximum weight of the composition fitness must lie'
                      ' in the interval [0,1].')
                print('Please change the value passed to the "max_weight" '
                      'keyword in the CompositionFitnessWeight block.')
                print('Quitting...')
                quit()
            else:
                comp_fitness_weight = general.CompositionFitnessWeight(
                    parameters['CompositionFitnessWeight'])
        else:
            comp_fitness_weight = general.CompositionFitnessWeight(
                parameters['CompositionFitnessWeight'])

    pool.selection = selection
    pool.comp_fitness_weight = comp_fitness_weight
    objects_dict['pool'] = pool

    job_specs = {}
    if 'JobSpecs' in parameters:
        job_specs = parameters['JobSpecs']

    if 'cores' in job_specs:
        if job_specs['cores'] > 8:
            print('Using max. default cpus_per_task: 8')
            job_specs['cores'] = 8
    else:
        # default cpus_per_task for a worker
        job_specs['cores'] = 1

    if 'walltime' not in job_specs:
        print('Using default wall time of 24:00:00')
        job_specs['walltime'] = '24:00:00'

    if 'memory' not in job_specs:
        # default job memory for a worker
        print('Using default total memory of 8GB per worker')
        job_specs['memory'] = '8GB'

    if 'project' not in job_specs:
        print('Please specify the project "#SBATCH -A" tag for worker.')
        quit()

    if 'queue' not in job_specs:
        print('Specify queue option in job_specs if the SLURM/PBS cluster'
              ' requires it. Otherwise ignore..')

    if 'interface' not in job_specs:
        print('Using default interface \'ib0\' (infiniband nodes)')
        job_specs['interface'] = 'ib0'

    objects_dict['job_specs'] = job_specs

    return objects_dict


def get_lat_match_params(parameters):
    """
    Returns the lattice matching constraints from input yaml files as a
    dictionary

    Sets default constraints if none present
    """

    lat_match_params = {}
    keys = ['max_area', 'max_mismatch', 'max_angle_diff', 'r1r2_tol',
            'separation', 'align_random', 'nlayers_substrate',
            'nlayers_2d']
    # all defaults
    match_constraints = [100, 0.05, 2, 0.06, 2, True, 1, 1, 1]

    # Set all default values to the keys
    for key, param in zip(keys, match_constraints):
        lat_match_params[key] = param

    # Check if any user specified constraints are available and replace
    if 'LatticeMatch' not in parameters:
        print('No lattice match constraints provided. Using defaults...')
    else:
        for key in keys:
            if key in parameters['LatticeMatch']:
                lat_match_params[key] = parameters['LatticeMatch'][key]

    return lat_match_params


def get_substrate_params(parameters):
    """
    Returns the total energy (enthalpy) of primitive substrate calculation
    and number of atoms in the primitive substrate cell

    This is mandatory and no defaults!
    """
    if 'Substrate' in parameters:
        sub_params = parameters['Substrate']
        return sub_params
    else:
        return None


def make_organism_creators(parameters, composition_space, constraints):
    """
    Returns a list containing organism creator objects.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        composition_space: the CompositionSpace of the search

        constraints: the Constraints of the search
    """

    if 'InitialPopulation' not in parameters:
        return make_default_organism_creator(composition_space, constraints)
    elif parameters['InitialPopulation'] in (None, 'default'):
        return make_default_organism_creator(composition_space, constraints)
    # make the specified creators
    else:
        # check that at least one valid option is used
        if 'random' not in parameters['InitialPopulation'] and 'from_files' \
                not in parameters['InitialPopulation']:
            print('At least one valid option for making structures for the '
                  'initial population must be provided.')
            print('Please use the "random" and/or "from_files" keywords in '
                  'the InitialPopulation block.')
            print('Quitting...')
            quit()

        initial_organism_creators = []

        # the random organism creator
        if 'random' in parameters['InitialPopulation']:
            random_organism_creator = organism_creators.RandomOrganismCreator(
                parameters['InitialPopulation']['random'], composition_space,
                constraints)
            initial_organism_creators.append(random_organism_creator)

        # the from files organism creator
        if 'from_files' not in parameters['InitialPopulation']:
            if composition_space.objective_function == 'pd':
                print('For phase diagram searches, reference structures at '
                      'each endpoint of the composition space must be '
                      'provided.')
                print('Please use the "from_files" keyword in the '
                      'InitialPopulation block to provide the reference '
                      'structures.')
                print('Quitting...')
                quit()
        elif parameters['InitialPopulation']['from_files'] is None:
            print('The path to the folder containing the files must be '
                  'provided. Please use the "path_to_folder" keyword.')
            print('Quitting...')
            quit()
        elif 'path_to_folder' not in parameters['InitialPopulation'][
                'from_files']:
            print('Incorrect keyword given after "from_files" in the '
                  'InitialPopulation block. Please use the "path_to_folder" '
                  'keyword.')
            print('Quitting...')
            quit()
        else:
            given_path = parameters['InitialPopulation']['from_files'][
                'path_to_folder']
            if given_path is None:
                print('The path to the folder containing the files for the '
                      'initial population must be provided. Please give the '
                      'path after the "path_to_folder" keyword.')
                print('Quitting...')
                quit()
            elif not os.path.exists(given_path):
                print('The given folder containing structures for the initial '
                      'population does not exist.')
                print('Quitting...')
                quit()
            elif len([f for f in os.listdir(given_path) if
                      os.path.isfile(os.path.join(given_path, f))]) == 0:
                print('The given folder containing structures for the initial '
                      'population does not contain any files.')
                print('Quitting...')
                quit()
            else:
                files_organism_creator = organism_creators.FileOrganismCreator(
                    given_path)

                # check that the files cover all composition space endpoints
                if composition_space.objective_function == 'pd':
                    cells = files_organism_creator.get_cells()

                    print("\n===== DEBUG: endpoint coverage check =====")
                    print("given_path =", given_path)
                    print("num parsed cells =", len(cells))
                    print("composition_space.objective_function =",
                          composition_space.objective_function)

                    print("\nFiles found in folder:")
                    try:
                        for fname in sorted(os.listdir(given_path)):
                            full = os.path.join(given_path, fname)
                            print("  {!r} | isfile={} | size={}".format(
                                fname,
                                os.path.isfile(full),
                                os.path.getsize(full) if os.path.isfile(full) else "NA"
                            ))
                    except Exception as e:
                        print("  Could not list folder contents:", repr(e))

                    print("\nExpected endpoints:")
                    for i, endpoint in enumerate(composition_space.endpoints):
                        try:
                            print("  endpoint[{}] = {} | reduced = {}".format(
                                i, endpoint, endpoint.reduced_composition))
                        except Exception as e:
                            print("  endpoint[{}] could not print: {}".format(
                                i, repr(e)))

                    print("\nParsed cells:")
                    for i, cell in enumerate(cells):
                        try:
                            print("  cell[{}] formula = {} | reduced = {}".format(
                                i,
                                cell.composition.formula,
                                cell.composition.reduced_composition))
                        except Exception as e:
                            print("  cell[{}] could not print: {}".format(
                                i, repr(e)))

                    provided_endpoints = []
                    for endpoint in composition_space.endpoints:
                        for cell in cells:
                            try:
                                is_match = cell.composition.reduced_composition.almost_equals(
                                    endpoint.reduced_composition)
                            except Exception as e:
                                print("  MATCH CHECK FAILED for cell {} and endpoint {}: {}".format(
                                    getattr(cell.composition, "formula", "UNKNOWN"),
                                    endpoint,
                                    repr(e)))
                                is_match = False

                            if is_match and endpoint not in provided_endpoints:
                                print("MATCH: cell {} matches endpoint {}".format(
                                    cell.composition.reduced_composition,
                                    endpoint.reduced_composition))
                                provided_endpoints.append(endpoint)

                    print("\nProvided endpoints:")
                    if len(provided_endpoints) == 0:
                        print("  NONE")
                    else:
                        for ep in provided_endpoints:
                            try:
                                print(" ", ep, "| reduced =", ep.reduced_composition)
                            except Exception as e:
                                print("  could not print provided endpoint:", repr(e))

                    print("===== END DEBUG =====\n")

                    for endpoint in composition_space.endpoints:
                        if endpoint not in provided_endpoints:
                            print('Error: valid structure files not provided '
                                  'to the initial population for all '
                                  'endpoints of the composition space.')
                            print('Quitting...')
                            quit()

                initial_organism_creators.append(files_organism_creator)

        return initial_organism_creators


def make_default_organism_creator(composition_space, constraints):
    """
    Returns a list containing a RandomOrganismCreator, or quits.

    Args:
        composition_space: the CompositionSpace of the search

        constraints: the Constraints of the search
    """

    if composition_space.objective_function == 'pd':
        print('For phase diagram searches, reference structures at each '
              'endpoint of the composition space must be provided in the '
              'initial population.')
        print('Please use the "from_files" keyword in the '
              'InitialPopulation block to provide the reference '
              'structures.')
        print('Quitting...')
        quit()
    else:
        random_organism_creator = organism_creators.RandomOrganismCreator(
            'default', composition_space, constraints)
        return [random_organism_creator]


def make_energy_calculator(parameters, geometry, composition_space):
    """
    Returns an EnergyCode object corresponding to which energy code was
    specified in the input file. Quits if an energy code object cannot be made.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        geometry: the Geometry for the search

        composition_space: the CompositionSpace of the search
    """

    if 'EnergyCode' not in parameters:
        print('A method for calculating energy must be provided. Please use '
              'the "EnergyCode" keyword.')
        print('Quitting...')
        quit()
    elif parameters['EnergyCode'] is None:
        print('An energy code must be specified after the "EnergyCode" '
              'keyword.')
        print('Quitting...')
        quit()
    elif 'gulp' in parameters['EnergyCode']:
        return make_gulp_energy_calculator(parameters, geometry)
    elif 'lammps' in parameters['EnergyCode']:
        return make_lammps_energy_calculator(parameters, geometry)
    elif 'vasp' in parameters['EnergyCode']:
        return make_vasp_energy_calculator(parameters, composition_space,
                                           geometry)
    else:
        print('The given energy code name is invalid.')
        print('Quitting...')
        quit()


def make_gulp_energy_calculator(parameters, geometry):
    """
    Returns a GulpEnergyCalculator object, or quits if one cannot be made.
    """

    if parameters['EnergyCode']['gulp'] is None:
        print('No GULP header or potential files given. Please use the '
              '"header_file" and "potential_file" keywords.')
        print('Quitting...')
        quit()
    else:
        if 'header_file' not in parameters['EnergyCode']['gulp']:
            print('A GULP header file must be provided. Please use the '
                  '"header_file" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['gulp']['header_file'] is None:
            print('No GULP header file given after the "header_file" '
                  'keyword. Please provide one.')
            print('Quitting...')
            quit()
        else:
            header_file_path = parameters['EnergyCode']['gulp'][
                'header_file']
            if not os.path.exists(header_file_path):
                print('The given GULP header file does not exist.')
                print('Quitting...')
                quit()

        if 'potential_file' not in parameters['EnergyCode']['gulp']:
            print('A GULP potential file must be provided. Please use the '
                  '"potential_file" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['gulp']['potential_file'] is None:
            print('No GULP potential file given after the '
                  '"potential_file" keyword. Please provide one.')
            print('Quitting...')
            quit()
        else:
            potential_file_path = parameters['EnergyCode']['gulp'][
                'potential_file']
            if not os.path.exists(potential_file_path):
                print('The given GULP potential file does not exist.')
                print('Quitting...')
                quit()

        return energy_calculators.GulpEnergyCalculator(
            header_file_path, potential_file_path, geometry)


def make_lammps_energy_calculator(parameters, geometry):
    """
    Returns a LammpsEnergyCalculator object, or quits if one cannot be made.
    """

    if parameters['EnergyCode']['lammps'] is None:
        print('No LAMMPS input script given. Please use the "input_script" '
              'keyword.')
        print('Quitting...')
        quit()
    else:
        if 'input_script' not in parameters['EnergyCode']['lammps']:
            print('A LAMMPS input script must be provided. Please use the '
                  '"header_file" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['lammps']['input_script'] is None:
            print('No LAMMPS input script given after the "input_script" '
                  'keyword. Please provide one.')
            print('Quitting...')
            quit()
        else:
            input_script_path = parameters['EnergyCode']['lammps'][
                'input_script']
            if not os.path.exists(input_script_path):
                print('The given LAMMPS input script does not exist.')
                print('Quitting...')
                quit()

        return energy_calculators.LammpsEnergyCalculator(
            input_script_path, geometry)


def make_vasp_energy_calculator(parameters, composition_space, geometry):
    """
    Returns a VaspEnergyCalculator object, or quits if one cannot be made.
    """

    if parameters['EnergyCode']['vasp'] is None:
        print('No VASP input files given.')
        print('Quitting...')
        quit()
    else:
        if 'incar' not in parameters['EnergyCode']['vasp']:
            print('An INCAR file must be provided. Please use the "incar" '
                  'keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['vasp']['incar'] is None:
            print('No INCAR file was given after the "incar" keyword. Please '
                  'provide one.')
            print('Quitting...')
            quit()
        else:
            incar_path = parameters['EnergyCode']['vasp']['incar']
            if not os.path.exists(incar_path):
                print('The given INCAR file does not exist.')
                print('Quitting...')
                quit()

        if 'kpoints' not in parameters['EnergyCode']['vasp']:
            print('A KPOINTS file must be provided. Please use the '
                  '"kpoints" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['vasp']['kpoints'] is None:
            print('No KPOINTS file was given after the "kpoints" keyword. '
                  'Please provide one.')
            print('Quitting...')
            quit()
        else:
            kpoints_path = parameters['EnergyCode']['vasp']['kpoints']
            if not os.path.exists(kpoints_path):
                print('The given KPOINTS file does not exist.')
                print('Quitting...')
                quit()

        if 'potcars' not in parameters['EnergyCode']['vasp']:
            print('POTCAR file(s) must be provided. Please use the '
                  '"potcars" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['vasp']['potcars'] is None:
            print('No POTCAR files were given after the "potcars" keyword. '
                  'Please provide them.')
            print('Quitting...')
            quit()
        else:
            potcar_paths = parameters['EnergyCode']['vasp']['potcars']
            elements_list = composition_space.get_all_elements()
            if len(potcar_paths) < len(elements_list):
                print('Not enough POTCAR files provided - one must be '
                      'given for each element in the composition space. '
                      'Please provide them.')
                print('Quitting...')
                quit()
            for element in elements_list:
                if element.symbol not in potcar_paths:
                    print('No POTCAR file given for {}. Please provide '
                          'one.'.format(element.symbol))
                    print('Quitting...')
                    quit()
            for key in potcar_paths:
                if potcar_paths[key] is None:
                    print('No POTCAR file given for {}. Please provide '
                          'one.'.format(key))
                    print('Quitting...')
                    quit()
                elif not os.path.exists(potcar_paths[key]):
                    print('The POTCAR file given for {} does not '
                          'exist.'.format(key))
                    print('Quitting...')
                    quit()

        if 'num_submits_to_converge' in \
                parameters['EnergyCode']['vasp']:
            num_submits_to_converge = \
                parameters['EnergyCode']['vasp']['num_submits_to_converge']
            if not isinstance(num_submits_to_converge, int):
                print('Error: Parameter num_submits_to_converge '
                      'should be an integer.')
                print('Quitting...')
                quit()
        else:
            print('No "num_submits_to_converge" option given. '
                  'Using default maximum value of 2 VASP'
                  ' calculations done on an organism to converge.')
            num_submits_to_converge = 2

        num_rerelax = 0
        if 'num_rerelax' in parameters['EnergyCode']['vasp']:
            num_rerelax = parameters['EnergyCode']['vasp']['num_rerelax']
            if not isinstance(num_rerelax, int):
                print('Error: Parameter num_rerelax should be an integer.')
                print('Quitting...')
                quit()

        print('VASP calculations on a structure: \n'
              'for energy convergence: {}  and\n'
              'for (re-)relaxation: {}'.format(
                  num_submits_to_converge, num_rerelax))

        return energy_calculators.VaspEnergyCalculator(
            incar_path, kpoints_path,
            potcar_paths, geometry,
            num_submits_to_converge=num_submits_to_converge,
            num_rerelax=num_rerelax)


def make_stopping_criteria(parameters, composition_space):
    """
    Returns a StoppingCriteria object.
    """

    if 'StoppingCriteria' not in parameters:
        return general.StoppingCriteria(None, composition_space)
    elif parameters['StoppingCriteria'] in (None, 'default'):
        return general.StoppingCriteria(None, composition_space)
    elif 'found_structure' in parameters['StoppingCriteria']:
        if parameters['StoppingCriteria']['found_structure'] in (None,
                                                                 'default'):
            return general.StoppingCriteria(parameters['StoppingCriteria'],
                                            composition_space)
        else:
            given_path = parameters['StoppingCriteria']['found_structure']
            if not os.path.exists(given_path):
                print('The file containing the structure to find does not '
                      'exist.')
                print('Quitting...')
                quit()
            elif not (os.path.basename(given_path).endswith('.cif') or
                      os.path.basename(given_path).startswith('POSCAR.')):
                print('File containing structure to find must be in POSCAR or '
                      'cif format and begin with POSCAR. or end with .cif, '
                      'respectively.')
                print('Quitting...')
                quit()
            else:
                try:
                    Structure.from_file(given_path)
                    return general.StoppingCriteria(
                        parameters['StoppingCriteria'], composition_space)
                except ValueError:
                    print('Error reading the structure to find from the given '
                          'file.')
                    print('Quitting...')
                    quit()
    else:
        return general.StoppingCriteria(parameters['StoppingCriteria'],
                                        composition_space)


def make_variations(parameters, default_fractions, composition_space):
    """
    Creates the variations, using default parameter values if needed.
    """

    if 'Variations' not in parameters:
        return make_default_variations(default_fractions, composition_space)
    elif parameters['Variations'] in (None, 'default'):
        return make_default_variations(default_fractions, composition_space)
    else:
        variations_list = []

        if 'Mating' not in parameters['Variations']:
            pass
        elif parameters['Variations']['Mating'] is None:
            print('If the "Mating" keyword is used, its "fraction" keyword '
                  'must also be set.')
            print('Quitting...')
            quit()
        else:
            if parameters['Variations']['Mating']['fraction'] in (None,
                                                                  'default'):
                print('The "fraction" kwyword is not optional and must '
                      'contain a valid entry (between 0 and 1) for the Mating '
                      'variation.')
                print('Quitting...')
                quit()
            else:
                mating = variations.Mating(parameters['Variations']['Mating'])
                variations_list.append(mating)

        if 'StructureMut' not in parameters['Variations']:
            pass
        elif parameters['Variations']['StructureMut'] is None:
            print('If the "StructureMut" keyword is used, its "fraction" '
                  'keyword must also be set.')
            print('Quitting...')
            quit()
        else:
            if parameters['Variations']['StructureMut']['fraction'] in (
                    None, 'default'):
                print('The "fraction" keyword is not optional and must '
                      'contain a valid entry (between 0 and 1) for the '
                      'StructureMut variation.')
                print('Quitting...')
                quit()
            else:
                structure_mut = variations.StructureMut(
                    parameters['Variations']['StructureMut'])
                variations_list.append(structure_mut)

        if 'NumAtomsMut' not in parameters['Variations']:
            pass
        elif parameters['Variations']['NumAtomsMut'] is None:
            print('If the "NumAtomsMut" keyword is used, its "fraction" '
                  'keyword must also be set.')
            print('Quitting...')
            quit()
        else:
            if parameters['Variations']['NumAtomsMut']['fraction'] in (
                    None, 'default'):
                print('The "fraction" keyword is not optional and must '
                      'contain a valid entry (between 0 and 1) for the '
                      'NumAtomsMut variation.')
                print('Quitting...')
                quit()
            else:
                num_atoms_mut = variations.NumAtomsMut(
                    parameters['Variations']['NumAtomsMut'])
                variations_list.append(num_atoms_mut)

        if 'Permutation' not in parameters['Variations']:
            pass
        elif parameters['Variations']['Permutation'] is None:
            print('If the "Permutation" keyword is used, its "fraction" '
                  'keyword must also be set.')
            print('Quitting...')
            quit()
        else:
            if parameters['Variations']['Permutation']['fraction'] in (
                    None, 'default'):
                print('The "fraction" keyword is not optional and must '
                      'contain a valid entry (between 0 and 1) for the '
                      'Permutation variation.')
                print('Quitting...')
            else:
                permutation = variations.Permutation(
                    parameters['Variations']['Permutation'], composition_space)
                variations_list.append(permutation)

        return variations_list


def make_default_variations(default_fractions, composition_space):
    """
    Creates the variations with default parameter values and the provided
    default fractions.
    """

    variations_list = []
    mating = variations.Mating({'fraction': default_fractions['mating']})
    structure_mut = variations.StructureMut(
        {'fraction': default_fractions['structure_mut']})
    num_atoms_mut = variations.NumAtomsMut(
        {'fraction': default_fractions['num_atoms_mut']})
    permutation = variations.Permutation(
        {'fraction': default_fractions['permutation']},
        composition_space)
    variations_list.append(mating)
    variations_list.append(structure_mut)
    variations_list.append(num_atoms_mut)
    variations_list.append(permutation)
    return variations_list