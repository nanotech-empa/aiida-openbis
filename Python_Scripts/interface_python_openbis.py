from aiida_openbis.utils import bisutils
import yaml
import utils
import os

# Connect to openBIS
session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

option = "0"
while option!="6":

    option = input("1: Create new space\n2: Create new project\n3: Create new experiment\n4: Create new sample\n5: Save openBIS to YML file\n6: Close\nIntroduce the option: ")

    # TODO: Option creating it from yml file
    if option == "1":
        print("\n----Create a new space----\n")

        new_space_name = input("Introduce the new space's name: ")

        # Convert the name to a code that openBIS API can understand
        new_space_code = utils.convert_name_to_code(new_space_name)

        # Check if the space it is being created already exists in the DB
        new_space_not_exists = session.get_spaces(code = new_space_code).df.empty

        if new_space_not_exists:
            new_space_description = input("Introduce the new space's description: ")
            values = {'description': new_space_description}

            session.new_space(code=new_space_code, **values).save()
            print(f"\nThe space {new_space_name} was successfully created!\n")
        else:
            print(f"\nThe space {new_space_name} already exists!\n")

    # TODO: Option creating it from yml file
    if option == "2":
        print("\n----Create a new project----\n")

        # Print all the spaces' codes available in the database
        all_spaces_database = session.get_spaces()
        all_spaces_database_codes = list(all_spaces_database.df.code)

        for idx,code in enumerate(all_spaces_database_codes):
            print(f"{idx + 1}: {code}")

        selected_space_code_idx = input("\nSelect a space from the list: ")
        selected_space_code_idx = int(selected_space_code_idx)
        selected_space_code = all_spaces_database_codes[selected_space_code_idx - 1]

        new_project_name = input("Introduce the new project's name: ")

        # Convert the name to a code that openBIS API can understand
        new_project_code = utils.convert_name_to_code(new_project_name)

        # Check if the project it is being created already exists in the DB
        new_project_not_exists = session.get_projects(space = selected_space_code,
                                                      code = new_project_code).df.empty

        if new_project_not_exists:
            new_project_description = input("Introduce the new project's description: ")
            values = {'space': selected_space_code, 'description': new_project_description}

            session.new_project(code=new_project_code, **values).save()
            print(f"The project {new_project_name} was successfully created!")
        else:
            print(f"The project {new_project_name} already exists!")
    
    # TODO: Option creating it from yml file
    if option == "3":
        print("\n----Create a new experiment----\n")

        # Print all the spaces' codes available in the database
        all_spaces_database = session.get_spaces()
        all_spaces_database_codes = list(all_spaces_database.df.code)

        for idx, code in enumerate(all_spaces_database_codes):
            print(f"{idx + 1}: {code}")

        selected_space_code_idx = input("\nSelect a space from the list: ")
        selected_space_code_idx = int(selected_space_code_idx)
        selected_space_code = all_spaces_database_codes[selected_space_code_idx - 1]

        # Print all the project' codes available in the database
        all_projects_space = session.get_projects(space = selected_space_code)
        all_projects_space_codes = list(all_projects_space.df.code)

        for idx, code in enumerate(all_projects_space_codes):
            print(f"{idx + 1}: {code}")

        selected_project_code_idx = input("\nSelect a project from the list: ")
        selected_project_code_idx = int(selected_project_code_idx)
        selected_project_code = all_projects_space_codes[selected_project_code_idx - 1]

        # Print all the experiment types available in the database
        all_experiment_types = session.get_experiment_types()
        all_experiment_types_codes = list(all_experiment_types.df.code)

        for idx, code in enumerate(all_experiment_types_codes):
            print(f"{idx + 1}: {code}")

        selected_experiment_type_code_idx = input("\nSelect an experiment type from the list: ")
        selected_experiment_type_code_idx = int(selected_experiment_type_code_idx)
        selected_experiment_type_code = all_experiment_types_codes[selected_experiment_type_code_idx - 1]

        new_experiment_name = input("Introduce the new experiment's name: ")

        # Convert the name to a code that openBIS API can understand
        new_experiment_code = utils.convert_name_to_code(new_experiment_name)

        # Check if the experiment it is being created already exists in the DB
        new_experiment_not_exists = session.get_experiments(space = selected_space_code,
                                                            project = selected_project_code,
                                                            code = new_experiment_code).df.empty

        if new_experiment_not_exists:
            new_experiment_description = input("Introduce the new experiment's description: ")
            values = {'project': f'/{selected_space_code}/{selected_project_code}', 
                      'type': selected_experiment_type_code, 
                      'description': new_experiment_description}

            session.new_experiment(code=new_experiment_code, **values).save()
            print(f"The experiment {new_experiment_name} was successfully created!")
        else:
            print(f"The experiment {new_experiment_name} already exists!")

    if option=="5":

        # Retrieve all spaces available in DB
        all_spaces = session.get_spaces()

        # Verify whether an YML file with the openBIS schema already exists
        openbis_schema_filepath = '/home/jovyan/work/aiida-openbis/Python_Scripts/openbis_schema.yml'
        openbis_schema_file_exists = os.path.isfile(openbis_schema_filepath)

        if openbis_schema_file_exists:
            all_openbis_dict = utils.load_yml_file(openbis_schema_filepath, 'r')
        else:
            all_openbis_dict = {'spaces': {}}

        for space in all_spaces:
            if space.code not in all_openbis_dict['spaces']:
                space_code = space.code
                space_description = space.description
                all_openbis_dict['spaces'][space_code] = {'description':space_description}

        # Retrieve all projects available in DB
        all_projects = session.get_projects()

        for project in all_projects:
            
            _, project_space, project_code = project.identifier.split("/")

            project_in_yml = False
            if 'projects' in all_openbis_dict:
                if project_code in all_openbis_dict['projects']:
                    if project_space == all_openbis_dict['projects'][project_code]['space']:
                        project_in_yml = True
            else:
                all_openbis_dict['projects'] = {}
            
            if project_in_yml == False:
                all_openbis_dict['projects'][project_code] = {'space': project_space, 'description': project.description}
        
        # Retrieve all experiments available in DB
        all_experiments = session.get_experiments()

        for experiment in all_experiments:
            
            _, experiment_space, experiment_project, experiment_code = experiment.identifier.split("/")
            experiment_space_project = f'/{experiment_space}/{experiment_project}'

            experiment_in_yml = False
            if 'experiments' in all_openbis_dict:
                if experiment_code in all_openbis_dict['experiments']:
                    if experiment_space_project == all_openbis_dict['experiments'][experiment_code]['project']:
                            experiment_in_yml = True
            else:
                all_openbis_dict['experiments'] = {}
            
            if experiment_in_yml == False:
                all_openbis_dict['experiments'][experiment_code] = {'space': experiment_space, 
                                                                    'project': experiment_project, 
                                                                    'type': experiment.type.code,
                                                                    'description': experiment.description}

        filepath = '/home/jovyan/work/aiida-openbis/Python_Scripts/openbis_schema.yml'
        mode = 'w'
        utils.save_yml_file(filepath, mode, all_openbis_dict)
        print(f"OpenBIS schema was successfully downloaded!")



# # Create new projects
# all_projects_database = session.get_projects()
# for key, value in items['projects'].items():
#     new_project_code = key.upper()
#     new_project_exists = utils.scanning_element_database(all_projects_database, new_project_code)
#     if new_project_exists == False:
#         print(new_project_code)
#         session.new_project(code=new_project_code, **value).save()

# # Create new experiments
# all_experiments_database = session.get_experiments()
# for key, value in items['experiments'].items():
#     new_experiment_code = key.upper()
#     new_experiment_exists = utils.scanning_element_database(all_experiments_database, new_experiment_code)
#     if new_experiment_exists == False:
#         print(new_experiment_code)
#         session.new_experiment(code=new_experiment_code, **value).save()