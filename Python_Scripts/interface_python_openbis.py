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
        all_spaces_database = session.get_spaces()
        new_space_exists = utils.scanning_element_database(all_spaces_database, new_space_code)

        if new_space_exists == False:
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
        all_projects_database = session.get_projects()
        new_project_exists = utils.scanning_element_database(all_projects_database, new_project_code)

        if new_project_exists == False:
            new_project_description = input("Introduce the new project's description: ")
            values = {'space': selected_space_code, 'description': new_project_description}

            session.new_project(code=new_project_code, **values).save()
            print(f"The project {new_project_name} was successfully created!")
        else:
            print(f"The project {new_project_name} already exists!")

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
            if space.registrator!='system' and space.code not in all_openbis_dict['spaces']:
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
            
            if project_in_yml == False and project.registrator!='system':
                all_openbis_dict['projects'][project_code] = {'space': project_space, 'description': project.description}

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