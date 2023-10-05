def scanning_element_database(all_elements, new_element_code):
    """Verifies if the element is already created in openBIS.

    Args:
        all_elements(pybis.things.Things): All elements that alread exist in the openBIS database.
        new_element_code(str): Element code that is being created.

    Returns:
        bool: Boolean representing if the element is already created.
    """

    for _, element in enumerate(all_elements):
        element_code = element.code
        if element_code == new_element_code:
            return True
    return False

def create_element(session,items,element_tag):
    if element_tag == 'spaces':
        all_elements_database = session.get_spaces()
    elif element_tag == 'projects':
        all_elements_database = session.get_projects()
    elif element_tag == 'experiments':
        all_elements_database = session.get_experiments()
    else:
        all_elements_database = session.get_samples()

    element_exists = False
    for key, value in items[element_tag].items():
        
        if element_tag in ['spaces', 'projects', 'experiments']:
            new_element_code = key.upper()
            element_exists = scanning_element_database(all_elements_database, new_element_code)

        if element_exists == False:
            if element_tag =='spaces':
                session.new_space(code=new_element_code,**value).save()
            elif element_tag == 'projects':
                session.new_project(code=new_element_code,**value).save()
            elif element_tag == 'experiments':
                session.new_experiment(code=new_element_code,**value).save()
            else:
                session.new_sample(**value).save()