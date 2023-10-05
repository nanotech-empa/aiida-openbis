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