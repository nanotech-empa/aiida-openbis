def scanning_object_types(all_object_types, new_object_code):
    """
    Verifies if the object type is already created in openBIS.

        Parameters:
                all_object_types(pybis.things.Things): All object types that alread exist in the openBIS architecture.
                new_object_code(string): Object code that is being created.
        
        Returns:
                object_exists(bool): Boolean representing if the object type is already created.
    """

    object_exists = False
    for _, object_type in enumerate(all_object_types):
        object_type_code = object_type.code
        if object_type_code == new_object_code:
            object_exists = True
            return object_exists
    return object_exists

def scanning_property_types(all_property_types, new_property_code):
    """
    Verifies if the property type is already created in openBIS.

        Parameters:
                all_property_types(pybis.things.Things): All property types that alread exist in the openBIS architecture.
                new_property_code(string): Property code that is being created.
        
        Returns:
                property_exists(bool): Boolean representing if the property type is already created.
    """

    property_exists = False
    for _, property_type in enumerate(all_property_types):
        property_type_code = property_type.code
        if property_type_code == new_property_code:
            property_exists = True
            return property_exists
    return property_exists
