def scanning_object_types(all_object_types,new_object_code):
    """
    Verifies if the object type is already created in openBIS.

        Parameters:
                all_object_types(pybis.things.Things): All object types that alread exist in the openBIS architecture.
                new_object_code(string): Object code that is being created.
        
        Returns:
                object_exists(bool): Boolean representing if the object type is already created.
    """

    for _, object_type in enumerate(all_object_types):
        
