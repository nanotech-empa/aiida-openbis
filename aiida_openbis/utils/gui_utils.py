def convert_name_to_code(name: str) -> str:
    """Convert the element name to a code that can be interpreted by the openBIS API

    Args:
        name (str): Attributed name

    Returns:
        str: openBIS code
    """
    code = name.upper()
    code = code.replace(" ","_")
    return code