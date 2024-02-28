from apples.OptionsBasic import OptionsBasic


def options_config():
    """
    A function to configure options.
    """
    parser = OptionsBasic('APPLES database')
    (options, args) = parser.parse()
    return options, args
