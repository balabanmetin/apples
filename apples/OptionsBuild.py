from apples.OptionsBasic import OptionsBasic


def options_config():
    parser = OptionsBasic("APPLES database")
    (options, args) = parser.parse()
    return options, args
