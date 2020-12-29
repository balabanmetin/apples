from apples.options_basic import options_basic


def options_config():
    parser = options_basic("applespkg")
    (options, args) = parser.parse()
    return (options, args)
