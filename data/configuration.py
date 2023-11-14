class Configuration:
    def __init__(self):
        # Initialize with configuration settings or parameters
        self.algorithm_type = ""
        self.settings = {}

    def set_algorithm_type(self, algorithm_type):
        self.algorithm_type = algorithm_type

    def get_algorithm_type(self):
        return self.algorithm_type

    def set_setting(self, key, value):
        self.settings[key] = value

    def get_setting(self, key):
        return self.settings.get(key)