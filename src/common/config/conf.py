import json
import os

from src.common.singleton import Singleton


class Conf(metaclass=Singleton):
    """This class holds all the methods that interact with the configuration file."""

    def __init__(self, file_path):
        file_path = os.path.abspath(file_path)
        with open(file_path) as config:
            self.conf = json.load(config)
