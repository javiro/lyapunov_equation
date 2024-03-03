import os

from src.common.conf import Conf
from src.common.constants import WORKSPACE, OUTPUTS
from src.common.singleton import Singleton


class Workspace(metaclass=Singleton):

    def __init__(self):
        """Workspace."""
        self.root = WORKSPACE
        if not os.path.exists(self.root):
            os.mkdir(self.root)
        if not os.path.exists(os.path.join(self.root, OUTPUTS)):
            os.mkdir(os.path.join(self.root, OUTPUTS))
        self.conf_path = os.path.join("common", "config", "default_conf.json")
        self.conf = Conf(self.conf_path)
