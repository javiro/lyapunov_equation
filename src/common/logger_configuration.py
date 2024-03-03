# Stdlib imports
import json
import logging
import logging.config
from os import makedirs, path

# Project imports
from src.common.singleton import Singleton

# Disable unneeded logging stuff, although they are not logged, they are calculated
# Disable logging current thread
logging.logThreads = 0
# Disable logging PID
logging.logProcesses = 0
# Disable logging module, function or file line
logging._srcfile = None


class LoggerManager(metaclass=Singleton):
    """Utility singleton class to configure loggers (only one for now):
    - one for the application traces

    This class should be initialized when launching the script. After that, just use the regular python logging system:

    To count on an application logger:

        >>> import logging
        >>> logger = logging.getLogger(__name__)

    """

    def __init__(self):
        """Logger initializer"""

        with open("logging_conf.json", "r") as logging_configuration_file:
            config_dict = json.load(logging_configuration_file)

        logger_handlers = ["file", "metrics_handler"]
        # Create default directories for log files
        for handler in logger_handlers:
            makedirs(path.dirname(config_dict["handlers"][handler]["filename"]), exist_ok=True)

        # Config logger
        logging.config.dictConfig(config_dict)
