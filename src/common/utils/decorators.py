import logging
import sys


def handle_config_parser_exception(message):
    """Decorator for logging errors accessing configuration parameters.

    :param message: beginning of the logging message.
    """
    logger = logging.getLogger()

    def decorate(func):
        def wrapper(self, *args, **kwargs):
            try:
                return func(self, *args, **kwargs)
            except KeyError as e:
                logger.error(message + str(e))
                sys.exit()
        return wrapper
    return decorate
