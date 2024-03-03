class ProjectException(Exception):
    """Generic Exception used by this project."""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)
